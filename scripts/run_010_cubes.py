"""Attempt 010, stage C ("cubes"): cube decompositions of the S3 splits nearest the 006 scheme.

Each cube is S3CubeCNF(split, cube, gran, use_w) (matmul.s3_cubes). CaDiCaL writes a DRAT proof while it
solves; an UNSAT answer is checked with drat-trim and the proof is deleted once VERIFIED. A SAT model is
decoded, verified mod 3, checked for S3 invariance and saved.

usage:
  python3 -I scripts/run_010_cubes.py calib  <cap> <cfg>...        rank 23, (1,5,0,2): the 006 scheme's own cube
  python3 -I scripts/run_010_cubes.py split  a b c d <cap> <cfg> [<max_cubes>]   every cube of a split
  python3 -I scripts/run_010_cubes.py one    a b c d <cap> <cfg> <cube index>...
  python3 -I scripts/run_010_cubes.py subw   a b c d <cap> <cfg> <cube index> <w index>...   depth-2 samples
cfg is gran[+w][+pin|+sk][+n1], e.g. r1+pin, r1+sk+n1 (see matmul.s3_cubes.parse_cfg). At most 3 solver processes run at once.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import expand, split_orbits  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.s3 import is_s3_invariant, s3_normal_form, s3_split, scheme_s3_split  # noqa: E402
from matmul.s3_cubes import cubes_cfg, encoder_cfg, normalize_cfg  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402
from matmul.verify import verify  # noqa: E402

WORK = ROOT / "scratch" / "010" / "cubes"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
LOG = ROOT / "results" / "logs" / "010-s3-f3" / "cubes.log"
ATTEMPT = "010-s3-f3"
S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"
WORKERS = 3
SEED0 = 120001


def sym_from_index(k):
    """The k-th symmetric 3x3 matrix over F3 (0 <= k < 729), entries (00, 01, 02, 11, 12, 22) in base 3."""
    e = [(k // 3 ** i) % 3 for i in range(6)]
    return [[e[0], e[1], e[2]], [e[1], e[3], e[4]], [e[2], e[4], e[5]]]


def run_cube(split, cube, idx, cfg, seed, cap, tag, w0=None):
    """w0: optional index of a symmetric matrix pinned as w of the first size-3 slot (depth-2 subcube)."""
    a, b, c, d = split
    r = 6 * a + 3 * b + 2 * c + d
    name = f"c_{a}{b}{c}{d}_{cfg.replace('+', '_')}_{tag}{idx}_s{seed}"
    t0 = time.time()
    if w0 is not None:
        name += f"_w{w0}"
    enc = encoder_cfg(split, cube, cfg)
    if w0 is not None:
        for z, x in zip(enc.orb[2 * a][2], (e for r in sym_from_index(w0) for e in r)):
            for lit in enc._is(z, x):
                enc.clauses.append([lit])
    WORK.mkdir(parents=True, exist_ok=True)
    cnf = WORK / f"{name}.cnf"
    cnf.write_text(enc.dimacs())
    proof = WORK / f"{name}.drat"
    out = dict(s3_split=list(split), s3_rank_max=r, cfg=cfg, cube_index=idx, w0_index=w0, cube=json.loads(json.dumps(cube)),
               solver="cadical", seed=seed, cap_seconds=cap, vars=enc.nv, clauses=len(enc.clauses),
               build_s=round(time.time() - t0, 2))
    t1 = time.time()
    try:
        res = subprocess.run(["cadical", "-q", f"--seed={seed}", str(cnf), str(proof)], capture_output=True,
                             text=True, timeout=cap)
        status = {10: "SAT", 20: "UNSAT"}.get(res.returncode, f"exit{res.returncode}")
    except subprocess.TimeoutExpired:
        status = "TIMEOUT"
    out["seconds"] = round(time.time() - t1, 2)
    out["status"] = status
    if status == "SAT":
        model = [int(x) for line in res.stdout.splitlines() if line.startswith("v") for x in line.split()[1:]]
        orbits, fixed = enc.decode(model)
        s = expand(3, orbits, fixed, {"attempt": ATTEMPT, "instance": name, "field": "GF(3)"})
        out.update(found_rank=s.rank, verified_mod3=verify(s, mod=3).ok, s3_invariant=is_s3_invariant(s))
        try:
            out["s3_split_found"] = list(scheme_s3_split(s))
        except Exception as e:  # noqa: BLE001
            out["s3_split_found"] = f"error: {e}"
        out["class_key"] = class_key(s, 3)
        FOUND.mkdir(parents=True, exist_ok=True)
        path = FOUND / f"010C_3x3x3_rank{s.rank}_{name}.json"
        s.save(path)
        out["scheme_file"] = str(path.relative_to(ROOT))
    if status == "UNSAT":
        out["proof_bytes"] = proof.stat().st_size
        t2 = time.time()
        chk = max(300, int(4 * out["seconds"]) + 60)
        try:
            r3 = subprocess.run([str(DRAT_TRIM), str(cnf), str(proof), "-t", str(chk)], capture_output=True, text=True,
                                timeout=chk + 30)
            out["drat_trim"] = "VERIFIED" if "s VERIFIED" in r3.stdout else (
                [l for l in r3.stdout.splitlines() if l.startswith("s ")] or ["no verdict"])[0]
        except subprocess.TimeoutExpired:
            out["drat_trim"] = "TIMEOUT"
        out["check_seconds"] = round(time.time() - t2, 2)
    if proof.exists() and (status != "UNSAT" or out.get("drat_trim") == "VERIFIED"):
        proof.unlink()
    cnf.unlink(missing_ok=True)  # reproducible from the split, cube, cfg
    return out


def record(res, stage_note):
    cap = res["cap_seconds"]
    record_run(RUNS, attempt=ATTEMPT, stage="C", format="3x3x3", field="GF(3)",
               method=f"SAT cadical on S3-invariant Brent equations mod 3 (native encoder), S3 split "
                      f"{tuple(res['s3_split'])}, cube {res['cube_index']} of decomposition {res['cfg']} ({stage_note})",
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3") and res.get("s3_invariant")) if res["status"] == "SAT"
               else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": WORKERS, "threads_each": 1}, result=res)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f:
        f.write(json.dumps(dict(res, note=stage_note)) + "\n")
    print(json.dumps(res), flush=True)


def run_jobs(jobs, note):
    """jobs: (split, cube, idx, cfg, seed, cap, tag). Stops early on a rank-22 SAT."""
    results = []
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(run_cube, *j): j for j in jobs}
        stop = False
        for f in cf.as_completed(futs):
            res = f.result()
            record(res, note)
            results.append(res)
            if res["status"] == "SAT" and res["s3_rank_max"] <= 22 and not stop:
                stop = True
                for g in futs:
                    g.cancel()
                print("RANK-22 SAT: stopping", flush=True)
    return results


def known_cube(cfg):
    orbits, fixed = s3_normal_form(*s3_split(*split_orbits(Scheme.load(S3_006))))
    _, _, cube, _ = normalize_cfg(orbits, fixed, (1, 5, 0, 2), cfg)
    return cube


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "calib":
        cap = float(sys.argv[2])
        cfgs = sys.argv[3:]
        jobs = []
        for k, cfg in enumerate(cfgs):
            cube = known_cube(cfg)
            idx = cubes_cfg((1, 5, 0, 2), cfg).index(cube)
            jobs.append(((1, 5, 0, 2), cube, idx, cfg, SEED0, cap, "calib"))
        run_jobs(jobs, "calibration: rank 23, the 006 scheme's own cube")
    elif mode in ("split", "one"):
        split = tuple(int(x) for x in sys.argv[2:6])
        cap = float(sys.argv[6])
        cfg = sys.argv[7]
        allc = cubes_cfg(split, cfg)
        if mode == "split":
            m = int(sys.argv[8]) if len(sys.argv) > 8 else len(allc)
            sel = list(range(min(m, len(allc))))
        else:
            sel = [int(x) for x in sys.argv[8:]]
        seed = SEED0 + 1
        jobs = [(split, allc[i], i, cfg, seed, cap, "k") for i in sel]
        run_jobs(jobs, f"rank-22 cubes, {len(allc)} in all")
    elif mode == "subw":
        # depth-2 sample: w of the first size-3 slot pinned to the given symmetric matrices (measurement only)
        split = tuple(int(x) for x in sys.argv[2:6])
        cap = float(sys.argv[6])
        cfg = sys.argv[7]
        ci = int(sys.argv[8])
        ws = [int(x) for x in sys.argv[9:]]
        allc = cubes_cfg(split, cfg)
        jobs = [(split, allc[ci], ci, cfg, SEED0 + 2, cap, "k", w) for w in ws]
        run_jobs(jobs, f"depth-2 sample (w of first size-3 slot pinned), cube {ci} of {len(allc)}; not a closure")
    else:
        raise SystemExit(__doc__)
