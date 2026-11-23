"""Re-check every proved statement about cyclic-invariant 3x3x3 schemes over F3 (attempt 009).

Each CNF is regenerated from the encoder (nothing is read from scratch/), solved with a DRAT proof by
CaDiCaL and by Kissat, and each proof is checked with drat-trim. One JSON line per (instance, solver)
goes to results/f3/recheck_009.jsonl. Proofs are deleted after checking.

  T1: no cyclic F3 scheme of rank <= 22 with <= 2 orbits
      (0,22) trivially UNSAT (no solver needed); (1,19); (2,16) class cubes 0..20 with the char-poly exclusion
  T2: for each seed scheme S, no cyclic F3 scheme of rank 22, split (7,1) or (6,4), shares >= 5 orbits with S
      168 cubes: 5 orbits of S frozen in slots 0..4, the rest free

usage: python3 -I scripts/recheck_cyclic_f3.py [--quick] [--budget CPU_SECONDS]
       --quick: T1 with CaDiCaL only (about 12 CPU-min)
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import subprocess
import sys
import threading
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import CyclicF3CNF, T, class_reps, split_orbits  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402

WORK = ROOT / "scratch" / "009"
OUT = ROOT / "results" / "f3" / "recheck_009.jsonl"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
SOLVE_CAP, CHECK_CAP = 600, 900


def seed_schemes():
    """The four seed schemes of 008, by name, as (orbits, fixed)."""
    f3 = ROOT / "results" / "f3"
    out = {}
    for f in sorted(f3.glob("006_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out["006-c0-s60001"] = split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {}))
    d = json.loads((ROOT / "tests" / "data" / "bilr_2fix_z3.json").read_text())
    out["bilr"] = ([(A, B, T(C)) for A, B, C in d["orbits"]], d["fixed"])
    for f in sorted(f3.glob("005_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out["005-" + f.stem.split("_")[-2]] = split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {}))
    assert len(out) == 4
    return out


def t1_instances():
    yield "T1", "a1_b19", lambda: CyclicF3CNF(3, 1, 19)
    for c in range(21):
        yield "T1", f"a2_b16_c{c}", (lambda c=c: CyclicF3CNF(3, 2, 16, pin_u1=class_reps()[c], min_class=c,
                                                             exclude_charpoly=True))


def t2_instances():
    for name, (orbits, _) in seed_schemes().items():
        for shape, (a, b) in (("reduce71", (7, 1)), ("reduce64", (6, 4))):
            for keep in combinations(range(7), 5):
                yield "T2", f"{shape}_{name}_keep{''.join(map(str, keep))}", (
                    lambda a=a, b=b, keep=keep, orbits=orbits: CyclicF3CNF(3, a, b, freeze=[orbits[i] for i in keep]))


def run(claim, name, build, solver, seed=90001):
    WORK.mkdir(parents=True, exist_ok=True)
    enc = build()
    text = enc.dimacs()
    tag = f"{name}_{solver}"
    cnf, proof = WORK / f"{tag}.cnf", WORK / f"{tag}.drat"
    cnf.write_text(text)
    out = dict(claim=claim, instance=name, solver=solver, seed=seed, vars=enc.nv, clauses=len(enc.clauses),
               cnf_sha256=hashlib.sha256(text.encode()).hexdigest())
    t0 = time.time()
    try:
        r = subprocess.run([solver, "-q", f"--seed={seed}", str(cnf), str(proof)], capture_output=True, text=True,
                           timeout=SOLVE_CAP)
        out["status"] = {10: "SAT", 20: "UNSAT"}.get(r.returncode, f"exit{r.returncode}")
    except subprocess.TimeoutExpired:
        out["status"] = "TIMEOUT"
    out["solve_s"] = round(time.time() - t0, 2)
    out["proof_bytes"] = proof.stat().st_size if proof.exists() else 0
    if out["status"] == "UNSAT":
        t1 = time.time()
        try:
            r = subprocess.run([str(DRAT_TRIM), str(cnf), str(proof), "-t", str(CHECK_CAP)], capture_output=True,
                               text=True, timeout=CHECK_CAP + 30)
            out["drat_trim"] = "VERIFIED" if "s VERIFIED" in r.stdout else (
                [l for l in r.stdout.splitlines() if l.startswith("s ")] or ["no verdict"])[0]
        except subprocess.TimeoutExpired:
            out["drat_trim"] = "TIMEOUT"
        out["check_s"] = round(time.time() - t1, 2)
    proof.unlink(missing_ok=True)
    cnf.unlink(missing_ok=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--budget", type=float, default=90 * 60)
    args = ap.parse_args()

    trivial = CyclicF3CNF(3, 0, 22).trivially_unsat
    print(json.dumps({"claim": "T1", "instance": "a0_b22", "trivially_unsat": trivial}), flush=True)
    if not trivial:
        raise SystemExit("a = 0 is not trivially UNSAT: stop")

    t1, t2 = list(t1_instances()), list(t2_instances())
    jobs = [(c, n, b, "cadical") for c, n, b in t1]
    if not args.quick:
        jobs += [(c, n, b, "kissat") for c, n, b in t1]
        jobs += [(c, n, b, "cadical") for c, n, b in t2]
        jobs += [(c, n, b, "kissat") for c, n, b in t2]

    lock, spent, skipped, results = threading.Lock(), [0.0], [], []

    def one(job):
        with lock:
            if spent[0] >= args.budget:
                skipped.append(job[1] + "_" + job[3])
                return
        res = run(*job)
        with lock:
            spent[0] += res["solve_s"] + res.get("check_s", 0)
            results.append(res)
            with OUT.open("a") as fh:
                fh.write(json.dumps(res) + "\n")
            print(json.dumps(res), flush=True)

    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        list(ex.map(one, jobs))

    bad = [r for r in results if r["status"] != "UNSAT" or r.get("drat_trim") != "VERIFIED"]
    summary = {"checked": len(results), "all_unsat_verified": not bad, "bad": [(r["instance"], r["solver"]) for r in bad],
               "not_run": len(skipped), "cpu_s": round(spent[0])}
    print(json.dumps(summary), flush=True)
    record_run(ROOT / "results" / "runs.jsonl", attempt="009-cyclic-f3-writeup", stage="recheck" + ("-quick" if args.quick else ""),
               format="3x3x3", field="GF(3)",
               method="re-check of T1 and T2: CNFs regenerated from the encoder, CaDiCaL and Kissat DRAT proofs checked by drat-trim",
               best_rank_found=None, verified=not bad and not skipped,
               budget={"max_seconds": args.budget, "processes": 2, "threads_each": 1}, result=summary)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
