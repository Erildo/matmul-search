"""Attempt 010, stage B ("sweep"): all 60 rank-22 S3 splits with the native encoder S3NativeF3CNF.

Order (hardness proxy): first by the number of sigma-orbits 2a+b (009's walls were in the orbit count: T1 closed
every split with <= 2 orbits, nothing with 3 closed), then by the number of free F3 values 27a+15b+9c+6d (log_3 of
the search space), smallest first.

usage: python3 -I scripts/run_010_sweep.py sizes
       python3 -I scripts/run_010_sweep.py run CAP START END        # CaDiCaL on ordered indices [START, END)
       python3 -I scripts/run_010_sweep.py kissat CAP IDX [IDX ...]  # Kissat on the given ordered indices

UNSAT: re-solved with a DRAT proof and checked by drat-trim (check cap 900 s). SAT: decoded, verified mod 3,
checked for S3 invariance, and everything still queued is cancelled.
Every result is printed and appended as one JSON line to results/logs/010-s3-f3/sweep.log, and recorded
(stage "B") in results/runs.jsonl.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.s3 import S3NativeF3CNF, is_s3_invariant, scheme_s3_split, s3_splits  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402
from matmul.verify import verify  # noqa: E402

WORK = ROOT / "scratch" / "010" / "sweep"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
LOG = ROOT / "results" / "logs" / "010-s3-f3" / "sweep.log"
ATTEMPT = "010-s3-f3"
SEED = 110001
CHECK_CAP = 900
CPU_BUDGET = 150 * 60  # CPU-seconds for this stream


def free_values(sp):
    a, b, c, d = sp
    return 27 * a + 15 * b + 9 * c + 6 * d


def ordered():
    return sorted(s3_splits(22, 3), key=lambda sp: (2 * sp[0] + sp[1], free_values(sp), sp))


def cpu_used():
    if not LOG.exists():
        return 0.0
    tot = 0.0
    for line in LOG.read_text().splitlines():
        r = json.loads(line)
        tot += r.get("seconds", 0) + r.get("proof_seconds", 0)
    return tot


def run_one(idx, sp, solver, cap):
    a, b, c, d = sp
    name = f"n3_r22_s3_{a}{b}{c}{d}_native_sweep_{solver}_s{SEED}"
    t = time.time()
    enc = S3NativeF3CNF(3, a, b, c, d)
    build = round(time.time() - t, 1)
    res = solve_instance(enc, name, solver, SEED, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="010", check_cap=CHECK_CAP,
                         meta=dict(order_index=idx, s3_split=list(sp), s3_rank_max=22, sigma_orbits=2 * a + b,
                                   free_values=free_values(sp), build_s=build, encoder="native"))
    if res["status"] == "SAT":
        s = Scheme.load(ROOT / res["scheme_file"])
        res.update(s3_invariant=is_s3_invariant(s), s3_split_found=list(scheme_s3_split(s)),
                   verify_mod3_recheck=verify(s, mod=3).ok, class_key=class_key(s, 3))
    return res


def record(res):
    sat = res["status"] == "SAT"
    record_run(RUNS, attempt=ATTEMPT, stage="B", format="3x3x3", field="GF(3)",
               method=f"SAT {res['solver']} on S3-invariant Brent equations mod 3 (native encoder), "
                      f"S3 split (a,b,c,d)={tuple(res['s3_split'])}, stage B sweep",
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3") and res.get("s3_invariant")) if sat
               else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": res["cap_seconds"], "processes": 2, "threads_each": 1}, result=res)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f:
        f.write(json.dumps(res) + "\n")
    print(json.dumps(res), flush=True)


def run_jobs(jobs, solver, cap):
    """jobs: list of (idx, split). Two processes; stop queueing on SAT or when the CPU budget is spent."""
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        pending = list(jobs)
        live = {}
        stop = False
        while pending or live:
            while not stop and pending and len(live) < 2:
                if cpu_used() + cap * len(live) + cap > CPU_BUDGET:
                    print(json.dumps(dict(note="CPU budget reached; not starting", remaining=[list(s) for _, s in pending])), flush=True)
                    stop = True
                    break
                idx, sp = pending.pop(0)
                live[ex.submit(run_one, idx, sp, solver, cap)] = sp
            if not live:
                break
            done, _ = cf.wait(live, return_when=cf.FIRST_COMPLETED)
            for f in done:
                live.pop(f)
                res = f.result()
                record(res)
                if res["status"] == "SAT":
                    print(json.dumps(dict(note="SAT AT RANK 22 - STOPPING", split=res["s3_split"])), flush=True)
                    stop = True
                    pending.clear()


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    order = ordered()
    if mode == "sizes":
        for i, sp in enumerate(order):
            e = S3NativeF3CNF(3, *sp)
            print(json.dumps(dict(idx=i, split=sp, sigma_orbits=2 * sp[0] + sp[1], free_values=free_values(sp),
                                  vars=e.nv, clauses=len(e.clauses))), flush=True)
    elif mode == "run":
        cap, lo, hi = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        run_jobs([(i, order[i]) for i in range(lo, hi)], "cadical", cap)
    elif mode == "runidx":
        cap = int(sys.argv[2])
        run_jobs([(int(i), order[int(i)]) for i in sys.argv[3:]], "cadical", cap)
    elif mode == "kissat":
        cap = int(sys.argv[2])
        run_jobs([(int(i), order[int(i)]) for i in sys.argv[3:]], "kissat", cap)
    else:
        raise SystemExit(__doc__)
    print(json.dumps(dict(note="cpu_used_total", seconds=round(cpu_used(), 1))), flush=True)
