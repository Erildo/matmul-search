"""Attempt 006: cyclic-invariant F3 search split into cubes by matrix class (see attempts/006-cyclic-f3-cubes.md).

A cube (n, a, b, c) is CyclicF3CNF(n, a, b, pin_u1=class_reps()[c], min_class=c): orbit slot 0's u is the
representative of class c, and (minimal class form) no orbit matrix lies in class 0 (c >= 1) or is of
rank one (c >= 2). c = None is the unpinned 005 instance. Solving, decoding, verifying, saving and the
drat-trim check work as in scripts/run_005.py. Every instance gets one line in results/runs.jsonl.

usage: python3 -I scripts/run_006.py STAGE        (STAGE = A, B, C, B2 or C2)
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import CyclicF3CNF, class_reps  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402

WORK = ROOT / "scratch" / "006"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
ATTEMPT = "006-cyclic-f3-cubes"


def run_instance(n, a, b, c, solver, seed, cap, check_cap=None, excl=False):
    name = f"n{n}_a{a}_b{b}_" + ("free" if c is None else f"c{c}") + ("x" if excl else "") + f"_{solver}_s{seed}"
    enc = CyclicF3CNF(n, a, b) if c is None else CyclicF3CNF(n, a, b, pin_u1=class_reps(n)[c], min_class=c,
                                                             exclude_charpoly=excl)
    return solve_instance(enc, name, solver, seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                          attempt=ATTEMPT, prefix="006", check_cap=check_cap,
                          meta=dict(cube=c, exclude_charpoly=excl))


def record(res, stage, cap):
    best = res.get("found_rank")
    cube = "unpinned" if res["cube"] is None else f"class cube {res['cube']}"
    if res.get("exclude_charpoly"):
        cube += " with char-poly exclusion"
    record_run(ROOT / "results" / "runs.jsonl", attempt=ATTEMPT, stage=stage,
               format=f"{res['n']}x{res['n']}x{res['n']}", field="GF(3)",
               method=f"SAT {res['solver']} on cyclic-invariant Brent equations mod 3, split a={res['a']} b={res['b']}, {cube}",
               best_rank_found=best,
               verified=bool(res.get("verified_mod3")) if res["status"] == "SAT" else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def run_many(jobs, stage, cap, check_cap, excl=False):
    results = []
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(run_instance, *j, cap, check_cap, excl) for j in jobs]
        for f in futs:
            res = f.result()
            record(res, stage, cap)
            print(json.dumps(res), flush=True)
            results.append(res)
    return results


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "A":
        run_many([(3, 7, 2, 0, "cadical", 60001), (3, 7, 2, 0, "cadical", 60002),
                  (3, 7, 2, None, "cadical", 60003), (3, 7, 2, None, "cadical", 60004)], "A", 300, None)
    elif stage == "B":
        run_many([(3, 2, 16, c, "cadical", 60011) for c in range(21)], "B", 120, 600)
    elif stage == "C":
        run_many([(3, 7, 1, c, "cadical", 60012) for c in range(21)], "C", 120, 600)
    elif stage == "B2":
        run_many([(3, 2, 16, c, "cadical", 60013) for c in (10, 19)], "B2", 300, 900, excl=True)
    elif stage == "C2":
        run_many([(3, 7, 1, c, "cadical", 60014) for c in range(10, 21)], "C2", 120, 600, excl=True)
    else:
        raise SystemExit(__doc__)
