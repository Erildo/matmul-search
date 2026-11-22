"""Attempt 007: cyclic F3 split (3,13), class cubes with every exclusion and all orbits nonzero, then
level-2 cubes by the characteristic polynomial of slot 1 (see attempts/007-cyclic-f3-three-orbits.md).

Level-1 cube c:      CyclicF3CNF(3, 3, 13, pin_u1=class_reps()[c], min_class=c, exclude_charpoly=True, all_nonzero=True)
Level-2 cube (c, p): the same with pin_cp1=p, for p in cp1_choices(c).
Every instance gets one line in results/runs.jsonl.

usage: python3 -I scripts/run_007.py sizes | A | B
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import CyclicF3CNF, class_reps, cp1_choices  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402

WORK = ROOT / "scratch" / "007"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
ATTEMPT = "007-cyclic-f3-three-orbits"
A, B = 3, 13
STAGE_B_BUDGET = 120 * 60  # CPU-seconds


def encoder(c, p=None):
    return CyclicF3CNF(3, A, B, pin_u1=class_reps()[c], min_class=c, exclude_charpoly=True, all_nonzero=True,
                       pin_cp1=p)


def run_instance(c, p, seed, cap, check_cap):
    name = f"n3_a{A}_b{B}_c{c}" + ("" if p is None else "_p" + "".join(map(str, p))) + f"_cadical_s{seed}"
    return solve_instance(encoder(c, p), name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT,
                          drat_trim=DRAT_TRIM, attempt=ATTEMPT, prefix="007", check_cap=check_cap,
                          meta=dict(cube=c, cp1=p, exclude_charpoly=True, all_nonzero=True))


def record(res, stage, cap):
    cube = f"class cube {res['cube']}" + ("" if res["cp1"] is None else f", slot-1 char-poly {tuple(res['cp1'])}")
    record_run(RUNS, attempt=ATTEMPT, stage=stage, format="3x3x3", field="GF(3)",
               method=f"SAT cadical on cyclic-invariant Brent equations mod 3, split a={A} b={B}, {cube}, "
                      "rank-1 bans + char-poly exclusion + all orbits nonzero",
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3")) if res["status"] == "SAT" else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def run_many(jobs, stage, cap, check_cap, budget=None):
    """jobs: (c, p, seed). With a budget (CPU-seconds), stop submitting once the spent time reaches it."""
    spent, results, pending = 0.0, [], list(jobs)
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        running = set()
        while pending or running:
            while pending and len(running) < 2 and (budget is None or spent < budget):
                running.add(ex.submit(run_instance, *pending.pop(0), cap, check_cap))
            if not running:
                break
            done, running = cf.wait(running, return_when=cf.FIRST_COMPLETED)
            for f in done:
                res = f.result()
                spent += res["seconds"] + res.get("proof_seconds", 0)
                record(res, stage, cap)
                print(json.dumps(res), flush=True)
                results.append(res)
    if pending:
        print(json.dumps({"not_run": [(c, p) for c, p, _ in pending], "spent_cpu_s": round(spent)}), flush=True)
    return results


def open_cubes():
    """Level-1 cubes that stage A left open (status other than UNSAT with a VERIFIED proof)."""
    closed, seen = set(), set()
    for line in RUNS.read_text().splitlines():
        r = json.loads(line)
        if r.get("attempt") == ATTEMPT and r.get("stage") == "A":
            c = r["result"]["cube"]
            seen.add(c)
            if r["result"]["status"] == "UNSAT" and r["result"].get("drat_trim") == "VERIFIED":
                closed.add(c)
    if seen != set(range(21)):
        raise SystemExit(f"stage A incomplete: cubes run {sorted(seen)}")
    return sorted(set(range(21)) - closed)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    if stage == "sizes":
        for c in range(21):
            t = time.time()
            e = encoder(c)
            print(json.dumps(dict(cube=c, vars=e.nv, clauses=len(e.clauses), level2_children=len(cp1_choices(c)),
                                  build_s=round(time.time() - t, 1))), flush=True)
    elif stage == "A":
        run_many([(c, None, 70001) for c in range(21)], "A", 120, 900)
    elif stage == "B":
        cubes = sorted(open_cubes(), key=lambda c: len(cp1_choices(c)))
        jobs = [(c, p, 70002) for c in cubes for p in cp1_choices(c)]
        print(json.dumps({"open_cubes": cubes, "level2_jobs": len(jobs)}), flush=True)
        run_many(jobs, "B", 120, 900, budget=STAGE_B_BUDGET)
    else:
        raise SystemExit(__doc__)
