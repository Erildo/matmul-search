"""Attempt 010: S3-invariant schemes over F3 (see attempts/010-s3-f3.md).

Instance (a, b, c, d): S3F3CNF(3, a, b, c, d), rank 6a + 3b + 2c + d, every slot may be zero.
Every instance gets one line in results/runs.jsonl; SAT models are decoded, verified mod 3, checked for S3
invariance and given their S3 split and class key.

usage: python3 -I scripts/run_010.py sizes | A
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
from matmul.s3 import S3F3CNF, is_s3_invariant, scheme_s3_split, s3_splits  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402

WORK = ROOT / "scratch" / "010"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
ATTEMPT = "010-s3-f3"

STAGE_A = [((1, 5, 0, 2), 300), ((3, 1, 0, 1), 120), ((2, 3, 0, 1), 120), ((1, 5, 0, 1), 120), ((0, 7, 0, 1), 120)]


def run_instance(split, seed, cap, check_cap):
    a, b, c, d = split
    r = 6 * a + 3 * b + 2 * c + d
    name = f"n3_r{r}_s3_{a}{b}{c}{d}_cadical_s{seed}"
    t = time.time()
    enc = S3F3CNF(3, a, b, c, d)
    build = round(time.time() - t, 1)
    res = solve_instance(enc, name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="010", check_cap=check_cap,
                         meta=dict(s3_split=list(split), s3_rank_max=r, build_s=build))
    if res["status"] == "SAT":
        s = Scheme.load(ROOT / res["scheme_file"])
        res.update(s3_invariant=is_s3_invariant(s), s3_split_found=list(scheme_s3_split(s)), class_key=class_key(s, 3))
    return res


def record(res, stage, cap):
    record_run(RUNS, attempt=ATTEMPT, stage=stage, format="3x3x3", field="GF(3)",
               method=f"SAT cadical on S3-invariant Brent equations mod 3, S3 split (a,b,c,d)={tuple(res['s3_split'])}",
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3") and res.get("s3_invariant")) if res["status"] == "SAT"
               else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def run_many(jobs, stage, check_cap):
    """jobs: (split, seed, cap)."""
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(run_instance, sp, seed, cap, check_cap) for sp, seed, cap in jobs]
        for f, (_, _, cap) in zip(futs, jobs):
            res = f.result()
            record(res, stage, cap)
            print(json.dumps(res), flush=True)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else ""
    if stage == "sizes":
        for sp in [s for s, _ in STAGE_A] + s3_splits(22, 3):
            e = S3F3CNF(3, *sp)
            print(json.dumps(dict(split=sp, vars=e.nv, clauses=len(e.clauses))), flush=True)
    elif stage == "A":
        run_many([(sp, 100001, cap) for sp, cap in STAGE_A], "A", 600)
    else:
        raise SystemExit(__doc__)
