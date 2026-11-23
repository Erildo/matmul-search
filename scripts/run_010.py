"""Attempt 010: S3-invariant schemes over F3 (see attempts/010-s3-f3.md).

Instance (a, b, c, d): S3F3CNF(3, a, b, c, d), rank 6a + 3b + 2c + d, every slot may be zero.
Every instance gets one line in results/runs.jsonl; SAT models are decoded, verified mod 3, checked for S3
invariance and given their S3 split and class key.

usage: python3 -I scripts/run_010.py sizes | A | A2 | F

A uses the tied encoder S3F3CNF; A2 onwards the native encoder S3NativeF3CNF (same instances, half the CNF).
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
from matmul.s3 import S3F3CNF, S3NativeF3CNF, is_s3_invariant, scheme_s3_split, s3_splits  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402

WORK = ROOT / "scratch" / "010"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
ATTEMPT = "010-s3-f3"

STAGE_A = [((1, 5, 0, 2), 300), ((3, 1, 0, 1), 120), ((2, 3, 0, 1), 120), ((1, 5, 0, 1), 120), ((0, 7, 0, 1), 120)]


def run_instance(split, seed, cap, check_cap, native=False):
    a, b, c, d = split
    r = 6 * a + 3 * b + 2 * c + d
    name = f"n3_r{r}_s3_{a}{b}{c}{d}" + ("_native" if native else "") + f"_cadical_s{seed}"
    t = time.time()
    enc = (S3NativeF3CNF if native else S3F3CNF)(3, a, b, c, d)
    build = round(time.time() - t, 1)
    res = solve_instance(enc, name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="010", check_cap=check_cap,
                         meta=dict(s3_split=list(split), s3_rank_max=r, build_s=build,
                                   encoder="native" if native else "tied"))
    if res["status"] == "SAT":
        s = Scheme.load(ROOT / res["scheme_file"])
        res.update(s3_invariant=is_s3_invariant(s), s3_split_found=list(scheme_s3_split(s)), class_key=class_key(s, 3))
    return res


S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"
# Freeze ladder (stage F): split (1,5,0,2), native encoder; orbit slots 0 (size 6; slot 1 is its partner) and
# 2..6 (size 3) hold the 006 scheme in S3 normal form except the listed free ones; the two size-1 slots are always free.
LADDER = [("3x1", [6]), ("3x2", [5, 6]), ("3x3", [4, 5, 6]), ("3x4", [3, 4, 5, 6]), ("3x5", [2, 3, 4, 5, 6]),
          ("6", [0]), ("6+3x1", [0, 6])]


def run_freeze(label, free, seed, cap):
    from matmul.cyclic import split_orbits
    from matmul.s3 import s3_normal_form, s3_split
    orbits, _ = s3_normal_form(*s3_split(*split_orbits(Scheme.load(S3_006))))
    t = time.time()
    enc = S3NativeF3CNF(3, 1, 5, 0, 2)
    c0 = len(enc.clauses)
    for k in [0] + list(range(2, 7)):
        if k in free:
            continue
        for vals, mat in zip(enc.orb[k], orbits[k]):
            for z, x in zip(vals, (e % 3 for r in mat for e in r)):
                enc.clauses += [[lit] for lit in enc._is(z, x)]
    frozen_units = len(enc.clauses) - c0
    build = round(time.time() - t, 1)
    name = f"n3_r23_s3_1502_native_free{label.replace('+', 'p')}_cadical_s{seed}"
    res = solve_instance(enc, name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="010", check_cap=None,
                         meta=dict(s3_split=[1, 5, 0, 2], s3_rank_max=23, build_s=build, encoder="native",
                                   free=label, free_slots=free, frozen_units=frozen_units))
    if res["status"] == "SAT":
        s = Scheme.load(ROOT / res["scheme_file"])
        res.update(s3_invariant=is_s3_invariant(s), s3_split_found=list(scheme_s3_split(s)), class_key=class_key(s, 3))
    return res


def record(res, stage, cap):
    record_run(RUNS, attempt=ATTEMPT, stage=stage, format="3x3x3", field="GF(3)",
               method=f"SAT cadical on S3-invariant Brent equations mod 3 ({res['encoder']} encoder), S3 split (a,b,c,d)={tuple(res['s3_split'])}"
                      + (f", 006 scheme frozen except slots {res['free']}" if "free" in res else ""),
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3") and res.get("s3_invariant")) if res["status"] == "SAT"
               else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def run_many(jobs, stage, check_cap, native=False):
    """jobs: (split, seed, cap)."""
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(run_instance, sp, seed, cap, check_cap, native) for sp, seed, cap in jobs]
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
    elif stage == "A2":
        # native encoder, rank-23 rediscovery only (420 CPU-s of stage A's allowance left)
        run_many([((1, 5, 0, 2), 100001, 300)], "A2", 600, native=True)
    elif stage == "F":
        with cf.ThreadPoolExecutor(max_workers=2) as ex:
            futs = [ex.submit(run_freeze, lab, free, 100003, 120) for lab, free in LADDER]
            for f in futs:
                res = f.result()
                record(res, "F", 120)
                print(json.dumps(res), flush=True)
    else:
        raise SystemExit(__doc__)
