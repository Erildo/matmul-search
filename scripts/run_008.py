"""Attempt 008: SAT neighbourhoods of the known cyclic F3 rank-23 schemes (see attempts/008-cyclic-f3-lns.md).

Shapes, for a seed scheme S (split (7,2)) and a distance k:
  reduce71  CyclicF3CNF(3, 7, 1, near=(S, 7-k, 7))  rank 22; the fixed term is forced to be I
  reduce64  CyclicF3CNF(3, 6, 4, near=(S, 7-k, 6))  rank 22
  walk      CyclicF3CNF(3, 7, 2, near=(S, 7-k, 6))  rank 23, differs from S in at least one orbit
Every instance gets one line in results/runs.jsonl.

usage: python3 -I scripts/run_008.py A
"""

from __future__ import annotations

import concurrent.futures as cf
import gzip
import json
import sys
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from matmul.cyclic import CyclicF3CNF, T, expand, split_orbits  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402

WORK = ROOT / "scratch" / "008"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
REPS = ROOT / "results" / "reps" / "f3_333_r23_reps.jsonl.gz"
ATTEMPT = "008-cyclic-f3-lns"
SHAPES = {"reduce71": (7, 1), "reduce64": (6, 4), "walk": (7, 2)}


def seed_schemes():
    """name -> (orbits, fixed), in the stage A order: 006 first, then BILR, then the two 005 schemes."""
    out = {}
    for f in sorted(FOUND.glob("006_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out["006-c0-s60001"] = split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {}))
    d = json.loads((ROOT / "tests" / "data" / "bilr_2fix_z3.json").read_text())
    out["bilr"] = ([(A, B, T(C)) for A, B, C in d["orbits"]], d["fixed"])
    for f in sorted(FOUND.glob("005_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out["005-" + f.stem.split("_")[-2]] = split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {}))
    return out


_known = None


def known_classes():
    global _known
    if _known is None:
        _known = set()
        if REPS.exists():
            with gzip.open(REPS, "rt") as fh:
                _known = {json.loads(l)["class"] for l in fh}
    return _known


def run_instance(seed_name, S, shape, k, seed, cap, check_cap):
    a, b = SHAPES[shape]
    m, M = 7 - k, (6 if shape != "reduce71" else 7)
    enc = CyclicF3CNF(3, a, b, near=(S[0], m, M))
    name = f"{shape}_{seed_name}_k{k}_cadical_s{seed}"
    res = solve_instance(enc, name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="008", check_cap=check_cap,
                         meta=dict(shape=shape, seed_scheme=seed_name, k=k, near_min=m, near_max=M))
    if res["status"] == "SAT":
        d = json.loads((ROOT / res["scheme_file"]).read_text())
        s = Scheme(3, 3, 3, d["U"], d["V"], d["W"], {})
        key = class_key(s, 3)
        seed_key = class_key(expand(3, *S), 3)
        res.update(class_key=key, same_class_as_seed=(key == seed_key), in_saved_reps=key in known_classes())
    return res


def record(res, stage):
    what = {"reduce71": "rank 22 split (7,1)", "reduce64": "rank 22 split (6,4)", "walk": "rank 23 split (7,2)"}[res["shape"]]
    record_run(RUNS, attempt=ATTEMPT, stage=stage, format="3x3x3", field="GF(3)",
               method=f"SAT cadical, cyclic-invariant Brent equations mod 3, {what}, sharing >= {res['near_min']} "
                      f"(<= {res['near_max']}) orbits with seed scheme {res['seed_scheme']}",
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3")) if res["status"] == "SAT" else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": res["cap_seconds"], "processes": 2, "threads_each": 1}, result=res)


class Budget:
    def __init__(self, seconds):
        self.left, self.lock, self.skipped = seconds, threading.Lock(), []

    def take(self, res):
        with self.lock:
            self.left -= res["seconds"] + res.get("proof_seconds", 0)

    def ok(self):
        with self.lock:
            return self.left > 0


def ladder(seed_name, S, shape, seed, budget, stage):
    """k = 1, 2, ... until the first non-UNSAT (or k = 7). Returns the results."""
    out = []
    for k in range(1, 8):
        if not budget.ok():
            budget.skipped.append((seed_name, shape, k))
            break
        res = run_instance(seed_name, S, shape, k, seed, 60 if k <= 3 else 180, 600)
        budget.take(res)
        record(res, stage)
        print(json.dumps(res), flush=True)
        out.append(res)
        if res["status"] != "UNSAT":
            break
    return out


def walks(seed_name, S, seed, budget, stage):
    out = []
    for k in (2, 3, 4):
        if not budget.ok():
            budget.skipped.append((seed_name, "walk", k))
            continue
        res = run_instance(seed_name, S, "walk", k, seed, 60, None)
        budget.take(res)
        record(res, stage)
        print(json.dumps(res), flush=True)
        out.append(res)
    return out


def stage_a():
    budget = Budget(40 * 60)
    tasks = []
    for name, S in seed_schemes().items():
        tasks += [(ladder, name, S, "reduce71"), (ladder, name, S, "reduce64"), (walks, name, S)]
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(t[0], *t[1:], 80001, budget, "A") for t in tasks]
        for f in futs:
            f.result()
    print(json.dumps({"stage": "A", "cpu_s_left": round(budget.left), "not_run": budget.skipped}), flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["A"]:
        stage_a()
    else:
        raise SystemExit(__doc__)
