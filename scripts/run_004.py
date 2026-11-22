"""Attempt 004: many short F3 walks for 3x3x3 rank 22, one pocket per seed.

Each seed runs flip_graph (ring Z3, target 22, 1 runner, 1 thread) from the naive
scheme under a hard per-seed cap. Seeds run in 2 worker processes until the
total budget is used. Any rank-22 scheme is verified mod 3 and copied to
results/f3/. One rank-23 scheme per seed is kept in scratch/004/reps/ with its
class key, to check that the seeds sampled different pockets.

usage: python3 -I scripts/run_004.py [total_seconds] [per_seed_seconds] [first_seed]
       (defaults 3600 CPU-seconds in total, 20 per seed, seed 1)
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402

FG = ROOT / "vendor" / "ternary_flip_graph" / "flip_graph"
WORK = ROOT / "scratch" / "004"
FOUND = ROOT / "results" / "f3"
RANK_RE = re.compile(r"_m(\d+)_")
WORKERS = 2


def one_seed(k: int, cap: float) -> dict:
    d = WORK / "runs" / str(k)
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    t0 = time.time()
    try:
        subprocess.run([str(FG), "-n1", "3", "-n2", "3", "-n3", "3", "--ring", "Z3", "--target-rank", "22",
                        "--count", "1", "--threads", "1", "--seed", str(k), "-o", str(d)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=cap)
        stop = "finished"
    except subprocess.TimeoutExpired:
        stop = "max_seconds"
    files = sorted(d.rglob("*.txt"), key=lambda f: f.stat().st_birthtime)
    ranks = [int(RANK_RE.search(f.name).group(1)) for f in files]
    best = min(ranks, default=27)
    t23 = next((round(f.stat().st_birthtime - t0, 2) for f, r in zip(files, ranks) if r <= 23), None)
    out = dict(seed=k, best_rank=best, seconds_to_23=t23, stop=stop)
    for f, r in zip(files, ranks):
        if r <= 22:  # the thing we are looking for: keep it safely, then verify
            FOUND.mkdir(parents=True, exist_ok=True)
            shutil.copy(f, FOUND / f"004_seed{k}_{f.name}")
            try:
                load_fg(f, mod=3)
                out["found_22"] = str(FOUND / f"004_seed{k}_{f.name}")
            except ValueError as e:
                out["found_22_failed_verify"] = f"{f.name}: {e}"
    last23 = [f for f, r in zip(files, ranks) if r == 23]
    if last23:
        s = load_fg(last23[-1], mod=3)
        out["class_23"] = class_key(s, 3)
        (WORK / "reps").mkdir(exist_ok=True)
        shutil.copy(last23[-1], WORK / "reps" / f"seed{k}.txt")
    shutil.rmtree(d)
    return out


def main(total: float, cap: float, first: int) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (ROOT / "scratch" / ".metadata_never_index").touch()
    n_seeds = int(total / cap)
    seeds = range(first, first + n_seeds)
    t0 = time.time()
    with cf.ProcessPoolExecutor(WORKERS) as ex:
        results = list(ex.map(one_seed, seeds, [cap] * n_seeds))
    wall = time.time() - t0
    best = min(r["best_rank"] for r in results)
    classes = Counter(r.get("class_23") for r in results if r.get("class_23"))
    found = [r for r in results if "found_22" in r]
    summary = dict(seeds=n_seeds, wall_seconds=round(wall, 1), best_rank=best,
                   reached_23=sum(1 for r in results if r["best_rank"] <= 23),
                   distinct_final_classes=len(classes), repeated_classes=sum(c - 1 for c in classes.values()),
                   found_22=[r["found_22"] for r in found])
    (WORK / "per_seed.jsonl").write_text("".join(json.dumps(r) + "\n" for r in results))
    record_run(ROOT / "results" / "runs.jsonl", attempt="004-f3-many-pockets", format="3x3x3",
               method="ternary_flip_graph Z3, many short independent walks", field="GF(3)",
               best_rank_found=best, verified=True if found or best == 23 else None,
               budget={"total_cpu_seconds": total, "per_seed_seconds": cap, "workers": WORKERS, "threads_each": 1,
                       "first_seed": first}, **summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(float(a[0]) if a else 3600.0, float(a[1]) if len(a) > 1 else 20.0, int(a[2]) if len(a) > 2 else 1)
