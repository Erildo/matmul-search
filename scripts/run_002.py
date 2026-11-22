"""Attempt 002, stage B: look for 3x3x3 rank 22 over F3 with ternary_flip_graph.

Runs two different search variants side by side, each a separate process with a
hard wall-clock cap (the tool has no time limit of its own). Afterwards it counts
what each variant reached, re-verifies the best scheme mod 3 with our verifier,
and appends one line per variant to results/runs.jsonl.

usage: python3 -I scripts/run_002.py [seconds]     (default 1800)
Work files go to scratch/002/ (git-ignored).
"""

from __future__ import annotations

import concurrent.futures as cf
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.record import record_run  # noqa: E402

FG = ROOT / "vendor" / "ternary_flip_graph" / "flip_graph"
WORK = ROOT / "scratch" / "002"
RANK_RE = re.compile(r"_m(\d+)_")

COMMON = ["-n1", "3", "-n2", "3", "-n3", "3", "--ring", "Z3", "--target-rank", "22",
          "--count", "16", "--threads", "2", "--save-metrics"]
VARIANTS = {
    # pool descent, keep only structurally different schemes, favour schemes with many flips
    "pool": ["--use-pool", "--pool-size", "2000", "--pool-unique-type", "structure",
             "--pool-select-strategy", "flips", "--reduce-probability", "0.01", "--seed", "2021"],
    # plain walks with plus transitions, plus some sandwiching to move between orbits
    "walk": ["--sandwiching-probability", "0.01", "--reduce-probability", "0.01", "--seed", "2022"],
}


def run(name: str, cap: float) -> dict:
    out = WORK / name
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(FG), *COMMON, *VARIANTS[name], "-o", str(out / "schemes"),
           "--metrics-path", str(out / "metrics.jsonl")]
    t0 = time.time()
    with (out / "stdout.txt").open("w") as log:
        try:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=cap)
            stop = "finished"
        except subprocess.TimeoutExpired:
            stop = "max_seconds"
    wall = time.time() - t0
    files = list((out / "schemes").rglob("*.txt"))  # pool mode writes into a 3x3x3/ subfolder
    ranks = [int(RANK_RE.search(f.name).group(1)) for f in files]
    best = min(ranks, default=27)
    verified = None
    if files:
        best_file = next(f for f, r in zip(files, ranks) if r == best)
        try:
            load_fg(best_file, mod=3)
            verified = True
        except ValueError:
            verified = False
    return dict(name=name, cmd=cmd, stop=stop, wall=round(wall, 1), best=best, verified=verified,
                files_by_rank={r: ranks.count(r) for r in sorted(set(ranks))})


def main(cap: float) -> None:
    with cf.ThreadPoolExecutor(len(VARIANTS)) as ex:
        results = list(ex.map(lambda n: run(n, cap), VARIANTS))
    for r in results:
        record_run(
            ROOT / "results" / "runs.jsonl",
            attempt="002-f3-rank22", format="3x3x3", method=f"ternary_flip_graph {r['name']}",
            field="GF(3)", best_rank_found=r["best"], verified=r["verified"],
            budget={"max_seconds": cap, "processes": len(VARIANTS), "threads_each": 2},
            stop_reason=r["stop"], wall_seconds=r["wall"], files_by_rank=r["files_by_rank"],
            command=" ".join(r["cmd"][1:]),
        )
        print(f"{r['name']}: best rank {r['best']} verified={r['verified']} stop={r['stop']} "
              f"wall={r['wall']}s files_by_rank={r['files_by_rank']}", flush=True)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 1800.0)
