"""Benchmark ternary_flip_graph: time to reach target ranks, with hard caps.

Each (case, seed) is one process with 1 thread and 1 runner, killed at its cap.
Seeds of a case run in parallel. Time-to-rank comes from the creation times of
the scheme files the tool saves on every improvement. Every final scheme is
re-checked by our own verifier. Results are appended to results/benchmarks.jsonl.

usage: python3 scripts/bench_fg.py [case ...]   (default: all cases)
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.record import git_commit  # noqa: E402

FG = ROOT / "vendor" / "ternary_flip_graph" / "flip_graph"
OUT = ROOT / "results" / "benchmarks.jsonl"

# name: (n, m, p, ring, target, cap_seconds)
CASES = {
    "222_ZT_7": (2, 2, 2, "ZT", 7, 30),
    "333_ZT_23": (3, 3, 3, "ZT", 23, 120),
    "333_Z2_23": (3, 3, 3, "Z2", 23, 120),
}
SEEDS = list(range(1, 11))
PARALLEL = 2  # keep the Mac usable; 4x4x4 and up were too heavy here
RANK_RE = re.compile(r"_m(\d+)_")


def run_one(name, n, m, p, ring, target, cap, seed):
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "out"
        cmd = [str(FG), "-n1", str(n), "-n2", str(m), "-n3", str(p), "--ring", ring,
               "--target-rank", str(target), "--seed", str(seed),
               "--threads", "1", "--count", "1", "-o", str(out)]
        t0 = time.time()
        try:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=cap)
            killed = False
        except subprocess.TimeoutExpired:  # run() kills the child on timeout
            killed = True
        wall = time.time() - t0
        files = sorted(out.glob("*.txt"), key=lambda f: f.stat().st_birthtime) if out.exists() else []
        curve = [(int(RANK_RE.search(f.name).group(1)), round(f.stat().st_birthtime - t0, 3)) for f in files]
        best = min((r for r, _ in curve), default=n * m * p)
        verified = None
        if files:
            best_file = min(files, key=lambda f: int(RANK_RE.search(f.name).group(1)))
            try:
                load_fg(best_file, mod=2 if ring == "Z2" else None)
                verified = True
            except ValueError:
                verified = False
        reached = best <= target
        t_target = next((t for r, t in curve if r <= target), None)
        return dict(case=name, format=f"{n}x{m}x{p}", ring=ring, target=target, seed=seed,
                    cap_seconds=cap, reached=reached, seconds_to_target=t_target,
                    best_rank=best, wall_seconds=round(wall, 2), killed=killed,
                    verified_by_us=verified, curve=curve)


def main(names):
    commit = git_commit(ROOT)
    for name in names:
        n, m, p, ring, target, cap = CASES[name]
        with cf.ThreadPoolExecutor(PARALLEL) as ex:
            runs = list(ex.map(lambda s: run_one(name, n, m, p, ring, target, cap, s), SEEDS))
        with OUT.open("a") as f:
            for r in runs:
                r.update(time=time.strftime("%Y-%m-%dT%H:%M:%S"), commit=commit, machine="Apple M5")
                f.write(json.dumps(r) + "\n")
        times = [r["seconds_to_target"] for r in runs if r["reached"]]
        med = f"{statistics.median(times):.2f}s" if times else "-"
        mx = f"{max(times):.2f}s" if times else "-"
        print(f"{name}: reached {len(times)}/{len(runs)}  median {med}  max {mx}  "
              f"best ranks {[r['best_rank'] for r in runs]}  verified {[r['verified_by_us'] for r in runs]}",
              flush=True)


if __name__ == "__main__":
    if not FG.exists():
        sys.exit(f"missing {FG}; build vendor/ternary_flip_graph first")
    main(sys.argv[1:] or list(CASES))
