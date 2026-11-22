"""Attempt 003: how many inequivalent 3x3x3 rank-23 schemes does a walk meet, over F3 vs GF(2)?

For each field, find_alternative_schemes walks from one rank-23 seed and saves
every new scheme (by exact coefficients), up to a count cap and a hard time cap.
We then verify every saved scheme and count distinct invariant classes
(matmul.invariants) in the order the schemes were found. A class count close
to the scheme count means almost every scheme is genuinely new.

usage: python3 -I scripts/run_003.py [max_count] [max_seconds]   (defaults 20000, 300)
Work files: scratch/003/ (git-ignored).
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.invariants import class_key, term_profile  # noqa: E402
from matmul.record import record_run  # noqa: E402

ALT = ROOT / "vendor" / "ternary_flip_graph" / "find_alternative_schemes"
WORK = ROOT / "scratch" / "003"
SEEDS = {"Z3": (ROOT / "tests" / "data" / "f3_333_r23_a.txt", 3),
         "Z2": (ROOT / "tests" / "data" / "gf2_333_r23_a.txt", 2)}
MARKS = (100, 300, 1000, 3000, 10000, 20000, 50000)


def generate(ring: str, max_count: int, cap: float) -> dict:
    seed_file, _ = SEEDS[ring]
    out = WORK / ring
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(ALT), "-i", str(seed_file), "--ring", ring, "--target-rank", "23",
           "--max-count", str(max_count), "--threads", "1", "--seed", "3", "-o", str(out / "schemes")]
    t0 = time.time()
    with (out / "stdout.txt").open("w") as log:
        try:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=cap)
            stop = "finished"
        except subprocess.TimeoutExpired:
            stop = "max_seconds"
    return dict(ring=ring, cmd=cmd, stop=stop, wall=round(time.time() - t0, 1))


def analyze(ring: str) -> dict:
    _, p = SEEDS[ring]
    files = sorted((WORK / ring / "schemes").rglob("*.txt"), key=lambda f: f.stat().st_birthtime)
    keys, profiles, curve, bad = set(), set(), [], 0
    for i, f in enumerate(files, 1):
        try:
            s = load_fg(f, mod=p)  # verifies; raises if wrong
        except ValueError:
            bad += 1
            continue
        keys.add(class_key(s, p))
        profiles.add(term_profile(s, p))
        if i in MARKS:
            curve.append([i, len(keys)])
    if not curve or curve[-1][0] != len(files):
        curve.append([len(files), len(keys)])
    return dict(schemes=len(files), failed_verify=bad, classes=len(keys),
                term_profiles=len(profiles), curve=curve)


def main(max_count: int, cap: float) -> None:
    with cf.ThreadPoolExecutor(2) as ex:
        gens = list(ex.map(lambda r: generate(r, max_count, cap), SEEDS))
    for g in gens:
        a = analyze(g["ring"])
        field = "GF(3)" if g["ring"] == "Z3" else "GF(2)"
        record_run(ROOT / "results" / "runs.jsonl", attempt="003-f3-landscape", format="3x3x3",
                   method="find_alternative_schemes + invariants", field=field, best_rank_found=23,
                   verified=a["failed_verify"] == 0,
                   budget={"max_count": max_count, "max_seconds": cap, "threads": 1},
                   stop_reason=g["stop"], wall_seconds=g["wall"], command=" ".join(g["cmd"][1:]), **a)
        print(field, json.dumps({k: v for k, v in a.items()}), "stop", g["stop"], "wall", g["wall"], flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20000, float(sys.argv[2]) if len(sys.argv) > 2 else 300.0)
