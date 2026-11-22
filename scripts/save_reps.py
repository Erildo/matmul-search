"""Rerun the 003 follow-up walks and keep one scheme per invariant class.

Same commands and seeds as scripts/run_003b.py (flip_graph --seed k, then
find_alternative_schemes --seed k, 20,000 schemes, 1 thread). For every scheme
that verifies, the first one seen in each class is kept. Output, per field:
results/reps/<field>_r23_reps.jsonl.gz, one JSON line per class:
{"class": key, "walk": k, "scheme": "<ternary_flip_graph txt>"}.
The class sets are then compared with the keys 003 saved (scratch/003b/keys/), as a
reproducibility check.

usage: python3 -I scripts/save_reps.py [K] [COUNT]     (defaults 5, 20000)
"""

from __future__ import annotations

import concurrent.futures as cf
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402

FG = ROOT / "vendor" / "ternary_flip_graph" / "flip_graph"
ALT = ROOT / "vendor" / "ternary_flip_graph" / "find_alternative_schemes"
WORK = ROOT / "scratch" / "reps"
OUT = ROOT / "results" / "reps"
MOD = {"Z3": 3, "Z2": 2}
NAME = {"Z3": "f3", "Z2": "gf2"}


def walk(ring: str, k: int, count: int) -> dict:
    p = MOD[ring]
    d = WORK / f"{ring}_{k}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    subprocess.run([str(FG), "-n1", "3", "-n2", "3", "-n3", "3", "--ring", ring, "--target-rank", "23",
                    "--seed", str(k), "--threads", "1", "--count", "1", "-o", str(d / "seed")],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    seed_file = next((d / "seed").rglob("*_m23_*.txt"))
    subprocess.run([str(ALT), "-i", str(seed_file), "--ring", ring, "--target-rank", "23",
                    "--max-count", str(count), "--threads", "1", "--seed", str(k), "-o", str(d / "alt")],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
    files = sorted((d / "alt").rglob("*.txt"), key=lambda f: f.stat().st_birthtime)
    reps, bad = {}, 0
    for f in files:
        try:
            key = class_key(load_fg(f, mod=p), p)
        except ValueError:
            bad += 1
            continue
        reps.setdefault(key, f.read_text())
    shutil.rmtree(d)
    return dict(ring=ring, k=k, schemes=len(files), failed_verify=bad, reps=reps)


def main(K: int, count: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(r, k) for r in MOD for k in range(1, K + 1)]
    with cf.ProcessPoolExecutor(2) as ex:
        results = list(ex.map(walk, *zip(*jobs), [count] * len(jobs)))
    for ring in MOD:
        ws = [w for w in results if w["ring"] == ring]
        keep: dict[str, tuple[int, str]] = {}
        repro = []
        for w in ws:
            for key, text in w["reps"].items():
                keep.setdefault(key, (w["k"], text))
            old = ROOT / "scratch" / "003b" / "keys" / f"{ring}_{w['k']}.json"
            if old.exists():
                repro.append(dict(walk=w["k"], same_classes=set(json.loads(old.read_text())) == set(w["reps"])))
        path = OUT / f"{NAME[ring]}_333_r23_reps.jsonl.gz"
        with gzip.open(path, "wt") as f:
            for key in sorted(keep):
                k, text = keep[key]
                f.write(json.dumps({"class": key, "walk": k, "scheme": text}) + "\n")
        field = "GF(3)" if ring == "Z3" else "GF(2)"
        summary = dict(classes=len(keep), schemes=sum(w["schemes"] for w in ws),
                       failed_verify=sum(w["failed_verify"] for w in ws), reproduces_003=repro,
                       file=str(path.relative_to(ROOT)), file_bytes=path.stat().st_size)
        record_run(ROOT / "results" / "runs.jsonl", attempt="003-f3-landscape", format="3x3x3",
                   method="save one representative per invariant class (rerun of 003 follow-up walks)",
                   field=field, best_rank_found=23, verified=summary["failed_verify"] == 0,
                   budget={"walks": K, "max_count": count, "processes": 2, "threads_each": 1}, **summary)
        print(field, json.dumps(summary), flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5, int(sys.argv[2]) if len(sys.argv) > 2 else 20000)
