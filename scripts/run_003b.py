"""Attempt 003 follow-up: do independent walks meet the same rank-23 classes? (F3 vs GF(2))

For each field and k = 1..K: make a fresh rank-23 scheme from the naive scheme
(flip_graph --seed k), walk from it with find_alternative_schemes (--seed k,
COUNT schemes), verify and classify every scheme. Then compare walks: classes
shared between walks, and capture-recapture estimates N = c_i * c_j / m_ij.
Class keys per walk are kept in scratch/003b/keys/; scheme files are deleted
after classification to save disk.

usage: python3 -I scripts/run_003b.py [K] [COUNT]     (defaults 5, 20000)
"""

from __future__ import annotations

import concurrent.futures as cf
import itertools
import json
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import load_fg  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402

FG = ROOT / "vendor" / "ternary_flip_graph" / "flip_graph"
ALT = ROOT / "vendor" / "ternary_flip_graph" / "find_alternative_schemes"
WORK = ROOT / "scratch" / "003b"
MOD = {"Z3": 3, "Z2": 2}


def walk(ring: str, k: int, count: int) -> dict:
    """Seed scheme + walk + classification for one (field, k). Runs in its own process."""
    p = MOD[ring]
    d = WORK / f"{ring}_{k}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    t0 = time.time()
    subprocess.run([str(FG), "-n1", "3", "-n2", "3", "-n3", "3", "--ring", ring, "--target-rank", "23",
                    "--seed", str(k), "--threads", "1", "--count", "1", "-o", str(d / "seed")],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    seed_files = list((d / "seed").rglob("*_m23_*.txt"))
    if not seed_files:
        return dict(ring=ring, k=k, error="no rank-23 seed within 120s")
    load_fg(seed_files[0], mod=p)  # the seed itself must verify
    subprocess.run([str(ALT), "-i", str(seed_files[0]), "--ring", ring, "--target-rank", "23",
                    "--max-count", str(count), "--threads", "1", "--seed", str(k), "-o", str(d / "alt")],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
    gen_s = time.time() - t0
    files = sorted((d / "alt").rglob("*.txt"), key=lambda f: f.stat().st_birthtime)
    keys, bad, late_new, cut = [], 0, 0, int(len(files) * 0.8)
    seen = set()
    for i, f in enumerate(files):
        try:
            key = class_key(load_fg(f, mod=p), p)
        except ValueError:
            bad += 1
            continue
        if key not in seen:
            seen.add(key)
            keys.append(key)
            if i >= cut:
                late_new += 1
    (WORK / "keys").mkdir(exist_ok=True)
    (WORK / "keys" / f"{ring}_{k}.json").write_text(json.dumps(keys))
    shutil.rmtree(d / "alt")  # keep only the seed scheme and the class keys
    late_n = len(files) - cut
    return dict(ring=ring, k=k, schemes=len(files), failed_verify=bad, classes=len(seen),
                late_rate=round(late_new / late_n, 4) if late_n else None,
                gen_seconds=round(gen_s, 1), total_seconds=round(time.time() - t0, 1))


def compare(ring: str, walks: list[dict]) -> dict:
    sets = {w["k"]: set(json.loads((WORK / "keys" / f"{ring}_{w['k']}.json").read_text())) for w in walks}
    pairs = []
    for i, j in itertools.combinations(sorted(sets), 2):
        m = len(sets[i] & sets[j])
        ci, cj = len(sets[i]), len(sets[j])
        pairs.append(dict(i=i, j=j, shared=m, share=round(m / min(ci, cj), 4),
                          estimate=round(ci * cj / m) if m else None))
    union = set().union(*sets.values())
    counts = {}
    for s in sets.values():
        for key in s:
            counts[key] = counts.get(key, 0) + 1
    ests = [p["estimate"] for p in pairs if p["estimate"]]
    return dict(union=len(union), in_two_or_more=sum(1 for c in counts.values() if c >= 2),
                mean_share=round(statistics.mean(p["share"] for p in pairs), 4),
                median_estimate=statistics.median(ests) if ests else None,
                pairs_with_no_overlap=sum(1 for p in pairs if not p["shared"]), pairs=pairs)


def main(K: int, count: int) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    (ROOT / "scratch" / ".metadata_never_index").touch()
    jobs = [(r, k) for r in MOD for k in range(1, K + 1)]
    with cf.ProcessPoolExecutor(2) as ex:  # 2 processes, each tool call uses 1 thread
        results = list(ex.map(walk, *zip(*jobs), [count] * len(jobs)))
    for ring in MOD:
        ws = [w for w in results if w["ring"] == ring and "error" not in w]
        errs = [w for w in results if w["ring"] == ring and "error" in w]
        cmp = compare(ring, ws) if len(ws) >= 2 else {}
        field = "GF(3)" if ring == "Z3" else "GF(2)"
        record_run(ROOT / "results" / "runs.jsonl", attempt="003-f3-landscape", format="3x3x3",
                   method="5 independent seeds: flip_graph -> find_alternative_schemes + invariants",
                   field=field, best_rank_found=23,
                   verified=all(w["failed_verify"] == 0 for w in ws) and not errs,
                   budget={"walks": K, "max_count": count, "processes": 2, "threads_each": 1},
                   walks=ws, errors=errs, comparison=cmp)
        print(field, "walks", json.dumps(ws), flush=True)
        print(field, "compare", json.dumps({k: v for k, v in cmp.items() if k != "pairs"}), flush=True)
        print(field, "pairs", json.dumps(cmp.get("pairs")), flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5, int(sys.argv[2]) if len(sys.argv) > 2 else 20000)
