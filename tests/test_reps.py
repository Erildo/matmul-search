"""Spot-check the saved representative files (the full check is scripts/check_reps.py).

Every file in results/reps/ must open, and its first schemes must verify,
match their stored class, and have distinct classes.
"""

import gzip
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from check_reps import check_line  # noqa: E402

FILES = sorted((ROOT / "results" / "reps").glob("*_333_r23_reps.jsonl.gz"))
SAMPLE = 40


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_rep_file_sample(path):
    p = 3 if path.name.startswith("f3") else 2
    with gzip.open(path, "rt") as f:
        lines = [next(f) for _ in range(SAMPLE)]
    keys = []
    for line in lines:
        key, err = check_line((line, p))
        assert err is None, err
        keys.append(key)
    assert len(set(keys)) == len(keys)
    assert all(json.loads(l)["class"] == k for l, k in zip(lines, keys))


def test_rep_files_exist():
    assert {p.name for p in FILES} >= {"f3_333_r23_reps.jsonl.gz", "gf2_333_r23_reps.jsonl.gz"}
