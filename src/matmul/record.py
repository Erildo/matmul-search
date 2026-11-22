"""Append-only, machine-readable record of every run: results/runs.jsonl.

One JSON object per line. Each run links to the attempt note in attempts/
that explains it. Never edit old lines; add a new run instead.
"""

from __future__ import annotations

import datetime as dt
import json
import subprocess
from pathlib import Path

REQUIRED = ("attempt", "format", "method", "field", "best_rank_found", "verified", "budget")


def git_commit(cwd: str | Path = ".") -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=cwd, capture_output=True, text=True, timeout=5
        )
        return out.stdout.strip() or None
    except Exception:
        return None


def record_run(path: str | Path, **run) -> dict:
    missing = [k for k in REQUIRED if k not in run]
    if missing:
        raise ValueError(f"run record missing fields: {missing}")
    run.setdefault("time", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    run.setdefault("commit", git_commit(Path(path).parent))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(json.dumps(run, sort_keys=True) + "\n")
    return run


def read_runs(path: str | Path) -> list[dict]:
    path = Path(path)
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
