"""Independently check a file of saved representatives (results/reps/*.jsonl.gz).

For every line: parse the scheme, verify it exactly over GF(p) (rank 23,
3x3x3), recompute its invariant class, and check that it matches the stored
class. Finally check that all classes are distinct. Since different invariants
prove inequivalence, a passing file is a checkable lower bound: "at least N
pairwise inequivalent 3x3x3 rank-23 schemes over GF(p)".

usage: python3 -I scripts/check_reps.py results/reps/f3_333_r23_reps.jsonl.gz [p]
       (p defaults to 3 for files starting with f3, 2 for gf2)
"""

from __future__ import annotations

import concurrent.futures as cf
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.import_fg import parse_fg  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.verify import verify  # noqa: E402


def check_line(args) -> tuple[str, str | None]:
    line, p = args
    rec = json.loads(line)
    for tw in (True, False):
        s = parse_fg(rec["scheme"], transpose_w=tw)
        if s.format == "3x3x3" and s.rank == 23 and verify(s, mod=p).ok:
            key = class_key(s, p)
            return key, None if key == rec["class"] else f"stored class {rec['class']} != recomputed {key}"
    return rec["class"], "does not verify"


def main(path: str, p: int | None) -> int:
    p = p or (3 if Path(path).name.startswith("f3") else 2)
    lines = gzip.open(path, "rt").read().splitlines()
    with cf.ProcessPoolExecutor(2) as ex:
        out = list(ex.map(check_line, [(l, p) for l in lines], chunksize=500))
    errors = [(k, e) for k, e in out if e]
    keys = [k for k, _ in out]
    dupes = len(keys) - len(set(keys))
    ok = not errors and not dupes
    print(json.dumps(dict(file=path, p=p, schemes=len(lines), errors=len(errors), duplicate_classes=dupes,
                          first_errors=errors[:5], ok=ok)))
    if ok:
        print(f"OK: at least {len(set(keys))} pairwise inequivalent 3x3x3 rank-23 schemes over GF({p})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None))
