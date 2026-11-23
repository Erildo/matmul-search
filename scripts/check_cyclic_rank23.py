"""Build or check the release file of cyclic-invariant 3x3x3 rank-23 schemes over F3 (attempt 009).

  build: collect BILR, the 005, 006 and 008 schemes, and write results/f3/cyclic_rank23.jsonl
  check: for every line, rebuild the scheme from its orbits and fixed terms, verify it mod 3, check cyclic invariance,
         recompute the class key and compare; then check the classes are pairwise distinct and (if the saved
         F3 representatives are present) which of them are among those. Prints a summary; exit 1 on any failure.

usage: python3 -I scripts/check_cyclic_rank23.py build | check [FILE]
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import T, expand, is_cyclic_invariant, normal_form, split_orbits  # noqa: E402
from matmul.invariants import class_key, rank_mod  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402
from matmul.verify import verify  # noqa: E402

F3 = ROOT / "results" / "f3"
RELEASE = F3 / "cyclic_rank23.jsonl"
REPS = ROOT / "results" / "reps" / "f3_333_r23_reps.jsonl.gz"


def _load(path):
    d = json.loads(path.read_text())
    return split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {}))


def build():
    d = json.loads((ROOT / "tests" / "data" / "bilr_2fix_z3.json").read_text())
    items = [("BILR S_2fix-Z3 (Ballard, Ikenmeyer, Landsberg, Ryder 2018), reduced mod 3", "tests/data/bilr_2fix_z3.json",
              ([(A, B, T(C)) for A, B, C in d["orbits"]], d["fixed"]))]
    for pat, src in (("005_3x3x3_rank23_*.json", "attempt 005 (SAT, unrestricted split (7,2))"),
                     ("006_3x3x3_rank23_*.json", "attempt 006 (SAT, class cube 0)"),
                     ("008_3x3x3_rank23_walk_*.json", "attempt 008 (SAT walk step)")):
        for f in sorted(F3.glob(pat)):
            items.append((src, str(f.relative_to(ROOT)), _load(f)))
    with RELEASE.open("w") as fh:
        for src, path, (orbits, fixed) in items:
            orbits, fixed = normal_form(orbits, fixed)
            s = expand(3, orbits, fixed)
            rec = dict(class_key=class_key(s, 3), source=src, file=path, split=[len(orbits), len(fixed)],
                       orbits_with_rank1_matrix=sum(any(rank_mod(m, 3) == 1 for m in t) for t in orbits),
                       orbits=[[m for m in t] for t in orbits], fixed=fixed)
            fh.write(json.dumps(rec) + "\n")
    print(f"wrote {len(items)} schemes to {RELEASE.relative_to(ROOT)}")


def check(path=RELEASE):
    recs = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    bad = []
    for r in recs:
        orbits = [tuple(t) for t in r["orbits"]]
        s = expand(3, orbits, r["fixed"])
        if s.rank != 23 or not verify(s, mod=3).ok:
            bad.append((r["class_key"], "does not verify mod 3 at rank 23"))
        if not is_cyclic_invariant(s):
            bad.append((r["class_key"], "not cyclic-invariant"))
        if class_key(s, 3) != r["class_key"]:
            bad.append((r["class_key"], "class key mismatch"))
    keys = [r["class_key"] for r in recs]
    if len(set(keys)) != len(keys):
        bad.append(("-", "classes not pairwise distinct"))
    in_reps = None
    if REPS.exists():
        with gzip.open(REPS, "rt") as fh:
            reps = {json.loads(l)["class"] for l in fh}
        in_reps = [r["class_key"] for r in recs if r["class_key"] in reps]
    print(json.dumps({"schemes": len(recs), "distinct_classes": len(set(keys)), "failures": bad,
                      "in_saved_f3_reps": in_reps}))
    return 0 if not bad else 1


if __name__ == "__main__":
    if sys.argv[1:2] == ["build"]:
        build()
    elif sys.argv[1:2] == ["check"]:
        sys.exit(check(*sys.argv[2:3]))
    else:
        raise SystemExit(__doc__)
