"""Attempt 010, stream "local": frozen-slot S3 encoder and the tau_M scan (src/matmul/s3_local.py)."""

import json
import random
import shutil
import subprocess
from itertools import combinations
from pathlib import Path

import pytest

from matmul.cyclic import expand
from matmul.invariants import class_key
from matmul.s3 import is_s3_invariant
from matmul.s3_local import (S3LocalF3CNF, group_slots, s3_scheme, s3_slots, slot_list, tau_M_list,
                             to_standard_s3)
from matmul.scheme import Scheme
from matmul.verify import verify

ROOT = Path(__file__).resolve().parent.parent
S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"


def known_items():
    return slot_list(s3_slots(Scheme.load(S3_006)))


def split_of(items):
    g = group_slots(items)
    return tuple(len(x) for x in g)


def test_tau_M_scan_and_change_of_basis_on_006():
    recs = [json.loads(l) for l in (ROOT / "results" / "f3" / "cyclic_rank23.jsonl").read_text().splitlines()]
    rec = next(r for r in recs if r["class_key"] == "eebc48753638175a")
    s = expand(3, [tuple(x) for x in rec["orbits"]], rec["fixed"])
    Ms = tau_M_list(s)
    assert Ms == [[[1, 0, 0], [0, 2, 1], [0, 1, 0]]]
    s2 = to_standard_s3(s, Ms[0])
    assert s2 is not None and verify(s2, mod=3).ok and is_s3_invariant(s2)
    assert class_key(s2, 3) == "eebc48753638175a"
    # the BILR class has no tau_M
    rec = next(r for r in recs if r["class_key"] == "094008dc304fad88")
    assert tau_M_list(expand(3, [tuple(x) for x in rec["orbits"]], rec["fixed"])) == []


def test_s3_scheme_round_trip():
    items = known_items()
    assert split_of(items) == (1, 5, 0, 2)
    s = s3_scheme(group_slots(items))
    assert s.rank == 23 and verify(s, mod=3).ok and is_s3_invariant(s)


@pytest.mark.timeout(120)
def test_known_scheme_satisfies_frozen_instances_at_rank_23():
    items = known_items()
    rng = random.Random(7)
    subsets = [()] + [c for k in (1, 2, 3, 4, 7, 8) for c in combinations(range(len(items)), k)]
    for idx in rng.sample(subsets[1:], 12) + [(), tuple(range(8))]:
        frozen = [items[i] for i in idx]
        free = [items[i] for i in range(len(items)) if i not in idx]
        enc = S3LocalF3CNF(3, group_slots(frozen), split_of(free))
        w = enc.witness(*enc.slot_contents(group_slots(frozen), group_slots(free)))
        assert enc.satisfied(w) == [], idx
        # with a larger free split (extra zero slots) it still holds
        sp = tuple(x + 1 for x in split_of(free))
        enc = S3LocalF3CNF(3, group_slots(frozen), sp)
        w = enc.witness(*enc.slot_contents(group_slots(frozen), group_slots(free)))
        assert enc.satisfied(w) == [], idx


def test_forbid_excludes_exactly_the_known_completion():
    items = known_items()
    frozen, free = items[:5], items[5:]
    enc = S3LocalF3CNF(3, group_slots(frozen), split_of(free), forbid=group_slots(free))
    w = enc.witness(*enc.slot_contents(group_slots(frozen), group_slots(free)))
    bad = enc.satisfied(w)
    assert bad == [enc.clauses[enc.forbid_clauses[1] - 1]]  # only the final "some slot differs" clause
    # forbidding other orbits keeps the known completion
    enc = S3LocalF3CNF(3, group_slots(frozen), split_of(free), forbid=group_slots(frozen))
    w = enc.witness(*enc.slot_contents(group_slots(frozen), group_slots(free)))
    assert enc.satisfied(w) == []


@pytest.mark.skipif(shutil.which("cadical") is None, reason="needs cadical")
@pytest.mark.timeout(120)
def test_solver_recovers_known_scheme_with_two_free_size3_slots_and_rank22_unsat_small():
    items = known_items()
    three_idx = [i for i, (ty, _) in enumerate(items) if ty == 3]
    free = [items[i] for i in three_idx[:2]]
    frozen = [x for x in items if not any(x is y for y in free)]
    enc = S3LocalF3CNF(3, group_slots(frozen), split_of(free))
    r = subprocess.run(["cadical", "-q", "/dev/stdin"], input=enc.dimacs(), capture_output=True, text=True)
    assert r.returncode == 10
    model = [int(x) for line in r.stdout.splitlines() if line.startswith("v") for x in line.split()[1:]]
    s = expand(3, *enc.decode(model))
    assert s.rank == 23 and verify(s, mod=3).ok and is_s3_invariant(s)
