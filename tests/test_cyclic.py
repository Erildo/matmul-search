"""Cyclic-invariant schemes and their CNF encoding over F3 (attempt 005)."""

import json
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from matmul.cyclic import CyclicF3CNF, T, expand, is_cyclic_invariant, normal_form, sigma
from matmul.known import strassen
from matmul.verify import verify

DATA = Path(__file__).resolve().parent / "data"


def bilr():
    d = json.loads((DATA / "bilr_2fix_z3.json").read_text())
    return [(A, B, T(C)) for A, B, C in d["orbits"]], d["fixed"]


def mat_sum(ms, n):
    return [[sum(m[i][j] for m in ms) % 3 for j in range(n)] for i in range(n)]


def test_sigma_has_order_three():
    rng = random.Random(0)
    t = tuple([[rng.randrange(3) for _ in range(3)] for _ in range(3)] for _ in range(3))
    assert sigma(sigma(sigma(t))) == t


def test_bilr_scheme_is_correct_and_cyclic():
    orbits, fixed = bilr()
    s = expand(3, orbits, fixed)
    assert s.rank == 23
    assert verify(s).ok and verify(s, mod=3).ok
    assert is_cyclic_invariant(s)


def test_fixed_matrices_sum_to_identity():
    # the necessary condition derived in attempts/005, on the two cyclic schemes we know
    I3 = [[int(i == j) for j in range(3)] for i in range(3)]
    assert mat_sum(bilr()[1], 3) == I3
    s = strassen()  # term 0 is the fixed term (I, I, I); the rest form two orbits
    assert is_cyclic_invariant(s)
    assert [[int(x) for x in r] for r in s.U[0]] == [[1, 0], [0, 1]]


def test_breaking_one_term_breaks_invariance():
    orbits, fixed = bilr()
    s = expand(3, orbits, fixed)
    s.W[0][0][0] += 1
    assert not is_cyclic_invariant(s)


def test_equation_count():
    assert len(CyclicF3CNF(3, 7, 2).equations) == 249  # (729 - 9) / 3 + 9
    assert len(CyclicF3CNF(2, 2, 1).equations) == 24   # (64 - 4) / 3 + 4


def test_all_fixed_is_trivially_unsat_but_orbits_are_not():
    assert CyclicF3CNF(2, 0, 6).trivially_unsat
    assert not CyclicF3CNF(2, 1, 3).trivially_unsat
    assert not CyclicF3CNF(3, 1, 1).trivially_unsat


@pytest.mark.timeout(60)
@pytest.mark.parametrize("break_sym", [False, True])
def test_known_scheme_satisfies_the_cnf(break_sym):
    orbits, fixed = bilr()
    if break_sym:
        orbits, fixed = normal_form(orbits, fixed)
        assert verify(expand(3, orbits, fixed), mod=3).ok
    cnf = CyclicF3CNF(3, 7, 2, break_sym=break_sym)
    assert cnf.satisfied(cnf.witness(orbits, fixed)) == []


@pytest.mark.timeout(30)
def test_known_scheme_with_spare_zero_slots_satisfies_larger_instance():
    orbits, fixed = normal_form(*bilr())
    zero = [[0] * 3 for _ in range(3)]
    cnf = CyclicF3CNF(3, 8, 2)
    cand = normal_form(orbits + [(zero, zero, zero)], fixed)
    assert cnf.satisfied(cnf.witness(*cand)) == []


def test_wrong_scheme_violates_the_cnf():
    orbits, fixed = bilr()
    orbits[0][0][0][1] = (orbits[0][0][0][1] + 1) % 3
    cnf = CyclicF3CNF(3, 7, 2, break_sym=False)
    assert cnf.satisfied(cnf.witness(orbits, fixed))


def test_decode_inverts_witness():
    orbits, fixed = normal_form(*bilr())
    cnf = CyclicF3CNF(3, 7, 2)
    w = cnf.witness(orbits, fixed)
    model = [v if w.get(v) else -v for v in range(1, cnf.n_primary + 1)]
    o2, f2 = cnf.decode(model)
    assert [tuple(map(list, o)) for o in o2] == [tuple(map(list, o)) for o in orbits] and f2 == fixed


@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_solver_finds_2x2_rank7_and_refutes_rank6(tmp_path):
    for (a, b), want in [((2, 1), 10), ((1, 3), 20)]:
        cnf = CyclicF3CNF(2, a, b)
        f = tmp_path / "x.cnf"
        f.write_text(cnf.dimacs())
        r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=8)
        assert r.returncode == want
        if want == 10:
            model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
            s = expand(2, *cnf.decode(model))
            assert s.rank == 7 and verify(s, mod=3).ok and is_cyclic_invariant(s)
