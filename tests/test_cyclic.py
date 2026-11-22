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


# ---- 006: classes of matrices and the class cubes ----------------------------------

from matmul.cyclic import (  # noqa: E402
    class_index, class_reps, conj, gl_generators, matrix_classes, rotate, split_orbits, to_class_form, to_rep,
)
from matmul.scheme import Scheme  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results" / "f3"


def known_cyclic():
    """The BILR scheme and the two 005 rank-23 schemes, as (orbits, fixed)."""
    out = [bilr()]
    for f in sorted(RESULTS.glob("005_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out.append(split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {})))
    return out


def test_split_orbits_inverts_expand():
    for orbits, fixed in known_cyclic():
        s = expand(3, orbits, fixed)
        assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s)
    assert len(known_cyclic()) == 3


@pytest.mark.timeout(30)
def test_class_counts():
    sim = matrix_classes(3, signed=False)
    signed = matrix_classes(3)
    assert len(sim) == 39  # q^3 + q^2 + q similarity classes of 3x3 matrices over F_q
    assert len(signed) == 22 and len(class_reps()) == 21
    assert sum(map(len, sim.values())) == sum(map(len, signed.values())) == 3 ** 9
    # each signed class is one similarity class or the union of x's and -x's
    for rep, members in signed.items():
        parts = [r for r, m in sim.items() if set(m) <= set(members)]
        assert sum(len(sim[r]) for r in parts) == len(members) and len(parts) in (1, 2)


def test_generators_are_invertible_and_conj_commutes_with_sigma():
    rng = random.Random(1)
    for g in gl_generators():
        t = tuple([[rng.randrange(3) for _ in range(3)] for _ in range(3)] for _ in range(3))
        assert conj(g, sigma(t)) == sigma(conj(g, t))


@pytest.mark.timeout(30)
def test_conjugation_keeps_known_schemes_correct_and_cyclic():
    rng = random.Random(2)
    gens = gl_generators()
    for orbits, fixed in known_cyclic():
        for _ in range(2):
            M = [[int(i == j) for j in range(3)] for i in range(3)]
            for _ in range(6):
                g = rng.choice(gens)
                M = [[sum(a * b for a, b in zip(r, c)) % 3 for c in zip(*g)] for r in M]
            o2 = [conj(M, t) for t in orbits]
            f2 = [conj(M, (x, x, T(x)))[0] for x in fixed]
            s = expand(3, o2, f2)
            assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s)


def test_to_rep_conjugates_to_the_representative():
    rng = random.Random(3)
    reps = {tuple(e for r in x for e in r) for x in class_reps()}
    for _ in range(20):
        x = [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
        if not any(any(r) for r in x):
            continue
        rep, s, M = to_rep(x)
        sx = [[e * s % 3 for e in r] for r in x]
        assert conj(M, (sx, sx, sx))[0] == rep and tuple(e for r in rep for e in r) in reps


@pytest.mark.timeout(120)
def test_known_schemes_in_class_form_satisfy_their_cube():
    for orbits, fixed in known_cyclic():
        for slot in range(len(orbits)):  # every orbit can serve as the pinned one
            o2, f2, rep = to_class_form(orbits, fixed, slot=slot)
            s = expand(3, o2, f2)
            assert verify(s, mod=3).ok and is_cyclic_invariant(s)
            cnf = CyclicF3CNF(3, 7, 2, pin_u1=rep)
            assert cnf.satisfied(cnf.witness(o2, f2)) == []
        # and the cube of another class rejects it
        other = next(r for r in class_reps() if r != rep)
        cnf = CyclicF3CNF(3, 7, 2, pin_u1=other)
        assert cnf.satisfied(cnf.witness(o2, f2))


def test_no_orbits_is_trivially_unsat_at_3x3():
    assert CyclicF3CNF(3, 0, 22).trivially_unsat


@pytest.mark.timeout(60)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_2x2_cubes_find_rank7_and_all_refute_rank6(tmp_path):
    f = tmp_path / "x.cnf"
    codes = {}
    for (a, b) in [(2, 1), (1, 3)]:
        for rep in class_reps(2):
            cnf = CyclicF3CNF(2, a, b, pin_u1=rep)
            f.write_text(cnf.dimacs())
            r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=8)
            codes.setdefault((a, b), []).append(r.returncode)
            if r.returncode == 10:
                model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
                s = expand(2, *cnf.decode(model))
                assert s.rank == 7 and verify(s, mod=3).ok and is_cyclic_invariant(s)
    assert 10 in codes[(2, 1)]
    assert set(codes[(1, 3)]) == {20}


def _rank(x):
    from matmul.invariants import rank_mod
    return rank_mod(x, 3)


@pytest.mark.timeout(60)
@pytest.mark.parametrize("min_class", [1, 2])
def test_min_class_gadget_matches_rank_and_trace(min_class):
    # only the class clauses are checked: random slot contents need not satisfy the Brent equations
    rng = random.Random(10 + min_class)
    cnf = CyclicF3CNF(3, 2, 1, pin_u1=class_reps()[min_class], min_class=min_class)
    lo, hi = cnf.class_clauses
    assert hi > lo
    rank1 = lambda x: _rank(x) == 1
    bad = (lambda x: rank1(x) and sum(x[i][i] for i in range(3)) % 3) if min_class == 1 else rank1
    seen = set()
    for trial in range(300):
        def mat():
            if rng.random() < 0.4:  # rank one: outer product
                a, b = [rng.randrange(3) for _ in range(3)], [rng.randrange(3) for _ in range(3)]
                return [[a[i] * b[j] % 3 for j in range(3)] for i in range(3)]
            return [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
        orbits = [(class_reps()[min_class], mat(), mat()), (mat(), mat(), mat())]
        fixed = [mat()]
        w = cnf.witness(orbits, fixed)
        ok = all(any(w.get(abs(l)) == (l > 0) for l in c) for c in cnf.clauses[lo:hi])
        want = not any(bad(m) for k, t in enumerate(orbits) for j, m in enumerate(t) if (k, j) != (0, 0))
        assert ok == want
        seen.add(want)
    assert seen == {True, False}


@pytest.mark.timeout(120)
def test_known_schemes_in_minimal_class_form_satisfy_their_cube():
    for orbits, fixed in known_cyclic():
        o2, f2, rep = to_class_form(orbits, fixed, minimal=True)
        c = class_reps().index(rep)
        assert verify(expand(3, o2, f2), mod=3).ok and is_cyclic_invariant(expand(3, o2, f2))
        cnf = CyclicF3CNF(3, 7, 2, pin_u1=rep, min_class=c)
        assert cnf.satisfied(cnf.witness(o2, f2)) == []
        # pinned to a class-1 member instead, the class-1 cube must reject it (other orbits hold class 0)
        if c == 0:
            for slot, t in enumerate(o2):
                for r in range(3):
                    if class_index(rotate(t, r)[0]) == 1:
                        o3 = list(o2)
                        o3[slot] = rotate(t, r)
                        o4, f4, rep1 = to_class_form(o3, f2, slot=slot)
                        cnf1 = CyclicF3CNF(3, 7, 2, pin_u1=rep1, min_class=1)
                        assert cnf1.satisfied(cnf1.witness(o4, f4))
                        break
