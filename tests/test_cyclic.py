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


@pytest.mark.timeout(60)
@pytest.mark.parametrize("c", [10, 19])
def test_charpoly_exclusion_gadget(c):
    from matmul.cyclic import charpoly, charpoly_exclusions
    rng = random.Random(20 + c)
    forb = set(charpoly_exclusions(c))
    reps = class_reps()
    # soundness of the table: no class >= c has a member with a forbidden characteristic polynomial
    for rep, members in matrix_classes(3).items():
        if any(rep) and reps.index([list(rep[i:i + 3]) for i in (0, 3, 6)]) >= c:
            assert not any(charpoly([list(m[i:i + 3]) for i in (0, 3, 6)]) in forb for m in members)
    cnf = CyclicF3CNF(3, 2, 1, pin_u1=reps[c], min_class=c, exclude_charpoly=True)
    lo, hi = cnf.class_clauses
    seen = set()
    for _ in range(300):
        mat = lambda: ([[0] * 3 for _ in range(3)] if rng.random() < 0.5
                       else [[rng.randrange(3) for _ in range(3)] for _ in range(3)])
        orbits = [(reps[c], mat(), mat()), (mat(), mat(), mat())]
        w = cnf.witness(orbits, [mat()])
        ok = all(any(w.get(abs(l)) == (l > 0) for l in cl) for cl in cnf.clauses[lo:hi])
        mats = [m for k, t in enumerate(orbits) for j, m in enumerate(t) if (k, j) != (0, 0)]
        want = not any(_rank(m) == 1 or charpoly(m) in forb for m in mats)
        assert ok == want
        seen.add(want)
    assert seen == {True, False}


# ---- 007: all orbits nonzero, and level-2 cubes by characteristic polynomial ------------

from itertools import product  # noqa: E402

from matmul.cyclic import (  # noqa: E402
    charpoly, charpoly_exclusions, cp1_choices, cp1_exclusions, cp_canon, cp_neg, to_level2_form,
)


def known_cyclic_all():
    """known_cyclic() plus the 006 rank-23 scheme."""
    out = known_cyclic()
    for f in sorted(RESULTS.glob("006_3x3x3_rank23_*.json")):
        d = json.loads(f.read_text())
        out.append(split_orbits(Scheme(3, 3, 3, d["U"], d["V"], d["W"], {})))
    assert len(out) == 4
    return out


def level2_cnf(c, p, a=7, b=2):
    return CyclicF3CNF(3, a, b, pin_u1=class_reps()[c], min_class=c, exclude_charpoly=True,
                       all_nonzero=True, pin_cp1=p)


def test_level1_cube_size_is_unchanged_from_006():
    cnf = CyclicF3CNF(3, 2, 16, pin_u1=class_reps()[0], min_class=0)
    assert (cnf.nv, len(cnf.clauses)) == (21339, 135996)


@pytest.mark.timeout(60)
def test_cp_tables_are_consistent():
    for q in product(range(3), repeat=3):
        assert cp_neg(cp_neg(q)) == q and cp_canon(q) == cp_canon(cp_neg(q))
    x = [[1, 2, 0], [0, 1, 1], [2, 0, 2]]
    assert charpoly([[-e % 3 for e in r] for r in x]) == cp_neg(charpoly(x))
    for c in range(21):
        ch = cp1_choices(c)
        assert ch and all(p == cp_canon(p) for p in ch)
        assert not set(ch) & set(charpoly_exclusions(c))  # carried by a class >= c, so allowed at level 1
        for p in ch:
            ex = set(cp1_exclusions(p))
            assert p not in ex and cp_neg(p) not in ex
            assert ex == {q for q in product(range(3), repeat=3) if cp_canon(q) < p}


@pytest.mark.timeout(120)
def test_known_schemes_in_level2_form_satisfy_their_cube():
    for orbits, fixed in known_cyclic_all():
        o2, f2, rep, p = to_level2_form(orbits, fixed)
        c = class_reps().index(rep)
        s = expand(3, o2, f2)
        assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s)
        assert p in cp1_choices(c) and charpoly(o2[1][0]) == p
        cnf = level2_cnf(c, p)
        assert cnf.satisfied(cnf.witness(o2, f2)) == []


def _slot1_from(o2, f2, i, r):
    """Put member r of orbit i (i >= 1) into slot 1, normalised as to_level2_form does, whatever its polynomial."""
    neg = lambda m: [[-e % 3 for e in row] for row in m]
    first = lambda m: next((e % 3 for row in m for e in row if e % 3), 1)
    u, v, w = rotate(o2[i], r)
    q = cp_canon(charpoly(u))
    if charpoly(u) != q or (q == cp_neg(q) and first(u) == 2):
        u, v = neg(u), neg(v)
    if first(v) == 2:
        v, w = neg(v), neg(w)
    rest, f3 = normal_form([t for k, t in enumerate(o2[1:], 1) if k != i], f2)
    return [o2[0], (u, v, w)] + rest, f3, q


@pytest.mark.timeout(240)
def test_non_minimal_slot1_is_rejected_by_the_level2_clauses():
    checked = 0
    for orbits, fixed in known_cyclic_all():
        o2, f2, rep, p = to_level2_form(orbits, fixed)
        c = class_reps().index(rep)
        for i in range(2, len(o2)):
            for r in range(3):
                o3, f3, q = _slot1_from(o2, f2, i, r)
                if q <= p:
                    continue
                assert verify(expand(3, o3, f3), mod=3).ok
                cnf = level2_cnf(c, q)
                lo, hi = cnf.cp1_clauses
                bad = cnf.satisfied(cnf.witness(o3, f3))
                assert bad and all(cl in cnf.clauses[lo:hi] for cl in bad)
                checked += 1
                break
            else:
                continue
            break
    assert checked >= 2


@pytest.mark.timeout(60)
def test_cp1_gadget_matches_brute_force():
    rng = random.Random(30)
    c = 2
    reps = class_reps()
    for p in (cp1_choices(c)[0], cp1_choices(c)[len(cp1_choices(c)) // 2], cp1_choices(c)[-1]):
        cnf = CyclicF3CNF(3, 3, 1, pin_u1=reps[c], pin_cp1=p)
        lo, hi = cnf.cp1_clauses
        ex = set(cp1_exclusions(p))
        seen = set()
        for _ in range(300):
            mat = lambda: [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
            u1 = mat()
            if rng.random() < 0.3:  # make the pin hold often enough to see both outcomes
                for _ in range(200):
                    if charpoly(u1) == p:
                        break
                    u1 = mat()
            # other matrices: mostly u1, -u1 or u1^T (polynomial p or its negative, never excluded)
            near = lambda: (rng.choice([u1, [[-e % 3 for e in r] for r in u1], T(u1)]) if rng.random() < 0.9
                            else mat())
            orbits = [(reps[c], mat(), mat()), (u1, near(), near()), (near(), near(), near())]
            w = cnf.witness(orbits, [mat()])
            ok = all(any(w.get(abs(l)) == (l > 0) for l in cl) for cl in cnf.clauses[lo:hi])
            want = charpoly(u1) == p and not any(charpoly(m) in ex for t in orbits[1:] for m in t)
            assert ok == want
            seen.add(want)
        assert seen == {True, False}


def test_all_nonzero_rejects_a_zero_matrix():
    orbits, fixed = known_cyclic_all()[0]
    o2, f2, rep, p = to_level2_form(orbits, fixed)
    cnf = level2_cnf(class_reps().index(rep), p)
    lo, hi = cnf.nonzero_clauses
    assert hi - lo == 7 * 3
    zero = [[0] * 3 for _ in range(3)]
    for k, m in [(0, 1), (3, 0), (6, 2)]:
        o3 = [list(t) for t in o2]
        o3[k][m] = zero
        w = cnf.witness([tuple(t) for t in o3], f2)
        assert not all(any(w.get(abs(l)) == (l > 0) for l in cl) for cl in cnf.clauses[lo:hi])


@pytest.mark.timeout(120)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_solver_recovers_slots_1_and_2_in_their_level2_cube(tmp_path):
    # fix every slot of a known scheme except orbit slots 1 and 2; the solver must complete them in the cube
    orbits, fixed = known_cyclic_all()[0]
    o2, f2, rep, p = to_level2_form(orbits, fixed)
    cnf = level2_cnf(class_reps().index(rep), p)
    w = cnf.witness(o2, f2)
    keep = [v for k in (0, 3, 4, 5, 6) for m in cnf.orb[k] for pq in m for v in pq]
    keep += [v for x in cnf.fix for pq in x for v in pq]
    cnf.clauses += [[v if w[v] else -v] for v in keep]
    f = tmp_path / "x.cnf"
    f.write_text(cnf.dimacs())
    r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=100)
    assert r.returncode == 10
    model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
    s = expand(3, *cnf.decode(model))
    assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s)


# ---- 008: neighbourhoods near(S, m, M) ------------------------------------------------

def _relabel(orbits, fixed, rng):
    """The same scheme, re-labelled: orbits permuted, members rotated, signs flipped, then put in normal form."""
    neg = lambda m: [[-e % 3 for e in r] for r in m]
    out = []
    for t in rng.sample(list(orbits), len(orbits)):
        u, v, w = rotate(t, rng.randrange(3))
        if rng.random() < 0.5:
            u, v = neg(u), neg(v)
        if rng.random() < 0.5:
            v, w = neg(v), neg(w)
        out.append((u, v, w))
    return normal_form(out, fixed)


@pytest.mark.timeout(120)
def test_seed_schemes_are_near_themselves_and_relabelled_copies_too():
    rng = random.Random(80)
    for orbits, fixed in known_cyclic_all():
        o, f = normal_form(orbits, fixed)
        assert verify(expand(3, o, f), mod=3).ok
        cnf = CyclicF3CNF(3, 7, 2, near=(orbits, 7, 7))
        assert cnf.satisfied(cnf.witness(o, f)) == []
        o2, f2 = _relabel(orbits, fixed, rng)
        assert verify(expand(3, o2, f2), mod=3).ok and o2 != o
        assert cnf.satisfied(cnf.witness(o2, f2)) == []
        cnf6 = CyclicF3CNF(3, 7, 2, near=(orbits, 0, 6))
        lo, hi = cnf6.near_clauses
        bad = cnf6.satisfied(cnf6.witness(o2, f2))
        assert bad and all(c in cnf6.clauses[lo:hi] for c in bad)


@pytest.mark.timeout(60)
def test_near_gadget_matches_a_direct_count():
    rng = random.Random(81)
    orbits, fixed = known_cyclic_all()[3]
    targets = set()
    for t in orbits:
        for r in range(3):
            (o,), _ = normal_form([rotate(t, r)], [])
            targets.add(tuple(e for m in o for row in m for e in row))
    for m, M in [(5, 7), (3, 6), (6, 6)]:
        cnf = CyclicF3CNF(3, 7, 2, near=(orbits, m, M))
        lo, hi = cnf.near_clauses
        seen = set()
        for _ in range(200):
            slots = []
            for _ in range(7):
                if rng.random() < 0.8:
                    (o,), _ = normal_form([rotate(rng.choice(orbits), rng.randrange(3))], [])
                else:
                    o = tuple([[rng.randrange(3) for _ in range(3)] for _ in range(3)] for _ in range(3))
                slots.append(o)
            w = cnf.witness(slots, fixed)
            ok = all(any(w.get(abs(l)) == (l > 0) for l in c) for c in cnf.clauses[lo:hi])
            cnt = sum(tuple(e % 3 for mm in t for row in mm for e in row) in targets for t in slots)
            assert ok == (m <= cnt <= M)
            seen.add(ok)
        assert seen == {True, False}


@pytest.mark.timeout(120)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_reduce71_at_distance_zero_is_unsat(tmp_path):
    orbits, fixed = known_cyclic_all()[0]
    cnf = CyclicF3CNF(3, 7, 1, near=(orbits, 7, 7))
    f = tmp_path / "x.cnf"
    f.write_text(cnf.dimacs())
    r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=100)
    assert r.returncode == 20


@pytest.mark.timeout(60)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_2x2_walk_step(tmp_path):
    # 2x2 rank 7 at split (2,1): from one solution, ask for a cyclic scheme sharing at most one orbit with it
    f = tmp_path / "x.cnf"
    cnf = CyclicF3CNF(2, 2, 1)
    f.write_text(cnf.dimacs())
    r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=30)
    model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
    S, fx = cnf.decode(model)
    cnf2 = CyclicF3CNF(2, 2, 1, near=(S, 0, 1))
    f.write_text(cnf2.dimacs())
    r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=30)
    assert r.returncode in (10, 20)
    if r.returncode == 10:
        model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
        o2, f2 = cnf2.decode(model)
        s = expand(2, o2, f2)
        assert s.rank == 7 and verify(s, mod=3).ok and is_cyclic_invariant(s)
        keys = lambda orbs: {tuple(e for m in normal_form([rotate(t, r)], [])[0][0] for row in m for e in row)
                             for t in orbs for r in range(3)}
        assert sum(tuple(e for m in t for row in m for e in row) in keys(S) for t in o2) <= 1
    print("2x2 walk step:", {10: "SAT", 20: "UNSAT"}[r.returncode])


@pytest.mark.timeout(60)
def test_freeze_keeps_chosen_orbits_and_breaks_symmetry_on_the_rest():
    rng = random.Random(82)
    orbits, fixed = known_cyclic_all()[1]
    for _ in range(5):
        keep = sorted(rng.sample(range(7), 5))
        frozen = [orbits[i] for i in keep]
        free, f2 = normal_form([orbits[i] for i in range(7) if i not in keep], fixed)
        cnf = CyclicF3CNF(3, 7, 2, freeze=frozen)
        assert cnf.satisfied(cnf.witness(frozen + free, f2)) == []
        # a different orbit in a frozen slot is rejected
        assert cnf.satisfied(cnf.witness([free[0]] + frozen[1:] + [frozen[0], free[1]], f2))


@pytest.mark.timeout(60)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_freeze_all_orbits_reduce71_is_unsat_and_walk_subset_is_sat(tmp_path):
    orbits, fixed = known_cyclic_all()[0]
    f = tmp_path / "x.cnf"
    f.write_text(CyclicF3CNF(3, 7, 1, freeze=orbits).dimacs())
    assert subprocess.run(["cadical", "-q", str(f)], capture_output=True, timeout=50).returncode == 20
    cnf = CyclicF3CNF(3, 7, 2, freeze=orbits[:5])
    f.write_text(cnf.dimacs())
    r = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=50)
    assert r.returncode == 10
    model = [int(x) for l in r.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
    s = expand(3, *cnf.decode(model))
    assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s)
