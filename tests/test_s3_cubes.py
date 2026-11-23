"""Tests for matmul.s3_cubes (attempt 010, stage C): the cube family covers every solution."""

import random
from collections import Counter
from pathlib import Path

import pytest

from matmul.cyclic import T, _key, _mm, _unkey, conj, expand, split_orbits
from matmul.s3 import S3NativeF3CNF, s3_normal_form, s3_split, tau, tensor_key
from matmul.s3_cubes import (RANK1, S3CubeCNF, cubes, cubes_cfg, det3, encoder_cfg, mclass, normalize, normalize_cfg,
                             orth_class_index, orth_classes, orth_group, orth_to_rep, rank3, sk, skkey, sk_to_rep)
from matmul.scheme import Scheme

ROOT = Path(__file__).resolve().parent.parent
S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"
CONFIGS = [("r1", False, False), ("r1", False, True), ("r1", True, False), ("r1", True, True),
           ("rank", False, False), ("rank", False, True)]


def known_slots():
    return s3_normal_form(*s3_split(*split_orbits(Scheme.load(S3_006))))


def violated(enc, orbits, fixed):
    assign = enc.witness(orbits, fixed)
    return [i for i, c in enumerate(enc.clauses) if not any(assign.get(abs(l)) == (l > 0) for l in c)]


def term_multiset(orbits, fixed, M=None):
    s = expand(3, orbits, fixed)
    terms = [tuple([[int(e) % 3 for e in r] for r in m] for m in t) for t in zip(s.U, s.V, s.W)]
    if M is not None:
        terms = [conj(M, t) for t in terms]
    return Counter(tensor_key(t) for t in terms)


def test_orthogonal_group_and_classes():
    I = [[int(i == j) for j in range(3)] for i in range(3)]
    assert len(orth_group()) == 48
    reps, index, _ = orth_classes()
    assert len(reps) == 480 and len(index) == 3 ** 9  # 479 nonzero classes up to sign, plus zero
    assert sum(1 for r in reps if any(r) and rank3(_unkey(r, 3)) == 1) == 17
    rng = random.Random(1)
    for _ in range(300):
        x = [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
        rep, s, M = orth_to_rep(x)
        assert _mm(T(M), M) == I
        Mi = T(M)
        assert _key(_mm(_mm(M, [[e * s % 3 for e in r] for r in x]), Mi)) == _key(rep)
        assert orth_class_index(rep) == orth_class_index(x)


def test_orthogonal_conjugation_commutes_with_tau():
    rng = random.Random(2)
    for k in orth_group():
        M = _unkey(k, 3)
        t = tuple([[rng.randrange(3) for _ in range(3)] for _ in range(3)] for _ in range(3))
        assert tensor_key(conj(M, tau(t))) == tensor_key(tau(conj(M, t)))


def _rand_mat(rng, sym=False):
    r = rng.random()
    if r < 0.15:
        m = [[0] * 3 for _ in range(3)]
    elif r < 0.6:  # rank <= 1
        a = [rng.randrange(3) for _ in range(3)]
        b = a if sym else [rng.randrange(3) for _ in range(3)]
        s = rng.choice((1, 2)) if sym else 1
        m = [[a[i] * b[j] * s % 3 for j in range(3)] for i in range(3)]
    else:
        m = [[rng.randrange(3) for _ in range(3)] for _ in range(3)]
    if sym:
        m = [[m[min(i, j)][max(i, j)] for j in range(3)] for i in range(3)]
    return m


def _rand_slots(rng, split):
    a, b, c, d = split
    six = [tuple(_rand_mat(rng) for _ in range(3)) for _ in range(a)]
    three = []
    for _ in range(b):
        u = _rand_mat(rng)
        three.append((u, T(u), _rand_mat(rng, sym=True)))
    two = [_rand_mat(rng) for _ in range(c)]
    one = [_rand_mat(rng, sym=True) for _ in range(d)]
    return six + [tau(t) for t in six] + three, two + [T(x) for x in two] + one


@pytest.mark.timeout(240)
@pytest.mark.parametrize("split", [(1, 5, 0, 2), (2, 3, 0, 1), (0, 7, 0, 1), (3, 1, 0, 1), (1, 2, 1, 1)])
def test_normal_form_on_random_contents(split):
    """Every slot assignment has an image (same tensor terms up to g_M) that meets one cube's extra clauses."""
    rng = random.Random(hash(split) % 1000)
    for gran, use_w, pin in CONFIGS:
        all_cubes = set(cubes(split, gran, use_w, pin))
        for _ in range(4):
            orbits, fixed = _rand_slots(rng, split)
            o2, f2, cube, M = normalize(orbits, fixed, split, gran, use_w, pin)
            assert cube in all_cubes
            assert term_multiset(o2, f2) == term_multiset(orbits, fixed, M)
            a, b = split[0], split[1]
            for j in range(b):
                u, v, w = o2[2 * a + j]
                assert v == T(u) and w == T(w)
            enc = S3CubeCNF(split, cube, gran, use_w)
            bad = violated(enc, o2, f2)
            assert all(i < enc.base_clauses for i in bad)  # only equations / fixed sum may fail


@pytest.mark.timeout(120)
def test_known_scheme_satisfies_its_cube():
    orbits, fixed = known_slots()
    split = (1, 5, 0, 2)
    for gran, use_w, pin in CONFIGS:
        o2, f2, cube, M = normalize(orbits, fixed, split, gran, use_w, pin)
        assert cube in cubes(split, gran, use_w, pin)
        enc = S3CubeCNF(split, cube, gran, use_w)
        assert violated(enc, o2, f2) == []
        # and it is excluded from a different cube with the same pin status
        other = next(cb for cb in cubes(split, gran, use_w, pin) if cb != cube and (cb[2] is None) == (cube[2] is None))
        enc2 = S3CubeCNF(split, other, gran, use_w)
        assert violated(enc2, o2, f2) != []


def test_cube_counts():
    assert len(cubes((0, 7, 0, 1), "r1")) == 8
    assert len(cubes((0, 7, 0, 1), "r1", pin=True)) == 7 * 17 + 1
    assert len(cubes((1, 5, 0, 2), "r1")) == 4 * 6


def test_mclass():
    z = [[0] * 3 for _ in range(3)]
    I = [[int(i == j) for j in range(3)] for i in range(3)]
    e = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
    assert mclass(e, "r1") == RANK1 and mclass(z, "r1") == 1 and mclass(I, "r1") == 1
    assert [mclass(m, "rank") for m in (e, z, I)] == [0, 3, 2]


SK_CONFIGS = ["r1+sk", "r1+sk+n1", "r1+w+sk", "rank+sk+n1", "r1+pin+n1"]


@pytest.mark.timeout(240)
@pytest.mark.parametrize("split", [(1, 5, 0, 1), (2, 3, 0, 1), (3, 1, 0, 1), (1, 4, 1, 3), (2, 2, 1, 1)])
def test_sk_normal_form_on_random_contents(split):
    rng = random.Random(7 + sum(split))
    for cfg in SK_CONFIGS:
        fam = set(cubes_cfg(split, cfg))
        for _ in range(4):
            orbits, fixed = _rand_slots(rng, split)
            o2, f2, cube, M = normalize_cfg(orbits, fixed, split, cfg)
            assert cube in fam
            assert term_multiset(o2, f2) == term_multiset(orbits, fixed, M)
            enc = encoder_cfg(split, cube, cfg)
            bad = violated(enc, o2, f2)
            # random contents break the equations, the fixed sum and possibly N1 (all before base_clauses_n1)
            assert all(i < enc.base_clauses_n1 for i in bad)


def test_sk_classes():
    # SO(3, F3) orbits of nonzero skew vectors up to sign = number of nonzero coordinates
    tr = sk_to_rep()
    assert len(tr) == 26
    for p, (Mk, e) in tr.items():
        M = _unkey(Mk, 3)
        K = [[0, p[0], p[1]], [-p[0] % 3, 0, p[2]], [-p[1] % 3, -p[2] % 3, 0]]
        q = sk(conj(M, (K, K, K))[0])
        assert skkey(q) == skkey(p)


def s3_schemes(tmp_path):
    import subprocess
    import sys
    files = sorted((ROOT / "results" / "f3").glob("010*_3x3x3_rank23_*.json"))
    lad = tmp_path / "laderman_s3.json"
    r = subprocess.run([sys.executable, "-I", str(ROOT / "scripts" / "check_010_theory.py"), "--save-laderman", str(lad)],
                       capture_output=True, text=True, timeout=120)
    if lad.exists():
        files.append(lad)
    out = []
    for f in files:
        s = Scheme.load(f)
        from matmul.s3 import is_s3_invariant, scheme_s3_split
        if is_s3_invariant(s):
            out.append((f.name, s, scheme_s3_split(s)))
    return out


@pytest.mark.timeout(300)
def test_known_s3_schemes_satisfy_their_cubes_with_n1(tmp_path):
    """Every saved S3 rank-23 scheme (and Laderman in S3 form) meets all clauses of its own cube, N1 included."""
    schemes = s3_schemes(tmp_path)
    assert any(name == "laderman_s3.json" for name, _, _ in schemes)
    assert len(schemes) >= 3
    for name, s, split in schemes:
        orbits, fixed = s3_normal_form(*s3_split(*split_orbits(s)))
        # N1 directly
        tot = sum(det3(sk(t[0]), sk(t[1]), sk(t[2])) for t in orbits[:split[0]]) % 3
        assert tot == 2, name
        for cfg in SK_CONFIGS + ["r1+w+pin+n1"]:
            o2, f2, cube, _ = normalize_cfg(orbits, fixed, split, cfg)
            assert cube in cubes_cfg(split, cfg)
            enc = encoder_cfg(split, cube, cfg)
            assert violated(enc, o2, f2) == [], (name, cfg)
