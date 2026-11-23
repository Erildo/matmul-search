"""Attempt 010: the transpose symmetry tau, the S3 ansatz over F3 and its CNF encoding."""

import json
import random
import shutil
import subprocess
from itertools import product
from pathlib import Path

import numpy as np
import pytest

from matmul.cyclic import _inv, _mm, expand, is_cyclic_invariant, sigma, split_orbits, T
from matmul.s3 import (S3F3CNF, is_s3_invariant, s3_normal_form, s3_split, s3_splits, scheme_s3_split, tau,
                       tensor_key)
from matmul.scheme import Scheme
from matmul.verify import verify

ROOT = Path(__file__).resolve().parent.parent
S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"


def load(path):
    d = json.loads(Path(path).read_text())
    return Scheme(3, 3, 3, d["U"], d["V"], d["W"], {})


def rand_mat(rng, n=3):
    return [[rng.randrange(3) for _ in range(n)] for _ in range(n)]


def test_tau_is_an_involution_inverting_sigma():
    rng = random.Random(1)
    for _ in range(50):
        t = tuple(rand_mat(rng) for _ in range(3))
        assert tau(tau(t)) == t
        assert tau(sigma(tau(t))) == sigma(sigma(t))


def test_tau_maps_correct_schemes_to_correct_schemes():
    recs = [json.loads(l) for l in (ROOT / "results" / "f3" / "cyclic_rank23.jsonl").read_text().splitlines()]
    for rec in recs[:3]:
        s = expand(3, [tuple(x) for x in rec["orbits"]], rec["fixed"])
        s2 = Scheme(3, 3, 3, [T(V) for V in s.V], [T(U) for U in s.U], [T(W) for W in s.W], {})
        assert verify(s2, mod=3).ok


@pytest.mark.timeout(60)
def test_nonsingular_symmetric_forms_are_congruent_to_plus_or_minus_identity():
    det = lambda M: round(np.linalg.det(np.array(M))) % 3
    syms = [[[k[0], k[1], k[2]], [k[1], k[3], k[4]], [k[2], k[4], k[5]]] for k in product(range(3), repeat=6)]
    nonsing = {tuple(map(tuple, M)) for M in syms if det(M)}
    assert len(nonsing) == 468
    reach = set()
    for k in product(range(3), repeat=9):
        N = [list(k[0:3]), list(k[3:6]), list(k[6:9])]
        if det(N) == 0:
            continue
        G = _mm(N, T(N))
        reach.add(tuple(map(tuple, G)))
        reach.add(tuple(tuple(2 * e % 3 for e in r) for r in G))
    assert reach == nonsing


@pytest.mark.timeout(60)
def test_tau_fixed_rank_one_tensors_have_the_size3_form_2x2():
    # every nonzero 2x2 rank-one tensor t with tau(t) = t (as a tensor) is (u, u^T, w) with w = w^T after rescaling
    mats = [np.array(k).reshape(2, 2) for k in product(range(3), repeat=4) if any(k)]
    fixed = 0
    for u in mats:
        for v in mats:
            for w in mats:
                t = np.einsum("a,b,c->abc", u.ravel(), v.ravel(), w.ravel()) % 3
                tt = np.einsum("a,b,c->abc", v.T.ravel(), u.T.ravel(), w.T.ravel()) % 3
                if not np.array_equal(t, tt):
                    continue
                fixed += 1
                ok = False
                for a, b in product((1, 2), repeat=2):
                    u2, v2, w2 = a * u % 3, b * v % 3, a * b * w % 3
                    ok |= np.array_equal(v2, u2.T) and np.array_equal(w2, w2.T)
                assert ok, (u, v, w)
    assert fixed > 0


@pytest.mark.timeout(60)
def test_s3_form_of_the_006_scheme():
    s = load(S3_006)
    assert s.rank == 23 and verify(s, mod=3).ok and is_cyclic_invariant(s) and is_s3_invariant(s)
    assert scheme_s3_split(s) == (1, 5, 0, 2)
    o, f = s3_normal_form(*s3_split(*split_orbits(s)))
    s2 = expand(3, o, f)
    assert verify(s2, mod=3).ok and is_s3_invariant(s2)
    cnf = S3F3CNF(3, 1, 5, 0, 2)
    assert cnf.satisfied(cnf.witness(o, f)) == []
    # breaking a tie: change the size-6 partner's w, keeping the scheme sigma-shaped but not tau-shaped
    o3 = [list(t) for t in o]
    o3[1][2] = [[(e + 1) % 3 for e in r] for r in o3[1][2]]
    lo, hi = cnf.tie_clauses
    bad = cnf.satisfied(cnf.witness([tuple(t) for t in o3], f))
    assert bad and any(c in cnf.clauses[lo:hi] for c in bad)
    # a larger instance (spare zero slots) also accepts it
    big = S3F3CNF(3, 2, 5, 1, 2)
    zero = tuple([[0] * 3 for _ in range(3)] for _ in range(3))
    orbs = [zero, o[0], zero, o[1]] + o[2:]
    fx = [[[0] * 3 for _ in range(3)]] * 2 + f
    assert big.satisfied(big.witness(orbs, fx)) == []


def test_not_tau_invariant_is_detected():
    recs = [json.loads(l) for l in (ROOT / "results" / "f3" / "cyclic_rank23.jsonl").read_text().splitlines()]
    s = expand(3, [tuple(x) for x in recs[0]["orbits"]], recs[0]["fixed"])  # BILR: cyclic, not tau-invariant
    assert is_cyclic_invariant(s) and not is_s3_invariant(s)


def test_split_counts():
    assert len(s3_splits(22, 3)) == 60 and len(s3_splits(23)) == 98
    assert all(6 * a + 3 * b + 2 * c + d == 22 and 2 * a + b >= 3 for a, b, c, d in s3_splits(22, 3))


def _solve(cnf, f, extra=()):
    text = cnf.dimacs()
    if extra:
        head = f"p cnf {cnf.nv} {len(cnf.clauses)}"
        text = text.replace(head, f"p cnf {cnf.nv} {len(cnf.clauses) + len(extra)}", 1) + "".join(
            " ".join(map(str, c)) + " 0\n" for c in extra)
    f.write_text(text)
    out = subprocess.run(["cadical", "-q", str(f)], capture_output=True, text=True, timeout=20)
    model = [int(x) for l in out.stdout.splitlines() if l.startswith("v") for x in l.split()[1:]]
    return out.returncode, model


@pytest.mark.timeout(120)
@pytest.mark.skipif(shutil.which("cadical") is None, reason="cadical not installed")
def test_2x2_s3_encoder_agrees_with_enumeration(tmp_path):
    """2x2 over F3, standard tau. The S3 encoder finds rank 8 (the naive scheme is S3-invariant) and refutes ranks 7 and 6.

    Independent check of the rank-7 refutation: enumerate every solution of the cyclic rank-7 encoder (split (2,1),
    normalised) and test each for invariance under tau_M for all nonsingular symmetric M. None is invariant under a
    tau_M with det M a square (the class of the standard tau = tau_I); all hits have det M = 2, a nonsquare. In even
    dimension -M has the same determinant class as M, so for 2x2 these are two inequivalent S3's (for 3x3 they merge).
    """
    from matmul.cyclic import CyclicF3CNF, conj
    f = tmp_path / "x.cnf"
    E = lambda i, j: [[int((r, c) == (i, j)) for c in range(2)] for r in range(2)]
    trip = list(product(range(2), repeat=3))
    naive = Scheme(2, 2, 2, [E(i, j) for i, j, l in trip], [E(j, l) for i, j, l in trip], [E(i, l) for i, j, l in trip], {})
    assert verify(naive, mod=3).ok and is_s3_invariant(naive)
    codes = {}
    for r in (8, 7, 6):
        for sp in s3_splits(r):
            cnf = S3F3CNF(2, *sp)
            code, model = (20, None) if cnf.trivially_unsat else _solve(cnf, f)
            codes.setdefault(r, []).append(code)
            if code == 10:
                s = expand(2, *cnf.decode(model))
                assert s.rank <= r and verify(s, mod=3).ok and is_s3_invariant(s)
    assert 10 in codes[8] and set(codes[7]) == {20} and set(codes[6]) == {20}
    # the enumeration
    syms = [[[a, b], [b, c]] for a, b, c in product(range(3), repeat=3) if (a * c - b * b) % 3]
    cnf = CyclicF3CNF(2, 2, 1)
    prim = [v for slot in cnf.orb for m in slot for pq in m for v in pq] + [v for x in cnf.fix for pq in x for v in pq]
    blocks, dets = [], set()
    while True:
        code, model = _solve(cnf, f, blocks)
        if code != 10:
            break
        assign = {abs(l): l > 0 for l in model}
        blocks.append([-v if assign.get(v) else v for v in prim])
        s = expand(2, *cnf.decode(model))
        terms = [tuple([[int(e) % 3 for e in row] for row in m] for m in (A, B, C)) for A, B, C in zip(s.U, s.V, s.W)]
        S = {tensor_key(t) for t in terms}
        dets |= {(M[0][0] * M[1][1] - M[0][1] ** 2) % 3 for M in syms if {tensor_key(conj(M, tau(t))) for t in terms} == S}
    assert code == 20 and len(blocks) == 36
    assert dets == {2}
