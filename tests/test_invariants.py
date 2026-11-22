"""The invariant must not change under any symmetry, and must separate known-different schemes.

Each symmetry move is first checked to map a correct scheme to a correct scheme
(so the moves themselves are right), then checked to leave the invariant unchanged.
"""

import random
from pathlib import Path

import pytest

from matmul.import_fg import load_fg
from matmul.invariants import class_key, invariant, rank_mod, term_profile
from matmul.known import naive, strassen
from matmul.scheme import Scheme
from matmul.verify import verify

DATA = Path(__file__).resolve().parent / "data"
P = 3


def f3_a():
    return load_fg(DATA / "f3_333_r23_a.txt", mod=P)


def f3_b():
    return load_fg(DATA / "f3_333_r23_b.txt", mod=P)


# ---- small linear algebra mod p --------------------------------------------

def matmul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % p for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A):
    return [list(r) for r in zip(*A)]


def inverse(A, p):
    n = len(A)
    M = [list(map(int, row)) + [int(i == j) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] % p)
        M[c], M[piv] = M[piv], M[c]
        inv = pow(M[c][c], -1, p)
        M[c] = [(x * inv) % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]


def random_invertible(n, p, rng):
    while True:
        A = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
        if rank_mod(A, p) == n:
            return A


def ints(s, p):
    from matmul.invariants import _factors
    return _factors(s, p)


def build(n, U, V, W, meta=None):
    return Scheme(n, n, n, U, V, W, meta or {})


# ---- the symmetry moves ------------------------------------------------------

def basis_change(s, p, rng):
    n = s.n
    Pm, Q, R = (random_invertible(n, p, rng) for _ in range(3))
    Qit, Rit, Pi = transpose(inverse(Q, p)), transpose(inverse(R, p)), inverse(Pm, p)
    U, V, W = ints(s, p)
    U2 = [matmul(matmul(transpose(Pm), u, p), Qit, p) for u in U]
    V2 = [matmul(matmul(transpose(Q), v, p), Rit, p) for v in V]
    W2 = [matmul(matmul(Pi, w, p), R, p) for w in W]
    return build(n, U2, V2, W2)


def rescale(s, p, rng):
    U, V, W = ints(s, p)
    out = ([], [], [])
    for u, v, w in zip(U, V, W):
        a, b = rng.randrange(1, p), rng.randrange(1, p)
        c = pow(a * b, -1, p)
        out[0].append([[x * a % p for x in r] for r in u])
        out[1].append([[x * b % p for x in r] for r in v])
        out[2].append([[x * c % p for x in r] for r in w])
    return build(s.n, *out)


def shuffle_terms(s, p, rng):
    U, V, W = ints(s, p)
    idx = list(range(s.rank))
    rng.shuffle(idx)
    return build(s.n, [U[i] for i in idx], [V[i] for i in idx], [W[i] for i in idx])


def transpose_move(s, p, rng=None):
    U, V, W = ints(s, p)
    return build(s.n, [transpose(v) for v in V], [transpose(u) for u in U], [transpose(w) for w in W])


def cyclic_move(s, p, rng=None):
    U, V, W = ints(s, p)
    return build(s.n, V, [transpose(w) for w in W], [transpose(u) for u in U])


MOVES = [basis_change, rescale, shuffle_terms, transpose_move, cyclic_move]


# ---- the moves are real symmetries --------------------------------------------

@pytest.mark.parametrize("move", MOVES, ids=lambda m: m.__name__)
@pytest.mark.parametrize("make", [f3_a, lambda: naive(3, 3, 3)], ids=["f3_r23", "naive"])
def test_move_keeps_scheme_correct(move, make):
    rng = random.Random(1)
    s2 = move(make(), P, rng)
    assert verify(s2, mod=P).ok


def test_moves_keep_strassen_correct_mod_3():
    rng = random.Random(2)
    s = strassen()
    for move in MOVES:
        s = move(s, P, rng)
        assert verify(s, mod=P).ok, move.__name__


# ---- the invariant does not change -----------------------------------------------

@pytest.mark.timeout(30)
@pytest.mark.parametrize("seed", range(4))
def test_invariant_unchanged_by_random_symmetries(seed):
    rng = random.Random(seed)
    s = f3_a()
    want = invariant(s, P)
    for _ in range(4):
        s = rng.choice(MOVES)(s, P, rng)
        assert invariant(s, P) == want
    assert verify(s, mod=P).ok


def test_each_move_alone_keeps_invariant():
    rng = random.Random(7)
    s = f3_b()
    want = class_key(s, P)
    for move in MOVES:
        assert class_key(move(s, P, rng), P) == want, move.__name__


# ---- the invariant separates schemes that differ -------------------------------------

def test_separates_two_known_different_schemes():
    a, b = f3_a(), f3_b()
    assert term_profile(a, P) != term_profile(b, P)  # coarse invariant already differs
    assert class_key(a, P) != class_key(b, P)


def test_separates_ranks():
    assert class_key(naive(3, 3, 3), P) != class_key(f3_a(), P)


def test_rank_mod():
    assert rank_mod([[1, 2], [2, 4]], 3) == 1
    assert rank_mod([[1, 2], [2, 1]], 3) == 1  # 1*1 - 2*2 = -3 = 0 mod 3
    assert rank_mod([[1, 2], [2, 1]], 5) == 2
    assert rank_mod([[0, 0], [0, 0]], 2) == 0
