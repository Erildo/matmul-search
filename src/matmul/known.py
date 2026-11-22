"""Reference schemes we trust, used as test fixtures and search starting points."""

from __future__ import annotations

from .scheme import Scheme


def _zeros(r, c):
    return [[0] * c for _ in range(r)]


def _unit(r, c, i, j):
    M = _zeros(r, c)
    M[i][j] = 1
    return M


def naive(n: int, m: int, p: int) -> Scheme:
    """The textbook algorithm: rank n*m*p."""
    U, V, W = [], [], []
    for i in range(n):
        for j in range(m):
            for l in range(p):
                U.append(_unit(n, m, i, j))
                V.append(_unit(m, p, j, l))
                W.append(_unit(n, p, i, l))
    return Scheme(n, m, p, U, V, W, {"name": "naive", "source": "textbook"})


def strassen() -> Scheme:
    """Strassen (1969): 2x2x2 in 7 multiplications."""
    # A = [[a11, a12], [a21, a22]], same for B, C.
    U = [
        [[1, 0], [0, 1]],   # M1 = (a11 + a22)(b11 + b22)
        [[0, 0], [1, 1]],   # M2 = (a21 + a22) b11
        [[1, 0], [0, 0]],   # M3 = a11 (b12 - b22)
        [[0, 0], [0, 1]],   # M4 = a22 (b21 - b11)
        [[1, 1], [0, 0]],   # M5 = (a11 + a12) b22
        [[-1, 0], [1, 0]],  # M6 = (a21 - a11)(b11 + b12)
        [[0, 1], [0, -1]],  # M7 = (a12 - a22)(b21 + b22)
    ]
    V = [
        [[1, 0], [0, 1]],
        [[1, 0], [0, 0]],
        [[0, 1], [0, -1]],
        [[-1, 0], [1, 0]],
        [[0, 0], [0, 1]],
        [[1, 1], [0, 0]],
        [[0, 0], [1, 1]],
    ]
    # c11 = M1+M4-M5+M7, c12 = M3+M5, c21 = M2+M4, c22 = M1-M2+M3+M6
    W = [
        [[1, 0], [0, 1]],
        [[0, 0], [1, -1]],
        [[0, 1], [0, 1]],
        [[1, 0], [1, 0]],
        [[-1, 1], [0, 0]],
        [[0, 0], [0, 1]],
        [[1, 0], [0, 0]],
    ]
    return Scheme(2, 2, 2, U, V, W, {"name": "strassen", "source": "Strassen 1969"})


def kron(s1: Scheme, s2: Scheme) -> Scheme:
    """Tensor (Kronecker) product: <n1,m1,p1> x <n2,m2,p2> -> <n1n2, m1m2, p1p2>, rank r1*r2."""

    def kmat(X, Y):
        r1, c1, r2, c2 = len(X), len(X[0]), len(Y), len(Y[0])
        return [
            [X[a // r2][b // c2] * Y[a % r2][b % c2] for b in range(c1 * c2)]
            for a in range(r1 * r2)
        ]

    U, V, W = [], [], []
    for k1 in range(s1.rank):
        for k2 in range(s2.rank):
            U.append(kmat(s1.U[k1], s2.U[k2]))
            V.append(kmat(s1.V[k1], s2.V[k2]))
            W.append(kmat(s1.W[k1], s2.W[k2]))
    return Scheme(
        s1.n * s2.n, s1.m * s2.m, s1.p * s2.p, U, V, W,
        {"name": f"kron({s1.meta.get('name')},{s2.meta.get('name')})"},
    )
