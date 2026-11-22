"""Exact verification of matrix multiplication schemes.

Two independent checks:

1. Tensor check (the real proof): the scheme is correct iff for every
   (i,j), (j',l), (i',l'):
       sum_k U_k[i][j] * V_k[j'][l] * W_k[i'][l'] = [j==j'][i==i'][l==l']
   computed exactly over Q, or modulo a prime p.

2. Random-matrix check (an independent sanity test): run the scheme on random
   integer matrices and compare with A @ B, exactly.

Both are deterministic and bounded: work is proportional to the number of
nonzero coefficients. Nothing here loops open-endedly.

CLI:  python -m matmul.verify schemes/foo.json [--mod 2]
"""

from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass, field
from fractions import Fraction

from .scheme import Scheme


@dataclass
class VerifyResult:
    ok: bool
    format: str
    rank: int
    field: str
    wrong_entries: int = 0
    examples: list = field(default_factory=list)  # [(a_index, b_index, c_index, got, want)]

    def summary(self) -> str:
        status = "CORRECT" if self.ok else f"WRONG ({self.wrong_entries} tensor entries differ)"
        return f"{self.format} rank {self.rank} over {self.field}: {status}"


def _reduce(x: Fraction, mod: int | None):
    if mod is None:
        return x
    den = x.denominator % mod
    if den == 0:
        raise ValueError(f"coefficient {x} has denominator divisible by {mod}; undefined mod {mod}")
    return (x.numerator * pow(den, -1, mod)) % mod


def _nonzeros(mat, mod):
    out = []
    for r, row in enumerate(mat):
        for c, x in enumerate(row):
            y = _reduce(x, mod)
            if y != 0:
                out.append(((r, c), y))
    return out


def verify(scheme: Scheme, mod: int | None = None, max_examples: int = 5) -> VerifyResult:
    """Exact tensor check. mod=None means over Q; mod=p means over GF(p)."""
    if mod is not None and mod < 2:
        raise ValueError("mod must be a prime >= 2")
    n, m, p = scheme.n, scheme.m, scheme.p
    acc: dict = {}
    for U, V, W in zip(scheme.U, scheme.V, scheme.W):
        u, v, w = _nonzeros(U, mod), _nonzeros(V, mod), _nonzeros(W, mod)
        for a, ua in u:
            for b, vb in v:
                uv = ua * vb
                for c, wc in w:
                    key = (a, b, c)
                    acc[key] = acc.get(key, 0) + uv * wc

    def target(a, b, c):
        (i, j), (j2, l), (i2, l2) = a, b, c
        return 1 if (j == j2 and i == i2 and l == l2) else 0

    wrong = []
    # entries that should be 1
    for i in range(n):
        for j in range(m):
            for l in range(p):
                key = ((i, j), (j, l), (i, l))
                got = acc.pop(key, 0)
                if mod is not None:
                    got %= mod
                if got != 1:
                    wrong.append((key[0], key[1], key[2], got, 1))
    # everything left over should be 0
    for key, got in acc.items():
        if mod is not None:
            got %= mod
        if got != 0:
            wrong.append((key[0], key[1], key[2], got, 0))

    return VerifyResult(
        ok=not wrong,
        format=scheme.format,
        rank=scheme.rank,
        field="Q" if mod is None else f"GF({mod})",
        wrong_entries=len(wrong),
        examples=wrong[:max_examples],
    )


def apply_scheme(scheme: Scheme, A, B):
    """Multiply A (n x m) by B (m x p) using the scheme, exactly."""
    n, m, p = scheme.n, scheme.m, scheme.p
    C = [[Fraction(0)] * p for _ in range(n)]
    for U, V, W in zip(scheme.U, scheme.V, scheme.W):
        left = sum(U[i][j] * A[i][j] for i in range(n) for j in range(m))
        right = sum(V[j][l] * B[j][l] for j in range(m) for l in range(p))
        prod = left * right
        if prod:
            for i in range(n):
                for l in range(p):
                    if W[i][l]:
                        C[i][l] += W[i][l] * prod
    return C


def random_check(scheme: Scheme, trials: int = 3, seed: int = 0, lo: int = -9, hi: int = 9) -> bool:
    """Independent check on random integer matrices (seeded, fixed number of trials)."""
    rng = random.Random(seed)
    n, m, p = scheme.n, scheme.m, scheme.p
    for _ in range(trials):
        A = [[rng.randint(lo, hi) for _ in range(m)] for _ in range(n)]
        B = [[rng.randint(lo, hi) for _ in range(p)] for _ in range(m)]
        want = [[sum(A[i][j] * B[j][l] for j in range(m)) for l in range(p)] for i in range(n)]
        if apply_scheme(scheme, A, B) != want:
            return False
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Exactly verify a matrix multiplication scheme")
    ap.add_argument("files", nargs="+")
    ap.add_argument("--mod", type=int, default=None, help="verify modulo this prime instead of over Q")
    args = ap.parse_args(argv)
    all_ok = True
    for f in args.files:
        s = Scheme.load(f)
        r = verify(s, mod=args.mod)
        print(f"{f}: {r.summary()}")
        for ex in r.examples:
            print(f"   entry A{ex[0]} B{ex[1]} C{ex[2]}: got {ex[3]}, want {ex[4]}")
        all_ok &= r.ok
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
