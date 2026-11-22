"""Invariants that tell inequivalent schemes apart (square formats, over GF(p)).

Two schemes for n x n x n are *equivalent* if one turns into the other by
  - a change of basis:  U -> P^T U Q^-T,  V -> Q^T V R^-T,  W -> P^-1 W R
    (P, Q, R invertible), applied to every term;
  - rescaling a term:  (U, V, W) -> (a U, b V, c W) with a b c = 1;
  - reordering the terms;
  - transposition:  (U, V, W) -> (V^T, U^T, W^T);
  - the cyclic symmetry:  (U, V, W) -> (V, W^T, U^T).
(The test suite checks that each of these really maps correct schemes to correct schemes.)

Everything computed here is unchanged by all of those moves, so different
invariants PROVE that two schemes are inequivalent. Equal invariants do not
prove equivalence: counting distinct invariants gives a LOWER BOUND on the
number of equivalence classes.

  - term profile: for each term, (rank U, rank V, rank W);
  - pair profile: for each pair of terms s, t and each factor X in {U, V, W},
    the multiset {rank(X_s + c X_t) : c != 0}. Rescaling a term only permutes c,
    and a change of basis keeps ranks;
  - each term carries its own triple plus the sorted pair entries it is in, so
    the invariant also sees how terms connect to each other.
The transposition and cyclic moves permute the three factor positions (all six
orders are reachable), so we take the smallest value over the six orders.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from itertools import combinations, permutations

from .scheme import Scheme


def _reduce(x: Fraction, p: int) -> int:
    den = x.denominator % p
    if den == 0:
        raise ValueError(f"coefficient {x} has denominator divisible by {p}")
    return (x.numerator * pow(den, -1, p)) % p


def rank_mod(mat: list[list[int]], p: int) -> int:
    """Rank of an integer matrix over GF(p) (Gaussian elimination)."""
    a = [[x % p for x in row] for row in mat]
    rows, cols = len(a), len(a[0]) if a else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if a[i][c]), None)
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [(x * inv) % p for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [(x - f * y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


def _factors(scheme: Scheme, p: int):
    red = lambda m: [[_reduce(x, p) for x in row] for row in m]  # noqa: E731
    return [[red(m) for m in X] for X in (scheme.U, scheme.V, scheme.W)]


def invariant(scheme: Scheme, p: int) -> tuple:
    """Canonical invariant of the scheme's equivalence class over GF(p)."""
    if not (scheme.n == scheme.m == scheme.p):
        raise ValueError("invariants are implemented for square formats only")
    F = _factors(scheme, p)
    r = scheme.rank
    term = [tuple(rank_mod(F[x][t], p) for x in range(3)) for t in range(r)]
    pair = {}
    for s, t in combinations(range(r), 2):
        entry = []
        for x in range(3):
            A, B = F[x][s], F[x][t]
            entry.append(tuple(sorted(
                rank_mod([[a + c * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)], p)
                for c in range(1, p))))
        pair[(s, t)] = tuple(entry)
    best = None
    for order in permutations(range(3)):
        per_term = []
        for t in range(r):
            mine = sorted(tuple(pair[tuple(sorted((t, u)))][x] for x in order) for u in range(r) if u != t)
            per_term.append((tuple(term[t][x] for x in order), tuple(mine)))
        cand = tuple(sorted(per_term))
        if best is None or cand < best:
            best = cand
    return best


def class_key(scheme: Scheme, p: int) -> str:
    """Short stable fingerprint of `invariant` (equal keys <=> equal invariants)."""
    return hashlib.sha1(repr(invariant(scheme, p)).encode()).hexdigest()[:16]


def term_profile(scheme: Scheme, p: int) -> tuple:
    """Coarse summary: sorted (rank U, rank V, rank W) triples, minimized over the six orders."""
    F = _factors(scheme, p)
    term = [tuple(rank_mod(F[x][t], p) for x in range(3)) for t in range(scheme.rank)]
    return min(tuple(sorted(tuple(tr[x] for x in o) for tr in term)) for o in permutations(range(3)))
