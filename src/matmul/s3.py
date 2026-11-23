"""S3-invariant schemes over F3: the cyclic symmetry sigma plus the transpose symmetry tau (attempt 010).

Convention as in matmul.cyclic: a term (U, V, W) contributes U[i][j] V[j][l] W[i][l] = tr(U V W^T).

    sigma: (U, V, W) -> (V, W^T, U^T)        tau: (U, V, W) -> (V^T, U^T, W^T)

tau keeps tr(U V W^T), tau^2 = id and tau sigma tau = sigma^-1, so <sigma, tau> = S3. Its orbits on the terms
of an invariant scheme have size 6, 3, 2 or 1; with a, b, c, d orbits of each size the rank is 6a + 3b + 2c + d.

Slot layout of S3F3CNF(n, a, b, c, d), on top of CyclicF3CNF(n, 2a + b, 2c + d):
  orbit slots  0..a-1        size 6, base t = (u, v, w)
               a..2a-1       size 6, partner: slot a + k = tau(slot k) = (v^T, u^T, w^T)
               2a..2a+b-1    size 3: (u, u^T, w) with w = w^T (the member fixed by tau, rescaled)
  fixed slots  0..c-1        size 2, base x
               c..2c-1       size 2, partner x^T
               2c..2c+d-1    size 1: x = x^T
Every slot may be zero, so an instance covers every smaller (a', b', c', d'). A size-6 slot whose t is itself
tau-fixed repeats a term; that only adds solutions (UNSAT stays sound), and SAT answers are verified.

Symmetry breaking (every S3-invariant scheme has a representative, see s3_normal_form):
  size 6: first nonzero of u and of v is 1 (rescaling (eu, dv, edw) acts consistently on the partner);
          bases in lex order;
  size 3: first nonzero of u is 1 ((eu, eu^T, w) keeps the tensor); slots in lex order;
  size 2: x <= x^T in lex order; bases in lex order;
  size 1: lex order.
"""

from __future__ import annotations

from .cyclic import P, CyclicF3CNF, T, _bitkey, expand, sigma, split_orbits
from .scheme import Scheme


def tau(t):
    u, v, w = t
    return (T(v), T(u), T(w))


def _flat(m):
    return tuple(e % P for r in m for e in r)


def tensor_key(t):
    """Canonical key of the rank-one tensor u (x) v (x) w over F3 (the representations are (au, bv, abw))."""
    return min((tuple(x * a % P for x in _flat(t[0])), tuple(x * b % P for x in _flat(t[1])),
                tuple(x * a * b % P for x in _flat(t[2]))) for a in (1, 2) for b in (1, 2))


def _terms(scheme: Scheme):
    return [tuple([[int(e) % P for e in r] for r in m] for m in (U, V, W)) for U, V, W in zip(scheme.U, scheme.V, scheme.W)]


def is_s3_invariant(scheme: Scheme) -> bool:
    keys = [tensor_key(t) for t in _terms(scheme)]
    S = set(keys)
    if len(S) != len(keys):
        return False
    return all(tensor_key(sigma(t)) in S and tensor_key(tau(t)) in S for t in _terms(scheme))


def _first(m):
    return next((e % P for r in m for e in r if e % P), 1)


def _sc(m, s):
    return [[e * s % P for e in r] for r in m]


def s3_split(orbits, fixed):
    """Classify a sigma-invariant scheme that is also tau-invariant into S3 slots.

    Returns (size6_bases, size3, size2_bases, size1) with size-3 terms written as (u, u^T, w), w = w^T.
    Raises ValueError if the scheme is not tau-invariant.
    """
    def nz(t):
        return all(any(e % P for r in m for e in r) for m in t)

    orbits = [t for t in orbits if nz(t)]
    fixed = [x for x in fixed if any(e % P for r in x for e in r)]
    okeys = [{tensor_key(t), tensor_key(sigma(t)), tensor_key(sigma(sigma(t)))} for t in orbits]
    used, six, three = set(), [], []
    for i, t in enumerate(orbits):
        if i in used:
            continue
        used.add(i)
        tk = {tensor_key(tau(m)) for m in (t, sigma(t), sigma(sigma(t)))}
        if tk == okeys[i]:
            m = next(m for m in (t, sigma(t), sigma(sigma(t))) if tensor_key(tau(m)) == tensor_key(m))
            u, v, w = m
            # v = alpha u^T as matrices up to the tensor scalings; write the term as (u, u^T, w')
            alpha = next(v[j][i2] * u[i2][j] % P for i2 in range(len(u)) for j in range(len(u)) if u[i2][j] % P)
            three.append((u, T(u), _sc(w, alpha)))
        else:
            j = next((j for j in range(len(orbits)) if j not in used and okeys[j] == tk), None)
            if j is None:
                raise ValueError("not tau-invariant")
            used.add(j)
            six.append(t)
    fk = {_flat(x) for x in fixed}
    two, one, seen = [], [], set()
    for x in fixed:
        if _flat(x) in seen:
            continue
        seen.add(_flat(x))
        if _flat(T(x)) == _flat(x):
            one.append(x)
        else:
            if _flat(T(x)) not in fk:
                raise ValueError("not tau-invariant (fixed part)")
            seen.add(_flat(T(x)))
            two.append(x)
    return six, three, two, one


def s3_normal_form(six, three, two, one):
    """Slot contents (orbits, fixed) for S3F3CNF that satisfy its symmetry breaking."""
    bases = []
    for u, v, w in six:
        s, r = _first(u), _first(v)
        bases.append((_sc(u, s), _sc(v, r), _sc(w, s * r)))
    bases.sort(key=_bitkey)
    thr = []
    for u, _, w in three:
        s = _first(u)
        thr.append((_sc(u, s), _sc(T(u), s), [[e % P for e in r] for r in w]))
    thr.sort(key=_bitkey)
    tb = sorted((min(x, T(x), key=lambda m: _bitkey([m])) for x in two), key=lambda m: _bitkey([m]))
    tb = [[[e % P for e in r] for r in x] for x in tb]
    on = sorted(([[e % P for e in r] for r in x] for x in one), key=lambda m: _bitkey([m]))
    orbits = bases + [tau(t) for t in bases] + thr
    fixed = tb + [T(x) for x in tb] + on
    return orbits, fixed


class S3F3CNF(CyclicF3CNF):
    """Brent equations mod 3 for an S3-invariant n x n x n scheme with S3 split (a, b, c, d); see the module doc."""

    def __init__(self, n: int, a: int, b: int, c: int, d: int, fixed_sum: bool = True, break_sym: bool = True):
        super().__init__(n, 2 * a + b, 2 * c + d, fixed_sum=fixed_sum, break_sym=False)
        self.s3 = (a, b, c, d)
        tr = lambda vals: [vals[j * n + i] for i in range(n) for j in range(n)]  # value list of the transpose
        c0 = len(self.clauses)
        for k in range(a):
            u, v, w = self.orb[k]
            u2, v2, w2 = self.orb[a + k]
            self._tie(u2, tr(v))
            self._tie(v2, tr(u))
            self._tie(w2, tr(w))
        for k in range(2 * a, 2 * a + b):
            u, v, w = self.orb[k]
            self._tie(v, tr(u))
            self._tie(w, tr(w))
        for k in range(c):
            self._tie(self.fix[c + k], tr(self.fix[k]))
        for k in range(2 * c, 2 * c + d):
            self._tie(self.fix[k], tr(self.fix[k]))
        self.tie_clauses = (c0, len(self.clauses))
        if break_sym:
            for k in range(a):
                self._first_nonzero_is_one(self.orb[k][0])
                self._first_nonzero_is_one(self.orb[k][1])
            for k in range(a - 1):
                self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
            for k in range(2 * a, 2 * a + b):
                self._first_nonzero_is_one(self.orb[k][0])
            for k in range(2 * a, 2 * a + b - 1):
                self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
            for k in range(c):
                self._lex_leq(self._bits(self.fix[k]), self._bits(tr(self.fix[k])))
            for k in range(c - 1):
                self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))
            for k in range(2 * c, 2 * c + d - 1):
                self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))

    def _tie(self, xs, ys):
        for (p, q), (p2, q2) in zip(xs, ys):
            if (p, q) == (p2, q2):
                continue
            self.clauses += [[-p, p2], [p, -p2], [-q, q2], [q, -q2]]


def s3_splits(r: int, min_sigma_orbits: int = 0):
    """All (a, b, c, d) with 6a + 3b + 2c + d = r and 2a + b >= max(1, min_sigma_orbits)."""
    out = []
    for a in range(r // 6 + 1):
        for b in range((r - 6 * a) // 3 + 1):
            for c in range((r - 6 * a - 3 * b) // 2 + 1):
                d = r - 6 * a - 3 * b - 2 * c
                if 2 * a + b >= max(1, min_sigma_orbits):
                    out.append((a, b, c, d))
    return out


def scheme_s3_split(scheme: Scheme):
    """S3 split (a, b, c, d) of an S3-invariant scheme."""
    six, three, two, one = s3_split(*split_orbits(scheme))
    return len(six), len(three), len(two), len(one)

