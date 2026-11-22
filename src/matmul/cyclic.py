"""Cyclic-invariant (Z3-symmetric) schemes for square matrix multiplication, and a CNF encoding over F3.

The cyclic symmetry of the matrix product tensor, in this repo's convention
(term (U, V, W) contributes U[i][j] * V[j][l] * W[i][l]), is

    sigma: (U, V, W) -> (V, W^T, U^T),       sigma^3 = id.

A sigma-invariant scheme is a union of
  - orbits of size 3:  (u, v, w), (v, w^T, u^T), (w^T, u, v^T), one representative each;
  - fixed terms:       a rank-one tensor equal to its own image. Up to rescaling this is
                       lambda * (x, x, x^T); over F3 every scalar is a cube (mu^3 = mu), so a
                       fixed term is exactly (x, x, x^T) for some nonzero matrix x.

Necessary condition over F3 (proof in attempts/005): the fixed matrices sum to the identity.
Setting all three arguments of the trilinear form to the same matrix kills every orbit
(it contributes 3 * (...) = 0 mod 3), and each fixed term contributes l_x(X)^3, which in
characteristic 3 equals the sum over the entries of x_e * X_e^3. The left side reduces to the
sum of X_ii^3.

The encoder writes the Brent equations of a sigma-invariant scheme with `a` orbit slots and
`b` fixed slots (any slot may be zero, so an (a, b) instance covers every split (a', b') with
a' <= a, b' <= b) as CNF. Each F3 value is two Boolean variables (p, n): value 1 if p, 2 if n,
0 if neither, never both.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product

from .scheme import Scheme

P = 3


# ---- the symmetry ----------------------------------------------------------------

def T(m):
    return [list(r) for r in zip(*m)]


def sigma(term):
    u, v, w = term
    return (v, T(w), T(u))


def expand(n: int, orbits, fixed, meta=None) -> Scheme:
    """Full scheme from orbit representatives (u, v, w) and fixed matrices x. Zero terms are dropped."""
    terms = []
    for t in orbits:
        t1 = sigma(t)
        t2 = sigma(t1)
        terms += [t, t1, t2]
    terms += [(x, x, T(x)) for x in fixed]
    terms = [t for t in terms if all(any(e % P for r in m for e in r) for m in t)]
    return Scheme(n, n, n, [t[0] for t in terms], [t[1] for t in terms], [t[2] for t in terms], meta or {})


def _tensor(U, V, W, p):
    fu = [x % p for r in U for x in r]
    fv = [x % p for r in V for x in r]
    fw = [x % p for r in W for x in r]
    return tuple(a * b * c % p for a in fu for b in fv for c in fw)


def is_cyclic_invariant(scheme: Scheme, p: int = P) -> bool:
    """True iff the multiset of rank-one tensors (mod p) is closed under sigma."""
    from collections import Counter

    def red(m):
        return [[int(x.numerator * pow(x.denominator, -1, p)) % p for x in r] for r in m]

    terms = [(red(u), red(v), red(w)) for u, v, w in zip(scheme.U, scheme.V, scheme.W)]
    a = Counter(_tensor(*t, p) for t in terms)
    b = Counter(_tensor(*sigma(t), p) for t in terms)
    return a == b


def split_orbits(scheme: Scheme, p: int = P):
    """Inverse of `expand` for a cyclic-invariant scheme mod p: (orbit representatives, fixed matrices)."""
    def red(m):
        return [[int(x.numerator * pow(x.denominator, -1, p)) % p if hasattr(x, "numerator") else x % p
                 for x in r] for r in m]

    terms = [(red(u), red(v), red(w)) for u, v, w in zip(scheme.U, scheme.V, scheme.W)]
    left = list(range(len(terms)))
    orbits, fixed = [], []
    while left:
        t = terms[left.pop(0)]
        k0 = _tensor(*t, p)
        if _tensor(*sigma(t), p) == k0:
            x = next(c for c in (t[0], [[-e % p for e in r] for r in t[0]])
                     if _tensor(c, c, T(c), p) == k0)
            fixed.append(x)
            continue
        for img in (sigma(t), sigma(sigma(t))):
            k = _tensor(*img, p)
            j = next(j for j in left if _tensor(*terms[j], p) == k)
            left.remove(j)
        orbits.append(t)
    return orbits, fixed


def _bitkey(mats):
    """The Boolean vector the encoder orders lexicographically: per value (p, n), False < True."""
    return tuple(b for m in mats for r in m for e in r for b in (e % 3 == 1, e % 3 == 2))


def normal_form(orbits, fixed):
    """Rescale and reorder slot contents so they satisfy the encoder's symmetry-breaking constraints.

    Each orbit (u, v, w) is rescaled to (u/s, v/t, w*s*t) with s, t the first nonzero entries of u, v;
    then orbit slots and fixed slots are sorted by `_bitkey`. The scheme is unchanged.
    """
    def first(m):
        return next((e % 3 for r in m for e in r if e % 3), 1)

    out = []
    for u, v, w in orbits:
        s, t = first(u), first(v)  # in F3, 1/s = s
        out.append(([[e * s % 3 for e in r] for r in u], [[e * t % 3 for e in r] for r in v],
                    [[e * s * t % 3 for e in r] for r in w]))
    out.sort(key=_bitkey)
    fx = sorted(([[e % 3 for e in r] for r in x] for x in fixed), key=lambda x: _bitkey([x]))
    return out, fx


# ---- conjugation symmetry and the classes of matrices ---------------------------------

def _mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % P for j in range(len(B[0]))] for i in range(len(A))]


def _inv(A):
    n = len(A)
    M = [[x % P for x in row] + [int(i == j) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c])
        M[c], M[piv] = M[piv], M[c]
        iv = M[c][c]  # in F3, 1/x = x
        M[c] = [x * iv % P for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % P for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]


def gl_generators(n: int = 3):
    """Generators of GL_n(F3): the transvections I + E_ij (i != j) and diag(2, 1, ..., 1)."""
    gens = []
    for i in range(n):
        for j in range(n):
            if i != j:
                g = [[int(r == c) for c in range(n)] for r in range(n)]
                g[i][j] = 1
                gens.append(g)
    d = [[int(r == c) for c in range(n)] for r in range(n)]
    d[0][0] = 2
    gens.append(d)
    return gens


def conj(M, term):
    """g_M: (u, v, w) -> (M u M^-1, M v M^-1, M^-T w M^T). Commutes with sigma and keeps the matrix product tensor."""
    Mi = _inv(M)
    u, v, w = term
    return (_mm(_mm(M, u), Mi), _mm(_mm(M, v), Mi), _mm(_mm(T(Mi), w), T(M)))


def _key(x):
    return tuple(e % P for r in x for e in r)


def _unkey(k, n):
    return [list(k[i * n:(i + 1) * n]) for i in range(n)]


def matrix_classes(n: int = 3, signed: bool = True):
    """Partition of all n x n matrices over F3 into conjugacy classes (signed=True: up to sign as well).

    Returns {rep_key: [member keys]}; the representative is the `_bitkey`-smallest member.
    """
    gens = [(g, _inv(g)) for g in gl_generators(n)]
    seen, out = set(), {}
    for k in product(range(P), repeat=n * n):
        if k in seen:
            continue
        cls, stack = [k], [k]
        seen.add(k)
        while stack:
            x = _unkey(stack.pop(), n)
            nb = [_key(_mm(_mm(g, x), gi)) for g, gi in gens]
            if signed:
                nb.append(_key([[-e for e in r] for r in x]))
            for y in nb:
                if y not in seen:
                    seen.add(y)
                    cls.append(y)
                    stack.append(y)
        rep = min(cls, key=lambda c: _bitkey([_unkey(c, n)]))
        out[rep] = cls
    return out


@lru_cache(maxsize=None)
def _class_rep_keys(n: int):
    reps = [r for r in matrix_classes(n) if any(r)]
    return tuple(sorted(reps, key=lambda c: _bitkey([_unkey(c, n)])))


def class_reps(n: int = 3):
    """Representatives of the nonzero matrices up to conjugation and sign, sorted by `_bitkey` (cached)."""
    return [_unkey(r, n) for r in _class_rep_keys(n)]


def to_rep(x):
    """(rep, s, M) with M (s x) M^-1 = rep, s in {1, 2}, rep the class representative of x != 0."""
    n = len(x)
    gens = gl_generators(n)
    start = (_key(x), 1)
    best = {start: [[int(i == j) for j in range(n)] for i in range(n)]}
    stack = [start]
    while stack:  # walk the class of x, recording a conjugating matrix for each member
        k, s = stack.pop()
        y, M = _unkey(k, n), best[(k, s)]
        nxt = [((_key(_mm(_mm(g, y), _inv(g))), s), _mm(g, M)) for g in gens]
        nxt.append(((_key([[-e for e in r] for r in y]), 3 - s), M))
        for st, M2 in nxt:
            if st not in best:
                best[st] = M2
                stack.append(st)
    (k, s), M = min(best.items(), key=lambda it: (_bitkey([_unkey(it[0][0], n)]), it[0][1]))
    return _unkey(k, n), s, M


def class_index(x, reps=None):
    """Index of the class of x != 0 in class_reps() (conjugation and sign)."""
    reps = reps or [_key(r) for r in class_reps(len(x))]
    return reps.index(_key(to_rep(x)[0]))


def charpoly(x):
    """(trace, sum of principal 2x2 minors, det) of a 3x3 matrix mod 3: the characteristic polynomial."""
    t = (x[0][0] + x[1][1] + x[2][2]) % P
    m = lambda a, b: x[a][a] * x[b][b] - x[a][b] * x[b][a]
    s = (m(0, 1) + m(0, 2) + m(1, 2)) % P
    d = (x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1]) - x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0])
         + x[0][2] * (x[1][0] * x[2][1] - x[1][1] * x[2][0])) % P
    return (t, s, d)


@lru_cache(maxsize=None)
def charpoly_exclusions(c: int):
    """Characteristic polynomials carried only by classes of index < c (the zero matrix counts as allowed).

    In cube c every orbit matrix has class >= c, so none of these can occur.
    """
    reps = list(_class_rep_keys(3))
    top = {}
    for rep, members in matrix_classes(3).items():
        i = reps.index(rep) if any(rep) else len(reps)
        for k in members:
            T = charpoly(_unkey(k, 3))
            top[T] = max(top.get(T, -1), i)
    return tuple(sorted(T for T, i in top.items() if i < c))


def rotate(t, r):
    """The r-th member of the sigma-orbit of t, as the orbit representative."""
    for _ in range(r % 3):
        t = sigma(t)
    return t


def to_class_form(orbits, fixed, slot=None, minimal=False):
    """Move one nonzero orbit to slot 0 with u equal to its class representative (the class-cube lemma).

    1. choose a nonzero orbit (the given `slot`, else the first nonzero one);
    2. rescale it by (u, v, w) -> (-u, -v, w) if needed and apply g_M to the whole scheme,
       so that u becomes the class representative;
    3. sign-normalise v of slot 0 by (u, v, w) -> (u, -v, -w), and put the other orbit slots and
       the fixed slots in `normal_form`.
    Returns (orbits, fixed, rep); the scheme is a symmetric image of the input.
    """
    def nz(t):
        return all(any(e % P for r in m for e in r) for m in t)

    if minimal:
        # the member (orbit, rotation) whose first matrix has the smallest class index;
        # the first matrices of the rotations of (u, v, w) are u, v and w^T
        orbits = list(orbits)
        reps = [_key(r) for r in class_reps(len(fixed[0]) if fixed else len(orbits[0][0]))]
        cands = [(class_index(rotate(t, r)[0], reps), i, r) for i, t in enumerate(orbits) if nz(t) for r in range(3)]
        _, slot, r = min(cands)
        orbits[slot] = rotate(orbits[slot], r)
    elif slot is None:
        slot = next(i for i, t in enumerate(orbits) if nz(t))
    u, v, w = orbits[slot]
    rep, s, M = to_rep(u)
    neg = (lambda m: [[-e % P for e in r] for r in m])
    first = [(neg(u), neg(v), w) if s == 2 else (u, v, w)]
    rest = [t for i, t in enumerate(orbits) if i != slot]
    imgs = [conj(M, t) for t in first + rest]
    fx = [conj(M, (x, x, T(x)))[0] for x in fixed]
    u0, v0, w0 = imgs[0]
    assert _key(u0) == _key(rep)
    t = next((e % P for r in v0 for e in r if e % P), 1)
    if t == 2:
        v0, w0 = neg(v0), neg(w0)
    o2, f2 = normal_form(imgs[1:], fx)
    return [(u0, v0, w0)] + o2, f2, rep


# ---- CNF encoding over F3 ------------------------------------------------------------

class CyclicF3CNF:
    """Brent equations mod 3 for a sigma-invariant n x n x n scheme with a orbit slots and b fixed slots.

    Options:
      fixed_sum: add the necessary condition "fixed matrices sum to I" (implied, but helps).
      break_sym: sign-normalise each orbit slot (first nonzero of u and of v is 1), and order the
                 orbit slots and the fixed slots lexicographically. Both keep at least one
                 representative of every solution (rescaling and reordering are symmetries).
    """

    def __init__(self, n: int, a: int, b: int, fixed_sum: bool = True, break_sym: bool = True,
                 pin_u1=None, min_class: int = 0, exclude_charpoly: bool = False):
        self.n, self.a, self.b = n, a, b
        self.pin_u1 = pin_u1
        self.min_class = min_class
        self.nv = 0
        self.clauses: list[list[int]] = []
        self._pair: dict = {}
        self._triple: dict = {}
        self.aux_defs: list = []  # (kind, out, ins) in creation order, for computing a witness
        nn = n * n
        # primary F3 values: orbit k has u, v, w (each n*n), fixed k has x (n*n)
        self.orb = [[[self._val() for _ in range(nn)] for _ in range(3)] for _ in range(a)]
        self.fix = [[self._val() for _ in range(nn)] for _ in range(b)]
        self.n_primary = self.nv
        self.equations = self._equations()
        for prods, target in self.equations:
            self._sum_equals([self._prod3(t) for t in prods], target)
        if fixed_sum and b > 0:
            for e in range(nn):
                i, j = divmod(e, n)
                self._sum_equals([self.fix[k][e] for k in range(b)], 1 if i == j else 0)
        k0 = 0
        if pin_u1 is not None:
            # cube of a class (see class_form): slot 0's u is fixed to the class representative,
            # slot 0 keeps the sign normalisation of v and leaves the orbit lex chain.
            if a < 1:
                raise ValueError("pin_u1 needs a >= 1")
            for vals, x in zip(self.orb[0][0], (e % 3 for r in pin_u1 for e in r)):
                for lit in self._is(vals, x):
                    self.clauses.append([lit])
            k0 = 1
        self.class_clauses = (len(self.clauses), len(self.clauses))
        if min_class > 0:
            # minimal-class cube (to_class_form(minimal=True)): no orbit matrix u, v, w lies in a class
            # below min_class. Classes 0 and 1 are the rank-1 matrices with trace != 0 and == 0; only
            # those two are encoded, so for min_class >= 2 the constraint is "no rank-1 orbit matrix".
            if n != 3:
                raise ValueError("min_class is for n = 3")
            c0 = len(self.clauses)
            for k in range(a):
                for m in range(3):
                    if k == 0 and m == 0 and pin_u1 is not None:
                        continue
                    self._not_rank1(self.orb[k][m], only_trace_nonzero=(min_class == 1))
                    if exclude_charpoly:
                        self._charpoly_not_in(self.orb[k][m], charpoly_exclusions(min_class))
            self.class_clauses = (c0, len(self.clauses))
        if break_sym:
            for k in range(a):
                if k > 0 or pin_u1 is None:
                    self._first_nonzero_is_one(self.orb[k][0])
                self._first_nonzero_is_one(self.orb[k][1])
            for k in range(k0, a - 1):
                self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
            for k in range(b - 1):
                self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))

    # -- variables
    def _new(self) -> int:
        self.nv += 1
        return self.nv

    def _val(self):
        p, q = self._new(), self._new()
        self.clauses.append([-p, -q])
        return (p, q)

    @staticmethod
    def _bits(vals):
        return [x for pq in vals for x in pq]

    # -- equations
    def _idx(self, i, j):
        return i * self.n + j

    def _terms(self):
        """Each term as three accessors (factor_values, transposed?) over primary values."""
        out = []
        for u, v, w in self.orb:
            out.append(((u, False), (v, False), (w, False)))  # (u, v, w)
            out.append(((v, False), (w, True), (u, True)))    # (v, w^T, u^T)
            out.append(((w, True), (u, False), (v, True)))    # (w^T, u, v^T)
        for x in self.fix:
            out.append(((x, False), (x, False), (x, True)))
        return out

    def _equations(self):
        n = self.n
        terms = self._terms()

        def at(acc, i, j):
            vals, tr = acc
            return vals[self._idx(j, i) if tr else self._idx(i, j)]

        seen = {}
        for i, j, j2, l, i2, l2 in product(range(n), repeat=6):
            prods = []
            for A, B, C in terms:
                prods.append(tuple(sorted((at(A, i, j), at(B, j2, l), at(C, i2, l2)))))
            target = 1 if (j == j2 and i == i2 and l == l2) else 0
            key = (tuple(sorted(prods)), target)
            seen.setdefault(key, None)
        # The same polynomial with two different targets: no solution at all. This happens for a = 0
        # (the tensor is not a sum of fixed terms); for a >= 1 it would mean the action is wrong (tested).
        byprods = {}
        self.trivially_unsat = False
        for prods, target in seen:
            if byprods.setdefault(prods, target) != target:
                self.trivially_unsat = True
        if self.trivially_unsat:
            self.clauses.append([])
        return [(list(p_), t) for p_, t in seen]

    # -- gates
    def _prod2(self, x, y):
        key = tuple(sorted((x, y)))
        if key in self._pair:
            return self._pair[key]
        (xp, xn), (yp, yn) = x, y
        z = self._val()
        zp, zn = z
        self.clauses += [
            [-xp, -yp, zp], [-xn, -yn, zp], [-xp, -yn, zn], [-xn, -yp, zn],
            [-zp, xp, xn], [-zp, yp, yn], [-zp, xp, yn], [-zp, xn, yp],
            [-zn, xp, xn], [-zn, yp, yn], [-zn, xp, yp], [-zn, xn, yn],
        ]
        self.aux_defs.append(("mul", z, (x, y)))
        self._pair[key] = z
        return z

    def _prod3(self, t):
        if t in self._triple:
            return self._triple[t]
        z = self._prod2(self._prod2(t[0], t[1]), t[2])
        self._triple[t] = z
        return z

    def _is(self, x, c):
        """Literals whose conjunction means x == c (negate them for a clause)."""
        p, q = x
        return [-p, -q] if c == 0 else ([p] if c == 1 else [q])

    def _add(self, x, y):
        z = self._val()
        for c1 in range(3):
            for c2 in range(3):
                pre = [-l for l in self._is(x, c1) + self._is(y, c2)]
                for lit in self._is(z, (c1 + c2) % 3):
                    self.clauses.append(pre + [lit])
        self.aux_defs.append(("add", z, (x, y)))
        return z

    def _sum_equals(self, vals, target):
        if not vals:
            if target % 3:
                self.clauses.append([])
            return
        while len(vals) > 1:  # balanced tree
            nxt = [self._add(vals[i], vals[i + 1]) for i in range(0, len(vals) - 1, 2)]
            if len(vals) % 2:
                nxt.append(vals[-1])
            vals = nxt
        for lit in self._is(vals[0], target % 3):
            self.clauses.append([lit])

    def _not_rank1(self, vals, only_trace_nonzero=False):
        """x is not of rank one (or, with only_trace_nonzero, not of rank one with nonzero trace)."""
        n = self.n
        x = lambda i, j: vals[i * n + j]
        neg = lambda z: (z[1], z[0])
        big = []
        for i1 in range(n):
            for i2 in range(i1 + 1, n):
                for j1 in range(n):
                    for j2 in range(j1 + 1, n):
                        m = self._add(self._prod2(x(i1, j1), x(i2, j2)), neg(self._prod2(x(i1, j2), x(i2, j1))))
                        big += list(m)  # minor != 0
        y = self._new()
        if only_trace_nonzero:
            t = self._add(self._add(x(0, 0), x(1, 1)), x(2, 2))
            self.clauses += [[-y, -t[0]], [-y, -t[1]]]  # y -> trace == 0 (this also covers x == 0)
            self.aux_defs.append(("is0", y, (t,)))
        else:
            self.clauses += [[-y, -l] for pq in vals for l in pq]  # y -> x == 0
            self.aux_defs.append(("is0", y, tuple(vals)))
        self.clauses.append(big + [y])

    def _charpoly_not_in(self, vals, forbidden):
        """The characteristic polynomial (trace, c2, det) of x is none of the forbidden triples."""
        if not forbidden:
            return
        x = lambda i, j: vals[i * 3 + j]
        neg = lambda z: (z[1], z[0])
        tr = self._add(self._add(x(0, 0), x(1, 1)), x(2, 2))
        minor = lambda a, b, c, d: self._add(self._prod2(x(a, c), x(b, d)), neg(self._prod2(x(a, d), x(b, c))))
        c2 = self._add(self._add(minor(0, 1, 0, 1), minor(0, 2, 0, 2)), minor(1, 2, 1, 2))
        det = self._add(self._add(self._prod2(x(0, 0), minor(1, 2, 1, 2)), neg(self._prod2(x(0, 1), minor(1, 2, 0, 2)))),
                        self._prod2(x(0, 2), minor(1, 2, 0, 1)))
        ne = lambda z, v: [z[0], z[1]] if v == 0 else ([-z[0]] if v == 1 else [-z[1]])
        for t, s2, d in forbidden:
            self.clauses.append(ne(tr, t) + ne(c2, s2) + ne(det, d))

    def _first_nonzero_is_one(self, vals):
        prefix = []
        for p, q in vals:
            self.clauses.append(prefix + [-q])
            prefix = prefix + [p, q]

    def _lex_leq(self, A, B):
        e = None  # "equal so far"; None means true
        for k, (x, y) in enumerate(zip(A, B)):
            pre = [] if e is None else [-e]
            self.clauses.append(pre + [-x, y])
            if k == len(A) - 1:
                break
            e2 = self._new()
            self.clauses.append(pre + [-x, -y, e2])
            self.clauses.append(pre + [x, y, e2])
            self.aux_defs.append(("eq", e2, (e, x, y)))
            e = e2

    # -- output and decoding
    def dimacs(self) -> str:
        lines = [f"p cnf {self.nv} {len(self.clauses)}"]
        lines += [" ".join(map(str, c)) + " 0" for c in self.clauses]
        return "\n".join(lines) + "\n"

    def _mat(self, vals, assign):
        n = self.n
        f = [(1 if assign.get(p) else 0) - (1 if assign.get(q) else 0) for p, q in vals]
        return [[f[i * n + j] % 3 for j in range(n)] for i in range(n)]

    def decode(self, model) -> tuple[list, list]:
        assign = {abs(l): l > 0 for l in model}
        orbits = [tuple(self._mat(m, assign) for m in o) for o in self.orb]
        fixed = [self._mat(x, assign) for x in self.fix]
        return orbits, fixed

    def witness(self, orbits, fixed) -> dict:
        """Full assignment (primary + auxiliary) for given slot contents, to test the encoding."""
        assign = {}

        def put(vals, mat):
            for (p, q), x in zip(vals, (e % 3 for r in mat for e in r)):
                assign[p], assign[q] = x == 1, x == 2

        for slot, o in zip(self.orb, orbits):
            for vals, mat in zip(slot, o):
                put(vals, mat)
        for vals, x in zip(self.fix, fixed):
            put(vals, x)

        def val(z):
            return (1 if assign[z[0]] else 0) + (2 if assign[z[1]] else 0)

        for kind, z, ins in self.aux_defs:
            if kind == "eq":
                e, x, y = ins
                assign[z] = (True if e is None else assign[e]) and assign[x] == assign[y]
                continue
            if kind == "is0":
                assign[z] = all(val(t) == 0 for t in ins)
                continue
            x, y = val(ins[0]), val(ins[1])
            r = x * y % 3 if kind == "mul" else (x + y) % 3
            assign[z[0]], assign[z[1]] = r == 1, r == 2
        return assign

    def satisfied(self, assign) -> list:
        """Clauses violated by a full assignment (empty list means all satisfied)."""
        return [c for c in self.clauses if not any(assign.get(abs(l)) == (l > 0) for l in c)]
