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


# ---- CNF encoding over F3 ------------------------------------------------------------

class CyclicF3CNF:
    """Brent equations mod 3 for a sigma-invariant n x n x n scheme with a orbit slots and b fixed slots.

    Options:
      fixed_sum: add the necessary condition "fixed matrices sum to I" (implied, but helps).
      break_sym: sign-normalise each orbit slot (first nonzero of u and of v is 1), and order the
                 orbit slots and the fixed slots lexicographically. Both keep at least one
                 representative of every solution (rescaling and reordering are symmetries).
    """

    def __init__(self, n: int, a: int, b: int, fixed_sum: bool = True, break_sym: bool = True):
        self.n, self.a, self.b = n, a, b
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
        if break_sym:
            for k in range(a):
                self._first_nonzero_is_one(self.orb[k][0])
                self._first_nonzero_is_one(self.orb[k][1])
            for k in range(a - 1):
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
            x, y = val(ins[0]), val(ins[1])
            r = x * y % 3 if kind == "mul" else (x + y) % 3
            assign[z[0]], assign[z[1]] = r == 1, r == 2
        return assign

    def satisfied(self, assign) -> list:
        """Clauses violated by a full assignment (empty list means all satisfied)."""
        return [c for c in self.clauses if not any(assign.get(abs(l)) == (l > 0) for l in c)]
