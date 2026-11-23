"""Attempt 010, stream "local": S3 local search and local non-existence theorems over F3.

- tau_M = g_M o tau for symmetric nonsingular M (attempt 010, fact 2) and a scan of cyclic schemes for
  tau_M-invariance, with the change of basis g_N that makes an invariant scheme invariant under the standard tau.
- S3LocalF3CNF: the native S3 encoder (S3NativeF3CNF) with some S3 slots frozen to given orbits and the
  rest free, with symmetry breaking on the free slots only. A scheme that contains the frozen orbits and
  whose other orbits fit the free split has a representative: normalise each free orbit as s3_normal_form
  does and sort the free slots of each type; the frozen slots are untouched. So UNSAT means: no S3-invariant
  scheme contains the frozen orbits plus at most the free split. Optionally `forbid` asks that at least one
  free slot differs from every listed orbit (local search away from a known scheme).

Slot layout as in matmul.s3, with frozen slots first inside each type:
  orbit slots: size-6 bases [frozen | free], their tau partners, size-3 slots [frozen | free]
  fixed slots: size-2 bases [frozen | free], their transposes, size-1 slots [frozen | free]
"""

from __future__ import annotations

from itertools import product

from .cyclic import P, CyclicF3CNF, T, _inv, _mm, conj, expand, sigma, split_orbits
from .s3 import S3F3CNF, S3NativeF3CNF, _first, _sc, s3_split, tau, tensor_key
from .scheme import Scheme


def _red(m):
    return [[int(e) % P for e in r] for r in m]


def _det3(M):
    (a, b, c), (d, e, f), (g, h, i) = M
    return (a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) % P


def gl3():
    for k in product(range(P), repeat=9):
        M = [list(k[0:3]), list(k[3:6]), list(k[6:9])]
        if _det3(M):
            yield M


def symmetric_nonsingular_up_to_sign():
    """The 234 symmetric nonsingular 3x3 matrices over F3 up to sign (g_{-M} = g_M)."""
    out, seen = [], set()
    for k in product(range(P), repeat=6):
        M = [[k[0], k[1], k[2]], [k[1], k[3], k[4]], [k[2], k[4], k[5]]]
        if not _det3(M):
            continue
        key = tuple(e for r in M for e in r)
        neg = tuple(2 * e % P for e in key)
        if neg in seen:
            continue
        seen.add(key)
        out.append(M)
    return out


def tau_M(M, t):
    return conj(M, tau(t))


def terms_of(scheme: Scheme):
    return [tuple(_red(m) for m in (U, V, W)) for U, V, W in zip(scheme.U, scheme.V, scheme.W)]


def invariant_under(terms, g) -> bool:
    keys = {tensor_key(t) for t in terms}
    return all(tensor_key(g(t)) in keys for t in terms)


def tau_M_list(scheme: Scheme):
    """All M (symmetric nonsingular, up to sign) with the scheme tau_M-invariant."""
    terms = terms_of(scheme)
    return [M for M in symmetric_nonsingular_up_to_sign() if invariant_under(terms, lambda t: tau_M(M, t))]


def to_standard_s3(scheme: Scheme, M):
    """g_N(scheme) for an N in GL(3, F3) that makes a tau_M-invariant cyclic scheme invariant under the standard tau."""
    terms = terms_of(scheme)
    for N in gl3():
        G = _mm(_mm(N, M), T(N))
        G = _red(G)
        if G not in ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[2, 0, 0], [0, 2, 0], [0, 0, 2]]):
            continue
        new = [tuple(_red(m) for m in conj(N, t)) for t in terms]
        if invariant_under(new, tau) and invariant_under(new, sigma):
            return Scheme(3, 3, 3, [t[0] for t in new], [t[1] for t in new], [t[2] for t in new],
                          dict(scheme.meta, change_of_basis_N=N, tau_M=M))
    return None


def s3_slots(scheme: Scheme):
    """(six, three, two, one) of an S3-invariant (standard tau) scheme, as from matmul.s3.s3_split."""
    return s3_split(*split_orbits(scheme))


def slot_list(slots):
    """Flat list of (type, content) for an S3 slot tuple; types 6, 3, 2, 1."""
    six, three, two, one = slots
    return [(6, t) for t in six] + [(3, t) for t in three] + [(2, x) for x in two] + [(1, x) for x in one]


def group_slots(items):
    six = [c for ty, c in items if ty == 6]
    three = [c for ty, c in items if ty == 3]
    two = [c for ty, c in items if ty == 2]
    one = [c for ty, c in items if ty == 1]
    return six, three, two, one


RANK = {6: 6, 3: 3, 2: 2, 1: 1}


def slots_rank(items):
    return sum(RANK[ty] for ty, _ in items)


def norm6(t):
    u, v, w = (_red(m) for m in t)
    s, r = _first(u), _first(v)
    return (_sc(u, s), _sc(v, r), _sc(w, s * r))


def reps6(t):
    """Base values the encoder can hold for the size-6 orbit of t (all 6 members, sign-normalised)."""
    out = []
    for m in (t, sigma(t), sigma(sigma(t)), tau(t), sigma(tau(t)), sigma(sigma(tau(t)))):
        x = norm6(m)
        if x not in out:
            out.append(x)
    return out


def rep3(t):
    u, _, w = t
    s = _first(u)
    return (_sc(u, s), _sc(T(u), s), _red(w))


class S3LocalF3CNF(S3NativeF3CNF):
    """S3-invariant Brent equations mod 3 with frozen S3 slots plus a free split (a, b, c, d).

    frozen: (six, three, two, one) as from s3_split (size-3 terms in the form (u, u^T, w), w = w^T).
    forbid: optional (six, three, two, one); then at least one free slot must hold an orbit not in it.
    """

    def __init__(self, n, frozen, free, fixed_sum=True, forbid=None):
        six, three, two, one = frozen
        fa, fb, fc, fd = free
        self.nfz = (len(six), len(three), len(two), len(one))
        self.free = tuple(free)
        a, b, c, d = (len(six) + fa, len(three) + fb, len(two) + fc, len(one) + fd)
        self.s3 = (a, b, c, d)
        CyclicF3CNF.__init__(self, n, 2 * a + b, 2 * c + d, fixed_sum=fixed_sum, break_sym=False)
        za, zb, zc, zd = self.nfz
        c0 = len(self.clauses)
        unit = lambda vals, mat: [[lit] for z, x in zip(vals, (e % P for r in mat for e in r)) for lit in self._is(z, x)]
        for k, t in enumerate(six):
            for vals, mat in zip(self.orb[k], t):
                self.clauses += unit(vals, mat)
        for k, t in enumerate(three):
            for vals, mat in zip(self.orb[2 * a + k], rep3(t)):
                self.clauses += unit(vals, mat)
        for k, x in enumerate(two):
            self.clauses += unit(self.fix[k], x)
        for k, x in enumerate(one):
            self.clauses += unit(self.fix[2 * c + k], x)
        self.frozen_clauses = (c0, len(self.clauses))
        c0 = len(self.clauses)
        tr = lambda vals: [vals[j * n + i] for i in range(n) for j in range(n)]
        for k in range(za, a):
            self._first_nonzero_is_one(self.orb[k][0])
            self._first_nonzero_is_one(self.orb[k][1])
        for k in range(za, a - 1):
            self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
        for k in range(2 * a + zb, 2 * a + b):
            self._first_nonzero_is_one(self.orb[k][0])
        for k in range(2 * a + zb, 2 * a + b - 1):
            self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
        for k in range(zc, c):
            self._lex_leq(self._bits(self.fix[k]), self._bits(tr(self.fix[k])))
        for k in range(zc, c - 1):
            self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))
        for k in range(2 * c + zd, 2 * c + d - 1):
            self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))
        self.break_clauses = (c0, len(self.clauses))
        c0 = len(self.clauses)
        if forbid is not None:
            f6, f3, f2, f1 = forbid
            groups = [
                ([sum(self.orb[k], []) for k in range(za, a)], [r for t in f6 for r in reps6(t)]),
                ([sum(self.orb[k], []) for k in range(2 * a + zb, 2 * a + b)], [rep3(t) for t in f3]),
                ([self.fix[k] for k in range(zc, c)], [m for x in f2 for m in (_red(x), _red(T(x)))]),
                ([self.fix[k] for k in range(2 * c + zd, 2 * c + d)], [_red(x) for x in f1]),
            ]
            flags = []
            for slots, targets in groups:
                keys = []
                for tg in targets:
                    mats = tg if isinstance(tg, tuple) else (tg,)
                    key = tuple(e for m in mats for r in m for e in r)
                    if key not in keys:
                        keys.append(key)
                for vals in slots:
                    es = []
                    for key in keys:
                        e = self._new()
                        lits = [l for z, cc in zip(vals, key) for l in self._is(z, cc)]
                        self.clauses += [[-e, l] for l in lits]
                        self.clauses.append([e] + [-l for l in lits])
                        self.aux_defs.append(("alleq", e, tuple(zip(vals, key))))
                        es.append(e)
                    f = self._new()
                    self.clauses.append([-f] + es)
                    self.clauses += [[-e, f] for e in es]
                    self.aux_defs.append(("or", f, tuple(es)))
                    flags.append(f)
            self.clauses.append([-f for f in flags])
        self.forbid_clauses = (c0, len(self.clauses))

    def slot_contents(self, slots_frozen, slots_free):
        """(orbits, fixed) slot contents for witness(): frozen slots plus free slots (each list per type, may be
        shorter than the free split; missing slots are zero). Free contents are normalised and sorted."""
        from .s3 import s3_normal_form
        a, b, c, d = self.s3
        n = self.n
        Z = [[0] * n for _ in range(n)]
        f6, f3, f2, f1 = slots_frozen
        g6, g3, g2, g1 = slots_free
        o_free, x_free = s3_normal_form(g6, g3, g2, g1)
        b6 = o_free[:len(g6)]
        t3 = o_free[2 * len(g6):]
        tb = x_free[:len(g2)]
        on = x_free[2 * len(g2):]
        pad = lambda lst, k, z: [z] * (k - len(lst)) + list(lst)  # zero slots first (zero is lex smallest)
        bases = [tuple(_red(m) for m in t) for t in f6] + pad(b6, a - len(f6), (Z, Z, Z))
        threes = [rep3(t) for t in f3] + pad(t3, b - len(f3), (Z, Z, Z))
        twos = [_red(x) for x in f2] + pad(tb, c - len(f2), Z)
        ones = [_red(x) for x in f1] + pad(on, d - len(f1), Z)
        orbits = bases + [tau(t) for t in bases] + threes
        fixed = twos + [T(x) for x in twos] + ones
        return orbits, fixed


def s3_scheme(slots, meta=None) -> Scheme:
    """Expand S3 slots (six, three, two, one) to a scheme."""
    six, three, two, one = slots
    orbits = [tuple(_red(m) for m in t) for t in six] + [tau(t) for t in six] + [rep3(t) for t in three]
    fixed = [_red(x) for x in two] + [_red(T(x)) for x in two] + [_red(x) for x in one]
    return expand(3, orbits, fixed, meta or {})
