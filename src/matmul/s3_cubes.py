"""Cube decompositions of the S3-invariant F3 instances (attempt 010, stage C).

An instance S3NativeF3CNF(3, a, b, c, d) is split into cubes. Each cube is an ordinary CNF (the native S3
encoding without its lex chains) plus constraints that fix a coarse, conjugation-invariant type for every
orbit slot, and a new symmetry breaking that is compatible with those types.

Matrix classes (`gran`):
  "r1":   0 = rank one, 1 = anything else (zero included);
  "rank": 0 = rank 1, 1 = rank 2, 2 = rank 3, 3 = zero.

Slot keys:
  size-3 slot (u, u^T, w), w = w^T:  (mc(u),) or, with use_w, (mc(u), mc(w));
  size-6 slot base t = (u, v, w):     (mc(u), mc(v), mc(w)), sorted. The six terms of the orbit have first
      factors u, v, w^T, v^T, u^T, w, and the six choices of base (t, sigma t, sigma^2 t, tau t, ...) realise
      all six permutations of the triple, so a base with sorted triple exists.

A cube is (keys6, keys3): non-decreasing tuples of keys for the a size-6 bases and the b size-3 slots.

Normal form (the soundness lemma, `normalize`). Every assignment of slot contents can be moved, by
tensor-preserving moves only, to one that satisfies the clauses of exactly one cube:
  1. size-6: choose the base whose class triple is sorted, rescale (e u, f v, e f w) so the first nonzero of u
     and of v is 1 (the tau partner follows);
  2. size-3: rescale (e u, e u^T, w) so the first nonzero of u is 1;
  3. sort the size-6 bases by (key, bitkey) and the size-3 slots by (key, bitkey);
  4. size-2 / size-1: as s3_normal_form (x <= x^T, sorted; sorted).
  Then the cube clauses are: each slot has its key; lex order between consecutive slots of the same type
  only when they have the same key (sorting by (key, bitkey) gives exactly this).
Optional orthogonal pin (`pin`), used only when keys3[0] has u of rank one: among the size-3 slots whose
key equals keys3[0], pick the one whose u has
the smallest orthogonal class (17 rank-one classes; conjugation by M with M^T M = +-I, and sign), conjugate the whole scheme by
g_M so that this u is the class representative, and move that slot to position 0 (no sign rule, no lex
link to slot 1). The other size-3 slots with the same key must then have u in a class >= the pinned one:
their u avoids every member of the smaller classes (explicit forbidden-matrix clauses). Classes depend on
mc only through the key, since mc is a conjugation invariant. g_M with M orthogonal commutes with sigma and
tau, so it preserves every slot form (checked by tests on random contents).
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement, product

from .cyclic import P, T, _bitkey, _inv, _key, _mm, _unkey, conj, sigma
from .s3 import S3NativeF3CNF, tau

# ---- matrix classes -----------------------------------------------------------------------


def rank3(m):
    A = [[e % P for e in r] for r in m]
    r = 0
    for c in range(3):
        p = next((i for i in range(r, 3) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        iv = A[r][c]
        A[r] = [x * iv % P for x in A[r]]
        for i in range(3):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % P for x, y in zip(A[i], A[r])]
        r += 1
    return r


def mclass(m, gran):
    r = rank3(m)
    if gran == "r1":
        return 0 if r == 1 else 1
    if gran == "rank":
        return 3 if r == 0 else r - 1
    raise ValueError(gran)


RANK1 = 0  # the class of the rank-one matrices, in both granularities


def n_classes(gran):
    return {"r1": 2, "rank": 4}[gran]


# ---- the orthogonal group ------------------------------------------------------------------


@lru_cache(maxsize=None)
def orth_group():
    """All M with M^T M = +-I. Over F3 M^T M = -I has no solution (det^2 = -1 is a nonsquare), so this is
    O(3, F3), 48 matrices; M and -M act the same way, so the acting group has order 24."""
    I = [[int(i == j) for j in range(3)] for i in range(3)]
    mI = [[(-int(i == j)) % P for j in range(3)] for i in range(3)]
    out = []
    for k in product(range(P), repeat=9):
        M = _unkey(k, 3)
        G = _mm(T(M), M)
        if G == I or G == mI:
            out.append(k)
    return tuple(out)


def _neg(m):
    return [[-e % P for e in r] for r in m]


@lru_cache(maxsize=None)
def orth_classes():
    """Classes of 3x3 matrices under x -> M x M^-1 (M in orth_group) and sign. Returns (reps, index, conj_to_rep):
    reps sorted by bitkey (class representative = bitkey-smallest member), index[key] = class index,
    to_rep[key] = (s, Mkey) with M (s x) M^-1 = rep."""
    Ms = [(_unkey(k, 3), _inv(_unkey(k, 3))) for k in orth_group()]
    index, to_rep, reps = {}, {}, []
    classes = []
    seen = set()
    for k in product(range(P), repeat=9):
        if k in seen:
            continue
        x = _unkey(k, 3)
        members = {}
        for M, Mi in Ms:
            for s in (1, 2):
                y = _key(_mm(_mm(M, [[e * s % P for e in r] for r in x]), Mi))
                members.setdefault(y, (s, M))
        seen.update(members)
        classes.append((k, members))
    # representative: bitkey-smallest member
    order = []
    for k, members in classes:
        rep = min(members, key=lambda c: _bitkey([_unkey(c, 3)]))
        order.append((rep, k, members))
    order.sort(key=lambda t: _bitkey([_unkey(t[0], 3)]))
    for i, (rep, k, members) in enumerate(order):
        reps.append(rep)
        # M0 (s0 x0) M0^-1 = rep for the base element x0 = k
        s0, M0 = members[rep]
        for y, (s, M) in members.items():
            # y = M (s x0) M^-1  ->  x0 = s M^-1 y M; rep = M0 s0 x0 M0^-1 = (s0 s) (M0 M^-1) y (M M0^-1)
            index[y] = i
            to_rep[y] = ((s0 * s) % P, _key(_mm(M0, _inv(M))))
    return tuple(reps), index, to_rep


def orth_class_index(x):
    return orth_classes()[1][_key(x)]


def orth_to_rep(x):
    """(rep, s, M) with M (s x) M^-1 = rep, M orthogonal up to sign."""
    reps, index, to_rep = orth_classes()
    s, Mk = to_rep[_key(x)]
    return _unkey(reps[index[_key(x)]], 3), s, _unkey(Mk, 3)


def pin_choices(gran, key_u):
    """Orthogonal class indices of the nonzero matrices whose mc is key_u."""
    reps = orth_classes()[0]
    return tuple(i for i, r in enumerate(reps) if any(r) and mclass(_unkey(r, 3), gran) == key_u)


# ---- cube enumeration -----------------------------------------------------------------------


def slot_keys(gran, use_w):
    k = n_classes(gran)
    k3 = [(x,) for x in range(k)] if not use_w else [(x, y) for x in range(k) for y in range(k)]
    k6 = list(combinations_with_replacement(range(k), 3))
    return k6, k3


def cubes(split, gran="r1", use_w=False, pin=False):
    """All cubes (keys6, keys3, pin_class) of an S3 split (a, b, c, d)."""
    a, b, _, _ = split
    k6, k3 = slot_keys(gran, use_w)
    out = []
    for keys6 in combinations_with_replacement(k6, a):
        for keys3 in combinations_with_replacement(k3, b):
            if pin and b > 0 and keys3[0][0] == RANK1:
                for ci in pin_choices(gran, RANK1):
                    out.append((keys6, keys3, ci))
            else:
                out.append((keys6, keys3, None))
    return out


# ---- skew parts of size-6 slots (decomposition "sk", from the KKK identity N1) ------------------
#
# sk(U) = (U12 - U21, U13 - U31, U23 - U32). Over the six representatives of a size-6 orbit, (sk u, sk v, sk w)
# is permuted up to signs, and g_M (M orthogonal) acts on sk vectors as the 24 maps of SO(3, F3), whose orbits
# on nonzero vectors up to sign are the vectors with 1, 2 and 3 nonzero coordinates (6, 12, 8 vectors).
# skkey: 0, 1, 2 for 1, 2, 3 nonzero coordinates, 3 for the zero vector. SK_REPS: the pinned representatives.
# N1 (theory stream, writeups/010-s3-theory-notes.md): sum over size-6 bases of det[sk u, sk v, sk w] = -1.

SK_REPS = {0: (0, 0, 1), 1: (0, 1, 1), 2: (1, 1, 1)}


def sk(U):
    return ((U[0][1] - U[1][0]) % P, (U[0][2] - U[2][0]) % P, (U[1][2] - U[2][1]) % P)


def skkey(p):
    n = sum(1 for x in p if x % P)
    return 3 if n == 0 else n - 1


def _reps6(t):
    reps, x = [], t
    for _ in range(3):
        reps.append(x)
        reps.append(tau(x))
        x = sigma(x)
    return reps


def det3(p, q, r):
    return (p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0])
            + p[2] * (q[0] * r[1] - q[1] * r[0])) % P


def _orient6_sk(t, signs=True):
    """The representative whose sk u has the smallest skkey (then bitkey-smallest); u, v sign-normalised."""
    reps = _reps6(t)
    m = min(skkey(sk(r[0])) for r in reps)
    best = min((r for r in reps if skkey(sk(r[0])) == m), key=lambda r: _bitkey(list(r)))
    u, v, w = best
    s, r = (_first(u) if signs else 1), _first(v)
    return (_sc(u, s), _sc(v, r), _sc(w, s * r))


def key6_sk(t):
    return (skkey(sk(t[0])),)


@lru_cache(maxsize=None)
def sk_to_rep():
    """For each nonzero sk vector p: (M key, e) with e * sk(g_M u) = rep of its class, for any u with sk u = p."""
    out = {}
    for k in orth_group():
        M = _unkey(k, 3)
        for p in product(range(P), repeat=3):
            if not any(p):
                continue
            K = [[0, p[0], p[1]], [-p[0] % P, 0, p[2]], [-p[1] % P, -p[2] % P, 0]]
            q = tuple(2 * x % P for x in sk(conj(M, (K, K, K))[0]))  # sk(K) = 2p, and sk is linear
            rep = SK_REPS[skkey(p)]
            for e in (1, 2):
                if tuple(x * e % P for x in q) == rep and p not in out:
                    out[p] = (k, e)
    return out


def cubes_sk(split, gran="r1", use_w=False):
    """Cubes of decomposition "sk": size-6 keys (skkey of sk u, sorted; slot 0 pinned when its key < 3),
    size-3 keys as in `cubes` (no orthogonal pin there: the group is spent on slot 0's sk u)."""
    a, b, _, _ = split
    _, k3 = slot_keys(gran, use_w)
    out = []
    for keys6 in combinations_with_replacement([(x,) for x in range(4)], a):
        for keys3 in combinations_with_replacement(k3, b):
            out.append((keys6, keys3, "sk" if a and keys6[0][0] < 3 else None))
    return out


def normalize_sk(orbits, fixed, split, gran="r1", use_w=False):
    """Normal form for decomposition "sk" (same contract as `normalize`)."""
    a, b, c, d = split
    orbits = [tuple(_md(m) for m in t) for t in orbits]
    fixed = [_md(x) for x in fixed]
    six = [_orient6_sk(t) for t in orbits[:a]]
    three = list(orbits[2 * a:2 * a + b])
    two = list(fixed[:c])
    one = list(fixed[2 * c:2 * c + d])
    Mused = [[int(i == j) for j in range(3)] for i in range(3)]
    pinned = None
    if a and min(key6_sk(t) for t in six)[0] < 3:
        i0 = min(range(a), key=lambda i: (key6_sk(six[i]), _bitkey(list(six[i]))))
        Mk, e = sk_to_rep()[sk(six[i0][0])]
        M = _unkey(Mk, 3)
        Mused = M
        six = [conj(M, t) for t in six]
        three = [conj(M, t) for t in three]
        two = [conj(M, (x, x, T(x)))[0] for x in two]
        one = [conj(M, (x, x, T(x)))[0] for x in one]
        u, v, w = six[i0]
        u, w = _sc(u, e), _sc(w, e)
        r = _first(v)
        pinned = (u, _sc(v, r), _sc(w, r))
        assert sk(pinned[0]) == SK_REPS[key6_sk(pinned)[0]]
        six = [_orient6_sk(t) for i, t in enumerate(six) if i != i0]
    six = sorted(six, key=lambda t: (key6_sk(t), _bitkey(list(t))))
    if pinned is not None:
        six = [pinned] + six
    three = sorted((_norm3(t) for t in three), key=lambda t: (key3(t, gran, use_w), _bitkey(list(t))))
    tb = sorted((min(x, T(x), key=lambda m: _bitkey([m])) for x in two), key=lambda m: _bitkey([m]))
    on = sorted(one, key=lambda m: _bitkey([m]))
    new_orbits = six + [tau(t) for t in six] + three
    new_fixed = tb + [T(x) for x in tb] + on
    cube = (tuple(key6_sk(t) for t in six), tuple(key3(t, gran, use_w) for t in three),
            "sk" if pinned is not None else None)
    return new_orbits, new_fixed, cube, Mused


# ---- normal form ----------------------------------------------------------------------------


def _first(m):
    return next((e % P for r in m for e in r if e % P), 1)


def _sc(m, s):
    return [[e * s % P for e in r] for r in m]


def _md(m):
    return [[e % P for e in r] for r in m]


def key6(t, gran):
    return tuple(sorted(mclass(m, gran) for m in t))


def key3(t, gran, use_w):
    u, _, w = t
    return (mclass(u, gran), mclass(w, gran)) if use_w else (mclass(u, gran),)


def _orient6(t, gran):
    reps = []
    x = t
    for _ in range(3):
        reps.append(x)
        reps.append(tau(x))
        x = sigma(x)
    cands = [r for r in reps if [mclass(m, gran) for m in r] == sorted(mclass(m, gran) for m in r)]
    best = min(cands, key=lambda r: _bitkey(list(r)))
    u, v, w = best
    s, r = _first(u), _first(v)
    return (_sc(u, s), _sc(v, r), _sc(w, s * r))


def _norm3(t):
    u, _, w = t
    s = _first(u)
    return (_sc(u, s), _sc(T(u), s), _md(w))


def normalize(orbits, fixed, split, gran="r1", use_w=False, pin=False):
    """Slot contents in the S3NativeF3CNF layout -> (orbits, fixed, cube, M) in normal form (see module doc);
    M is the orthogonal matrix applied (identity without pin).

    orbits: 2a + b terms (size-6 bases, their tau partners, size-3 slots (u, u^T, w)); fixed: 2c + d matrices.
    The expanded tensor is unchanged.
    """
    a, b, c, d = split
    orbits = [tuple(_md(m) for m in t) for t in orbits]
    fixed = [_md(x) for x in fixed]
    six = list(orbits[:a])
    three = list(orbits[2 * a:2 * a + b])
    two = list(fixed[:c])
    one = list(fixed[2 * c:2 * c + d])
    pin_class = None
    Mused = [[int(i == j) for j in range(3)] for i in range(3)]
    if pin and b > 0 and min(key3(t, gran, use_w) for t in three)[0] == RANK1:
        k0 = min(key3(t, gran, use_w) for t in three)
        cand = [(orth_class_index(t[0]), i) for i, t in enumerate(three) if key3(t, gran, use_w) == k0]
        pin_class, i0 = min(cand)
        rep, s, M = orth_to_rep(three[i0][0])
        Mused = M
        # apply g_M to everything (slot forms are preserved; sign s fixed on the pinned slot below)
        six = [conj(M, t) for t in six]
        three = [conj(M, t) for t in three]
        two = [conj(M, (x, x, T(x)))[0] for x in two]
        one = [conj(M, (x, x, T(x)))[0] for x in one]
        u, _, w = three[i0]
        u = _sc(u, s)
        assert _key(u) == _key(rep)
        three[i0] = (u, T(u), _md(w))
        p0 = three.pop(i0)
    six = sorted((_orient6(t, gran) for t in six), key=lambda t: (key6(t, gran), _bitkey(list(t))))
    three = sorted((_norm3(t) for t in three), key=lambda t: (key3(t, gran, use_w), _bitkey(list(t))))
    if pin_class is not None:
        three = [p0] + three
    tb = sorted((min(x, T(x), key=lambda m: _bitkey([m])) for x in two), key=lambda m: _bitkey([m]))
    on = sorted(one, key=lambda m: _bitkey([m]))
    new_orbits = six + [tau(t) for t in six] + three
    new_fixed = tb + [T(x) for x in tb] + on
    cube = (tuple(key6(t, gran) for t in six), tuple(key3(t, gran, use_w) for t in three), pin_class)
    return new_orbits, new_fixed, cube, Mused


# ---- the cube encoder -----------------------------------------------------------------------


class S3CubeCNF(S3NativeF3CNF):
    """S3NativeF3CNF(3, a, b, c, d) restricted to one cube; see the module doc."""

    def __init__(self, split, cube, gran="r1", use_w=False, extra_units=None, n1=False):
        a, b, c, d = split
        super().__init__(3, a, b, c, d, fixed_sum=True, break_sym=False)
        self.split, self.cube, self.gran, self.use_w = tuple(split), cube, gran, use_w
        keys6, keys3, pin_class = cube
        assert len(keys6) == a and len(keys3) == b
        self.base_clauses = len(self.clauses)
        tr = lambda vals: [vals[j * 3 + i] for i in range(3) for j in range(3)]
        skmode = len(keys6) > 0 and len(keys6[0]) == 1  # decomposition "sk" (keys6 entries are 1-tuples)
        if pin_class == "sk":
            pin_class = None  # the size-3 slots are not pinned in this decomposition
            sk_pinned = True
        else:
            sk_pinned = False
        if n1:
            # (N1) sum over size-6 bases of det[sk u, sk v, sk w] = -1 (implied; redundant)
            dets = []
            for k in range(a):
                p, q, r = (self._skvals(m) for m in self.orb[k])
                dets.append(self._det(p, q, r))
            terms = [t for dd in dets for t in dd]
            self._sum_equals(terms, 2)
        self.base_clauses_n1 = len(self.clauses)
        if skmode:
            for k in range(a):
                u, v, w = self.orb[k]
                ck = keys6[k][0]
                self._skkey_in(self._skvals(u), {ck})
                self._skkey_in(self._skvals(v), set(range(ck, 4)))
                self._skkey_in(self._skvals(w), set(range(ck, 4)))
                if k == 0 and sk_pinned:
                    for z, x in zip(self._skvals(u), SK_REPS[ck]):
                        for lit in self._is(z, x):
                            self.clauses.append([lit])
                else:
                    self._first_nonzero_is_one(u)
                self._first_nonzero_is_one(v)
            for k in range(1 if sk_pinned else 0, a - 1):
                if keys6[k] == keys6[k + 1]:
                    self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
        for k in (range(0) if skmode else range(a)):
            u, v, w = self.orb[k]
            for m, cl in zip((u, v, w), keys6[k]):
                self._in_class(m, cl)
            self._first_nonzero_is_one(u)
            self._first_nonzero_is_one(v)
        for k in (range(0) if skmode else range(a - 1)):
            if keys6[k] == keys6[k + 1]:
                self._lex_leq(self._bits(sum(self.orb[k], [])), self._bits(sum(self.orb[k + 1], [])))
        # size 3
        for j in range(b):
            u, _, w = self.orb[2 * a + j]
            self._in_class(u, keys3[j][0])
            if use_w:
                self._in_class(w, keys3[j][1])
            if not (j == 0 and pin_class is not None):
                self._first_nonzero_is_one(u)
        for j in range(b - 1):
            if j == 0 and pin_class is not None:
                continue
            if keys3[j] == keys3[j + 1]:
                self._lex_leq(self._bits(sum(self.orb[2 * a + j], [])), self._bits(sum(self.orb[2 * a + j + 1], [])))
        if pin_class is not None:
            reps, index, _ = orth_classes()
            u0 = self.orb[2 * a][0]
            target = reps[pin_class]
            for z, x in zip(u0, target):
                for lit in self._is(z, x):
                    self.clauses.append([lit])
            # other size-3 slots with the same key: u not in a smaller orthogonal class (all rank one)
            bad = [k for k, i in index.items() if i < pin_class and mclass(_unkey(k, 3), self.gran) == RANK1]
            for j in range(1, b):
                if keys3[j] != keys3[0]:
                    continue
                u = self.orb[2 * a + j][0]
                for k in bad:
                    self.clauses.append([-l for z, x in zip(u, k) for l in self._is(z, x)])
        # size 2 / size 1 (as S3F3CNF._s3_break)
        for k in range(c):
            self._lex_leq(self._bits(self.fix[k]), self._bits(tr(self.fix[k])))
        for k in range(c - 1):
            self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))
        for k in range(2 * c, 2 * c + d - 1):
            self._lex_leq(self._bits(self.fix[k]), self._bits(self.fix[k + 1]))
        for z, x in (extra_units or []):
            for lit in self._is(z, x):
                self.clauses.append([lit])

    def _skvals(self, vals):
        key = ("sk",) + tuple(vals)
        if key not in self._cp:
            x = lambda i, j: vals[i * 3 + j]
            neg = lambda z: (z[1], z[0])
            self._cp[key] = [self._add(x(0, 1), neg(x(1, 0))), self._add(x(0, 2), neg(x(2, 0))),
                             self._add(x(1, 2), neg(x(2, 1)))]
        return self._cp[key]

    def _det(self, p, q, r):
        """The six signed products of det[p, q, r] (F3 values; a minus sign swaps the pair)."""
        neg = lambda z: (z[1], z[0])
        out = []
        for (i, j, k), sg in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((1, 0, 2), -1),
                              ((2, 1, 0), -1)):
            z = self._prod2(self._prod2(p[i], q[j]), r[k])
            out.append(z if sg == 1 else neg(z))
        return out

    def _skkey_in(self, pv, allowed):
        """skkey of the sk vector pv lies in `allowed` (0, 1, 2: 1, 2, 3 nonzero coordinates; 3: zero)."""
        if allowed >= {0, 1, 2, 3}:
            return
        nz = []
        for z in pv:
            f = self._new()
            self.clauses += [[-f, z[0], z[1]], [f, -z[0]], [f, -z[1]]]
            self.aux_defs.append(("or", f, z))
            nz.append(f)
        for pat in product((0, 1), repeat=3):
            cnt = sum(pat)
            key = 3 if cnt == 0 else cnt - 1
            if key not in allowed:
                self.clauses.append([-f if bit else f for f, bit in zip(nz, pat)])

    def _minors(self, vals):
        x = lambda i, j: vals[i * 3 + j]
        neg = lambda z: (z[1], z[0])
        return [self._add(self._prod2(x(i1, j1), x(i2, j2)), neg(self._prod2(x(i1, j2), x(i2, j1))))
                for i1 in range(3) for i2 in range(i1 + 1, 3) for j1 in range(3) for j2 in range(j1 + 1, 3)]

    def _in_class(self, vals, cl):
        g = self.gran
        if g == "r1":
            if cl == 0:
                self._rank1(vals)
            else:
                self._not_rank1(vals)
            return
        # "rank": 0 = rank 1, 1 = rank 2, 2 = rank 3, 3 = zero
        if cl == 3:
            for pq in vals:
                for l in pq:
                    self.clauses.append([-l])
        elif cl == 0:
            self._rank1(vals)
        elif cl == 1:
            self._not_rank1(vals)
            self.clauses.append(self._bits(vals))  # nonzero
            det = self._charpoly(vals)[2]
            self.clauses += [[-det[0]], [-det[1]]]
        else:
            det = self._charpoly(vals)[2]
            self.clauses.append([det[0], det[1]])

    def _rank1(self, vals):
        self.clauses.append(self._bits(vals))
        for m in self._minors(vals):
            self.clauses += [[-m[0]], [-m[1]]]


def tensor_of(orbits, fixed):
    """The expanded tensor mod 3 (27 x 27 x 27 dict of nonzero entries) of slot contents, for tests."""
    from .cyclic import expand
    s = expand(3, orbits, fixed)
    acc = {}
    for U, V, W in zip(s.U, s.V, s.W):
        fu = [int(x) % P for r in U for x in r]
        fv = [int(x) % P for r in V for x in r]
        fw = [int(x) % P for r in W for x in r]
        for i, x in enumerate(fu):
            if not x:
                continue
            for j, y in enumerate(fv):
                if not y:
                    continue
                for k, z in enumerate(fw):
                    if z:
                        acc[(i, j, k)] = (acc.get((i, j, k), 0) + x * y * z) % P
    return {k: v for k, v in acc.items() if v}


# ---- decompositions by name ------------------------------------------------------------------


def parse_cfg(cfg):
    """'gran[+w][+pin|+sk][+n1]' -> dict. n1 adds the KKK identity as a redundant constraint (no new cubes)."""
    parts = cfg.split("+")
    opts = dict(gran=parts[0], use_w="w" in parts[1:], pin="pin" in parts[1:], sk="sk" in parts[1:],
                n1="n1" in parts[1:])
    if opts["pin"] and opts["sk"]:
        raise ValueError("pin and sk both spend the orthogonal group")
    return opts


def cubes_cfg(split, cfg):
    o = parse_cfg(cfg)
    if o["sk"]:
        return cubes_sk(split, o["gran"], o["use_w"])
    return cubes(split, o["gran"], o["use_w"], o["pin"])


def normalize_cfg(orbits, fixed, split, cfg):
    o = parse_cfg(cfg)
    if o["sk"]:
        return normalize_sk(orbits, fixed, split, o["gran"], o["use_w"])
    return normalize(orbits, fixed, split, o["gran"], o["use_w"], o["pin"])


def encoder_cfg(split, cube, cfg):
    o = parse_cfg(cfg)
    return S3CubeCNF(split, cube, o["gran"], o["use_w"], n1=o["n1"])
