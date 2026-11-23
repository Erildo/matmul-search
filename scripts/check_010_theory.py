"""Numerical checks for writeups/010-s3-theory-notes.md (S3-invariant 3x3x3 schemes over F3).

Convention as in matmul.cyclic / matmul.s3: a term (U, V, W) is the trilinear form
l_U(A) l_V(B) l_W(C) with l_U(A) = sum U_ij A_ij, and the scheme computes tr(A B C^T).

Checks (each prints PASS/FAIL; exit 1 on any failure):
  1. block contributions of every S3 orbit type, on random S3-shaped orbits:
       KKK (A, B, C skew):     size 3, 2, 1 orbits restrict to 0; a size-6 orbit restricts to c0 * det[sk u, sk v, sk w] * eps
       KKS (A, B skew, C sym): C-flattening rank <= 3, 1, 1, 0 for sizes 6, 3, 2, 1
       SSK (A, B sym, C skew): size-6 orbit = sum of 3 wedge terms (alpha^beta) (x) gamma, size 3 = 1 wedge term, fixed = 0
  2. the blocks of the matrix multiplication tensor itself: T_KKK = kappa * eps (kappa != 0), C-flattening rank of T_KKS = 6
  3. the necessary conditions on the known S3 schemes (006 in S3 form; Laderman conjugated into standard S3):
       fixed-sum, sum over size-6 orbits of det[sk u, sk v, sk w] = kappa / c0, 3a + b + c >= 6
  4. Laderman's algorithm (Burichenko's tensor form) is invariant under Burichenko's S3 = <Phi3, Phi4> over F3 and over F2,
     and is equivalent to a scheme invariant under the standard <sigma, tau> (explicit change of basis, verified);
     its S3 split is (1,4,1,3); its class key differs from the 006 scheme's
  5. the SSK wedge rank: the minimum number of decomposable 2-forms whose span contains the 3 slices of T_SSK
     (exhaustive over subspaces), which bounds 3a + b from below
  6. which of the 60 rank-22 splits each condition kills

Also checks every other S3-invariant rank-23 scheme found in results/f3/010*.json (whatever is present).

usage: python3 -I scripts/check_010_theory.py [--save-laderman PATH]
"""

from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import is_cyclic_invariant, sigma, split_orbits  # noqa: E402
from matmul.invariants import class_key  # noqa: E402
from matmul.s3 import is_s3_invariant, s3_split, s3_splits, tau  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402
from matmul.verify import verify  # noqa: E402

P = 3
FAIL = []


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{info}]" if info else ""))
    if not ok:
        FAIL.append(name)


# ---- linear algebra mod p ------------------------------------------------------------

def rank_p(M, p=P):
    a = np.array(M, dtype=np.int64) % p
    a = a.reshape(a.shape[0], -1) if a.ndim > 1 else a.reshape(1, -1)
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if a[i, c]), None)
        if piv is None:
            continue
        a[[r, piv]] = a[[piv, r]]
        a[r] = a[r] * pow(int(a[r, c]), -1, p) % p
        for i in range(rows):
            if i != r and a[i, c]:
                a[i] = (a[i] - a[i, c] * a[r]) % p
        r += 1
        if r == rows:
            break
    return r


def inv_p(m, p=P):
    m = np.array(m, dtype=np.int64) % p
    n = len(m)
    a = np.concatenate([m, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        piv = next(i for i in range(c, n) if a[i, c] % p)
        a[[c, piv]] = a[[piv, c]]
        a[c] = a[c] * pow(int(a[c, c]), -1, p) % p
        for i in range(n):
            if i != c and a[i, c]:
                a[i] = (a[i] - a[i, c] * a[c]) % p
    return a[:, n:]


# ---- Sym / Skew bases and block restriction ---------------------------------------

def E(i, j):
    m = np.zeros((3, 3), dtype=np.int64)
    m[i, j] = 1
    return m


PAIRS = [(0, 1), (0, 2), (1, 2)]
SYM = [E(i, i) for i in range(3)] + [E(i, j) + E(j, i) for i, j in PAIRS]
SKEW = [E(i, j) - E(j, i) for i, j in PAIRS]
BASIS = {"S": SYM, "K": SKEW}


def ell(U, X):
    return int((np.array(U) * X).sum())


def restrict(terms, kinds, p=P):
    """Restriction of sum of rank-one trilinear forms to (A, B, C) in the subspaces named by kinds, e.g. 'KKS'."""
    bs = [BASIS[k] for k in kinds]
    out = np.zeros(tuple(len(b) for b in bs), dtype=np.int64)
    for U, V, W in terms:
        x = np.array([ell(U, X) for X in bs[0]])
        y = np.array([ell(V, X) for X in bs[1]])
        z = np.array([ell(W, X) for X in bs[2]])
        out += np.einsum("i,j,k->ijk", x, y, z)
    return out % p


def mm_block(kinds, p=P):
    bs = [BASIS[k] for k in kinds]
    return np.array([[[int(np.trace(A @ B @ C.T)) for C in bs[2]] for B in bs[1]] for A in bs[0]]) % p


EPS = np.zeros((3, 3, 3), dtype=np.int64)
for perm in itertools.permutations(range(3)):
    EPS[perm] = round(np.linalg.det(np.eye(3)[list(perm)]))
EPS %= P


def sk(U):
    U = np.array(U)
    return np.array([U[i, j] - U[j, i] for i, j in PAIRS]) % P


def det3(a, b, c):
    return int(round(np.linalg.det(np.array([a, b, c], dtype=float)))) % P


# ---- S3 orbits -------------------------------------------------------------------------

def L(m):
    return [[int(e) % P for e in r] for r in np.array(m)]


def Tm(m):
    return L(np.array(m).T)


def orbit_terms(kind, data):
    if kind == 6:
        t = tuple(L(m) for m in data)
        sig = [t, sigma(t), sigma(sigma(t))]
        return sig + [tau(s) for s in sig]
    if kind == 3:
        u, w = data
        t = (L(u), Tm(u), L(w))
        return [t, sigma(t), sigma(sigma(t))]
    if kind == 2:
        x = L(data)
        return [(x, x, Tm(x)), (Tm(x), Tm(x), x)]
    x = L(data)
    return [(x, x, x)]


def rnd(rng, sym=False):
    m = np.array([[rng.randrange(3) for _ in range(3)] for _ in range(3)])
    return (m + m.T) % P if sym else m


def random_orbit(kind, rng):
    if kind == 6:
        return (rnd(rng), rnd(rng), rnd(rng))
    if kind == 3:
        return (rnd(rng), rnd(rng, sym=True))
    if kind == 2:
        return rnd(rng)
    return rnd(rng, sym=True)


def is_wedge_sum(blk, k):
    """Is the 6x6x3 SSK block a sum of k terms (alpha ^ beta) (x) gamma? Only for k in {0, 1}; k = 3 uses the pair form."""
    if k == 0:
        return not blk.any()
    alt = all(((blk[:, :, c] + blk[:, :, c].T) % P == 0).all() for c in range(3))
    flat_c = rank_p(blk.reshape(36, 3).T)
    return alt and flat_c <= 1 and all(rank_p(blk[:, :, c]) <= 2 for c in range(3))


def pair_wedges(terms6):
    """The SSK block of a size-6 orbit, split into its three tau-pairs {s, tau s}; each must be a single wedge term."""
    sig = terms6[:3]
    return [restrict([s, tau(s)], "SSK") for s in sig]


# ---- 1, 2: orbit-type contributions and the tensor's own blocks -----------------------

def check_blocks(samples=400, seed=140001):
    rng = random.Random(seed)
    T_kkk = mm_block("KKK")
    nz = np.argwhere(EPS)[0]
    kappa = int(T_kkk[tuple(nz)]) * int(EPS[tuple(nz)]) % P
    check("T_KKK = kappa * eps with kappa != 0", kappa != 0 and ((T_kkk - kappa * EPS) % P == 0).all(), f"kappa = {kappa}")
    T_kks = mm_block("KKS")
    rk = rank_p(T_kks.reshape(9, 6).T)
    check("C-flattening rank of T_KKS (Skew x Skew x Sym) is 6", rk == 6, f"rank {rk}")
    c0s = set()
    bad = {k: 0 for k in ("kkk", "kks", "ssk")}
    maxr = {6: 0, 3: 0, 2: 0, 1: 0}
    for _ in range(samples):
        for kind in (6, 3, 2, 1):
            d = random_orbit(kind, rng)
            terms = orbit_terms(kind, d)
            kkk = restrict(terms, "KKK")
            if kind == 6:
                det = det3(sk(d[0]), sk(d[1]), sk(d[2]))
                coef = int(kkk[tuple(nz)]) * int(EPS[tuple(nz)]) % P
                if not ((kkk - coef * EPS) % P == 0).all():
                    bad["kkk"] += 1
                if det:
                    c0s.add(coef * pow(det, -1, P) % P)
                elif coef:
                    bad["kkk"] += 1
            elif kkk.any():
                bad["kkk"] += 1
            kks = restrict(terms, "KKS")
            maxr[kind] = max(maxr[kind], rank_p(kks.reshape(9, 6).T))
            ssk = restrict(terms, "SSK")
            if kind == 6:
                pw = pair_wedges(terms)
                if not all(is_wedge_sum(b, 1) for b in pw) or ((sum(pw) - ssk) % P).any():
                    bad["ssk"] += 1
            elif kind == 3:
                if not is_wedge_sum(ssk, 1):
                    bad["ssk"] += 1
            elif ssk.any():
                bad["ssk"] += 1
    check("KKK: size 3/2/1 orbits vanish, size-6 orbit = c0 * det[sk u, sk v, sk w] * eps with one constant c0",
          bad["kkk"] == 0 and len(c0s) == 1, f"c0 = {sorted(c0s)}, {samples} samples per type")
    check("KKS: C-flattening rank per orbit <= 3, 1, 1, 0 (sizes 6, 3, 2, 1)",
          maxr[6] <= 3 and maxr[3] <= 1 and maxr[2] <= 1 and maxr[1] == 0, f"max ranks {maxr}")
    check("SSK: size-6 orbit = 3 wedge terms (one per tau-pair), size-3 = 1 wedge term, fixed terms = 0", bad["ssk"] == 0)
    return kappa, c0s.pop()


# ---- 3: the known S3 schemes ------------------------------------------------------------

def s3_parts(scheme):
    o, f = split_orbits(scheme)
    return s3_split(o, f)


def check_scheme(name, scheme, kappa, c0):
    ok = verify(scheme, mod=P).ok and is_cyclic_invariant(scheme) and is_s3_invariant(scheme)
    six, three, two, one = s3_parts(scheme)
    a, b, c, d = len(six), len(three), len(two), len(one)
    check(f"{name}: verifies mod 3, sigma- and tau-invariant", ok, f"split ({a},{b},{c},{d}), key {class_key(scheme, P)}")
    fs = (sum((np.array(x) for x in one), np.zeros((3, 3), dtype=np.int64))
          + sum((np.array(x) + np.array(x).T for x in two), np.zeros((3, 3), dtype=np.int64))) % P
    check(f"{name}: fixed-sum  sum_1 x + sum_2 (x + x^T) = I", (fs == np.eye(3, dtype=np.int64)).all())
    s = sum(det3(sk(u), sk(v), sk(w)) for u, v, w in six) % P
    check(f"{name}: KKK identity  sum_6 det[sk u, sk v, sk w] = kappa / c0", s == kappa * pow(c0, -1, P) % P,
          f"sum = {s}, kappa/c0 = {kappa * pow(c0, -1, P) % P}")
    check(f"{name}: 3a + b + c >= 6", 3 * a + b + c >= 6, f"3a+b+c = {3 * a + b + c}")
    # every block of the scheme's tensor equals the block of the matrix multiplication tensor
    terms = list(zip(scheme.U, scheme.V, scheme.W))
    terms = [tuple(L([[int(e) for e in r] for r in m]) for m in t) for t in terms]
    allb = all((restrict(terms, k) == mm_block(k)).all() for k in ("KKK", "KKS", "KSK", "SKK", "SSK", "SKS", "KSS", "SSS"))
    check(f"{name}: all 8 Sym/Skew blocks match the tensor", allb)
    return (a, b, c, d)


# ---- 4: Laderman --------------------------------------------------------------------------

LADERMAN = [  # Burichenko (arXiv:1508.01110, Prop. 5.1): t = x (x) y (x) z, tensor sum e_ij (x) e_jk (x) e_ki
    ("+11 +12 +13 -21 -22 -32 -33", "+22", "+21"), ("+11 -21", "-12 +22", "+12 +22"),
    ("+22", "-11 +12 +21 -22 -23 -31 +33", "+12"), ("-11 +21 +22", "+11 -12 +22", "+21 +12 +22"),
    ("+21 +22", "-11 +12", "+21 +22"), ("+11", "+11", "+11 +21 +31 +12 +22 +13 +33"),
    ("-11 +31 +32", "+11 -13 +23", "+31 +13 +33"), ("-11 +31", "+13 -23", "+13 +33"),
    ("+31 +32", "-11 +13", "+31 +33"), ("+11 +12 +13 -22 -23 -31 -32", "+23", "+31"),
    ("+32", "-11 +13 +21 -22 -23 -31 +32", "+13"), ("-13 +32 +33", "+22 +31 -32", "+21 +13 +23"),
    ("+13 -33", "+22 -32", "+13 +23"), ("+13", "+31", "+11 +21 +31 +12 +32 +13 +23"),
    ("+32 +33", "-31 +32", "+21 +23"), ("-13 +22 +23", "+23 +31 -33", "+31 +12 +32"),
    ("+13 -23", "+23 -33", "+12 +32"), ("+22 +23", "-31 +33", "+31 +32"),
    ("+12", "+21", "+11"), ("+23", "+32", "+22"), ("+21", "+13", "+32"), ("+31", "+12", "+23"), ("+33", "+33", "+33")]


def _mat(s):
    m = np.zeros((3, 3), dtype=np.int64)
    for tok in s.split():
        m[int(tok[1]) - 1, int(tok[2]) - 1] += 1 if tok[0] == "+" else -1
    return m


def laderman_b():
    return [tuple(_mat(s) for s in t) for t in LADERMAN]


def to_scheme(btriples, p, meta=None):
    """Burichenko (x, y, z) -> this repo's (U, V, W) = (x, y, z^T)."""
    return Scheme(3, 3, 3, [L2(x, p) for x, _, _ in btriples], [L2(y, p) for _, y, _ in btriples],
                  [L2(z.T, p) for _, _, z in btriples], meta or {})


def L2(m, p):
    return [[int(e) % p for e in r] for r in np.array(m)]


def tkey(t, p):
    x, y, z = [np.array([[int(e) for e in r] for r in np.array(m)], dtype=np.int64) % p for m in t]
    return (np.einsum("a,b,c->abc", x.ravel(), y.ravel(), z.ravel()) % p).tobytes()


PI12 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
EP1 = np.diag([-1, 1, 1])
EP2 = np.diag([1, -1, 1])


def phi3(t):
    x, y, z = t
    return (y.T @ EP2, EP2 @ x.T, z.T)


def phi4(t):
    x, y, z = t
    return (EP1 @ z @ PI12, PI12 @ x @ PI12 @ EP1, EP1 @ PI12 @ y @ EP1)


def check_laderman(p):
    B = laderman_b()
    lad = to_scheme(B, p)
    labz = Scheme(3, 3, 3, [x.tolist() for x, _, _ in B], [y.tolist() for _, y, _ in B], [z.T.tolist() for _, _, z in B], {})
    check(f"Laderman (Burichenko's tensor form) verifies over Z and mod {p}", verify(labz).ok and verify(lad, mod=p).ok)
    keys = {tkey(t, p) for t in B}
    check(f"Laderman is invariant under Phi3 and Phi4 mod {p}",
          all(tkey(phi3(t), p) in keys and tkey(phi4(t), p) in keys for t in B))
    # Phi4 = T(I, pi12, eps1) s^-1 T(I, pi12, eps1)^-1 with s: x(x)y(x)z -> y(x)z(x)x (this repo's sigma);
    # conjugate it back, then search the g_N (they commute with sigma) for one that makes tau a symmetry.
    # T(a, b, c): (x, y, z) -> (a x b^-1, b y c^-1, c z a^-1); apply it with a = I, b = pi12^-1, c = eps1^-1
    Bc = [((x @ PI12) % p, (inv_p(PI12, p) @ y @ EP1) % p, (inv_p(EP1 % p, p) @ z) % p) for x, y, z in B]
    cyc = to_scheme(Bc, p)
    check(f"T(I, pi12, eps1)^-1 . Laderman is sigma-invariant mod {p}", verify(cyc, mod=p).ok and is_cyclic_invariant(cyc, p))
    terms = [tuple(np.array([[int(e) for e in r] for r in m], dtype=np.int64) for m in t) for t in zip(cyc.U, cyc.V, cyc.W)]
    hit = None
    for flat in itertools.product(range(p), repeat=9):
        N = np.array(flat).reshape(3, 3)
        if round(np.linalg.det(N)) % p == 0:
            continue
        Ni = inv_p(N, p)
        tt = [((N @ u @ Ni) % p, (N @ v @ Ni) % p, (Ni.T @ w @ N.T) % p) for u, v, w in terms]
        s = Scheme(3, 3, 3, [L2(t[0], p) for t in tt], [L2(t[1], p) for t in tt], [L2(t[2], p) for t in tt], {})
        if p == 3 and is_s3_invariant(s):
            hit = (N, s)
            break
        if p == 2 and _s3_inv_mod(tt, p):
            hit = (N, s)
            break
    check(f"Laderman is equivalent mod {p} to a standard <sigma, tau>-invariant scheme",
          hit is not None and verify(hit[1], mod=p).ok, f"g_N with N = {hit[0].tolist() if hit else None}")
    return hit[1] if hit else None


def _s3_inv_mod(tt, p):
    K = {tkey(t, p) for t in tt}
    def tau_np(t):
        u, v, w = t
        return (v.T, u.T, w.T)
    def sig_np(t):
        u, v, w = t
        return (v, w.T, u.T)
    return len(K) == len(tt) and all(tkey(sig_np(t), p) in K and tkey(tau_np(t), p) in K for t in tt)


# ---- 5: SSK wedge rank ------------------------------------------------------------------

def ssk_wedge_rank(max_r):
    """Minimum r such that the 3 slices of T_SSK (restricted to a complement of I in Sym) lie in the span of r
    decomposable 2-forms. Exhaustive: r = 3 + k, over all k-dimensional subspaces of Lambda^2 / slices."""
    T = mm_block("SSK")[1:, 1:, :]           # drop E11: Sym = <E11> + rest, and I = E11 + E22 + E33 is not in the rest
    # the restriction to any complement of I gives the same wedge rank (I is in the radical of every slice)
    rad = all((mm_block("SSK")[:, :, c] @ np.array([1, 1, 1, 0, 0, 0]) % P == 0).all() for c in range(3))
    check("I is in the radical of every SSK slice", rad)
    n = 5
    idx = [(i, j) for i in range(n) for j in range(i + 1, n)]       # coordinates of Lambda^2 F3^5 (10)
    def vec(S):
        return np.array([S[i, j] % P for i, j in idx])
    slices = [vec(T[:, :, c]) for c in range(3)]
    Lmat = np.array(slices)
    check("the three SSK slices are linearly independent", rank_p(Lmat) == 3)
    # decomposable <=> the 5x5 skew matrix has rank <= 2 <=> all 4x4 Pfaffians vanish
    allv = np.array(list(itertools.product(range(P), repeat=10)), dtype=np.int64)
    def pf4(v, a, b, c, d):
        g = lambda i, j: v[:, idx.index((i, j))]  # noqa: E731
        return (g(a, b) * g(c, d) - g(a, c) * g(b, d) + g(a, d) * g(b, c)) % P
    dec = np.ones(len(allv), dtype=bool)
    for q in itertools.combinations(range(n), 4):
        dec &= pf4(allv, *q) == 0
    pw = 3 ** np.arange(9, -1, -1)
    decset = np.zeros(3 ** 10, dtype=bool)
    decset[(allv @ pw)] = dec
    # complement of span(slices): pick 7 unit vectors extending the slices to a basis
    basis = [s for s in slices]
    for e in range(10):
        u = np.zeros(10, dtype=np.int64)
        u[e] = 1
        if rank_p(np.array(basis + [u])) > len(basis):
            basis.append(u)
    comp = np.array(basis[3:])
    span_L = np.array([(np.array(c) @ Lmat) % P for c in itertools.product(range(P), repeat=3)])
    result = None
    for k in range(0, max_r - 2):
        found = False
        # all k-dim subspaces of the complement (represented by reduced row echelon forms over F3^7)
        for Z in _rref_subspaces(7, k):
            Zv = (np.array(Z, dtype=np.int64) @ comp) % P if k else np.zeros((0, 10), dtype=np.int64)
            gens = np.concatenate([Lmat, Zv]) if k else Lmat
            coeffs = np.array(list(itertools.product(range(P), repeat=3 + k)), dtype=np.int64)
            W = (coeffs @ gens) % P
            D = W[decset[W @ pw]]
            if len(D) and rank_p(D) == 3 + k:
                found = True
                break
        print(f"     wedge rank <= {3 + k}: {'yes' if found else 'no'}")
        if found:
            result = 3 + k
            break
    return result


def _rref_subspaces(n, k):
    if k == 0:
        yield []
        return
    for pivots in itertools.combinations(range(n), k):
        free = [(r, c) for r in range(k) for c in range(n) if c > pivots[r] and c not in pivots]
        for vals in itertools.product(range(P), repeat=len(free)):
            M = [[0] * n for _ in range(k)]
            for r, pc in enumerate(pivots):
                M[r][pc] = 1
            for (r, c), v in zip(free, vals):
                M[r][c] = v
            yield M


# ---- 6: splits ------------------------------------------------------------------------------

def check_splits(wedge):
    S = [s for s in s3_splits(22) if 2 * s[0] + s[1] >= 3]
    check("60 rank-22 splits with 2a + b >= 3", len(S) == 60, f"{len(S)}")
    k_a = [s for s in S if s[0] == 0]
    k_kks = [s for s in S if s[0] >= 1 and 3 * s[0] + s[1] + s[2] < 6]
    k_w = [s for s in S if s[0] >= 1 and 3 * s[0] + s[1] + s[2] >= 6 and wedge and 3 * s[0] + s[1] < wedge]
    left = [s for s in S if s not in k_a + k_kks + k_w]
    print(f"     killed by a >= 1 (KKK): {len(k_a)}: {k_a}")
    print(f"     killed by 3a + b + c >= 6 (KKS), among a >= 1: {len(k_kks)}: {k_kks}")
    print(f"     killed by 3a + b >= {wedge} (SSK wedge rank), among the rest: {len(k_w)}: {k_w}")
    print(f"     left: {len(left)}: {left}")
    return left


# ---- 7: orthogonal action on skew parts (symmetry breaking for the size-6 orbit) ---------

def check_orthogonal_skew(schemes):
    """g_M with M orthogonal up to sign commutes with sigma and tau. It acts on skew parts by k -> det(M) M k.
    Its orbits on nonzero k in F3^3 are the three norm classes k.k = 0, 1, 2."""
    Os = []
    for flat in itertools.product(range(P), repeat=9):
        M = np.array(flat).reshape(3, 3)
        if ((M @ M.T) % P == np.eye(3, dtype=np.int64)).all() or ((M @ M.T) % P == (2 * np.eye(3, dtype=np.int64))).all():
            Os.append(M)
    rng = random.Random(140002)
    J = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]]) % P
    ok = True
    for M in Os:
        Mi = inv_p(M)
        d = round(np.linalg.det(M)) % P
        for _ in range(5):
            U = rnd(rng)
            img = sk((M @ U @ Mi) % P)
            # with K = [k]_x (cross-product coordinates), M K M^T = [det(M) M k]_x for M M^T = I;
            # sk(K) = (K12, K13, K23) = J k with J k = (-k3, k2, -k1)
            if not any(((img - c * (J @ M @ inv_p(J) @ sk(U))) % P == 0).all() for c in (1, 2)):
                ok = False
    check("g_M (M orthogonal up to sign) acts on skew parts by k -> +-M k", ok, f"{len(Os)} matrices M (group of order {len(Os) // 2} up to sign)")
    vecs = [np.array(v) for v in itertools.product(range(P), repeat=3) if any(v)]
    seen, orbs = set(), []
    for v in vecs:
        if tuple(v) in seen:
            continue
        orb = {tuple((s * (M @ v)) % P) for M in Os for s in (1, 2)}
        seen |= orb
        orbs.append((int(v @ v) % P, len(orb)))
    check("orbits on nonzero skew parts are the norm classes", sorted(orbs) == [(0, 8), (1, 6), (2, 12)], f"(norm, size): {sorted(orbs)}")
    for name, sch in schemes:
        six = s3_parts(sch)[0]
        print(f"     {name}: size-6 base norms (sk u, sk v, sk w).(same) = "
              f"{[tuple(int(sk(m) @ sk(m)) % P for m in t) for t in six]}")


def main(argv):
    kappa, c0 = check_blocks()
    s006 = Scheme.load(ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json")
    sp006 = check_scheme("006 S3 scheme", s006, kappa, c0)
    lad3 = check_laderman(3)
    if lad3 is not None:
        spl = check_scheme("Laderman in standard S3 form", lad3, kappa, c0)
        check("Laderman S3 form is inequivalent to the 006 S3 scheme (class keys differ)",
              class_key(lad3, P) != class_key(s006, P), f"splits {sp006} vs {spl}")
    if lad3 is not None and "--save-laderman" in argv:
        lad3.meta = {"source": "Laderman 1976 (Burichenko 2015 tensor form) conjugated by T(I, pi12, eps1)^-1 then g_N, "
                     "N = [[0,1,1],[0,1,2],[1,0,0]]", "field": "GF(3)", "s3_split": [1, 4, 1, 3]}
        lad3.save(argv[argv.index("--save-laderman") + 1])
    check_laderman(2)
    named = [("006", s006)] + ([("Laderman", lad3)] if lad3 is not None else [])
    seen = {class_key(s006, P)} | ({class_key(lad3, P)} if lad3 is not None else set())
    for f in sorted((ROOT / "results" / "f3").glob("010*.json")):
        try:
            sch = Scheme.load(f)
        except Exception:
            continue
        if sch.rank > 23 or not verify(sch, mod=P).ok or not is_s3_invariant(sch):
            continue
        k = class_key(sch, P)
        if k in seen:
            continue
        seen.add(k)
        check_scheme(f.name, sch, kappa, c0)
        named.append((f.name, sch))
    print(f"     distinct S3-invariant rank-23 classes checked: {len(seen)}")
    check_orthogonal_skew(named)
    wedge = ssk_wedge_rank(6)
    print(f"     SSK wedge rank = {wedge if wedge else '> max tried'}")
    check_splits(wedge if wedge else 7)
    print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
