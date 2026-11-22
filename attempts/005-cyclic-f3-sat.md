---
id: 005-cyclic-f3-sat
date: 2026-11-22
target: 3x3x3 rank 22 over F3, restricted to schemes invariant under the cyclic symmetry
method: SAT (CaDiCaL 3.1.0, Kissat 4.0.4) on the Brent equations mod 3 under the cyclic ansatz; UNSAT proofs checked with drat-trim
status: partial
best_rank_found: 23 over F3 (cyclic-invariant, found by SAT from scratch); no 22; split (1,19) proved impossible
budget: 2 processes × 1 thread; stage 0 ≤ 60 s per instance; stage A 2 × 20 min; stage B 7 instances × 10 min (≤ 110 CPU-minutes in all)
---

# 005: cyclic-invariant schemes over F3, by SAT

## Why
002–004 show that plain flip-graph walks over F3 stall at 23, both when they run long (002) and when they run wide (004).
The F3 rank-23 landscape is split into hundreds or more pockets (003, 004), so random walks can't cover it.
This attempt restricts the search to schemes that are invariant under the cyclic symmetry of the tensor. That space
is about three times smaller and splits into a few finite cases, so a SAT solver can search it *completely*. A complete
search can end in either outcome: a scheme, or a checkable proof that no scheme of that shape exists.

## The symmetric ansatz (worked out before running)
Convention: a term (U, V, W) contributes U[i][j]·V[j][l]·W[i][l]. The cyclic move is σ: (U, V, W) → (V, Wᵀ, Uᵀ),
σ³ = id, and σ maps correct schemes to correct schemes (tested in `tests/test_invariants.py`).

**Orbits.** A σ-invariant scheme is a union of size-3 orbits {(u,v,w), (v,wᵀ,uᵀ), (wᵀ,u,vᵀ)} and fixed terms.
With a orbits and b fixed terms the rank is 3a + b.

**Fixed terms.** If σ(t) = t as a rank-one tensor, then (v, wᵀ, uᵀ) = (αu, βv, γw) with αβγ = 1, so v = αu, w = αβ·uᵀ,
and t = λ·(x ⊗ x ⊗ xᵀ) for x = u. Over F3 every scalar is a cube (μ³ = μ), so λ can be absorbed:
**every fixed term is exactly (x, x, xᵀ) with x a nonzero 3x3 matrix over F3.**
(This matches Ballard et al., whose fixed terms are cubes a⊗a⊗a in their convention.)

**A necessary condition: the fixed matrices sum to I.** Take the trilinear form of the scheme and set all three
arguments to the same matrix X (in the convention where σ permutes factors plainly). The left side is tr(X³). An orbit contributes
3·ℓ_u(X)ℓ_v(X)ℓ_w(X) ≡ 0 mod 3. A fixed term contributes ℓ_x(X)³, and in characteristic 3
(Σ x_e X_e)³ = Σ x_e³ X_e³ = Σ x_e X_e³. In tr(X³) = Σ X_ij X_jk X_ki every monomial other than X_ii³ comes with a
multiple of 3 (its cyclic rotations), so tr(X³) ≡ Σ X_ii³ mod 3. Comparing coefficients of X_e³: **Σ_fixed x = I over F3.**
Checked on Ballard et al.'s Z3-invariant rank-23 scheme (its two fixed matrices sum to I). Consequences:
- b ≥ 1. An all-orbit scheme is impossible over F3.
- If b = 1, the fixed term is forced: (I, I, I).
- a ≥ 1 too, since the tensor is not symmetric, so it can't be a sum of fixed terms alone.

**Splits for rank 22.** 22 = 3a + b with a ≥ 1, b ≥ 1, b ≡ 1 mod 3:
(a, b) ∈ {(7,1), (6,4), (5,7), (4,10), (3,13), (2,16), (1,19)}.
Each slot in the encoding may be zero, so the instance (a, b) also covers every (a', b') with a' ≤ a and b' ≤ b
(ranks ≤ 22). The seven instances together cover **every σ-invariant scheme of rank ≤ 22 over F3**.

**Size.** Unknowns: 27a + 9b = 198 F3 values for *every* rank-22 split (207 at rank 23), or 396 Boolean variables.
Equations: of the 729 Brent equations only 249 are distinct once the ansatz is substituted (computed, and
matches (729 − 9)/3 + 9).

## Method
`src/matmul/cyclic.py` (tested in `tests/test_cyclic.py`):
- each F3 value is 2 Boolean variables (p, n) with ¬(p ∧ n); products are built from shared 2-input product gates;
  each equation is a balanced tree of mod-3 adders whose root is fixed to the target;
- the implied condition Σ fixed = I is added as 9 more sums;
- symmetry breaking that keeps a representative of every solution: in each orbit slot, the first nonzero entry of u and of v is 1
  (rescaling), and orbit slots and fixed slots are each in lexicographic order (reordering).
  Rotation within an orbit and simultaneous conjugation x → MxM⁻¹ are *not* broken (possible improvement).
- every SAT model is decoded, expanded and verified with `matmul.verify --mod 3`, and checked for σ-invariance;
- every UNSAT is re-run with a DRAT proof and checked with drat-trim (an independent checker).

Seeds: solver seeds 50001, 50002, ... (earlier attempts used 1–180 and 9001–9008).

## Plan, budget and decision rule (fixed before running)
Check `uptime` first. At most 2 solver processes, 1 thread each, each under a hard timeout.

- **Stage 0, method check on 2x2 (≤ 60 s each).** 2x2 rank 7 = 3·2 + 1 must be SAT (Strassen's scheme is cyclic-invariant),
  and 2x2 rank 6 (splits (1,3) and (0,6)) must be UNSAT with a proof that drat-trim accepts. If either fails, stop and fix the encoding.
- **Stage A, sanity at 3x3 rank 23 (2 × 20 min).** Split (7, 2), the split of Ballard et al.'s integer scheme, which reduces mod 3 to a valid F3 scheme.
  CaDiCaL and Kissat run side by side from scratch.
  - Found and verified: the encoding can find 3x3 schemes. Use the faster solver for stage B.
  - Not found in 20 min: finding power at 3x3 is unproven. Stage B still runs, because UNSAT for small a might come fast,
    but a timeout in stage B then means nothing for finding.
- **Stage B, rank 22 (7 instances × 10 min, 2 at a time, smallest a first).**
  - **SAT:** verify mod 3 and check σ-invariance. If it holds, this is a rank-22 scheme over F3. Stop, save it outside `schemes/`
    (that folder holds schemes valid over Q), check the catalogue, then try lifting.
  - **UNSAT with a drat-trim-checked proof:** a theorem: *no σ-invariant scheme over F3 with at most a orbits and at most b fixed terms*.
    If all seven close, there is no σ-invariant rank-≤22 scheme over F3 at all.
  - **Timeout:** no conclusion for that split.
- Expectation, stated in advance: Yang (2024) could not get past rank 15 for cyclic-only symmetry over GF(2) with SAT in
  more than 900,000 s. So a full UNSAT at rank 22 is unlikely in this budget. The cases most likely to close are those with few orbits.

## Result
Load before each stage: 1–3 (`uptime`). Every instance is one line in `results/runs.jsonl` (attempt `005-cyclic-f3-sat`,
fields `stage` and `result`). Found schemes are in `results/f3/005_*.json`, and every one passes `matmul.verify --mod 3` and is cyclic-invariant.
Commands: `python3 -I scripts/run_005.py 0`, then `A`, then `B cadical`.

**Stage 0 (2x2, method check): passed.**

| Instance | CaDiCaL | Kissat |
|---|---|---|
| rank ≤ 7, split (2,1) | SAT, 0.08 s, rank 7, verified | SAT, 0.03 s, rank 7, verified |
| rank ≤ 6, split (1,3) | UNSAT, proof VERIFIED by drat-trim | UNSAT, proof VERIFIED |
| rank ≤ 6, split (0,6) | UNSAT (trivially: two identical equations with different targets), proof VERIFIED | same |

**Stage A (3x3 rank 23, split (7,2), from scratch): passed.**

| Solver | Time | Result |
|---|---|---|
| CaDiCaL, seed 50002 | 257 s | rank 23, 7 orbits + 2 fixed, verified mod 3, cyclic-invariant |
| Kissat, seed 50002 | 907 s | rank 23, 7 orbits + 2 fixed, verified mod 3, cyclic-invariant |

The two schemes are inequivalent to each other and to Ballard et al.'s scheme (`matmul.invariants.class_key`: b2c74bf5…, a04dac5d… vs 094008dc…).
The encoding can find 3x3 schemes. By the rule, stage B used CaDiCaL (the faster solver).

**Stage B (3x3 rank ≤ 22, CaDiCaL seed 50003, 10 min per split).** 21,084–22,932 variables and 134k–146k clauses per instance.

| Split (a orbits, b fixed) | Result |
|---|---|
| (1, 19) | **UNSAT in 1.1 s, proof (3.3 MB DRAT) VERIFIED by drat-trim** |
| (2, 16) | timeout (600 s) |
| (3, 13) | timeout |
| (4, 10) | timeout |
| (5, 7) | timeout |
| (6, 4) | timeout |
| (7, 1) | timeout |

**What is proved (and under which assumptions).** Over F3 there is no scheme for 3x3x3 of rank ≤ 22 that is invariant
under σ and has at most one orbit of size 3 (and at most 19 fixed terms). This includes every all-fixed scheme.
The proof rests on:
1. the reduction from the claim to the CNF: fixed terms have the form (x, x, xᵀ) over F3, the fixed matrices sum to I
   (both derived above), and the symmetry breaking keeps a representative of every solution (rescaling and reordering).
   The CNF generator is tested: a known scheme satisfies it with and without symmetry breaking, a corrupted one doesn't,
   and the 249 equations match the count we derived. But the generator itself is not formally verified;
2. drat-trim's check of CaDiCaL's proof against that CNF file.
It is a small, nearly easy case (one orbit can carry little of the tensor), and the only rank-22 split we closed.
**The other six splits are open**, so this attempt does not decide whether a cyclic-invariant rank-22 scheme over F3 exists.

**By the rule fixed beforehand:** SAT never happened. (1,19) is UNSAT with a checked proof, a theorem for that split only.
The six timeouts give no conclusion. Stage A succeeded, so those timeouts at least come from a solver that *can* find
3x3 schemes of this shape at rank 23 in minutes. That is weak evidence (not proof) that a rank-22 scheme in these splits
is not easy to find.

Compute: stage 0 under 2 s; stage A 1,164 CPU-s; stage B 3,601 CPU-s; about 80 CPU-minutes in all, 2 processes at most.

## What went wrong / surprised us
- The first stage-0 run crashed on split (0,6). The generator treated "same polynomial, two targets" as a bug, but with a = 0
  it is a genuine contradiction: the tensor is not a sum of fixed terms. Now it emits an empty clause. The crashed run had already
  written 3 lines to `results/runs.jsonl`, which duplicate the first three stage-0 lines.
- Finding rank 23 took minutes, not hours. Splits with 2+ orbits at rank 22 didn't close in 10 minutes, and we don't know
  how far each is from closing.
- Stage A had to be run first. Without it, the six timeouts would have meant nothing about finding power.

## What we would do better
- **Break more symmetry.** Simultaneous conjugation x → MxM⁻¹ (|PGL₃(F3)| = 5,616) and the choice of rotation inside each orbit
  are not broken. For b = 1 the fixed term is forced to be I (conjugation-invariant). For b ≥ 2 the first fixed matrix could be fixed to
  a conjugacy-class representative, which gives a case split that also parallelises well.
- **Cube and conquer:** split each open instance into many small cubes (for example by the rank or support of the first orbit's u),
  so progress can be measured and the proof comes in pieces.
- Try **local search** (Heule et al. found it beat CDCL on finding 3x3 schemes over GF(2)) for the finding side, and a
  cyclic flip-graph walk seeded from the two new rank-23 schemes.
- Calibrate first: time several rank-23 splits ((7,2), (6,5), (4,11)) and several seeds, to know what a "normal" time-to-SAT is
  before choosing the rank-22 cap.

## Next
006 candidates, best first:
1. The split (2,16) is the next smallest. Add conjugation symmetry breaking and run it longer; it may close, extending the theorem to "≤ 2 orbits".
2. Split (7,1) is the closest relative of the rank-23 split (7,2) that SAT solves in minutes, and its fixed term is forced to be I. A longer or
   cube-and-conquer run there is the best bet for *finding* a 22.
3. Cyclic flip graph over F3 from the stage-A schemes (finding side).
