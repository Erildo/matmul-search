# Cyclic-invariant 3×3 matrix multiplication over F3: what is and isn't possible at rank 22

*matmul-search, attempts 005–009 (November 2026). Every claim below links to the test, log line or attempt note it rests on.*

## Summary
We searched for rank-22 schemes for 3×3 matrix multiplication over F3 that are invariant under the cyclic symmetry of the matrix multiplication tensor.
We found none. Two kinds of rank 22 are now ruled out, each with machine-checked proofs from two independent SAT solvers:
- **T1.** There is no cyclic-invariant scheme of rank ≤ 22 with at most 2 orbits.
- **T2.** There is no cyclic-invariant scheme of rank 22 (split (7,1) or (6,4)) that shares 5 or more of its 7 orbits with any of four known cyclic rank-23 schemes.

Along the way we found **10 new cyclic-invariant rank-23 schemes over F3**. Together with the one of Ballard, Ikenmeyer, Landsberg and Ryder, that makes **11 pairwise
inequivalent** schemes, released with a check script.

All 11 have the same shape:
- every one has split (7,2);
- exactly one orbit has all three matrices of rank 2;
- the two fixed matrices have ranks 1 and 2;
- no matrix anywhere has rank 3.

We also report where the method stops: the next cases are out of reach at the caps we used (section 7).

## 1. Statements
Convention: a term (U, V, W) contributes U[i][j]·V[j][l]·W[i][l], and the cyclic symmetry is σ: (U, V, W) → (V, Wᵀ, Uᵀ).
A scheme is **cyclic-invariant** if σ permutes its terms. It is then a union of **orbits** {t, σt, σ²t} of size 3 and **fixed terms**, and its rank is 3a + b,
with a orbits and b fixed terms (the **split** (a, b)). Two schemes **share an orbit** if both contain the same three rank-one tensors.

- **T1.** Over F3, no cyclic-invariant scheme for ⟨3,3,3⟩ of rank ≤ 22 has a ≤ 2.
- **T2.** Let S be one of BILR, 005-cadical, 005-kissat or 006 (section 6). Over F3, no cyclic-invariant scheme for ⟨3,3,3⟩ of rank 22 with split (7,1) or (6,4)
  shares ≥ 5 orbits with S.
- **D1.** The 11 schemes in [`results/f3/cyclic_rank23.jsonl`](../results/f3/cyclic_rank23.jsonl) are correct over F3, cyclic-invariant, of rank 23, and pairwise inequivalent.

T1 and T2 are proofs. Everything in section 7 is evidence only.

## 2. Two facts about fixed terms over F3
**F1. A fixed term is (x, x, xᵀ).** If σt = t for t = u ⊗ v ⊗ w, then v = αu and w = βuᵀ, so t = λ·(x, x, xᵀ). Every element of F3 is a cube, so λ can be absorbed into x.

**F2. The fixed matrices sum to I.** Put the same matrix X into all three arguments of the scheme's trilinear form. The product side is tr(X³). In characteristic 3, every monomial of tr(X³)
except Σ X_ii³ occurs a multiple of 3 times. An orbit contributes 3·(…) ≡ 0. A fixed term contributes (Σ x_e X_e)³ = Σ x_e X_e³. Comparing coefficients gives Σ x = I.
So b ≥ 1. If b = 1, the fixed term is (I, I, Iᵀ). And a ≥ 1, because the tensor is not a sum of fixed terms ([`trivially_unsat`](../src/matmul/cyclic.py) for (0, 22)).
Details in [attempt 005](../attempts/005-cyclic-f3-sat.md).

## 3. The encoding
`CyclicF3CNF(n, a, b)` ([`src/matmul/cyclic.py`](../src/matmul/cyclic.py)):
- Variables: each orbit slot holds (u, v, w) and each fixed slot holds x, with each F3 entry stored as two Boolean variables.
- Equations: the Brent equations of the expanded scheme, with products and sums built as small gate networks.
- Any slot may be zero, so an instance (a, b) covers every split (a′, b′) with a′ ≤ a and b′ ≤ b.

**Checked, not assumed.** The test `test_equations_match_the_expanded_tensor` covers every split used here and checks two things:
- the encoder's equations account for all n⁶ index tuples;
- for random slot contents, each equation's value equals the tensor entry computed independently from the expanded list of terms, and the clauses are satisfied exactly when every equation holds.

**Symmetry breaking.**
- Each orbit is rescaled so that the first nonzero entries of u and v are 1, using (u, v, w) → (εu, δv, εδw).
- Orbit slots and fixed slots are each sorted lexicographically.

Both moves keep the set of tensors, so every scheme has a representative that satisfies them (`normal_form`, tested on every known scheme).

## 4. T1: class cubes
- (0, 22) is trivially UNSAT (F2).
- (1, 19) is UNSAT as one instance.
- (2, 16) is split into 21 cubes.

**Classes.** The nonzero 3×3 matrices over F3 fall into 21 classes up to conjugation and sign. Order them by their smallest member. Classes 0 and 1 are the rank-1 matrices with nonzero and zero trace.

**Cube c.** Slot 0's u is pinned to the representative of class c. No orbit matrix is in a class below c, enforced in three ways:
- for c ≥ 1, no rank-1 matrix with nonzero trace;
- for c ≥ 2, no rank-1 matrix at all;
- every characteristic polynomial carried only by classes below c is excluded, counting the zero matrix as allowed.

**Lemma (coverage).** Every cyclic-invariant scheme with a nonzero orbit has a symmetric image in cube c, where c is the smallest class among the first matrices u, v, wᵀ of its orbit members.

*Proof.* Rotate the minimizing member into slot 0. Rescale by (−u, −v, w) if needed. Conjugate the whole scheme by g_M: (u, v, w) → (MuM⁻¹, MvM⁻¹, M⁻ᵀwMᵀ). This keeps the tensor, commutes with σ, and maps fixed terms to fixed terms. Choose M so that u becomes the representative of class c.
Classes are invariant under conjugation, sign and transposition, so every orbit matrix keeps a class ≥ c. Then normalise the other slots.
Full proof and tests in [attempt 006](../attempts/006-cyclic-f3-cubes.md). ∎

All 21 cubes, (1,19) and (0,22) are UNSAT, so T1 follows.

## 5. T2: neighbourhoods of the known schemes
For S with orbits S₁..S₇ and a 5-subset K of {1..7}, the **frozen cube** (S, split, K) fixes orbit slots 0–4 to S's stored representatives of the orbits in K.
The other orbit slots and the fixed slots are free, under their usual symmetry breaking.

**Coverage.** Suppose a scheme shares the orbits in K with S. Rotation and (u, v, w) → (εu, δv, εδw) keep each orbit's tensors, so each shared orbit can be written as S's stored representative.
Put the shared orbits in slots 0–4 and put the rest in normal form. That gives a solution of the cube.
There are 2 splits × C(7, 5) = 21 subsets × 4 seeds, so **168 cubes, all UNSAT**. Sharing ≥ 6 orbits is included in sharing ≥ 5.

## 6. The 11 cyclic rank-23 schemes (D1)
| class key | found by | orbits whose three matrices have ranks (1,1,1) / (1,1,2) / (2,2,2) | fixed-matrix ranks |
|---|---|---|---|
| `094008dc304fad88` | BILR S_2fix-Z3, reduced mod 3 | 5 / 1 / 1 | 1, 2 |
| `b2c74bf59a566dc7` | 005, SAT (CaDiCaL), split (7,2) from scratch | 4 / 2 / 1 | 1, 2 |
| `a04dac5d40238167` | 005, SAT (Kissat), split (7,2) from scratch | 4 / 2 / 1 | 1, 2 |
| `eebc48753638175a` | 006, SAT, class cube 0 | 5 / 1 / 1 | 1, 2 |
| `ac0136d1e3ad372d` | 008, walk from 005-cadical, k = 2 | 4 / 2 / 1 | 1, 2 |
| `6d732e700b4ee825` | 008, walk from 005-cadical, k = 3 | 4 / 2 / 1 | 1, 2 |
| `ad5563bc5f3b2158` | 008, walk from 005-kissat, k = 2 | 3 / 3 / 1 | 1, 2 |
| `346a978816ebcbb6` | 008, walk from 005-kissat, k = 3 | 5 / 1 / 1 | 1, 2 |
| `6f888e94663cd913` | 008, walk from 006, k = 2 | 5 / 1 / 1 | 1, 2 |
| `ec8ce52996e9a191` | 008, walk from 006, k = 3 | 5 / 1 / 1 | 1, 2 |
| `3f98b9cc15d7710b` | 008, walk from BILR, k = 2 | 5 / 1 / 1 | 1, 2 |

- "Class key" is the invariant fingerprint of [`matmul.invariants`](../src/matmul/invariants.py). Different keys prove inequivalence.
- None of the 11 keys is among the 51,555 F3 rank-23 classes we saved from unrestricted flip-graph walks (attempt 003).
- A **walk step** at distance k asks SAT for a cyclic rank-23 scheme that shares ≥ 7 − k orbits with S and differs from it in at least one ([attempt 008](../attempts/008-cyclic-f3-lns.md)).

**Observation (not proved).** Every known cyclic F3 rank-23 scheme has these properties:
- split (7,2);
- one orbit with all three matrices of rank 2;
- fixed matrices of ranks 1 and 2;
- no matrix of rank 3.

A rank-22 scheme of split (7,1) would need I (rank 3) as its fixed term, which breaks this pattern.

## 7. Where it stops (evidence, not proof)
| question | result | cap |
|---|---|---|
| (3,13): 3 orbits, 21 class cubes | 0/21 closed | 120 s each ([007](../attempts/007-cyclic-f3-three-orbits.md)) |
| (7,1): 7 orbits, 21 class cubes | 2/21 closed | 120 s each ([006](../attempts/006-cyclic-f3-cubes.md)) |
| distance 2 as one SAT instance | TIMEOUT, though its 21 frozen cubes close in about 1 min in total | 60 s (008) |
| distance 3 from 006, split (7,1), frozen cubes | 0/31 closed (4 not run) | 20 s each (008) |
| cyclic flip moves available in the known schemes | 1–3 per scheme; flips never change the number of fixed terms | – (008) |

Each wall is sudden. At 2 orbits, (2,16)'s 21 cubes all closed (19 in 1–115 s, the last two once the polynomial exclusion was added), and at 3 orbits none did. Every distance-2 cube closed in under 13 s, and at distance 3 none closed in 20 s.
Rank 22 over F3 without symmetry is open. Over F2 it is now known to be impossible (Tomar 2026, below).

## 8. Re-checking
```bash
python3 -I scripts/recheck_cyclic_f3.py --quick    # T1, CaDiCaL only: about 12 CPU-min
python3 -I scripts/recheck_cyclic_f3.py            # T1 and T2, CaDiCaL and Kissat, every proof checked by drat-trim
python3 -I scripts/check_cyclic_rank23.py check    # D1
./run_tests.sh                                     # includes the encoder test and the lemma tests
```
The re-check rebuilds every CNF from the code and records each one's SHA-256 in [`results/f3/recheck_009.jsonl`](../results/f3/recheck_009.jsonl).
It needs `cadical`, `kissat` and drat-trim (`vendor/drat-trim`; build instructions in attempt 006).

**Re-check result:** RECHECK_RESULT

**Trusted base.** What is proved:
- each CNF is unsatisfiable (two solvers, each proof checked by drat-trim);
- each CNF encodes the Brent equations of the expanded scheme (the complete equation test above).

What is argued in the text and tested:
- F1 and F2;
- soundness of the symmetry breaking;
- the coverage lemma for class cubes and the class/polynomial tables (`tests/test_cyclic.py`: the class tests, `test_charpoly_exclusion_gadget`, the minimal-class-form tests);
- coverage of the frozen cubes.

What is trusted: drat-trim, and Python's arithmetic.

## 9. Related work
- **Ballard, Ikenmeyer, Landsberg, Ryder** ([arXiv:1801.00843](https://arxiv.org/abs/1801.00843), JPAA 2019): the first cyclic-invariant rank-23 decompositions,
  with integer coefficients. Their S_2fix-Z3, reduced mod 3, is our BILR.
- **Yang** ([arXiv:2402.01011](https://arxiv.org/abs/2402.01011)): SAT under cyclic symmetry over F2. It rules out rank ≤ 15 with the cyclic symmetry alone, and ≤ 21 with more symmetry.
- **Moosbauer & Poole** ([arXiv:2502.04514](https://arxiv.org/abs/2502.04514)): flip graphs with symmetry over F2, with Hensel lifting.
- **Perminov** ([arXiv:2511.20317](https://arxiv.org/abs/2511.20317)): ternary-coefficient flip graphs, with no symmetry constraint.
- **Rudich and Rousseau** ([arXiv:2610.01639](https://arxiv.org/abs/2610.01639)): rank ≥ 22 over Z and F2. **Tomar** ([arXiv:2610.11823](https://arxiv.org/abs/2610.11823)): rank = 23 over F2.

As of 2026-11-23 we found no prior search for cyclic-invariant schemes over F3 (sources checked in [PRIOR_WORK.md](../PRIOR_WORK.md) and in attempt 009).
We claim only that: we found none.
