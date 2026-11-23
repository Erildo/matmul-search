---
id: 009-cyclic-f3-writeup
date: 2026-11-23
target: a checkable write-up of what attempts 005–008 established about cyclic-invariant 3x3x3 schemes over F3
method: independent re-check of every proved statement (two solvers, two checked proofs, CNFs regenerated from code), an encoder-vs-tensor test, a data release of the cyclic rank-23 schemes, and a self-contained document with every claim tied to its evidence
status: running
best_rank_found: – (no search)
budget: 2 processes × 1 thread; re-check ≤ 90 CPU-min; no new searches
---

# 009: write-up of the cyclic F3 results

## Why
Attempts 005–008 produced results that are spread over four notes, several scripts and `runs.jsonl`. Some of those results are proved, and some are only evidence.
They are only worth something if a reader can (1) see exactly what is claimed, (2) see what each claim rests on, and (3) re-check it with one command.
Three things stand in the way today:
- **The DRAT proofs are gone.** They were deleted after drat-trim verified them, as designed, because they are large. The CNFs live only in git-ignored `scratch/`.
- **Each UNSAT has one anchor.** Each one rests on one CaDiCaL run and one drat-trim check.
- **The trusted base is implicit.** A drat-trim check proves that the *CNF* is unsatisfiable, but the mathematical claim also needs the encoder and the coverage lemmas to be right.

## What is claimed (the write-up's spine)
Definition: two cyclic schemes **share an orbit** if both contain the same three rank-one tensors {t, σt, σ²t}.

Proved statements, each with checked proofs:
- **T1 (005, 006).** No cyclic-invariant 3x3x3 scheme over F3 of rank ≤ 22 has ≤ 2 orbits. It rests on:
  - a = 0 is impossible;
  - instance (a, b) covers every split with fewer slots, because a slot may be zero;
  - (1,19) is UNSAT;
  - all 21 class cubes of (2,16) are UNSAT, which needs the cube-coverage lemma and the soundness of the char-poly exclusion.
- **T2 (008).** For each of the four seed schemes S (BILR, the two 005 schemes and the 006 scheme), no cyclic-invariant F3 scheme of rank 22 with split (7,1) or (6,4)
  shares ≥ 5 orbits with S. The proof is the 168 frozen-subset cubes alone (see revision R1).
- **Supporting facts.** Over F3 a fixed term is (x, x, xᵀ), and the fixed matrices sum to I.

Data:
- **D1.** 11 pairwise inequivalent cyclic-invariant rank-23 schemes over F3, all of split (7,2): BILR, plus 10 found by us in 005, 006 and 008.
  None of the 10 is among the 51,555 saved F3 rank-23 representatives.

Evidence only, not proof, to be reported as such:
- **E1.** The walls, each measured at a stated cap:
  - (3,13): 0/21 cubes at 120 s;
  - (7,1): 2/21 cubes at 120 s;
  - distance 2 as one near() instance: TIMEOUT at 60 s;
  - distance 3 on the 006 seed: 31/31 TIMEOUT at 20 s.
- **E2.** The cyclic flip neighbourhood of the known schemes is nearly empty: 1–3 shared-factor pairs. Orbit flips also cannot change the number of fixed terms.

## Revision before running (2026-11-23)
- **R1. T2 needs only the k = 2 frozen cubes.** "Shares ≥ 5 orbits" includes "shares ≥ 6", so the 8 near() instances at k = 1 add nothing.
  That takes `near()` out of T2's trusted base: it stays a walk tool only. The coverage argument is then a single one:
  if a scheme shares orbits O₁..O₅ with S, the rotation and the three sign moves (u, v, w) → (εu, δv, εδw) with ε, δ = ±1 keep each orbit's three tensors.
  So each shared orbit can be written as S's stored representative and moved to the frozen slots 0–4. The free orbit slots and the fixed slots then still admit
  the normal form, which is all their symmetry breaking asks for. The frozen slots are outside the lex chain.
- **R2. Run the compute first.** The re-check is pure compute and needs nothing from stages 0, 2 or 3. It starts first, in the background, and the rest is done while it runs.
- **R3. Kissat proofs are binary DRAT.** drat-trim detects the format by itself; checked on two toy CNFs. No conversion needed.
- **R4. Make the equation test complete, not sampled.**
  - The encoder records, for each equation it writes, *every* index tuple (i, j, j′, l, i′, l′) that produced it. Today the de-duplication discards them.
  - The test checks two things:
    - the index lists partition all n⁶ tuples;
    - for random slot contents, each equation's value (computed from the decoded matrices, the way the CNF's products and sums do) equals
      Σ U[i][j] V[j′][l] W[i′][l′] − δ at every one of its tuples, computed independently from `expand(...)`'s terms.

  That shows that the CNF's equations are exactly the Brent equations of the expanded scheme. It covers every split used in T1 and T2.

## Plan
- **Stage 1: one-command re-check (≤ 90 CPU-min, started first).**
  - `scripts/recheck_cyclic_f3.py` regenerates every CNF behind T1 and T2 from the encoder, with no files from `scratch/`. For each one it:
    - solves with CaDiCaL and checks the proof with drat-trim;
    - solves with Kissat and checks that proof with drat-trim;
    - writes the CNF's SHA-256, both statuses, both verdicts and the times to `results/f3/recheck_009.jsonl`.

    Proofs are deleted after the check, as before.
  - Instances:
    - **T1:** the a = 0 trivial check; (1,19); (2,16) cubes 0–20 with the char-poly exclusion on in all of them. This is a uniform configuration, valid under the lemma.
    - **T2:** 168 frozen cubes.
  - Caps: solve 600 s, check 900 s. Expected cost: about 12 CPU-min for T1 and about 14 for T2 with CaDiCaL including proofs, plus Kissat at an unknown factor.
    Order: T1 with CaDiCaL, T1 with Kissat, T2 with CaDiCaL, then T2 with Kissat until the budget is spent.
  - `--quick` for readers: T1 only, CaDiCaL only.
- **Stage 2: encoder check (R4).** Add the index lists to the encoder and the test. Then list the remaining trusted base with a pointer for each:
  - symmetry-breaking soundness;
  - the cube-coverage lemma;
  - the char-poly table;
  - freeze coverage (R1);
  - drat-trim.
- **Stage 0: literature freshness (while stage 1 runs).**
  - Re-check the cyclic section of PRIOR_WORK.md against the Sedoglavic catalogue, and against arXiv for 3x3 cyclic/Z3/symmetric decompositions over F3 or GF(3), and work citing Yang 2024 or Tomar 2026.
  - Outcome: each novelty claim ("first Z3-invariant search over F3", "first F3 cyclic classes beyond BILR") is kept, softened or dropped.
- **Stage 3: data release.**
  - `results/f3/cyclic_rank23.jsonl`: one line per class, giving orbits, fixed terms, source, class key, and the number of orbits containing a rank-1 matrix.
  - `scripts/check_cyclic_rank23.py` verifies each scheme mod 3, checks cyclic invariance, recomputes class keys, checks they are pairwise distinct and not in the saved reps,
    and prints a summary. A test runs it.
- **Stage 4: the document `writeups/cyclic-f3.md`**, with the outline as before:
  1. results;
  2. setup, with the two F3 facts proved;
  3. the encoding;
  4. class cubes and the lemma;
  5. the local theorem (R1);
  6. the census;
  7. where it stops (E1, E2);
  8. re-checking and the trusted base;
  9. related work.

  Every claim links to its evidence: a test, a line of `recheck_009.jsonl`, or an attempt note.
- **Stage 5: review and release.**
  - Read the document cold, starting with its weakest link.
  - Update LOG.md (007, 008 and 009 are missing), TARGETS.md and README.md.
  - Commit to `main` and push. As of 2026-11-23 there are no branches: everything goes straight to `main`, and each run's log is kept in `results/logs/<attempt>/`.
  - Offer a shareable web page as a separate step.

## Decision rule (fixed before running)
- **Any re-check disagreement:** a SAT answer, a TIMEOUT, or a failed proof in either checker. Stop. Withdraw the affected claim until the cause is found, and record it here.
- **The equation test fails:** every UNSAT claim is suspended until the encoder is fixed and stage 1 is re-run.
- **Stage 0 finds prior F3 cyclic work:** rewrite the novelty claims, and compare our 11 classes with any data it has.
- **Kissat over budget:** T1 with Kissat in full, T2 with Kissat on the cubes reached, and say how many in the document.

## Out of scope
New searches. Gaps the write-up exposes go under Next.

## Stage results so far
- **Stage 0 (literature, 2026-11-23).** Two web searches:
  - (1) cyclic/Z3-invariant 3x3 decompositions over GF(3) with SAT;
  - (2) 2025–26 arXiv work on 3x3 rank 22 with symmetry or flip graphs in characteristic 3.

  Neither found a cyclic-invariant search over F3. New or confirmed items: Perminov's ternary meta flip graphs ([arXiv:2511.20317](https://arxiv.org/abs/2511.20317)) work with
  {−1, 0, 1} coefficients and sign symmetry only, with no Z3 constraint. Rudich and Rousseau ([arXiv:2610.01639](https://arxiv.org/abs/2610.01639)) prove ≥ 22 over Z and F2, and are already in PRIOR_WORK.
  **Novelty claims softened** to "as of 2026-11-23 we found no prior cyclic-invariant search over F3". No claim of being first.
- **Stage 2.** The equation-map test passes for (2,2,1), (2,1,3), (3,1,19), (3,2,16), (3,7,1) and (3,6,4). Every equation is checked at every one of its index tuples against the tensor
  computed from `expand`, and the gate network agrees with the equations. Suite: 137 passed.
- **Stage 3.** `results/f3/cyclic_rank23.jsonl` holds 11 schemes, and `scripts/check_cyclic_rank23.py check` reports 11 distinct classes and no failures. **None of the 11, BILR included, is among
  the 51,555 saved F3 rank-23 classes**, so the cyclic schemes are not in the part of the F3 locus the unrestricted walks visited.
  The invariants turned up a pattern shared by all 11 (recorded as an observation in the write-up):
  - split (7,2);
  - one orbit whose three matrices all have rank 2;
  - the other six orbits have rank profile (1,1,1) or (1,1,2);
  - fixed matrices of ranks 1 and 2;
  - no rank-3 matrix anywhere.

  Split (7,1) needs I as its fixed term, which breaks this pattern.
- **Stage 4.** Drafted in `writeups/cyclic-f3.md`. The re-check result is filled in when stage 1 ends.
