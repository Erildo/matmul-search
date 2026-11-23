---
id: 009-cyclic-f3-writeup
date: 2026-11-23
target: a checkable write-up of what attempts 005–008 established about cyclic-invariant 3x3x3 schemes over F3
method: independent re-check of every proved statement (two solvers, two proof checks, regenerated CNFs), a data release of the cyclic rank-23 schemes, and a self-contained document with every claim tied to its evidence
status: planned
best_rank_found: – (no search)
budget: 2 processes × 1 thread; re-check ≤ 90 CPU-min; no new searches
---

# 009: write-up of the cyclic F3 results

## Why
Attempts 005–008 produced results that are spread over four notes, several scripts and `runs.jsonl`. Some of those results are proved, and some are only evidence.
They are only worth something if a reader can (1) see exactly what is claimed, (2) see what each claim rests on, and (3) re-check it with one command.
Three things stand in the way today:
- **The DRAT proofs are gone.** They were deleted after drat-trim verified them, as designed, because they are large. The CNFs live only in git-ignored `scratch/`.
  A reader has to regenerate the instances from the code, and no single script does that today.
- **Each UNSAT has one anchor.** Each one rests on one CaDiCaL run and one drat-trim check. A second solver's checked proof would remove any doubt about the solver or the proof pipeline.
- **The trusted base is implicit.** A drat-trim check proves that the *CNF* is unsatisfiable, but the mathematical claim also needs the encoder to be right and the coverage lemmas to be right.
  The write-up must name every link in that chain and point to its proof or test.

## What is claimed (the write-up's spine)
Proved statements, each with checked proofs:
- **T1 (005, 006).** No cyclic-invariant 3x3x3 scheme over F3 of rank ≤ 22 has ≤ 2 orbits. It rests on:
  - a = 0 is impossible;
  - instance (a, b) covers every split with fewer slots;
  - (1,19) is UNSAT;
  - all 21 class cubes of (2,16) are UNSAT, which needs the cube-coverage lemma and the soundness of the char-poly exclusion.
- **T2 (008).** For each of the four seed schemes S (BILR, the two 005 schemes and the 006 scheme), no cyclic-invariant F3 scheme of rank 22 with split (7,1) or (6,4)
  shares ≥ 5 orbits with S. It rests on the k = 1 near() instances, the k = 2 frozen-subset cubes, and their coverage arguments.
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

## Plan
- **Stage 0: literature freshness (no compute).**
  - Re-check PRIOR_WORK.md's cyclic section against new sources since 2026-11-22: the Sedoglavic catalogue, an arXiv search for 3x3 cyclic/Z3/symmetric
    decompositions over F3 or GF(3), the Moosbauer–Poole code, and anything citing Yang 2024 or Tomar 2026.
  - Outcome: a list of the novelty claims we may make ("first Z3-invariant search over F3", "first F3 cyclic classes beyond BILR"),
    each kept, softened or dropped.
- **Stage 1: one-command re-check (≤ 90 CPU-min).**
  - `scripts/recheck_cyclic_f3.py` regenerates every CNF behind T1 and T2 from the encoder (no files from `scratch/`). For each one it:
    - solves with CaDiCaL and checks the proof with drat-trim;
    - solves with Kissat and checks *that* proof with drat-trim, as the second anchor;
    - records the CNF's SHA-256, both statuses and both verdicts in `results/f3/recheck_009.jsonl`.

    The instance set:
    - T1: (1,19); (2,16) cubes 0–20, all with the char-poly exclusion on. This is a uniform configuration, and it is just as valid under the lemma, even though 19 of them were
      first proved without the exclusion. Also the a = 0 trivial check;
    - T2: 8 near() k = 1 instances and 168 frozen k = 2 cubes.
  - Expected cost: about 12 CPU-min for T1 and about 15 for T2 with CaDiCaL including proofs, plus Kissat at an unknown factor (assume ≤ 3×). If over budget, run Kissat on T1 in full
    and on a stated random half of T2's cubes.
  - `--quick` mode for readers: T1 only, CaDiCaL only, about 12 min on one core.
- **Stage 2: the encoder check that drat-trim cannot give.**
  - A new test compares the encoder's equation map with the residual tensor computed directly. For random slot contents, each Brent equation's left side
    (evaluated from the CNF's product/sum structure) must equal the matching entry of Σ(terms) computed by `expand`. This covers every split used in T1 and T2.
    With it, "the CNF says the Brent equations" is tested, not assumed.
  - List the remaining trusted base with a pointer for each:
    - the symmetry-breaking soundness (sign and lex);
    - the cube-coverage lemma;
    - the char-poly table;
    - near() and freeze coverage;
    - drat-trim itself.
- **Stage 3: data release.**
  - `results/f3/cyclic_rank23.jsonl`: one line per class, giving orbits, fixed terms, source attempt, class key, and the number of orbits with a rank-1 matrix.
  - `scripts/check_cyclic_rank23.py` verifies each scheme mod 3, checks cyclic invariance, recomputes class keys, checks they are pairwise distinct and not in the saved reps,
    and prints a one-line summary. A test runs it.
  - Optional, if cheap: a lifting check. Does any of the 10 lift to Z with coefficients in {−1, 0, 1}? Report it, don't chase it.
- **Stage 4: the document `writeups/cyclic-f3.md`.**
  1. Results (T1, T2, D1) in three sentences, then the precise statements.
  2. Setup: the convention, σ, orbits and fixed terms, and the two F3 facts with proofs.
  3. The encoding: variables, equations, symmetry breaking, with soundness arguments.
  4. Class cubes: classes up to conjugation and sign, minimal-class pinning, the char-poly exclusion, the coverage lemma with proof.
  5. Local theorems: near() and frozen cubes, and why their union covers the neighbourhood.
  6. The rank-23 census (D1): how each scheme was found, and a table of invariants.
  7. Where it stops (E1, E2): every timeout with its cap, and why we stopped.
  8. Re-checking: commands, runtimes, the trusted base.
  9. Related work, from stage 0.

  Every claim links to its evidence: a test name, a line of `recheck_009.jsonl`, or an attempt note.
- **Stage 5: review and release.**
  - Read the document cold, starting with its weakest link (the coverage lemma, then the encoder) and checking each proof step against the code.
  - Then update LOG.md (entries for 007, 008 and 009 are missing), TARGETS.md (the F3 cyclic row: T2 and the class count), and README.md ("Attempts 000–009").
  - Commit, merge into `main`, delete the branch, push.
  - Offer a shareable web page of the write-up as a separate step.

## Decision rule (fixed before running)
- **Any re-check disagreement:** a SAT answer, a TIMEOUT where a proof was claimed, or a proof that fails in either checker. Stop. Withdraw the affected claim from the
  write-up until the cause is found, and record it here. Don't publish a theorem with one anchor missing without saying so.
- **The equation-map test fails:** every UNSAT claim is suspended until the encoder is fixed and stage 1 is re-run.
- **Stage 0 finds prior F3 cyclic work:** rewrite the novelty claims, and if the prior work has data, compare our 11 classes with theirs.
- **Kissat is far over budget:** re-check T1 in full with Kissat, re-check T2 with Kissat on a random half, and say so in the document.

## Out of scope
New searches of any kind. If the write-up exposes a gap (for example, a cheap extension of T2 to another seed), note it under Next instead of running it.
