---
id: 007-cyclic-f3-three-orbits
date: 2026-11-22
target: 3x3x3 rank 22 over F3, cyclic-invariant schemes with 3 orbits (split (3,13)); continues 006
method: SAT (CaDiCaL) on the 006 class cubes with every exclusion on from the start, all orbits forced nonzero, and open cubes split again by the characteristic polynomial of a second orbit matrix; drat-trim proof per cube
status: failed
best_rank_found: –
budget: 2 processes × 1 thread; stage A ≤ 21 cubes × 120 s (+ proofs) ≤ 60 CPU-min; stage B ≤ 120 CPU-min (skipped); 42 CPU-min used
---

# 007: cyclic-invariant F3, three orbits

## Hand-over (read this first in a new session)
- 006 is on `main` (commit `72819e7`). Branch 007 from `main`.
- Tools: `cadical` 3.1.0 and `kissat` 4.0.4 (Homebrew). drat-trim is in the git-ignored `vendor/drat-trim/` of the main checkout.
  Rebuild it in a fresh worktree:
  `git clone --depth 1 https://github.com/marijnheule/drat-trim vendor/drat-trim && cc -O2 -o vendor/drat-trim/drat-trim vendor/drat-trim/drat-trim.c`
- Code: `src/matmul/cyclic.py` (`CyclicF3CNF(n, a, b, pin_u1=, min_class=, exclude_charpoly=)`, `to_class_form(..., minimal=True)`,
  `charpoly`, `charpoly_exclusions`, `class_reps`, `class_index`). Runner: `scripts/run_006.py` (`run_instance`). Move `run_instance`
  into `src/matmul/` now that two attempts use it.
- House rules: check `uptime` first. Use at most 2 processes with 1 thread each, every one under a timeout. **Seeds start at 70001.**
  Verify every scheme with `matmul.verify --mod 3`. F3 schemes go in `results/f3/`. Write one `record_run` line per instance.

## Why
006 proved that no cyclic F3 scheme of rank ≤ 22 has ≤ 2 orbits. The next instance is (3,13): 3 orbits + 13 fixed terms = 22.
It also covers every split with 3 orbits and fewer fixed terms, so closing it extends the theorem to **≤ 3 orbits**.
In 006, (2,16) closed 19 of its 21 cubes at 120 s, and the characteristic-polynomial exclusion closed the last two. (7,1) closed only 2 of 21. (3,13) is
the next step, and the smallest one that tells us how hardness grows with the number of orbits.

Two things are new compared with 006:
1. **All orbits nonzero.** In (3,13), a zero orbit slot leaves a scheme of split (2,b′) with b′ ≤ 13, and 006 refuted all of those.
   So we may require every orbit slot to be nonzero. That is free pruning that 006 could not use.
2. **A second split level.** Open cubes are split by the characteristic polynomial of a second orbit matrix. This reuses the clause
   gadget from 006 (`_charpoly_not_in`) and needs no new class-membership encoding.

## Lemma to prove in the note before running (level-2 cubes)
Fix the sign pairing p ~ p̄, where p = (t, s, d) and p̄ = (−t, s, −d) is the characteristic polynomial of −x. Fix a total order on the
characteristic polynomials.

*Claim.* Every cyclic F3 scheme of split (3,13) has a symmetric image that satisfies the level-1 cube c (006 lemma, all three orbits nonzero),
and in addition the level-2 cube (c, p):
- slot 1's u has characteristic polynomial p;
- every matrix u, v, w of orbit slots 1 and 2 has a characteristic polynomial ≥ p or p̄ in the fixed order. Equivalently, the polynomials q with
  q < p and q̄ < p are excluded on those slots.

Sketch:
1. Apply the 006 transformation, which puts the minimal-class member in slot 0 with u = rep_c. The conjugation is now used up.
2. Among the 6 members of orbits 1 and 2, pick the one whose first matrix has the smallest characteristic polynomial up to sign. Put it in slot 1,
   rotated so that this matrix is u. The first matrices of a member's rotations are u, v and wᵀ, and wᵀ has the same polynomial as w.
3. Use (u, v, w) → (−u, −v, w) to make u's polynomial the chosen one of {p, p̄}. (Fix the choice: the smaller of the two.) Then use
   (u, v, w) → (u, −v, −w) to sign-normalise v. These two moves are independent, so both normalisations hold.
4. Slot 2 keeps the ordinary sign normalisation. There is no lex constraint between slots 1 and 2, since slot 1 is now chosen by a rule. Fixed slots keep their lex order.
Classes stay ≥ c, so the level-1 constraints still hold. The cubes (c, p) over every p that some class ≥ c carries cover cube c.

## Plan
- **Stage 0: code and tests (no budget needed).**
  - Add `all_nonzero=True` to `CyclicF3CNF`: every orbit slot has u, v, w ≠ 0 (one clause per matrix). Forcing u alone is not enough,
    because (u, 0, w) is a zero tensor that the equations allow.
  - Add `pin_cp1=p`: slot 1's u has polynomial p, and orbit slots 1..a−1 exclude every q with q, q̄ < p. Slot 1 leaves the orbit lex chain.
  - Add `to_class_form(..., level2=True)`, the transformation in the lemma.
  - Tests:
    - the 3 known 7-orbit schemes and the new 006 scheme, put in level-2 form, satisfy their own (7,2) cube (c, p);
    - an image that is not minimal is rejected by its cube;
    - on 2x2, the level-2 cubes of (2,1) find rank 7 and those of (1,3) refute rank 6;
    - the order and pairing are consistent: every q excluded in (c, p) is carried by no allowed matrix;
    - `all_nonzero` rejects a slot in which any one of u, v, w is zero.
  - Record the cube sizes of (3,13) at level 1 and the number of level-2 children per cube.
- **Stage A: (3,13), all 21 level-1 cubes, cap 120 s, seed 70001.** Rank-1 bans, the characteristic-polynomial exclusion and `all_nonzero` are on from the start.
  Each UNSAT cube gets a proof run and a drat-trim check (check cap 900 s). Two cubes run at a time.
- **Stage B: level-2 cubes of every cube still open, cap 120 s, seed 70002.** Order: cubes with the fewest children first, so whole level-1 cubes close
  early. Stop when the stage B budget (120 CPU-min) is spent, and list what was not run.
- Progress measure: level-1 cubes closed out of 21, and level-2 cubes closed out of those tried.

## Decision rule (fixed before running)
- A SAT model that passes `matmul.verify --mod 3` and `is_cyclic_invariant`: stop. Save it, check the catalogue, try lifting (see 002).
- All 21 cubes closed, every proof VERIFIED: theorem **"no cyclic F3 scheme of rank ≤ 22 with ≤ 3 orbits"**.
  - If the median cube time is < 30 s, 008 is (4,10) the same way.
  - Otherwise, 008 is the finding side (below).
- Some level-2 cubes still open: record them. Don't add a third level or longer caps without a new idea. The candidate idea: break the remaining
  symmetry, the centralizer of rep_c acting on slot 1. It is largest for the rank-1 classes 0 and 1.
- Stage A closes fewer than half the cubes and stage B closes less than half its children: the proof side has reached its limit at 3 orbits. Write down how hardness grows
  with a (2: 21/21; 3: ?; 7: 2/21) and move to the finding side.

## Backup / finding side (008)
A cyclic flip-graph walk over F3 for (7,1) and (7,2). It flips orbit representatives and applies the flip to all three members. Seed it from the four known cyclic rank-23
schemes (BILR, both 005 schemes, `results/f3/006_3x3x3_rank23_n3_a7_b2_c0_cadical_s60001.json`). Our own code is needed, because
`ternary_flip_graph` has no symmetric mode. This is the decision rule from 006 for the finding side. It runs as its own attempt, because it needs new
code and is independent of the cube work.

## Stage 0 result
- `matmul.cyclic` gains:
  - the options `all_nonzero` and `pin_cp1`, with clause ranges `nonzero_clauses` and `cp1_clauses` for tests;
  - `to_level2_form(orbits, fixed)`, the lemma's transformation. It is a separate function rather than the planned `to_class_form(..., level2=True)`, and it returns `(orbits, fixed, rep, p)`;
  - `cp_neg`, `cp_canon`, `cp1_choices(c)` and `cp1_exclusions(p)`.

  The characteristic-polynomial gadget is now built once per matrix and shared by the level-1 exclusion and the level-2 pin.
  Level-1 instances are unchanged: (2,16) cube 0 is still 21,339 vars and 135,996 clauses, and a test checks this.
- `run_instance` moved to `matmul.satrun.solve_instance`. `scripts/run_006.py` now calls it, and a smoke test reproduces 2x2 rank 7 (SAT) and rank 6 (UNSAT, VERIFIED).
  The runner for this attempt is `scripts/run_007.py` (`sizes`, `A`, `B`). Stage B reads the cubes still open from stage A's lines in `runs.jsonl`.
- The order on polynomials is tuple order on (t, c₂, det). p is canonical when p ≤ p̄.
- `cp1_exclusions(p)` contains (0,0,0) for every p ≠ (0,0,0), so the zero matrix is excluded from orbit slots ≥ 1. This is sound only because all orbits are nonzero.
  It also means the level-2 cubes cover only schemes whose orbits are all nonzero, which `all_nonzero` requires anyway.
- When p = p̄, the sign of slot 1's u is not fixed by p. Slot 1 then keeps the first-nonzero normalisation of u as well.
- Tests: 124 passed. New tests:
  - the four known 7-orbit schemes, BILR, both 005 schemes and the 006 scheme, satisfy their own (7,2) level-2 cube after `to_level2_form`;
  - a slot 1 with a larger polynomial violates exactly the level-2 clauses;
  - the level-2 gadget agrees with brute force for three values of p;
  - the polynomial tables are consistent (pairing, canonical forms, choices not excluded at level 1);
  - `all_nonzero` rejects a zero u, v or w;
  - with slots 0 and 3–6 and the fixed slots of BILR fixed, CaDiCaL completes slots 1 and 2 inside the level-2 cube to a verified cyclic rank-23 scheme.
- **Change from the plan:** the 2x2 level-2 test is replaced by the partial-scheme solver test above, because the characteristic-polynomial gadget is for 3x3 only.
- Sizes for (3,13):
  - level-1 cube 0 has 21,647 vars and 138,016 clauses; cubes 3–20 have 22,327 vars and 142.5–142.7k clauses;
  - the number of level-2 children per cube is 15, 15, 15, 14, 13, 12, 12, 11, 10, 9, 9, 8, 7, 6, 5, 5, 5, 4, 3, 2, 1, which is **181 in all**;
  - at 120 s each, the worst case for stage B is 362 CPU-min, three times its 120 CPU-min budget. The fewest-children-first order is there so whole cubes close inside the budget.

## Stage A result ((3,13), 21 level-1 cubes, cap 120 s, CaDiCaL seed 70001)
**All 21 cubes timed out.** Nothing was refuted and nothing was found. Compute used: 42 CPU-min, all spent in the solves; no proof runs were needed.
Every cube had rank-1 bans, the characteristic-polynomial exclusion and `all_nonzero` on.

For comparison:

| split | orbits | level-1 cubes closed at 120 s |
|---|---|---|
| (2,16) | 2 | 19/21 (21/21 with the exclusion) |
| (3,13) | 3 | 0/21 (exclusion on) |
| (7,1) | 7 | 0/21 (2/21 with the exclusion) |

The third orbit is where the cubes stop closing at this cap. The free pruning from `all_nonzero` was not enough.

Notes for stage B, written before it runs:
- Cube 20 has a single level-2 child. That child pins slot 1's polynomial to the only value it can take, and it drops the lex order between slots 1 and 2.
  It is the stage A cube, slightly weaker, and running it would only repeat a timeout. The same holds, nearly, for cubes with 2–3 children.
  The fewest-children-first order would spend stage B's first runs on exactly these.

## Decision: stage B skipped
Stage B was skipped, against the plan, and the decision was made after stage A. The reasons:
- 0/21 at level 1;
- the 120 CPU-min budget covers about a third of the 181 children;
- the fewest-children-first order would have opened with near-copies of stage A cubes, as noted above.

The decision rule's last branch (proof side at its limit) is taken on stage A alone.

## Result
- (3,13): 0 of 21 level-1 cubes closed at 120 s. The theorem stays at **≤ 2 orbits** (006).
- How the cubes close at 120 s as the number of orbits a grows: a = 2, 21/21; a = 3, 0/21; a = 7, 2/21. Conjugation breaking plus the class exclusions are enough at a = 2, but not at a = 3.
- Stage 0 code (`all_nonzero`, level-2 cubes, `matmul.satrun`) is tested and kept for later use.

## What went wrong / surprised us
- The jump from 2 to 3 orbits is abrupt. Most (2,16) cubes took under 5 s, and no (3,13) cube closed in 120 s.
- The level-2 child counts are lopsided: cubes 0–2 have 15 children each, cube 20 has 1. A child of a cube with few children is barely stronger
  than its parent. Splitting by polynomial helps most where the classes are small, which is the opposite of where the open cubes are.

## What we would do better
- Before planning a stage, check that its first jobs differ from runs that already failed. Here the run order would have repeated stage A.

## Next
008: the finding side. Before planning it, we checked whether a cyclic flip walk could move at all. In the four known cyclic rank-23 schemes, the number of
(orbit, orbit) member pairs sharing a factor is 1 (BILR), 1, 1 (both 005 schemes) and 3 (006). None share a factor with a fixed term.
A cyclic flip walk would barely move, and orbit flips cannot change the number of fixed terms. So they can never reach a split of rank 22
from (7,2): (7,1) needs one fixed term fewer, and an orbit vanishing gives rank 20. 008 therefore uses SAT re-solving of a few freed orbits
(large-neighbourhood search) in place of flips; see `attempts/008-cyclic-f3-lns.md`.
