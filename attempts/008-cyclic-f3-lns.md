---
id: 008-cyclic-f3-lns
date: 2026-11-22
target: 3x3x3 rank 22 over F3, cyclic-invariant schemes near the known cyclic rank-23 schemes (finding side; follows 007)
method: SAT neighbourhoods — "a cyclic scheme of a given split sharing at least m orbits with scheme S", one instance per (S, split, m); a walk at rank 23 and reduction attempts to rank 22, with DRAT proofs for the UNSAT ones
status: running
best_rank_found: –
budget: 2 processes × 1 thread; stage A ≤ 40 CPU-min; stage B ≤ 120 CPU-min; ≤ 160 CPU-min in all
---

# 008: cyclic F3 neighbourhood search

## Hand-over
- 007 is on `main` (commit `3b35377`). 008 is on branch `claude/008-cyclic-lns`, which is merged into `main` and deleted when 008 ends.
- Tools and house rules as in 007:
  - `cadical`, plus drat-trim in `vendor/`;
  - check `uptime` first, and use at most 2 processes × 1 thread, each under a timeout;
  - verify every scheme with `matmul.verify --mod 3`, and save F3 schemes in `results/f3/`;
  - write one `record_run` line per instance;
  - **seeds start at 80001**.
- Seed schemes, all of split (7,2): BILR (`tests/data/bilr_2fix_z3.json`), `results/f3/005_3x3x3_rank23_*.json` (2 files), and
  `results/f3/006_3x3x3_rank23_n3_a7_b2_c0_cadical_s60001.json`. They are pairwise inequivalent (006).

## Why not a cyclic flip walk (the 006 backup)
Checked before planning:
- In the four seed schemes, the number of (orbit, orbit-member) pairs that share a factor (same position, up to sign) is 1, 1, 1 and 3.
  No orbit shares a factor with a fixed term. A walk would have one to three moves available.
- An orbit flip changes two orbits and keeps the number of fixed terms. Starting from (7,2), it can reach a lower rank only if an orbit vanishes, which gives (6,2) = rank 20.
  It can never reach a split of rank 22: (7,1) needs one fixed term fewer, and (6,4) two more.

## Idea: neighbourhoods as single SAT instances
For a cyclic scheme S with orbits S₁..S₇, add to `CyclicF3CNF` the constraint **near(S, m, M)**: between m and M orbit slots hold an orbit of S.
"Hold" means equal to some member rotation of some Sⱼ, in sign normal form.
- **Reduction to rank 22, distance k:** split (7,1) (its fixed term is forced to be I) or (6,4), with near(S, 7 − k, ·).
  UNSAT means there is **no cyclic rank-22 scheme of that split sharing ≥ 7 − k orbits with S**. That is a local theorem, and its DRAT proof is checked.
- **Walk step, distance k:** split (7,2) with near(S, 7 − k, 6). Any SAT answer differs from S in at least one orbit, so it is a different rank-23 scheme.
- k = 7 is the whole split again, which timed out in 006. **The largest k that closes is the measure of progress**, between the trivial k = 0 and the full problem.

What this changes from the first draft of this plan:
1. **One instance per k, not one per subset.** The draft froze a chosen subset of 7 − k orbits, which takes C(7, k) instances per k: 35 at k = 3 or 4.
   Here the solver chooses which orbits to keep, learned clauses are shared across all choices, and UNSAT covers every subset at once. That makes it a statement about
   the neighbourhood, not about the subsets we happened to try.
2. **The draft's walk was not sound.** It blocked the old freed orbits in their 3^k rotations, but a free slot could take an old orbit from a *kept*
   slot's list, or a permutation of the freed ones, and return the same scheme. near(S, ·, 6) counts every slot against every orbit of S in every rotation, so a SAT answer is a new scheme.
3. **The symmetry breaking stays.** Sign and lex normalisation apply to all slots as usual. That is sound, because near() compares against S's orbits in every rotation and in sign normal form,
   so the normal form of any qualifying scheme still qualifies.
4. **A filter that does not work, checked and dropped:** "the residual must have small flattening rank". The flattenings of a 3x3x3 tensor are 9 × 81, so
   their rank is at most 9 = 3k already at k = 3. Not worth coding.

Encoding:
- e_{i,j,r} ⇔ "slot i equals N(rot_r(Sⱼ))". These are 7 × 7 × 3 indicators per instance, each the conjunction of 27 value literals.
- s_i ⇔ OR of the e_{i,j,r}.
- At least m of the s_i: for every set of 8 − m slots, a clause saying one of them is true. At most M: for every set of M + 1 slots, a clause saying one of them is false.
  With 7 slots this is at most 35 clauses.
- If one orbit of S occurs in two slots, it is counted twice. That only weakens the constraint, so UNSAT stays sound for "sharing ≥ m distinct orbits".

## Plan
- **Stage 0: code and tests (no budget needed).**
  - Add `near=(orbits_S, m, M)` to `CyclicF3CNF`, with clause range `near_clauses`, and witness support for the new indicators.
  - Add `scripts/run_008.py` with shapes `reduce71`, `reduce64` and `walk`, on top of `matmul.satrun.solve_instance`. SAT answers record `class_key` and whether
    the scheme is among the 51,555 saved F3 rank-23 representatives (`scratch/003`, if present).
  - Tests:
    - each seed scheme satisfies near(itself, 7, 7) at split (7,2);
    - so does a re-labelled copy: orbits permuted, members rotated, signs flipped, then put in `normal_form`;
    - near(itself, 7, 6) rejects it;
    - the gadget agrees with brute force (random slot contents against a direct count);
    - reduce71 at k = 0 is UNSAT, because the two fixed terms differ as a tensor from (I, I, I), and the solver confirms it;
    - a walk on 2x2 (split (2,1), near(S, 1, 1)) finds a verified cyclic rank-7 scheme different from S, or is UNSAT. Record which.
- **Stage A: the k ladder (≤ 40 CPU-min, seed 80001).**
  - For each seed scheme and each of reduce71 and reduce64: k = 1, 2, 3, … (reduce64 starts at k = 1, because 6 slots hold at most 6 orbits of S).
    Cap 60 s for k ≤ 3 and 180 s for k ≥ 4. **Stop a ladder at its first timeout.** Every UNSAT gets a proof run and a drat-trim check (cap 600 s).
  - Walk viability: for each seed scheme, a walk step at k = 2, 3 and 4, cap 60 s. Record SAT/UNSAT and whether the class is new.
  - Order: the 006 scheme first, then BILR and the two 005 schemes, ladders before walks. Hard stop at 40 CPU-min, and list what was not run.
- **Stage B: walk with reduction ladders (≤ 120 CPU-min, 2 walkers × 60 min, seeds 80002–80003).**
  - Walker 1 starts from the 006 scheme, walker 2 from BILR.
  - Each step:
    - move by a walk step at k_w, the smallest k whose stage A walk step was SAT in ≤ 30 s, raising k when a step is UNSAT;
    - from the new scheme, run both reductions at k_r, the largest k that closed on every seed in stage A, cap 60 s.
  - Proofs only for a sample, one step in ten, to keep the budget for search.
  - Log every step to `scratch/008/`. Write one `record_run` line per instance.
- **Measures:**
  - the largest k closed per seed and split (stage A);
  - distinct cyclic rank-23 classes visited (we know 4);
  - local theorems proved;
  - new classes compared with the 51,555 saved representatives.

## Decision rule (fixed before running)
- A rank-22 model that passes `matmul.verify --mod 3` and `is_cyclic_invariant`: stop. Save it, check the catalogue, try lifting (see 002).
- **Stage A:**
  - Ladders reach k ≥ 4 on every seed, and walk steps are SAT: go to stage B.
  - Ladders stall at k ≤ 2: even small neighbourhoods are too hard for SAT. Record that and skip stage B. The proof side (006–007) already met the same wall.
  - Walk steps UNSAT up to k = 4 on every seed: the known schemes are isolated at that distance. That is a result in itself.
    Stage B then only runs ladders, at the highest k reached, with no walk.
- **Stage B ends with no 22:** record the census and the local theorems. Reductions UNSAT at k_r across ≥ 50 classes is evidence (not proof)
  that no cyclic 22 lies within k_r orbits of the part of the cyclic rank-23 locus we can reach. That is the decision point named after 007:
  write up the cyclic results (the ≤ 2-orbit theorem, the new rank-23 schemes, the census and the local theorems) before choosing a new symmetry or method.
- **The walk stalls** (fewer than 20 classes in 60 min): don't add seeds. Record it as for isolation.

## Backup
If near() instances are slow even at k = 2, fall back to the first draft's frozen subsets at k ≤ 3 (C(7,k) small instances each). Keep the walk
sound by adding near(S, ·, 6) to it, which costs little with most slots frozen.

## Stage 0 result
- `CyclicF3CNF(..., near=(S, m, M))`, with clause range `near_clauses`, slot flags `near_flags` and witness support (`alleq`, `or`).
  For the four seed schemes it adds 21 targets per slot (7 orbits × 3 rotations, deduplicated) and at most 35 cardinality clauses.
- Instance sizes, near(BILR, 4, ·): reduce71 has 23,086 vars and 153,532 clauses; reduce64 has 22,756 / 150,464; walk has 24,057 / 159,678.
- `scripts/run_008.py A` runs the stage A ladders and walk steps. A ladder stops at its first non-UNSAT. The budget is checked between instances, so the
  hard stop can be overrun by at most one instance and its proof.
- Tests: 128 passed. New tests:
  - the seed schemes satisfy near(itself, 7, 7), and so does a copy re-labelled by permutation, rotation and sign;
  - near(·, 0, 6) rejects the copy through the near clauses only;
  - the gadget agrees with a direct count;
  - reduce71 at k = 0 is UNSAT;
  - a 2x2 walk step (split (2,1), sharing ≤ 1 orbit) is SAT and gives a verified cyclic rank-7 scheme.
- Smoke run, not recorded (seed 1, scheme discarded): a walk step at k = 2 from BILR was SAT in 9.5 s. Its scheme verifies mod 3, is cyclic,
  and has class `798722b6a87a9844`. That differs from BILR's class and from all 51,555 saved F3 rank-23 classes. So the walk is not frozen.

## Stage A result (seed 80001, 15.6 of 40 CPU-min)
| seed scheme | reduce71 k=1 | reduce71 k=2 | reduce64 k=1 | reduce64 k=2 | walk k=2 | walk k=3 | walk k=4 |
|---|---|---|---|---|---|---|---|
| 006 | UNSAT 5.4 s ✓ | TIMEOUT 60 s | UNSAT 1.9 s ✓ | TIMEOUT | SAT 1.6 s | SAT 10.9 s | TIMEOUT |
| BILR | UNSAT 4.1 s ✓ | TIMEOUT | UNSAT 2.1 s ✓ | TIMEOUT | SAT 16.2 s | TIMEOUT | TIMEOUT |
| 005-cadical | UNSAT 7.4 s ✓ | TIMEOUT | UNSAT 3.7 s ✓ | TIMEOUT | SAT 12.2 s | SAT 11.8 s | TIMEOUT |
| 005-kissat | UNSAT 4.8 s ✓ | TIMEOUT | UNSAT 2.7 s ✓ | TIMEOUT | SAT 22.7 s | SAT 4.5 s | TIMEOUT |

✓ = DRAT proof VERIFIED by drat-trim.

- **Local theorems (k = 1), 8 of 8:** no cyclic F3 rank-22 scheme of split (7,1) or (6,4) shares 6 orbits with any of the four seed schemes.
- **k = 2 times out on all 8 ladders.** As in 007 (2 → 3 orbits), the step from k = 1 to k = 2 is abrupt.
- **Walk: 7 SAT answers, 7 new classes.** All are verified mod 3 and cyclic, all differ from their seed's class and from each other, and none is among the 51,555 saved F3
  rank-23 representatives: `6f888e94663cd913`, `ec8ce52996e9a191`, `3f98b9cc15d7710b`, `ac0136d1e3ad372d`, `6d732e700b4ee825`,
  `ad5563bc5f3b2158`, `346a978816ebcbb6`. Saved as `results/f3/008_*_walk_*.json`. Known cyclic F3 rank-23 classes: 4 → 11.
- **Decision rule:** the ladders stall at k ≤ 2, so stage B is skipped. The plan's backup applies, because near() is slow at k = 2. Freeze each subset of 7 − k shared orbits separately:
  21 small instances per (seed, shape) at k = 2. The union of those cubes is the same local theorem as the near() instance, because every scheme sharing ≥ 5 orbits with S
  shares some specific 5. The backup runs within the 24 CPU-min left of stage A's budget.
