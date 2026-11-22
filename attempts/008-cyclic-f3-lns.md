---
id: 008-cyclic-f3-lns
date: 2026-11-22
target: 3x3x3 rank 22 over F3, cyclic-invariant schemes, by local search around the known cyclic rank-23 schemes (finding side; follows 007)
method: large-neighbourhood search with SAT — free k orbits of a cyclic rank-23 scheme and re-solve them, either at the same shape (a walk) or at a rank-22 shape (a reduction attempt)
status: planned
best_rank_found: –
budget: 2 processes × 1 thread; stage A ≤ 40 CPU-min; stage B ≤ 120 CPU-min; ≤ 160 CPU-min in all
---

# 008: cyclic F3 large-neighbourhood search

## Hand-over
- Branch `claude/007-three-orbits` holds 007 (commits `0c74e3f`, `a6f7f62` and the 007 close-out). Branch 008 from it, or from `main` once it is merged.
- Tools and house rules as in 007:
  - `cadical`, plus drat-trim in `vendor/`;
  - check `uptime` first, and use at most 2 processes × 1 thread, each under a timeout;
  - verify every scheme with `matmul.verify --mod 3`, and save F3 schemes in `results/f3/`;
  - write one `record_run` line per instance or per walker;
  - **seeds start at 80001**.
- Seed schemes, all of split (7,2): BILR (`tests/data/bilr_2fix_z3.json`), `results/f3/005_3x3x3_rank23_*.json` (2 files), and
  `results/f3/006_3x3x3_rank23_n3_a7_b2_c0_cadical_s60001.json`. Loaders: `split_orbits` and the `bilr()` helper in `tests/test_cyclic.py`.

## Why not a cyclic flip walk (the 006 backup)
Checked before planning:
- In the four seed schemes, the number of (orbit, orbit-member) pairs that share a factor (same position, up to sign) is 1, 1, 1 and 3.
  No orbit shares a factor with a fixed term. A walk would have one to three moves available and would hardly leave its start.
- An orbit flip changes two orbits and keeps the number of fixed terms. Starting from (7,2), it can reach a lower rank only if an orbit vanishes,
  which gives (6,2) = rank 20. It can never reach a split of rank 22: (7,1) needs one fixed term fewer, and (6,4) two more.
  Flips alone cannot find a cyclic 22 from these schemes.

## Idea
Re-solve a part of a known scheme. Keep 7 − k orbits of a cyclic rank-23 scheme and free the other k orbits (and the fixed terms). What is left
is a SAT instance that is small for small k, because the kept orbits enter as unit clauses. The freed part is solved in one of three shapes:
- **walk**: the same shape (k orbits + 2 fixed), with blocking clauses against the old freed part. A SAT answer is a different rank-23 scheme;
- **reduce to (7,1)**: k orbits + 1 fixed term. The fixed term is forced to be I, since the fixed terms sum to I. A SAT answer is a rank-22 scheme;
- **reduce to (6,4)**: k − 1 orbits + 4 fixed terms. A SAT answer is a rank-22 scheme.

This is the cyclic version of "repair a scheme with one fewer term". It also covers the move a flip cannot make: changing the number of fixed terms.

## Plan
- **Stage 0: code and tests (no budget needed).**
  - Add `keep=[orbits]` to `CyclicF3CNF`. The first len(keep) orbit slots are fixed by unit clauses to the given orbits. Sign normalisation and lex order
    apply only to the free orbit slots and to the fixed slots.
  - Add `block=[(orbits, fixed)]`: a clause excluding each given content of the free slots. For the walk, block the old freed part in all
    3^k member rotations, in normal form, so a mere re-labelling is not accepted.
  - Add `scripts/run_008.py` with the three shapes on top of `matmul.satrun.solve_instance`. Each walk step records the scheme's
    `class_key` (`matmul.invariants`).
  - Tests:
    - every seed scheme, with k orbits freed (k = 1..4), satisfies its own same-shape instance without `block`;
    - with `block`, the original freed part violates the block clauses;
    - for k = 0, the (7,1) reduction is UNSAT, because the two fixed terms are not equal as a tensor to (I, I, I);
    - a 2x2 analogue: from cyclic rank 7 (2,1), freeing both orbits and reducing to (1,3) gives UNSAT, consistent with 005.
- **Stage A: calibration (≤ 40 CPU-min, seed 80001).**
  - Reductions: for each seed scheme, all 35 subsets with k = 3, in both shapes (cap 20 s), and then all 35 with k = 4 in both shapes (cap 60 s),
    until the stage budget is spent. Record SAT, UNSAT or TIMEOUT, and the solve time.
  - Walk viability: for each seed scheme, 10 random subsets at k = 3 and 10 at k = 4, same shape with `block` (cap 30 s). Record how many are
    SAT, and how many of those have a `class_key` different from the seed's.
- **Stage B: walk with reduction attempts (≤ 120 CPU-min, 2 walkers × 60 min, seeds 80002–80003).**
  - Walker 1 starts from the 006 scheme, which has the most shared factors, and walker 2 from BILR.
  - Each step: choose a random k-subset (k from stage A: the largest k whose same-shape median is ≤ 5 s) and move to the SAT answer.
    From the new scheme, try both reductions on 3 random k′-subsets (k′ from stage A, the largest whose reductions finish within the cap).
  - Log every step (classes, times, statuses) to `scratch/008/`. Write one summary line per walker in `runs.jsonl`.
- **Measures:** distinct cyclic rank-23 classes visited (we now know 4); reduction attempts and their outcomes; solve times as the walk moves.

## Decision rule (fixed before running)
- A rank-22 model that passes `matmul.verify --mod 3` and `is_cyclic_invariant`: stop. Save it, check the catalogue, try lifting (see 002).
- **Stage A:**
  - if same-shape re-solves at k = 3 or 4 give a new class in at least 30% of SAT answers, go to stage B;
  - if they almost always return the seed's class, the cyclic locus is rigid near these schemes. Try k = 5 once (10 subsets, cap 60 s);
    if it is still rigid, record the rigidity as the result and skip stage B.
- **Stage B ends with no 22:** record the census of classes and the number of reduction attempts. Reductions UNSAT at every k′ tried, across ≥ 100 classes,
  is evidence (not proof) that no cyclic 22 lies within k′ orbits of the cyclic rank-23 locus we can reach. That is the decision point named after 007:
  the cyclic route has run its course. Write up what we have (the ≤ 2-orbit theorem, the new rank-23 schemes and the census) before choosing a new symmetry or method.
- **Stage B walk stalls** (fewer than 20 classes in 60 min): don't add seeds. Record the stall and treat it like rigidity.

## Backup
If SAT re-solving is too slow even at k = 3, use a cheaper neighbourhood: free k = 2 orbits but also the fixed terms, and enumerate all solutions (CaDiCaL with
blocking clauses, cap 60 s per subset). Then treat the set of rank-23 schemes reachable that way as a graph and search it breadth-first.
