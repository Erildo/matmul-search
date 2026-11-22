---
id: 004-f3-many-pockets
date: 2026-11-22
target: 3x3x3 rank 22 over F3
method: ternary_flip_graph, ring Z3: many short independent walks from the naive scheme, one per seed
status: failed
best_rank_found: 23 over F3 (no 22 in 180 pockets)
budget: 2 processes × 1 thread × 30 min (3,600 CPU-seconds), 20 s per seed → about 180 seeds
---

# 004: search for rank 22 over F3 in many pockets

## Why
003 showed the F3 rank-23 landscape is split into pockets. Five independent walks shared no class,
and some walks used up their pocket (one ran dry at about 4,900 classes). 002 spent an hour of
compute inside a few pockets. If rank 22 is reachable from only some pockets, many short
walks from different seeds have a much better chance than a few long ones.

## Design
- Seed k = 1..N: `flip_graph -n1 3 -n2 3 -n3 3 --ring Z3 --target-rank 22 --count 1 --threads 1 --seed k`,
  killed after 20 s. Each walk starts from the naive scheme, so it lands in a pocket chosen by its seed.
- Per seed, record: time to reach 23 (from file creation times), best rank, and the class key
  (`matmul.invariants`) of its final rank-23 scheme, to confirm the seeds really sample different pockets.
- Keep one rank-23 scheme per seed (`scratch/004/reps/`, git-ignored). Any rank-22 scheme is verified
  mod 3 at once and saved to `results/f3/`, outside `schemes/`, because `schemes/` holds schemes valid over Q.

Why 20 s: one walk reaches 23 in under 1 s (002 stage A, max 0.65 s), so 20 s gives about 19 s of rank-22
search per pocket, about 10⁷ flips. That's the scale at which flip-graph walks usually find what's in reach.

## Decision rule (fixed before running)
- **Any rank-22 scheme that passes `matmul.verify --mod 3`:** stop. Check it against the catalogue, then try lifting to Q
  (it would need even denominators; see 002).
- **No 22, and the final rank-23 schemes of different seeds are pairwise inequivalent:** the pockets were really
  sampled. About 180 pockets without a 22 is evidence (not proof) that 22 is not one short walk away in
  typical pockets. 005 then goes deeper (pool search) in the pockets that look most promising, or changes method.
- **No 22, and many seeds end in the same class:** the seeds are not sampling different pockets. Fix that before concluding.

## Result
`python3 -I scripts/run_004.py 3600 20 1` (seeds 1–180), 1,801 s wall time. One summary line in `results/runs.jsonl`; per-seed lines in
`scratch/004/per_seed.jsonl`, and one rank-23 scheme per seed in `scratch/004/reps/` (both git-ignored).

| Seeds | Reached 23 | Time to 23 (median / max) | Found 22 | Distinct final classes |
|---|---|---|---|---|
| 180 | 180 | 0.21 s / 2.22 s | **0** | **180** (no repeats) |

- No rank-22 scheme in any pocket. All 180 walks ran the full 20 s.
- The 180 final rank-23 schemes are pairwise inequivalent, so the seeds really sampled different pockets
  (as far as the final classes show).
- **By the rule fixed beforehand:** branch 2. About 180 pockets × about 19 s each turned up no 22. That's evidence (not proof)
  that 22 is not one short walk away in typical pockets.

**Overlap with 003.** 6 of the 180 final classes appear among 003's 51,555 F3 classes. Five of them are a design slip:
seeds 1–5 here used the same random seeds as 003's walks 1–5, so they started in the same pockets by construction.
The only independent hit: seed 107 ended in the pocket of 003's walk 3.
If pockets were equally likely, 1 hit among 175 independent seeds against 5 known pockets would suggest
on the order of 5 × 175 ≈ 900 pockets. With a single hit this is a very rough estimate (it could easily be a few hundred
or several thousand), and pockets are surely not equally likely.

## What went wrong / surprised us
- Seed ranges overlapped with 003 (seeds 1–5 in both), which made 5 of the 6 overlaps trivial. Use disjoint seed ranges
  across attempts (for example, 004 seeds from 10,001).
- The rank-22 save path (copy before verifying) was never exercised, because no 22 appeared. It is untested.

## What we would do better
- Decide what makes a pocket "promising" before going deep. Right now we have no signal that tells pockets apart
  except their invariants.

## Next
Options for 005:
1. **Go deeper in a few pockets** (pool search from a seed's scheme). But we lack a signal for which pockets are promising.
2. **Change method:** restrict to schemes with a symmetry (for example, invariant under the cyclic symmetry), so the F3 search space
   gets small enough to search systematically, using SAT or exhaustive flips within the restricted space.
3. **Make the counts checkable:** save one representative per class from 003 (the "save representatives" option).
