---
id: 006-cyclic-f3-cubes
date: 2026-11-22 (planned)
target: 3x3x3 rank 22 over F3, cyclic-invariant schemes (continues 005)
method: SAT (CaDiCaL) with conjugation symmetry broken by a case split on the similarity class of one matrix; one cube per class; drat-trim proof per cube
status: planned
best_rank_found: –
budget: 2 processes × 1 thread; stage A ≤ 30 CPU-min; stage B ≤ 21 cubes × 120 s; stage C ≤ 21 cubes × 120 s (≤ 115 CPU-min in all)
---

# 006: cyclic-invariant F3 search, split into cubes by conjugacy class

## Hand-over (read this first in a new session)
State after 005 (`attempts/005-cyclic-f3-sat.md`, commit `2cb137e` on branch `claude/005-symmetric-f3-723f13`):
- **Check the branch first.** `git log --oneline main | grep 2cb137e`. If 005 is not on `main` yet, branch 006 from `claude/005-symmetric-f3-723f13`, not from `main`.
- **Tools.** `cadical` (3.1.0) and `kissat` (4.0.4) are installed system-wide by Homebrew. **drat-trim exists only in the 005 worktree's git-ignored `vendor/`.** Rebuild it in the checkout you work in:
  `git clone --depth 1 https://github.com/marijnheule/drat-trim vendor/drat-trim && cc -O2 -o vendor/drat-trim/drat-trim vendor/drat-trim/drat-trim.c`
- **Code.** `src/matmul/cyclic.py`:
  - `CyclicF3CNF(n, a, b, fixed_sum=True, break_sym=True)` builds the CNF. It has `.dimacs()`, `.decode(model)`, `.witness(orbits, fixed)`, `.satisfied(assign)` and `.trivially_unsat`.
  - The module also has `expand`, `is_cyclic_invariant`, `normal_form` and `sigma`.
  - `scripts/run_005.py` runs a single instance: solve, decode, verify, save to `results/f3/`, re-solve with a proof, run drat-trim, and record with `record_run`. Reuse it (move `run_instance` into `src/matmul/` if 006 needs it too).
  - Tests: `tests/test_cyclic.py`. Known cyclic schemes: `tests/data/bilr_2fix_z3.json` and `results/f3/005_3x3x3_rank23_*.json`.
- **Facts already established.** Over F3 a fixed term is (x, x, xᵀ). The fixed matrices sum to I. A scheme of fixed terms alone is impossible, so a ≥ 1, b ≥ 1. The rank-22 splits are (7,1), (6,4), (5,7), (4,10), (3,13), (2,16) and (1,19). (1,19) is proved UNSAT. The other six timed out after 10 min each with CaDiCaL. Rank 23 at (7,2) took 257 s with CaDiCaL, 907 s with Kissat.
- **House rules.** Check `uptime` first. Use at most 2 processes with 1 thread each, every one under a timeout. **Seeds start at 60001** (used so far: 1–180, 9001–9008, 50001–50003). Verify every scheme with `matmul.verify --mod 3`. Save F3 schemes in `results/f3/`, never in `schemes/`. Write one `record_run` line per instance. Commit on a branch; the user merges to `main`.

## Why
In 005, six of the seven rank-22 splits stayed open after 10 minutes. The encoding breaks rescaling and reordering, but not the
largest remaining symmetry. Simultaneous conjugation g_M: (U, V, W) → (MUM⁻¹, MVM⁻¹, M⁻ᵀWMᵀ) commutes with σ, so it
maps cyclic-invariant schemes to cyclic-invariant schemes. Modulo scalars, it acts by a group of order |PGL₃(F3)| = 5,616. So every
solution currently has up to about 5,600 copies in the search space. Breaking this symmetry should shrink each split a lot, and the case split
also cuts each split into independent cubes, so progress can be measured (cubes closed out of 21).

## The case split (to prove in the note before running)
Claim: every cyclic-invariant scheme with a ≥ 1 orbits can be brought, by symmetries, into a form where **orbit slot 1 is
nonzero and its u is the chosen representative of one of the 21 classes** of nonzero 3x3 matrices over F3. A class here is a
conjugacy class up to sign (computed in this session: 39 similarity classes including 0, 22 up to sign including 0, so 21 nonzero ones).
Sketch:
1. A nonzero orbit has u ≠ 0, otherwise its terms are zero tensors. Move a nonzero orbit to slot 1; orbit slots can be reordered.
2. Rescale (u, v, w) → (−u, −v, w) if needed, then apply g_M to bring u to its class representative. g_M changes the other slots as well,
   but it keeps the scheme correct and cyclic-invariant (to be checked by a test, as for the moves in `tests/test_invariants.py`).
3. Then restore the constraints on the *other* slots: sign normalisation of slots 2..a, lex order on slots 2..a, and lex order on the fixed slots.
   Slot 1 gets no sign normalisation of u (fixed by the pin), keeps it for v, and is left out of the orbit lex chain.
The 21 cubes together cover every solution. So the split is UNSAT if and only if all 21 cubes are UNSAT, given (i) this lemma and (ii) a test that the
known schemes, transformed this way, satisfy the cube of their own class.

## Plan
- **Stage 0, code and tests (no budget needed).**
  - Add `class_reps()` (union-find over the 3^9 matrices with conjugation by generators of GL₃(F3), plus x ~ −x; the representative
    is the `_bitkey`-smallest member) and `to_class_form(orbits, fixed)` (the transformation in the lemma) to `matmul.cyclic`.
  - Add `pin_u1=rep` to `CyclicF3CNF`.
  - Tests: 39 and 22 class counts; class sizes sum to 3^9; g_M keeps the BILR scheme and both 005 rank-23 schemes correct and cyclic-invariant;
    each of those schemes, put in class form, satisfies the cube of its own class.
  - Check that `CyclicF3CNF(3, 0, 22).trivially_unsat` holds, which confirms a = 0 is impossible at 3x3.
- **Stage A, calibration at rank 23, split (7,2) (≤ 30 CPU-min).** Find the classes of u among the orbit members of the
  two 005 schemes. Run those cubes (cap 300 s, seeds 60001–60002), plus unpinned (7,2) with seeds 60003–60004 (cap 300 s),
  to see how much seed-to-seed spread 005's 257 s hides.
  - Pinned cubes find 23 at least as fast as unpinned: pinning helps finding, go on.
  - Pinned is clearly slower (all pinned runs time out while unpinned finds): the pinning hurts the solver. Go to stage B only for proofs (UNSAT side), and say so.
- **Stage B, split (2,16), all 21 cubes × 120 s.** Smallest open split first, the most likely to close.
  - Every cube UNSAT and every proof VERIFIED: theorem **"no cyclic-invariant F3 scheme of rank ≤ 22 with ≤ 2 orbits"**.
  - Some cubes open: record which classes stay open. The number closed out of 21 is the progress measure. Don't rerun longer without a new idea.
- **Stage C, split (7,1), all 21 cubes × 120 s.** This is the rank-22 split closest to (7,2), which solves fast at rank 23. Its single fixed term is forced to be I.
  This is the most likely place to *find* a 22.
  - SAT: verify mod 3 and check cyclic invariance. If both pass, it is a rank-22 scheme over F3. Stop, save it, check the catalogue, try lifting (it would need even denominators, see 002).
  - UNSAT cubes: record them with proofs. Timeouts: list them.
- Stages B and C run 2 cubes at a time. Compute bound ≤ 115 CPU-min, at most about 60 min of wall time.

## Decision rule (fixed before running)
- Any rank-22 scheme that passes `matmul.verify --mod 3` and `is_cyclic_invariant`: stop, then do the follow-up above.
- (2,16) fully closed: record the theorem, and plan 007 to close (3,13) the same way (cubes split further by the class of the second u if needed).
- Few or no cubes close in B and C: conjugation was not the bottleneck. 007 switches the finding side to a cyclic flip-graph walk over F3
  (flip orbit representatives, applied to all three members) seeded from the two 005 rank-23 schemes, or to local search.

## Backup
A cyclic flip-graph walk over F3 (Moosbauer–Poole style, our own code, since `ternary_flip_graph` has no symmetric mode), seeded from
`results/f3/005_3x3x3_rank23_*.json`.

## Result
(to fill in)

## What went wrong / surprised us

## What we would do better

## Next
