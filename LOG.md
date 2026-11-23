# Log

Newest first. One entry per working session: what we did, and where to pick up.

## 2026-11-23: attempt 010, S3-invariant F3 schemes, stage A

- Stage 0: S3 encoder (σ plus transpose τ), the 006 scheme in S3 form, tests. For 3x3 the S3 is unique up to equivalence; 2x2 has two.
- **Stage A: all 5 instances timed out** (rank 23 at split (1,5,0,2) in 300 s, and the four nearest rank-22 splits in 120 s). 780 CPU-s.
- By the decision rule, stage B (all 60 rank-22 splits) is on hold: the encoding is too slow even where a solution is known.
- Orthogonal-class cubes rejected before building: the symmetry left after S3 has order 24, which gives 65–479 cubes per split.
- Native S3 encoder (ties built into the variables): 249 → 138 equations, CNF half the size, tests pass. **Rank 23 is still not rediscovered in 300 s.**
- **Pick up here:** freeze part of the 006 scheme and shrink the frozen part, to find the distance at which search stops working; then decide on a pin.

## 2026-11-23: attempts 007–009, end of the cyclic F3 route

- **007:** split (3,13) with every exclusion on, plus "all orbits nonzero" (allowed by 006's theorem). **0 of 21 cubes closed at 120 s.** Stage B was skipped by decision. The proof side stops at ≤ 2 orbits.
- **008:** a cyclic flip walk is pointless: the known schemes have 1–3 flip moves, and flips cannot change the number of fixed terms. Replaced by SAT neighbourhoods.
  - **Local theorem:** no cyclic rank-22 scheme, split (7,1) or (6,4), shares ≥ 5 orbits with any of the four known schemes (168 frozen cubes). Distance 3 times out.
  - **7 new cyclic rank-23 classes** from walk steps.
- **009:** write-up `writeups/cyclic-f3.md`.
  - Re-check: 380/380 runs UNSAT with proofs verified (CaDiCaL and Kissat, drat-trim). Encoder tested against the tensor directly.
  - 11 cyclic F3 rank-23 classes released (`results/f3/cyclic_rank23.jsonl`). All share one structure (no rank-3 matrix, one (2,2,2) orbit).
- Workflow change: no branches. Everything goes to `main`, and run logs go in `results/logs/`.
- **Pick up here:** choose 010. The candidates are listed in `attempts/009-cyclic-f3-writeup.md` → Next: a larger symmetry group (S₃) over F3, or a (7,1) search built around I.

## 2026-11-22: attempt 006, class cubes for cyclic F3 schemes

- Split each case of 005 into 21 SAT cubes, one per class of nonzero 3x3 matrices up to conjugation and sign. Cube c pins the simplest orbit
  matrix to the class representative and forbids simpler matrices elsewhere (no rank-1 matrices, and no characteristic polynomial of a lower
  class). There is a written proof that the cubes cover every solution, plus tests.
- **Theorem, with drat-trim-checked proofs: no cyclic-invariant F3 scheme of rank ≤ 22 with ≤ 2 orbits.** All 21 cubes of (2,16) are UNSAT, and (1,19) was refuted in 005.
- Split (7,1): 2 of 21 cubes refuted, the rest time out at 120 s. No rank 22 found.
- A new cyclic F3 rank-23 scheme, inequivalent to all known ones (`results/f3/006_*`).
- **Pick up here:** 007, close (3,13) with the same cubes (split further by a second matrix if needed), and start a cyclic flip-graph walk for the finding side.

## 2026-11-22: attempt 005, cyclic-invariant schemes over F3 by SAT

- Literature (PRIOR_WORK.md): Z3-invariant rank-23 schemes exist with integer coefficients (Ballard et al. 2018). SAT under cyclic symmetry has only been run over GF(2) (Yang 2024, stuck past rank 15). Moosbauer–Poole's symmetric flip graphs are GF(2) only. ternary_flip_graph has no symmetric mode. We found no Z3-invariant search over F3.
- Derived: over F3, fixed terms are exactly (x, x, xᵀ), and the fixed matrices must sum to I. So rank 22 splits as 3a + b with (a, b) in {(7,1), …, (1,19)}, 198 unknowns in every case.
- Built `matmul.cyclic` (CNF encoder, tested), installed CaDiCaL and Kissat, and built drat-trim in vendor/.
- Sanity checks passed: 2x2 rank 7 found and rank 6 refuted with checked proofs. 3x3 rank 23 found from scratch in 257 s (CaDiCaL) and 907 s (Kissat): two new inequivalent cyclic F3 schemes.
- Rank 22: split (1,19) **UNSAT, proof checked by drat-trim**. The other six splits timed out at 10 min each. So: no cyclic-invariant F3 rank-≤22 scheme with ≤ 1 orbit. Everything else is open.
- Planned 006 (attempts/006-cyclic-f3-cubes.md): break the conjugation symmetry by splitting each case on the class of one matrix (21 classes), with one SAT cube per class and a checked proof per cube.
- **Pick up here:** read the hand-over section of attempts/006-cyclic-f3-cubes.md (branch, tools, drat-trim rebuild, seeds from 60001), then start stage 0.

## 2026-11-22: planning 002

- Found that 3x3x3 rank over F2 and Z is now proved to be 23 (arXiv:2610.11823). Our GF(2) plan was dead on arrival.
- Consequence: a rank-22 rational scheme would need even denominators; such schemes reduce to F3, where the rank is open.
- Planned attempt 002: flip-graph search over F3 (attempts/002-f3-rank22.md).
- Ran 002 step 0 (no published F3 search found) and stage A (F3 reaches 23 in under 1s, 10/10, verified mod 3).
- Ran 002 stage B: 2×30 min over F3, no rank 22. Rank 23 reached easily, but only 164 structurally different rank-23 schemes turned up (vs 2000 at ranks 24–26).
- Found that the tool's "structure" uniqueness is a flip fingerprint, not equivalence. The 164 isn't a real count.
- Planned 003: count F3 rank-23 classes with our own invariants, with GF(2) as the control (attempts/003-f3-landscape.md).
- Ran 003: built and tested orbit invariants. Both F3 and GF(2) have thousands of proven-inequivalent rank-23 schemes (8,176 and 15,808 classes in 20k schemes each). F3 slows down (14% new in the last 10k vs 77% for GF(2)), but from one seed only.
- 003 follow-up, 5 seeds per field: F3 walks share **no** classes (51,555 in the union; some walks get trapped in pockets of about 4k), while GF(2) walks overlap heavily (37,537 in the union, capture-recapture about 37k). F3 is large and fragmented.
- Possible side result: at least 37,537 inequivalent GF(2) rank-23 schemes (published: "more than 17,000"). Needs saved representatives before any claim.
- Ran 004: 180 independent F3 walks × 20 s, each in a different pocket (180 distinct final classes). No rank 22.
- Only 1 independent repeat pocket against 003 (seed 107), a hint of several hundred pockets or more (very rough).
- Saved representatives: 37,537 GF(2) and 51,555 F3 rank-23 schemes, pairwise inequivalent, in `results/reps/`. Independently re-checked with `scripts/check_reps.py`. The rerun reproduced 003 exactly.
- 005 (cyclic-symmetric rank-22 search over F3) started in a separate session on its own branch.
- **Pick up here:** 005's results; then check the literature's equivalence notion before calling the GF(2) count new.

## 2026-11-22: benchmarking the best method

- Collected reported run times (BENCHMARKS.md). Papers rarely state them.
- Built Perminov's ternary_flip_graph and timed it on this Mac. 3x3x3 reaches rank 23
  in 0.06–4.3s, 10/10 runs, all re-verified by us. See attempts/001-benchmark-flip-graph.md.
- 4x4x4 was too heavy for the Mac, so the scope is now 3x3x3 and smaller.
- **Pick up here:** choose a method to attack 3x3x3 rank 22 (or a secondary 3x3x3 goal).

## 2026-11-22: setup

- Chose the problem: fewer multiplications for small matrix products.
- Built the exact verifier, budgets and test suite (see attempts/000-test-harness.md).
- Surveyed prior work (PRIOR_WORK.md). The field is active, but nobody publishes a journey log.
- **Pick up here:** choose the first target in TARGETS.md and confirm its record
  against the Sedoglavic catalogue.
