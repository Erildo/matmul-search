# Log

Newest first. One entry per working session: what we did, and where to pick up.

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
- **Pick up here:** choose 005: a deeper pocket search, a symmetry-restricted search, or saving representatives.

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
