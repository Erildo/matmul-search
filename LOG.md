# Log

Newest first. One entry per working session: what we did, and where to pick up.

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
