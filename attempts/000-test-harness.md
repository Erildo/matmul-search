---
id: 000-test-harness
date: 2026-11-22
target: infrastructure
method: exact verifier + test suite
status: success
best_rank_found: –
budget: tests capped at 10s each, suite at 120s
---

# 000: test harness

## Hypothesis
Before searching, we need a judge we trust: an exact verifier that accepts
every correct scheme and rejects every broken one.

## Setup
- `src/matmul/verify.py`: exact tensor check over Q or GF(p), plus an independent
  check on random matrices.
- `src/matmul/budget.py`: a hard time/step cap that every search must use.
- Command: `./run_tests.sh`

## Result
57 tests pass in 0.3s. Naive, Strassen and Strassen⊗Strassen (4x4x4 rank 49)
verify. All 21 single-sign corruptions of Strassen are rejected. A scheme that is
correct only mod 2 is told apart from one correct over Q. A deliberately infinite
test was killed at its 1s limit.

## What went wrong / surprised us
Nothing yet. Note: Python's SIGALRM cannot interrupt long-running native (C)
code, which is why run_tests.sh also has an outer `timeout --signal=KILL`.

## What we would do better
The verifier is pure Python. That's fine up to about 6x6x6. For large formats we may need a
faster integer path.

## Next
001: first real search (target to be chosen).
