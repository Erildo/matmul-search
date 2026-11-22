---
id: 001-benchmark-flip-graph
date: 2026-11-22
target: 2x2x2, 3x3x3 (4x4x4 dropped)
method: ternary_flip_graph (Perminov), plain flip-graph random walk
status: success
best_rank_found: 3x3x3 rank 23 (= best known)
budget: 30s / 120s cap per run, 10 seeds, 2 in parallel
---

# 001: benchmark the best existing method

## Hypothesis
Rediscovering known ranks should take seconds. We need real numbers so our
budgets are not too short.

## Setup
- Tool: vendor/ternary_flip_graph @ 1ab814b. Built with
  `make flip_graph CXX=clang++ FLAGS="-Wall -O3 -std=c++17 -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp"`
- Command: `python3 -I scripts/bench_fg.py`

## Result
| Case | Success | Median | Max |
|---|---|---|---|
| 2x2x2 → 7 | 10/10 | 0.08s | 0.08s |
| 3x3x3 → 23 over {-1,0,1} | 10/10 | 0.78s | 4.33s |
| 3x3x3 → 23 over GF(2) | 10/10 | 0.06s | 0.06s |

All 30 final schemes passed our own verifier (new importer: src/matmul/import_fg.py).

## What went wrong / surprised us
- 4x4x4 with 4 parallel runs made the Mac too sluggish. The user stopped it, and no 4x4x4 numbers were recorded.
- The tool has no time limit of its own; it runs until the target or forever. We always wrap it in a timeout.
- GF(2) is about 10x faster than {-1,0,1} for 3x3x3.

## What we would do better
- Run at most 2 processes at once.
- Check machine load before longer runs.

## Next
002: first attempt at something new for 3x3x3.
