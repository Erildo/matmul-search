# matmul-search

An open research log on finding **faster algorithms for multiplying small matrices**.

## The problem

Multiplying two n×n matrices the textbook way needs n³ scalar multiplications.
In 1969 Strassen showed that 2×2 matrices need only 7 instead of 8. Applied
recursively, this beats n³ for every size. For most small sizes (3×3, 4×4,
5×5, rectangular cases) the minimum number of multiplications (the *tensor
rank*) is still unknown. Any new upper bound is a real result, because matrix
multiplication sits underneath AI, simulation and graphics.

## Rules

1. **Exact verification only.** A scheme counts only if `matmul.verify` proves it
   correct over Q or GF(p), never "numerically close".
2. **Everything is bounded.** Every search takes a `Budget` (max seconds and/or max steps)
   and stops when the budget runs out. Tests are capped at 10s each and the suite at 120s.
3. **Every attempt is recorded**, including failures.

## Layout

| Path | What |
|---|---|
| `src/matmul/` | scheme format, verifier, budgets, run recorder, known schemes |
| `tests/` | test suite; `./run_tests.sh` runs it |
| `schemes/` | verified schemes as JSON (`<format>_rank<r>_<name>.json`); tests verify every file here |
| `attempts/` | one markdown note per attempt (copy `TEMPLATE.md`) |
| `results/runs.jsonl` | machine-readable record of every run, append-only |
| `TARGETS.md` | formats we target, best known ranks, our best |
| `PRIOR_WORK.md` | what others have done, and where the gaps are |
| `LOG.md` | session journal, newest first; "pick up here" pointer |

## Workflow for an attempt

1. Copy `attempts/TEMPLATE.md` to `attempts/NNN-name.md`, then write the hypothesis and the budget.
2. Run the search with its budget. It appends a line to `results/runs.jsonl`.
3. Any scheme it finds is saved to `schemes/` and must pass `./run_tests.sh`.
4. Fill in the result, what went wrong, and what to do better.
5. Add a `LOG.md` entry, update `TARGETS.md`, and commit.

## Usage

```bash
./run_tests.sh                                # all tests, hard time limits
PYTHONPATH=src python3 -m matmul.verify schemes/*.json
PYTHONPATH=src python3 -m matmul.verify my.json --mod 2
```

## Status

Attempts 000–009 done. See LOG.md for the latest state, and `writeups/cyclic-f3.md` for the cyclic F3 results.
