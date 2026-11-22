# Benchmarks: how long the best methods take

Why this file exists: our time budgets must be long enough for a real method to
succeed. Otherwise a "failure" only means we stopped too early.

## Our machine

Apple M5, 10 cores (4 performance + 6 efficiency), 32 GB RAM, macOS.

## What the papers report (2026-11-22)

The best current method is **flip-graph random walks** (Kauers & Moosbauer 2022 and
later versions). Papers report run times very rarely. Everything we could find:

| Format | Field | Rank reached | Reported cost | Hardware | Source |
|---|---|---|---|---|---|
| 3x3x3 | GF(2) | 23 (best known) | >95% of walks reach 23 within 10^7 flips | not stated | Kauers & Moosbauer 2022 (from a summary of the text, not checked word for word) |
| any | – | – | paths of 10^8 flips "within minutes", almost no memory | not stated | Kauers & Moosbauer 2022 (summary) |
| 5x5x5 | GF(2) | 96 → 95 | "within seconds", starting from AlphaTensor's rank-96 scheme | not stated | Kauers & Moosbauer 2022 (summary) |
| 5x5x5 | Q | 93 | "within few minutes on a laptop" (symmetric search, flip limit 10^8) | laptop | Moosbauer & Poole 2025 |
| 6x6x6 | Q | 153 | ~30 min on 160 cores ≈ **80 CPU-hours** (flip limit 10^9) | 160-core machine | Moosbauer & Poole 2025 |
| speed | – | – | 1–6 billion flips per minute | varies | Moosbauer & Poole 2025 (summary) |
| 4x4x4 | Q | "competitive rank" | "relatively short search" (no number) | – | Perminov 2026 |
| 3x3x6, 2x4x5 | Z, GF(2), GF(3) | stalls at 42 and 33 (best known 40 and 32) | does not reach the record | – | Perminov 2026 |
| new records | many | 1–5 improvements per day | meta flip graph running for **weeks** | 20-core Xeon | Perminov 2026 |

Sources: [arXiv 2212.01175](https://arxiv.org/abs/2212.01175), [arXiv 2502.04514](https://arxiv.org/html/2502.04514),
[arXiv 2603.02398](https://arxiv.org/html/2603.02398), [arXiv 2510.19787](https://arxiv.org/html/2510.19787) (gives no timings at all).

## What this tells us so far

- **Rediscovering known small ranks takes seconds to minutes.** For 3x3x3 → 23 and
  4x4x4 → 47 over GF(2), a budget of minutes on one core should be enough.
- **New records take days to weeks** of multi-core time, even for the experts.
  A 10-minute budget is far too short for that. Record-hunting budgets should be
  hours per run, repeated across many runs.
- No paper gives seconds-to-rank for small formats on stated hardware.
  **We need to measure it ourselves** (next section).

## Local measurements (2026-11-22)

Tool: Perminov's `ternary_flip_graph`, commit 1ab814b (MIT license), built with Apple clang
and Homebrew libomp. Each run: 1 thread, 1 runner, starting from the textbook scheme,
stopping at the target rank, killed at its cap. 10 seeds, 2 runs at a time.
**Every final scheme was re-verified by our own verifier.** Raw data: `results/benchmarks.jsonl`.
Rerun with `python3 -I scripts/bench_fg.py`.

| Format | Ring | Target | Runs | Success | Median | Max | Cap |
|---|---|---|---|---|---|---|---|
| 2x2x2 | {-1,0,1} | 7 | 10 | 10/10 | 0.08s | 0.08s | 30s |
| 3x3x3 | {-1,0,1} | 23 | 10 | 10/10 | 0.78s | 4.33s | 120s |
| 3x3x3 | GF(2) | 23 | 10 | 10/10 | 0.06s | 0.06s | 120s |
| 3x3x3 | F3 (mod 3) | 23 | 10 | 10/10 | 0.19s | 0.65s | 120s |

**4x4x4 and larger: not measured.** Four parallel 4x4x4 runs made the Mac too
sluggish, so the user stopped them. We are scoping the project to **3x3x3 and smaller**.

## Budgets we will use

- Reaching the known record (3x3x3 → 23) takes **under 5s on one core**.
  A search that can't reach 23 within **60s** is broken, not unlucky.
  This is our sanity check for any new method.
- The real open question for 3x3x3 is whether **rank 22** exists (best known 23,
  lower bound 19). Hundreds of thousands of rank-23 schemes have been found and none
  reduce to 22, so a plain flip-graph search is very unlikely to find it. A new idea is needed.
  Runs aimed at 22 should be short, bounded, and repeated (e.g. 60–300s each),
  with at most 2 running at once.
