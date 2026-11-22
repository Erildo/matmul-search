# Targets

The records we measure ourselves against. Before trusting a "best known" number,
confirm it against the [Sedoglavic catalogue](https://fmm.univ-lille.fr) and set Confirmed to yes.

| Format | Field | Best known rank | Lower bound | Source | Confirmed | Our best |
|---|---|---|---|---|---|---|
| 2x2x2 | any | 7 | 7 (optimal) | Strassen 1969; Winograd 1971 | yes | 7 (reference) |
| 3x3x3 | Q | 23 | 19 | Laderman 1976; Bläser 2003 | no | 23 (reproduced with ternary_flip_graph, <5s) |
| 3x3x3 | F2 and Z | 23 (optimal) | 23 | Tomar 2026, arXiv:2610.11823 (preprint, Lean-checked) | no | 23 |
| 3x3x3 | F3 | 23 | ? (open) | – | no | 23 (002: no 22 in 2×30 min; 004: no 22 in 180 pockets × 20 s) |
| 3x3x3 | F3, cyclic-invariant | 23 (Ballard et al. 2018, mod 3; 005: two more found by SAT; 006: a third) | 22 ruled out with ≤ 2 orbits (005, 006; checked proofs); otherwise open | 005, 006 | – | 23 |
| 4x4x4 | GF(2) | 47 | ? | AlphaTensor 2022 | no | – |
| 4x4x4 | C | 48 | ? | AlphaEvolve 2025 | no | – |
| 4x4x4 | Q | 48 or 49? | ? | see arXiv 2602.13171 | no | 49 (Strassen², in tests) |
| 5x5x5 | Q | ~93? | ? | Moosbauer & Poole 2025 | no | – |

**Current scope: 3x3x3 and smaller** (larger sizes are too heavy for this Mac; see BENCHMARKS.md).

Add a row whenever we pick up a new target.
