# Prior work

What others have done, so we don't repeat it and know where the gaps are.
Searched on 2026-11-22. Claims marked *(unverified)* come from abstracts or
secondary sources we haven't checked yet.

## The field is very active

- **Flip graphs** (Kauers & Moosbauer, 2022 onward): random walks over the
  space of correct schemes, with "flips" that keep a scheme correct and
  "reductions" that lower its rank. This is the main method behind most recent records.
  - Symmetry-aware flip graphs, Moosbauer & Poole 2025: big improvements for 5x5x5 and 6x6x6. [arXiv 2502.04514](https://arxiv.org/pdf/2502.04514)
  - Meta flip graph, Kauers 2025: improved ranks for about 24–30 formats, starting partly from AlphaEvolve schemes. [arXiv 2510.19787](https://arxiv.org/html/2510.19787)
  - Flip graphs for structured products. [arXiv 2511.10786](https://arxiv.org/pdf/2511.10786)
- **AlphaTensor** (DeepMind 2022): RL search; found 4x4x4 rank 47 over GF(2).
- **AlphaEvolve** (DeepMind 2025): LLM-guided search; 4x4x4 rank 48 with complex
  coefficients, plus about a dozen other formats.
  - Turning complex schemes into rational ones. [arXiv 2602.13171](https://arxiv.org/pdf/2602.13171)
- **Open-source flip graph framework** (Perminov, 2026): C++ code, searched every format up to
  16x16x16, improved 79 formats. [arXiv 2603.02398](https://arxiv.org/html/2603.02398),
  [code: dronperminov/FastMatrixMultiplication](https://github.com/dronperminov/FastMatrixMultiplication)
- **Finite fields**: new CP decompositions over GF(2) and others. [arXiv 2505.09824](https://arxiv.org/pdf/2505.09824)
- **3x3x3 rank over F2 and Z is exactly 23** (Tomar, Oct 2026, [arXiv:2610.11823](https://arxiv.org/abs/2610.11823),
  preprint, certificates checked in Lean). Builds on "at least 22 over Z" ([arXiv:2610.01639](https://arxiv.org/abs/2610.01639)).
  Over characteristic 0 the rank is still between 19 and 23; F3 is also open. Rank-22 searches over GF(2) or {-1,0,1} are now pointless.
- **Over 13,000 inequivalent rank-23 schemes over GF(2)** found with SAT (Heule, Kauers, Seidl, "New ways to multiply 3x3 matrices",
  [arXiv:1905.10192](https://arxiv.org/abs/1905.10192)).
- **Fewer additions at rank 23** is a crowded secondary goal: 58, 55 and 52 additions reported in 2025–26
  ([2512.21980](https://arxiv.org/abs/2512.21980), [2610.27451](https://arxiv.org/abs/2610.27451), [2609.06588](https://arxiv.org/abs/2609.06588)) *(unverified)*.
- **Catalogue of best known ranks**: Sedoglavic, "Yet another catalogue of fast matrix
  multiplication algorithms", <https://fmm.univ-lille.fr>. This is the reference for whether a result is new.

## What we didn't find

Nobody keeps a public, honest, step-by-step **journey log**: attempts, dead ends,
and what was learned. Published work shows only the final schemes. That is our gap:
a reproducible open notebook, where every run is recorded and every result is verified exactly.

## What this means for strategy

- Small square formats (2–6) are heavily searched by strong tools. A new
  record there is unlikely early on.
- Better odds: **rectangular formats**, **specific fields** (GF(2), GF(3), integers only),
  and secondary goals such as **fewer additions** or **sparser schemes** for a known rank.
- Before calling anything "new", check it against the Sedoglavic catalogue and Perminov's tables.
