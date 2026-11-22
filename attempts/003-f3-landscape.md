---
id: 003-f3-landscape
date: 2026-11-22
target: 3x3x3 rank 23 over F3, with GF(2) as the control
method: find_alternative_schemes (ternary_flip_graph) + our own orbit invariants
status: partial
best_rank_found: 23 (no rank change was sought)
budget: 2 processes × 1 thread, at most 20000 schemes or 300s per field
---

# 003: how many genuinely different rank-23 schemes exist over F3?

## Why
In 002 the pool filled to 2000 "different" schemes at ranks 24–26 but stopped at **164** at rank 23.
"Different" there means a different *flip fingerprint* (the tool's `structure` hash, built from
which flips a scheme has available). Two schemes related by a change of basis can have different
fingerprints, and two genuinely different schemes could share one. So 164 isn't a real count yet.

The real count decides what to do next:
- **Few classes that stop growing:** the F3 rank-23 world is small enough to study scheme by scheme.
  Next, look for rank 22 next to each class, or search only schemes with a symmetry.
- **Many classes that keep growing:** the 002 stall was the search circling, not the landscape being
  small. Next, change how we explore (longer pools, the tool's meta flip graph).

## The control
GF(2) is a case where the answer is roughly known: Heule–Kauers–Seidl found more than 13,000
inequivalent rank-23 schemes ([arXiv:1905.10192](https://arxiv.org/abs/1905.10192)). Running exactly the same
procedure over GF(2) tells us what our method sees in a big landscape. If GF(2) also levels off at a
few hundred, the leveling off is our method, not F3.

## Plan
1. **Our invariant** (`src/matmul/invariants.py` + tests). Two schemes are equivalent if one maps to
   the other by (a) changes of basis (U, V, W) → (P U Q⁻¹, Q V R⁻¹, R W P⁻¹), (b) rescaling a term,
   (c) reordering terms, (d) the cyclic and transpose symmetries of the product. We compute
   quantities that can't change under those moves, so a different value proves the schemes are different:
   - the multiset of rank triples (rank U_t, rank V_t, rank W_t), taken over the symmetric versions;
   - a finer one: for each pair of terms, the ranks of U_s ± U_t (likewise for V and W), as a multiset.
   - Tests: the invariant must not change under random basis changes, rescalings, reorderings and the
     cyclic/transpose maps; it must separate a few schemes we know differ.
   - It gives a **lower bound** on the number of classes (equal invariants don't prove equivalence).
2. **Build** `find_alternative_schemes` and `analyze_schemes` with the same flags as `flip_graph`.
3. **Run** `find_alternative_schemes --ring Z3 --target-rank 23` from one 002 scheme and, in parallel,
   `--ring Z2` from one GF(2) rank-23 scheme. 10 minutes each, 1 thread each, under a hard timeout.
   Record the time of every saved scheme.
4. **Analyze:** verify every scheme (mod 3 or mod 2), compute invariants, and plot the number of
   distinct classes against time for both fields. Also run the 165 schemes from 002 through it.
   Compare our classes with the tool's `type invariant` as a cross-check.

## Decision rule (fixed before running)
- F3 levels off (no new class in the last 30% of the run) while GF(2) keeps growing: **small** → 004
  studies the classes one by one and searches for 22 next to them, or with a symmetry restriction.
- Both keep growing: **large** → 004 changes the exploration method.
- Both level off: our invariant or the generator is too weak. Improve them before drawing conclusions.

## Result
**Invariants** (`src/matmul/invariants.py`, tests in `tests/test_invariants.py`, 80 tests pass). The tests check
that every symmetry move maps correct schemes to correct schemes (on naive, Strassen and an F3 rank-23 scheme)
and leaves the invariant unchanged, alone and in random combinations.

**002's schemes:** all 165 are pairwise inequivalent (165 classes, 96 term profiles).

**Smoke test** (F3, seed 1, 15s, not recorded in runs.jsonl): 40,223 schemes. In the first 10,000: 98/100,
975/1,000, 9,608/10,000 distinct classes.

**Main run** `python3 -I scripts/run_003.py 20000 300` (generator seed 3, 1 thread each; lines in `results/runs.jsonl`):

| Field | Schemes | Failed verify | Classes (lower bound) | Term profiles | Classes after 1k / 3k / 10k / 20k |
|---|---|---|---|---|---|
| F3 | 20,000 (5.8s) | 0 | **8,176** | 89 | 982 / 2,928 / 6,746 / 8,176 |
| GF(2) | 20,000 (12.6s) | 0 | **15,808** | 303 | 997 / 2,783 / 8,126 / 15,808 |

Reading it:
- Both landscapes are large: thousands of proven-inequivalent rank-23 schemes within seconds of one seed.
  The 164 in 002 came from the pool search, not from a small landscape.
- F3 slows down: between 10k and 20k schemes only 14% were new classes (1,430 of 10,000), against 77% for
  GF(2). GF(2) also has about 3 times as many term profiles.
- This isn't the invariant being weaker over F3. The F3 invariant uses c ∈ {1,2}, GF(2)'s only c = 1, so
  it's finer over F3 if anything.
- But a single walk explores near its seed, and the seed matters (seed 1 gave 96% new at 10k, seed 3 only 67%).
  So the slowdown may be one walk's neighbourhood running out, not F3 as a whole.

**Against the decision rule:** neither fixed outcome. F3 doesn't level off (new classes still appear at the
end), but it slows clearly, while GF(2) keeps growing. Intermediate, so no branch is chosen yet.

## What went wrong / surprised us
- The generator is extremely fast (40k schemes in 15s). The cap has to be a scheme count, not a time.
- Spotlight indexed the 40k smoke files and pushed the load average to 15. `scratch/.metadata_never_index` now marks the folder.
- The decision rule didn't cover "slows but doesn't stop". Next time state a threshold for the slope too.

## What we would do better
- Several walks from several unrelated seeds, then count classes in the *union*. One walk can't tell its own neighbourhood from the whole landscape.
- Fix the slope threshold in advance (for example, "new-class rate below 5% over the last 20%").

## Follow-up: 5 unrelated seeds per field (rule fixed before running, 2026-11-22)
Setup: for each field and each k = 1..5, a fresh rank-23 scheme from the naive scheme (`flip_graph --seed k`),
then `find_alternative_schemes --seed k`, 20,000 schemes. Every scheme is verified and classified. Per field we report:
- per walk: classes, and the **late rate** (share of new classes among the walk's last 20% of schemes);
- **overlap**: classes found by two or more walks, and for each pair of walks i, j the capture-recapture
  estimate N̂ = cᵢ·cⱼ / mᵢⱼ (cᵢ classes in walk i, mᵢⱼ shared). This estimates the size of the region the walks sample.

Rule:
- **Small**: F3 walks share at least 10% of their classes on average (mean mᵢⱼ / min(cᵢ, cⱼ) ≥ 0.10), and the
  GF(2) share is at least 5 times lower → F3 region of about N̂ classes. 004 studies it directly.
- **Large**: F3 walks share less than 1% → treat F3 as large. 004 changes the exploration method.
- **In between** (1–10%): report N̂ as a rough size and choose 004 on cost.

(Results in "Follow-up results" below.)

## Next
Settle the intermediate result cheaply: F3 and GF(2), 5 unrelated seeds each, 20k schemes per seed, count classes
in the union and the per-field slope. Then apply the decision rule (with the slope threshold) and choose 004.
