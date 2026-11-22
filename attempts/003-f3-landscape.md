---
id: 003-f3-landscape
date: 2026-11-22
target: 3x3x3 rank 23 over F3, with GF(2) as the control
method: find_alternative_schemes (ternary_flip_graph) + our own orbit invariants
status: planned
best_rank_found: –
budget: 2 processes × 1 thread × 10 min per stage (quieter than 002)
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
(to fill in)

## What went wrong / surprised us

## What we would do better

## Next
004, chosen by the decision rule above.
