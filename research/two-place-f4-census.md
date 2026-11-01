# Registered full F4 coefficient census of the coupled two-place pencil

## Target and scope

This run enumerates every polynomial

```math
 q(t)=(1+t)\sum_{j=0}^{10}c_jt^j,\qquad c_j\in F4.
```

There are exactly **`4^11=4,194,304`** coefficient vectors. Multiplication
by `1+t` is injective and its image is exactly the degree-at-most-eleven
polynomials vanishing at one. Thus the registered run covers the whole
finite F4 pencil, including zero, at its inherited physical heights.
It is not an all-field census or an all-height optimization.

The predecessor is [PR57](https://github.com/1337isnot1337/integer-factoring/pull/57)
at exact head `8a3ea3e6aeacba7d0a53764de5579e3707bb2ca7`. Its retained
198-case catalogue and all-field **balanced** `8/11` obstruction are
credited and preserved. The new census targets arbitrary coefficient
polynomials and asymmetric supporting lines. The
[determinant classification](two-place-determinant-target.md) separately
closes one sufficient asymmetric escape over all fields; it does not
settle the full asymmetric criterion.

The full census is **in progress at this checkpoint**. Its
[manifest](../evidence/two-place-f4-census-manifest.json) fixes the
complete candidate family, OS-random seed, source and binary hashes,
mathematical input hash, and enumeration order before execution.
No completed coverage claim is made here. The factoring score remains
**1/2**.

## Exact enumeration and coverage

The [input constructor](../tools/two_place_f4_census.py) verifies the
square-zero generator, all nineteen sandwich identities, the exact
commutator factorization and the physical-energy compatibility. It
reconstructs the retained [integer input](../evidence/two-place-f4-census-input.txt)
without importing a module from the unmerged branch.

The [native scanner](../tools/scan_two_place_f4_pencil.cpp) uses a 22-bit
code for the eleven F4 coefficients. For enumeration position `i`, it
uses `(A*i+B) mod 2^22`, with odd `A`. This is a bijection, because an
odd integer is invertible modulo a power of two. The
[runner](../tools/run_two_place_f4_census.py) validates a contiguous
chunk ledger, the inverse permutation of every retained representative,
all outcome counts, and every representative's actual polynomial.

Every candidate has an exact strict endpoint check and both complete
25-order profiles at integer weights. The native polynomial field
arithmetic is exact; these weights require no floating approximation.
An unresolved promised pivot is recorded separately. Every zero-cut
order sum must equal the conjugation invariant 444. A chunk error
cannot be reported as successful coverage.

Each newly encountered complete profile class is independently replayed
by packed-F4 fraction-free polynomial elimination, without a truncated
series window. Only then is its asymmetric pressure fitted numerically.
Those fits use the same intrinsic scale, normalized supporting weights,
and both physical pair charges. An apparent `R<1` is a candidate for
separate rational certification, not an accepted bound. Positive fitted
values do not prove an all-angle or all-scale exclusion.

The runner checkpoints every chunk and holds a process lock. A resume
rejects changed sources, binaries or mathematical input. Compile and
replay the focused controls with `make verify-two-place-f4-census`.
The default proof gate also includes `make verify-two-place-determinant`.

## Independent and predecessor controls

The [checker](../tools/check_two_place_f4_census.py) matches all 198
retained source cases, comprising 9,900 complete orders. Twelve
additional coefficient vectors have independent full fraction-free
replays of 600 orders. A nonidentity affine permutation checks all
512 visits in a separate coverage-control chunk.

The [receipt](../evidence/two-place-f4-census-controls.json) retains
actual polynomials and profiles. The
[source packet](../evidence/two-place-pencil-source-8a3ea3e.json.xz)
preserves the predecessor receipt and exact head. The first native
benchmark covered all 4,096 vectors with `deg r<=5`, found 25 complete
profile classes, and had no unresolved pivots or endpoint failures.
This benchmark is a scoped finite control, not completion of the full run.

## Earlier registered searches

The [Pascal search packet](../evidence/pascal-registered-search-round1.json.xz)
retains all 96 registered permutation/flag attempts from frozen source
`566b596`, their original logs, parameters, and OS-random selection.
Ninety attempts produced search records; six stopped at their startup
gradient guard and are explicitly retained as failures. The recorded
searches made 741,832 objective evaluations. No stronger contraction
certificate was found. The same archive retains the separate full
16-letter symbol search: 4,563 queries over its 900-second budget,
two completed restarts, and no improvement over its initial reference.
Neither statement is a global optimality result.

The [earlier asymmetric scan](../evidence/two-place-asymmetric-retained-profiles.json.xz)
checks all 66 distinct retained profiles, including the reference, over
its declared numerical search region. It reproduces the symmetric
normalization before optimizing. Its best fitted supporting value is
above one. The full coefficient census extends that finite catalogue;
it does not infer anything about untested fields or real gauges.

Substantial OpenAI assistant support was used for the construction,
proof, implementations, replays and text. No new Lean formalization or
improved uniform factoring exponent is claimed.
