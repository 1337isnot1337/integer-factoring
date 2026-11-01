# Certified sequential low-rank Smith updates

The structured Smith engine now evaluates a sum of up to 32 rank-one
updates, while propagating every remaining update through the exact row
and column operations. It retains absolute series-error bounds between
phases. This makes multi-term physical correction searches practical
without treating a truncated series prediction as a complete profile.

The construction supplies an exact computational interface and finite
controls. It does not improve the factoring exponent. Its first use is
the [nine-term flag family](nine-term-flag-templates.md), beyond the
[eight-term fixed-gauge obstruction](low-rank-two-copy-flag-floors.md).

## 1. Keep the diagonal units and every future vector

Suppose a local reference is diagonal in regular frames. Write the
desired matrix as

\[
 D+\sum_{j=1}^{r}u_jv_j^T.
\]

The [rank-one engine](certified-structured-smith.md) diagonalizes the
first update. If its operations are `L` and `R`, the result is

\[
 D_1+\sum_{j=2}^{r}(Lu_j)(v_j^TR).
\]

The new diagonal `D_1` contains the actual pivot series, including
their units and precision bounds. Keeping only their valuations would
give the wrong input for the next update. The engine repeats this
operation until it has incorporated every correction.

For a pivot `p=A[i,j]`, its regular elimination uses the multipliers
`A[k,j]/p` and `A[i,l]/p`. Each future column probe transforms as

\[
 u_k\longmapsto u_k-\frac{A_{kj}}p u_i,
\]

and each future row probe transforms as

\[
 v_l\longmapsto v_l-v_j\frac{A_{il}}p.
\]

In characteristic two the subtractions are additions. The removed
pivot coordinates are appended to the corresponding output probes.
Row and column permutations are tracked separately when the pivot is
off the diagonal. The six residual-state branches use the same literal
pivot row and column to compute these maps; they do not infer a probe
map from the pivot valuations alone.

The next phase receives all pivot coordinates and all remaining ones.
Future corrections can therefore couple coordinates that were separated
in a previous phase. No block is discarded merely because an earlier
update has already diagonalized it.

## 2. Precision is retained across phases

The arithmetic uses the representation

\[
 z^v p(z)+O(z^E),
\]

with a fixed relative storage window. Addition, multiplication,
division, cancellation, and storage truncation propagate the absolute
error exponent `E` as in the certified rank-one engine. A zero visible
polynomial with finite `E` is an unknown tail, not an exact zero.

Every pivot must have a known leading coefficient and valuation no
larger than the certified lower bound of every active entry. Hence the
row and column eliminations are regular. The minimum candidate search
uses exact coefficient arithmetic and integer order comparisons.

After a phase, the retained diagonal and transformed probes keep their
existing error bounds. They are never reimported as exact polynomials.
The unit-inverse cache only reuses the polynomial inverse computation;
each division still computes its own error bound from its two operands.
An unresolved pivot, an intermediate singular proxy, or a lost precision
window causes an explicit unresolved result. It cannot become a claimed
negative result or a reported complete profile.

For rank `r` and matrix size `n`, the structured algebra uses
`O(r^2 n^2)` series operations, including propagation of all remaining
probes. This count charges `r` by value and is not an algorithm for
arbitrary large matrices at a constant cost. Polynomial multiplication
and precision storage remain additional costs.

## 3. Transfer to the complete physical crossing

The [two-copy adapter](../tools/low_rank_two_copy_updates.py) starts with
the retained five-letter reference and its tensor square. The physical
alphabet has 25 letters and each cut matrix has size 625. Each correction
is a literal four-vector tensor times an endpoint-zero polynomial of
degree at most 36. Its four physical vectors may correlate the two
copies arbitrarily.

The reference frames satisfy the exact finite-jet identities through
order 128, with invertible constant terms. The native sandwich routine
computes the transformed physical vectors in those same jets. Each
structured proxy is thus congruent to the full original weighted cut
matrix, after regular frame changes, modulo `z^128`.

If the certified engine exposes all 625 orders below 128, those orders
are also the complete orders of the original matrix. The ranks on every
truncated quotient up to that precision agree, which determines all
orders below the cutoff. The reported orders restore the common local
monomial shifts, including the infinity degree shift. Increasing work
precision inside the elimination does not silently increase the frame
precision or relax this final cutoff condition.

## 4. Independent finite controls

The [small-matrix checker](../tools/check_certified_low_rank_smith_f4.py)
constructs 672 complete random polynomial matrices of dimensions
1, 2, 3, 4, 8, 12, and 16, using 1, 2, 3, or 6 updates. It compares
the structured answer with a separately assembled complete matrix at
work precisions 256 and 512. The retained evidence records:

- 1,344 complete-matrix profile comparisons;
- 120 independent binary block-Toeplitz rank checks;
- 63 reversed-update-order comparisons;
- 47,424 literal endpoint entries; and
- nonzero coverage of all six branches while future probes are active.

Its minimum certified pivot error margin is 249. In particular, the
controls exercise off-diagonal and arrow transitions with nontrivial
probe propagation, rather than testing only diagonal updates.

The [physical checker](../tools/check_low_rank_two_copy_updates.py)
separately inserts every four-leg correction into the complete original
625-by-625 matrices. Its four cases use 2, 6, 9, and 12 terms. The
twelve-term case includes three copies of each of four updates with
coefficients `1, omega, omega^2`; they cancel exactly and restore the
reference. This checks cancellation across consecutive phases.

That checker verifies 3,125,000 literal endpoint entries and 16 complete
profile comparisons at work precisions 256 and 512. The minimum
certified pivot error margin is 116. Every final order lies below the
frame cutoff. These are independent finite controls, not a formal
verification of the native implementation.

The records are
[small-matrix controls](../evidence/certified-low-rank-smith-f4-controls.json)
and [physical controls](../evidence/low-rank-two-copy-update-controls.json).
Both retain the exact implementation source hashes.

```sh
python3 tools/check_certified_low_rank_smith_f4.py > /tmp/low-rank-controls.json
python3 tools/check_low_rank_two_copy_updates.py > /tmp/low-rank-physical-controls.json
cmp /tmp/low-rank-controls.json evidence/certified-low-rank-smith-f4-controls.json
cmp /tmp/low-rank-physical-controls.json evidence/low-rank-two-copy-update-controls.json
```

The C++17 implementation is
[`native_certified_low_rank_smith_f4.cpp`](../tools/native_certified_low_rank_smith_f4.cpp),
with its [Python interface](../tools/native_certified_low_rank_smith_f4.py).
The separate [nine-term search](../tools/search_nine_term_flag_updates.py)
retains exact profiles and integer all-scale prefix comparisons. Any
reported numerical pressure improvement still requires its own rigorous
pressure-sign certificate before a stronger factoring bound is claimed.

No Lean source is added or changed.
