# Exact interval-4 codec image and the interval-8 prefix proxy

This note continues the [local-solver boundary](local-solver-boundary.md).
The remaining object there is the joint image of the two cut codec bases,
not the raw tensor support.  The native weighted oracle exposes complete
orders and attaining nested minors, but not the Smith transformation
matrices.  A finite exact unweighted specialization is now available.

## Interval 4: exact Smith codecs

For the interval-4 reflection Gaussian, the two cut matrices are
`16 x 16` over `F4[x]`.  The checker computes exact Smith decompositions
`A = U S V`, verifies both reconstructions, and records the invariant
degrees.  For cut 0 they are

```text
0, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 12, 12, 12, 12, 16,
```

and for cut 1 they are all zero: the cut-1 matrix is unimodular over
`F4[x]` in this specialization.  The transformation degrees are at most
`7` for cut 0 and `15` for cut 1.

Over the tensor support, the exact unweighted row-codec image and
column-codec image each contain `56` distinct pairs, out of `16^2 = 256`
possible pairs.  The raw tensor support is `104`, so the codec image is
smaller than the support.

The exposed nested-minor prefix partition is not merely the same
cardinality: it is identical to the exact unweighted row- and column-codec
partition on all `104` support entries.  There are `48` two-element
classes and `8` singleton classes.  Thus at interval 4 the oracle's
prefix data already determines the unweighted Smith codec image exactly,
not just its size.

Raw tensor values are shared by the two cuts at each four-index entry, so
only `16` distinct raw value pairs occur.  This value collision is not a
codec image and is not used below.

## Interval 8: the cheaper prefix proxy

Exact Smith transformations for the `256 x 256` interval-8 cut matrices
would be expensive in the current Python exact arithmetic.  The exposed
nested-minor prefix image is cheap, so the checker records it separately.
The interval-8 support is `22328` of `65536` entries.  The prefix image
is `3656` row pairs and `3656` column pairs, a fraction

```text
3656 / 65536 = 0.0557861328125
```

of the full product.  The tensor has `429` distinct polynomial values.

The interval-4 exact codec and prefix counts agree at `56`, but the
interval-8 prefix number is not a Smith codec image.  It is a finite
proxy for the still-missing weighted codec computation.

The prefix count is also stable under deterministic weight changes.  The
checker repeats the interval-8 count at nine gauges: all heights zero,
positive heights on either end, alternating heights, three small random
height vectors, and two wider random vectors.  Every gauge gives exactly
`3656` row pairs and `3656` column pairs.  The selected row/column
permutations do change; only the joint image size is invariant.  This
reduces the concern that `3656` is an artifact of one gauge, but it is
still a finite observation about the exposed prefix image, not a proof
about the weighted Smith codec.

The row and column marginals are not the joint four-tuple.  For the full
basis, the checker also counts

```text
(row code cut 0, column code cut 0, row code cut 1, column code cut 1).
```

At interval 4 the exact transformed four-tuple has `104` distinct values,
equal to the support, and the nested-minor four-tuple also has `104`.  At
interval 8 the nested-minor four-tuple has `22328`, again equal to the
support.  This is forced for any injective row/column basis: the four
codes determine the two row labels and two column labels, hence the
four-index tensor entry.  Therefore a full-rank Smith basis cannot
compress the two-cut joint image beyond the raw support factor.

A codec compression must come from a strict rank truncation or a
non-injective scalar code, not from the invertible Smith transformations.
If "the two cut codec bases" means the full Smith bases, the local-solver
route is closed at the support factor.  The next bounded object is the
prefix-truncated image at `k < n`, which the existing nested-minor data
can already expose.

## Prefix truncation

For a prefix length `k`, retain the first `k` selected nested minors and
send every later row or column to the code `-1`.  The interval-8 row
image then grows subquadratically:

| `k` | row pairs | column pairs | joint four-tuple |
| ---: | ---: | ---: | ---: |
| 32 | 211 | 202 | 1208 |
| 64 | 520 | 503 | 3351 |
| 128 | 1422 | 1398 | 10925 |
| 192 | 2513 | 2524 | 19405 |
| 256 | 3656 | 3656 | 22328 |

At `k=128` the row image is `1422/16384 = 0.0868` of the square of the
code range; at `k=256` it is `3656/65536 = 0.0558`.  This is a genuine
submultiplicative row-codec marginal.  The full four-tuple still reaches
the support at `k=256`, so the row marginal alone does not compress the
complete row-and-column pad.  A composition theorem would need to combine
the row marginal with the column code and the independent digit blocks;
the component rule in the bundle notes is the relevant template.

A one-sided truncation, keeping `k` row codes and the full column code,
has the following interval-8 joint image sizes.  The last column is the
factor by which the code count is smaller than the support `22328`.

| `k` | row-truncated, full column | support/code factor |
| ---: | ---: | ---: |
| 16 | 5897 | 3.786 |
| 32 | 8152 | 2.739 |
| 64 | 11860 | 1.882 |
| 128 | 17750 | 1.258 |
| 256 | 22328 | 1.000 |

## Lossless lower bound and the value-class diagnostic

The truncation table is not a lossless compression.  If a code is a
truncated pair `(r_k, c_k)` and two support entries share that code, a
decoder must receive additional residual information to separate them.
More generally, for any deterministic lossless codebook `C` and residual
alphabet `R`, injectivity on the support gives

```text
|C| * |R| >= support.
```

Thus shrinking the codebook by a factor `F = support/|C|` requires a
residual alphabet of size at least `F` on average.  The full-basis
four-tuple has `F=1`; the `k=32` one-sided table has `F~2.74` and the
two-sided `k=32` table has `F~18.5`, but those factors are paid back as
residual resolution.  The truncation data therefore locate where
resolution work is deferred, not a free reduction of the local-solver
image.

The raw value classes give the same warning in a scalar form.  At
interval 4 there are `16` distinct raw values over `104` support entries,
with largest class `24`.  At interval 8 there are `429` distinct values
over `22328` entries, with largest class `522` and average class
`22328/429 ~= 52.0`.  A value-only hash has `429` codes and must resolve
classes as large as `522`; it is not a lossless codec by itself.  These
counts are diagnostics for a possible non-injective scalar code, not an
instance of one.

## Research question

Let `B_0` and `B_1` be the weighted two-cut codec bases.  Does the exact
joint image of `(B_0 entry, B_1 entry)` over the tensor support stay near
the interval-8 prefix proxy, or does it grow to the full product?  A
stable submultiplicative image would keep the two-cut local-solver route
open; a full-product image would close that route for the current source.
The present finite diagnostic does not answer this question.

## Scope

The interval-4 codec is exact only for the unweighted `F4[x]`
specialization.  It is not the native weighted oracle's codec.  The
interval-8 number is only the exposed prefix proxy.  No transfer,
algorithm, or factoring exponent is claimed, and the complete score
remains `1/2`.

## Reproduction

```sh
make verify-two-cut-codec-joint-image
```

The receipt is `evidence/two-cut-codec-joint-image.json`.
