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
