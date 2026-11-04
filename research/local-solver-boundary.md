# Local-solver boundary: local compressions are real, the two-cut product is not yet compressed

## Result and scope

The local-solver route has been split into six checked diagnostics:

1. The source identity minor is an equality join, not a Cartesian product
   ([note](hash-join-identity-minor.md)): at `s=1`, `8,192` message-plus-output
   operations replace `1,048,576` pair comparisons.
2. The source fan-out is a forward map
   ([note](fanout-forward-join.md)): at `h=2`, `65,536` messages replace
   `65,536*(2^32+1)^4` input-message pairs.
3. A proper Fourier-coset bundle compresses jointly
   ([note](bundle-joint-compression.md)): at `h=2`, `65,536` joint pairs replace
   the `2^32` code-space product.
4. Two bundles compress only when their input digits overlap
   ([note](two-bundle-joint-image.md)): shared inputs give factor `65,536`;
   disjoint inputs give factor `1`.
5. A shared pair plus an independent bundle follows the component rule
   ([note](three-bundle-joint-image.md)): the joint image is the shared graph
   times the independent rank.
6. The Gaussian two-cut support is about one third full
   ([note](two-cut-support-density.md)): factors `2.4615` at interval `4` and
   `2.9351` at interval `8`.

The local compressions are real, but the two-cut product is not compressed by
them. The component rule explains why: independent digit blocks multiply. The
Gaussian support saves only a small constant. This note records the boundary
and the next object; it does not change the factoring exponent.

## The remaining object

The raw-support route is bounded by the support density. The larger local
factors live in the **codec**: the row and column codes are functions of the
same digits, and a joint graph table can be much smaller than the product of
the two code ranges. To test the same effect across the two cuts, we need the
two cut codec bases, not just their Smith orders or attaining minors.

The current native oracle exposes the complete orders and the attaining nested
minors (`RectangularOracle.prefixes`), but not the Smith transformation
matrices or the cut codec bases. A bounded implementation can add an output of
the two cut transformation matrices to
[`native_rectangular_smith_f4.cpp`](../tools/native_rectangular_smith_f4.cpp)
and compute the joint image of the two bases.

A finite exact unweighted start is recorded in
[the codec-image note](two-cut-codec-joint-image.md): the interval-4 Smith
transformations give `56` row-codec and `56` column-codec pairs out of `256`,
while the exposed interval-8 prefix proxy is `3656/65536`. The weighted codec
and its transfer remain open.

The same check records the full four-tuple `(row0,column0,row1,column1)`.
At interval 4 the exact transformed four-tuple has `104` values, equal to the
support; the interval-8 prefix four-tuple has `22328`, also equal to the
support. This is forced for an injective full basis: the four codes determine
the four tensor indices. A full-rank codec therefore cannot beat the raw
support factor; only a strict rank truncation or a non-injective scalar code
could. If the intended two-cut bases are the full Smith bases, this closes
that version of the route at the support factor. The codec-image note also
tabulates the prefix-truncated image at `k<n` and records a lossless lower
bound: if a codebook `C` is paired with a residual alphabet `R`, injectivity
on the support forces `|C||R| >= support`, so truncation defers collision
resolution rather than removing it. A changed source or a certified
sublinear collision resolver is required for this route.

## Research question

Let `B_0` and `B_1` be the two cut codec bases, with dimensions `P_0` and
`P_1`. Is the joint image of `(B_0 entry, B_1 entry)` over the tensor support
submultiplicative in `P_0 P_1`? If yes, the two-cut pad can be compressed and
the local-solver route remains open. If the image is the full product for
independent cut digit sets, the route is closed for the current source
construction and a changed source family is required.

This is the next concrete theorem/diagnostic; it is not a parameter sweep.

## Reproduction

```sh
make verify-local-solver-boundary
```

This runs the six cited finite checkers. Their scopes and limitations are
recorded in the individual notes.
