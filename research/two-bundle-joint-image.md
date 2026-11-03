# Two bundles compress only when their input digits overlap

## Result and scope

The `p=2` Fourier-coset bundle compresses its row/column code product
because both codecs read the same input digits. Across two bundles the
same test separates the cases:

- **Shared inputs** (two cosets on the same fan): the joint image has
  one pair per input. At `h=2` it has `65,536` pairs instead of the
  `2^32` rank product; at `h=1` the factor is `4`.
- **Disjoint inputs** (two independent fans): at `h=1` the joint image
  is the full `16`-pair rank product, factor `1`.

Thus the joint compression is local to shared input digits. It composes
when bundles share digits, and does not compose across independent
inputs. This is a finite codec diagnostic; it does not prove a full
two-cut pad theorem or change the factoring exponent.

## Mechanism

For two bundles `(R_1,C_1)` and `(R_2,C_2)`, the joint image is

```math
\{(R_1(x_1),C_1(x_1),R_2(x_2),C_2(x_2))\}.
```

When `x_1=x_2`, the pair is the graph of a single input, so the image
has at most the input count. When `x_1` and `x_2` are independent, the
two graphs multiply and the image is the full rank product (the `h=1`
case confirms this exactly). A composition theorem must therefore track
which bundles share digits; a uniform claim that all bundle boundaries
compress jointly would be false in the disjoint case.

This is consistent with the joint-boundary-compression discussion in
[the half-top structure note](half-top-structure.md) and with
[the single-bundle diagnostic](bundle-joint-compression.md). It also
explains why the fan-out forward map and the bundle graph do not by
themselves reduce a pad built from independent cut codebooks.

## What remains

The source's two cuts are built from bundles on the two operands. Some
digits are shared through the common input; others are independent.
A full compression needs a digit-sharing graph and a bound on the joint
image of each connected component. This note supplies the two extreme
cases, not the general composition.

The next bounded diagnostic is a three-bundle chain with one shared and
one independent digit block, measuring whether the image is the product
of the shared-block graph and the independent rank. That would give the
composition rule a nonuniform pad analysis needs.

## Reproduction

```sh
python3 tools/check_two_bundle_joint_image.py \
  --output evidence/two-bundle-joint-image.json
make verify-two-bundle-joint-image
```

The checker enumerates the shared-input cases at `h=1,2` and the
disjoint-input case at `h=1` exactly.
