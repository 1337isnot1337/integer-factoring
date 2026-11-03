# A proper Fourier-coset bundle compresses jointly, not as a code product

## Result and scope

For the source `p=2` Fourier-coset bundle, the row code and column code
are both functions of the same input digits. The set of joint
`(row, column)` codes is therefore the image of one map, not the
Cartesian product of the two code spaces.

At `h=2` the bundle rank is `65,536 = 256^2`. The joint image has
exactly `65,536` pairs, while the product of the two code spaces has
`65,536^2 = 2^32` pairs. A joint table indexed by the shared input
therefore has `65,536` entries instead of `2^32`, a factor of `65,536`.
At `h=1` the factor is `4`.

This is a finite source-codec diagnostic for the `p=2` bundle. It does
not prove a full two-cut pad theorem, a composition over all bundles, or
a factoring rate. The complete score remains `1/2`.

## Mechanism

The source codec assigns to input digits `x` a row code `R(x)` and a
column code `C(x)`. For the proper coset the joint map

```math
x\mapsto (R(x),C(x))
```

is injective on the source digit domain. The two codebooks are
correlated because they are evaluated on the same input; the joint
boundary is the graph of the pair, not the full product of the code
ranges. A table of the graph has one entry per input.

This is a concrete instance of the joint boundary compression named as
an unexcluded route in
[the half-top structure note](half-top-structure.md). It is stronger
than the separate forward-map observation in
[the fan-out note](fanout-forward-join.md): here the compression is
between the two sides of the codec, not only between input and message.

## What remains

The full evaluator pad is a product over cut codebooks and later gate
tensors. The finite bundle graph shows that one proper coset does not
require the product on its own. A complete compression must handle a
composition of bundles, the middle and inverse transforms, and the
carry constraints. Those are not addressed here.

The next bounded diagnostic is a two-bundle joint image: compose two
`p=2` codecs with disjoint and with overlapping input digit sets, and
measure the joint image against the product of their ranks. If the image
grows submultiplicatively for overlapping inputs, the graph table can be
composed; if it is the full product for disjoint inputs, the compression
is local to shared inputs.

## Reproduction

```sh
python3 tools/check_bundle_joint_compression.py \
  --output evidence/bundle-joint-compression.json
make verify-bundle-joint-compression
```

The checker enumerates the source inputs at `h=1,2`, computes both codecs,
and asserts that the joint image equals the bundle rank and the source
input count.
