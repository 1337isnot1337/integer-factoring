# Component rule for a shared pair plus an independent bundle

## Result and scope

A two-bundle shared-input pair has a graph image of one pair per input.
Adding one independent bundle multiplies that graph by the independent
bundle rank. At `h=2` the joint image has `65,536 * 65,536 = 2^32`
entries instead of the full three-bundle rank product `2^48`; the
compression factor is `65,536`. At `h=1` the factor is `4`, checked by
exact enumeration.

This is a finite codec diagnostic. It supplies the composition rule for
one shared component and one independent component, not a general
digit-sharing graph, a full two-cut pad theorem, or a factoring rate.
The complete score remains `1/2`.

## Mechanism

Let the shared pair on digits `x` have graph image `G`, and let the
third bundle on independent digits `y` have rank `R_3`. The joint map is

```math
(x,y)\mapsto (R_1(x),C_1(x),R_2(x),C_2(x),R_3(y),C_3(y)).
```

For fixed `x`, the third coordinate ranges over the independent rank;
for fixed `y`, the first two coordinates range over the shared graph.
Because the digit sets are disjoint, the image is the Cartesian product
`G x R_3`, not the full rank product. This is the component rule.

Combined with [the two-bundle split](two-bundle-joint-image.md), it says
a pad analysis can group bundles into connected components of the
digit-sharing graph: each component contributes its graph size, and
independent components multiply. A uniform claim that every boundary
compresses jointly would be false, as the disjoint case shows.

## What remains

The source's actual operand/cut graph must be partitioned into connected
components. Within a component, the graph may be larger than one input
per bundle; the two-bundle and three-bundle tests only cover the
identity-sharing extreme. A general component bound is the next
structural target.

## Reproduction

```sh
python3 tools/check_three_bundle_joint_image.py \
  --output evidence/three-bundle-joint-image.json
make verify-three-bundle-joint-image
```

The `h=1` row is exact by enumeration; the `h=2` row uses the component
product, since enumerating `2^32` pairs is unnecessary.
