# A colour-asymmetric order-twelve floor for the spatial-mode pair kernel

## Result and scope

The [spatial-mode fusion construction](spatial-mode-colour-fusion.md)
uses distinct-colour XOR-zero quartic gates and attains the order-twelve
prefix `0^1,6^(2d),12^(5d)` for `d=q-1`.  This note tests the natural
escape of using **colour-asymmetric** quartic weights, including weights
that depend on the assigned colours and on the support.  They do not
lower the order-twelve block.

**Theorem.**  Fix the spatial-mode pair kernel and the fixed four-parent
leg grouping used in the q8 construction.  Complete it by any creation
gates supported on four **distinct** nonzero colours with XOR zero, of
valuation at least twelve at both cuts, with arbitrary coefficients in
the characteristic-two field.  Then at each cut the new order-twelve
Smith multiplicity is at least

```math
(q-1)+4(q-1)=5(q-1).                                  \tag{1}
```

The canonical spatial-mode completion attains equality.  At `q=8` the
floor is `35`.

Combining this with the [existing q8 leading-prefix
floor](q8-leading-prefix-floor.md) gives a scoped corollary: every such
completion at `q=8` has common-scale, unweighted paired-type pressure at
least `5/7>1/sqrt(2)`.  Thus recolouring the quartic gates cannot rescue
this fixed q8 pair kernel.  The result does not cover a different pair
kernel, a different leg grouping, gates below valuation twelve, other
alphabets, weighted scales or unequal scales.  It proves no factoring
lower bound; the complete score remains `1/2`.

## 1. The charge blocks at order twelve

Write `q=2^s>=4` and fix one cut.  Charge is total colour XOR.  The
two-occupied sector has two colours on each side.  For a nonzero charge
`c`, its unordered pairs are

```math
\mathcal P_c=\{\{a,b\}:a,b\ne0,\ a\ne b,\ a\oplus b=c\},
\qquad |\mathcal P_c|=(q-2)/2.                        \tag{2}
```

For a pair `P` on the row side and `Q` on the column side, let `M_c(P,Q)`
be the `12 by 12` spatial block of the exact order-twelve coefficient.
For `P=Q`, a quartic gate would need four distinct colours on its support.
The row and column vertices together use only the two colours in `P`, so
the quartic contribution is zero.  The block is therefore the Gaussian
block `B` computed in the spatial-mode note.  That block has rank four:
it is the tensor product of the two rank-two spatial leadings.

We do not need the off-diagonal blocks.  Charge conservation makes `M_c`
block diagonal over `c`, and the `P=Q` block is a literal submatrix of
`M_c`.  Hence

```math
\operatorname{rank} M_c\ge \operatorname{rank} B=4,\qquad c\ne0.
                                                               \tag{3}
```

The one-occupied order-six sector has total occupation two, while `M_c`
has total occupation four.  The total-occupation grading separates them,
so these four channels are new at order twelve rather than rows already
spanned at order six.

Charge zero is the repeated-colour sector: an unordered pair with XOR zero
has equal colours.  A distinct-colour quartic gate cannot act on a support
carrying only one repeated colour, so the `q-1` same-colour Gaussian
channels are unchanged.  Therefore the new order-twelve total is at least

```math
(q-1)+\sum_{c\ne0}4=5(q-1).                            \tag{4}
```

For the canonical completion every off-diagonal block is the same `B`,
so `M_c=J_{(q-2)/2}\otimes B` has rank four.  Charge zero contributes
`q-1`, and equality holds.

## 2. The q8 consequence

At `q=8`, `d=7`.  The one-particle sector is fixed by the Gaussian pair
kernel and contributes exactly `2d=14` order-six channels, two on each
nonzero charge.  The theorem gives at least `5d=35` further channels at
order twelve.  The [q8 leading-prefix
theorem](q8-leading-prefix-floor.md) assumes precisely a zero order, at
least fourteen orders at most six, and at least thirty-five further
orders at most twelve; it does not require the selected orders to be the
entire low prefix.  It therefore applies with its `5/7` floor.

The colour-asymmetric escape from the [colour-alternation
bound](colour-alternation-joint-cold-curve.md) is thus closed for this
fixed q8 pair kernel, not merely for the symmetric choice.  To reduce the
order-twelve floor one must reduce the rank-four diagonal Gaussian block
`B` or change the order-six one-particle sector.  Both require changing
the pair kernel's leading spatial maps, not reweighting the quartic gates.
Moreover, an arbitrary pair-kernel term of valuation at least twelve
cannot lower (1): in total occupation four it appears only in products of
valuation at least eighteen, so it cannot change the order-twelve
occupation-four block.  The finite order-twelve deformation controls are
in the [pair-kernel note](pair-kernel-order12-gauge.md).  Terms below
valuation twelve or a different order-six leading pair map are the next
assumptions to test.

## 3. Finite controls

```sh
make verify-colour-asymmetric-order12-floor
```

The [checker](../tools/check_colour_asymmetric_order12_floor.py) rebuilds
the q8 charge blocks for both cuts and all eight charges.  For each case
it computes the Smith prefix under four completions: the canonical
completion, the Gaussian-only block, and two deterministic asymmetric
subsets of the distinct-colour gates.  The canonical new order-twelve
increment is `4` on each nonzero charge and `7` on charge zero.  The
Gaussian-only value is `12` on each nonzero charge, and every sampled
asymmetric subset is at least `4`.  The checker asserts the floor and the
canonical equality; the [receipt](../evidence/colour-asymmetric-order12-floor.json)
records all sixteen cases.  The proof covers arbitrary field
coefficients; the finite subsets only exercise the block inequality.

## 4. Reproduction

```sh
make verify-colour-asymmetric-order12-floor
```

The target runs the checker and `make verify-docs`.
