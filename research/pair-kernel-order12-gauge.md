# Pair-kernel order-twelve deformations at q8

## Result and scope

The [colour-asymmetric order-twelve floor](colour-asymmetric-order12-floor.md)
closes reweighting of the distinct-colour quartic gates.  The next local
escape is to change the pair kernel itself at order twelve while keeping
its order-six leading block.  The valuation grading below proves that
**arbitrary** order-twelve pair-kernel deformations cannot lower the
protected order-twelve floor.  The note also records finite q8 controls
for seven rank-one supports; these are not a uniform gauge theorem for
the exact low prefix.

Take the canonical spatial-mode pair matrix `P` and its exact leading
order-six block `P_6`.  The nonzero row space of `P_6` is

```math
R=\operatorname{span}\{(1,1,0,0),(0,0,1,1)\}.          \tag{1}
```

With `p_6=P_{0,2}` and `x` the cut variable, set

```math
q=p_6^2(1+x).                                         \tag{2}
```

The polynomial `q` has exact valuation twelve at the cut and endpoint
zero, so adding it to selected pair entries preserves the endpoint and
the order-six block.  Seven deterministic rank-one support patterns
`u v^T` are tested: the empty pattern, the three nonzero supports with
both `u` and `v` in `R`, one pattern with neither factor in `R`, and the
two patterns with exactly one factor in `R`.  For both cuts and all eight
charges the canonical low prefix is

```math
0^1,\ 6^{14},\ 12^{35},\ 15^{84},\ 16^{70},\ 17^{14},\ 18^{208}.
                                                               \tag{3}
```

None of the seven deformations lowers the order-twelve increment below
`35`.  The three rank-one supports with both factors in `R` leave the
entire prefix (3) unchanged.  The two patterns with exactly one factor in
`R` keep order-twelve at `35` but move fourteen channels from order
sixteen to orders fifteen and seventeen.  The support with neither
factor in `R` raises order-twelve to `42`.

## Valuation-grading stability

Keep the hypotheses of the colour-asymmetric theorem: the spatial-mode
pair kernel, the fixed four-parent leg grouping, distinct-colour XOR-zero
quartic gates of valuation at least twelve, and the endpoint and cold
ranks.  Add an arbitrary pair-kernel term of valuation at least twelve at
each cut.

The order-twelve `P=Q` block in total occupation four is a sum of
products of two pair entries.  The unperturbed pair kernel has leading
valuation six, so a product using one new order-twelve entry has
valuation at least `12+6=18`; a product using two has valuation at least
twenty-four.  It therefore does not change the order-twelve
occupation-four block.  A distinct-colour quartic gate cannot act on the
`P=Q` block, which uses only two colours.  Hence the rank-four Gaussian
block `B` is unchanged, and the same argument leaves the `q-1`
same-colour occupation-four channels unchanged.

The new pair terms can act inside the lower occupation-two sector at
order twelve.  Those channels are disjoint from the occupation-four
block and can only increase the total low prefix.  They cannot cancel
`B`.  Consequently the floor `5(q-1)` of the colour-asymmetric theorem
persists under arbitrary pair-kernel terms of valuation at least twelve.
At q8 the protected-prefix floor `5/7` also persists: the extra channels
have order at most twelve, and the existing profile certificate may leave
them in its conservative tail.  This stability does not cover
pair-kernel terms below valuation twelve, which can change the order-six
leading maps or the `P=Q` block directly.

## Pair-block gauge observation

On the pair-generated block alone, a deformation of the form

```math
P\longmapsto P+q\,u v^T,\qquad v=w^T P_6,\quad v\in R, \tag{4}
```

is a left multiplication

```math
P+q\,u v^T=(I+q\,u w^T)P.                             \tag{5}
```

The matrix `I+q\,u w^T` is unipotent and increases the filtration by
twelve.  It is therefore a filtration-preserving change of the
one-particle basis, so it cannot change the pair-generated Smith new
ranks at orders at most twelve.  This explains the pair-block part of the invariance of the
three all-in-`R` patterns.

The quartic contribution is not covered by the pair-block calculation.
The finite q8 control shows that it remains unchanged for the three
rank-one supports with both `u` and `v` in `R`.  A support with exactly
one factor in `R` shifts fourteen quartic channels from order sixteen to
orders fifteen and seventeen; a support with neither factor in `R`
raises order twelve.  The control does not establish a uniform
invariance for every quartic completion or every support matrix.

## Finite control

| pattern | support | order 6 | order 12 | low prefix |
| --- | --- | ---: | ---: | --- |
| canonical | none | 14 | 35 | (3) |
| all cross | `u=v=(1,1,1,1)` | 14 | 35 | (3) |
| rows 0,1 all columns | `u=(1,1,0,0)`, `v=(1,1,1,1)` | 14 | 35 | (3) |
| all rows, columns 2,3 | `u=(1,1,1,1)`, `v=(0,0,1,1)` | 14 | 35 | (3) |
| single entry | `u=e_0`, `v=e_2` | 14 | 42 | `0^1,6^14,12^42,15^84,16^56,17^14,18^248` |
| row 0, columns 2,3 | `u=e_0`, `v=(0,0,1,1)` | 14 | 35 | `0^1,6^14,12^35,15^91,16^56,17^21,18^208` |
| rows 0,1, column 2 | `u=(1,1,0,0)`, `v=e_2` | 14 | 35 | same as previous row |

The [checker](../tools/check_pair_kernel_order12_deformations.py) rebuilds
both q8 cuts and all eight charges, asserts the canonical prefix (3), and
asserts the tabulated result for every pattern.  The
[receipt](../evidence/pair-kernel-order12-deformations.json) records all
fourteen axis-pattern cases.

## What is and is not proved

- The valuation-grading argument excludes arbitrary pair-kernel
  deformations of valuation at least twelve from lowering the protected
  order-twelve floor.  It does not claim the exact finite low prefix is
  gauge-invariant.
- The seven rank-one support patterns are a finite q8 illustration for
  one endpoint-zero order-twelve polynomial `q`.  The pair-block gauge
  calculation (5) is exact under its hypotheses, but the quartic block
  and the exact prefix variation require the finite check.
- Pair-kernel terms below valuation twelve, a different leading pair
  kernel, a different leg grouping, and other alphabets remain open.  No
  pressure or factoring exponent is changed; the complete score remains
  `1/2`.

## Reproduction

```sh
make verify-pair-kernel-order12-deformations
```

The target runs the checker and `make verify-docs`.
