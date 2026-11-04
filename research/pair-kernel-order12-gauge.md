# Pair-kernel order-twelve deformations at q8

## Result and scope

The [colour-asymmetric order-twelve floor](colour-asymmetric-order12-floor.md)
closes reweighting of the distinct-colour quartic gates.  The next local
escape is to change the pair kernel itself at order twelve while keeping
its order-six leading block.  This note records a bounded q8 control for a
specific family of such deformations.  It is finite evidence, not a
uniform gauge theorem.

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

- This is a finite q8 control for one order-twelve endpoint-zero
  polynomial `q` and seven support patterns.  It is not an all-kernel or
  all-deformation theorem.
- The pair-block gauge calculation (5) is exact under its hypotheses, but
  the quartic block and the full Smith prefix require the finite check.
- A non-gauge pair-kernel deformation, a different leading pair kernel,
  lower-valuation deformations, and other alphabets remain open.  No
  pressure or factoring exponent is changed; the complete score remains
  `1/2`.

## Reproduction

```sh
make verify-pair-kernel-order12-deformations
```

The target runs the checker and `make verify-docs`.
