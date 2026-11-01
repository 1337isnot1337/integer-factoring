# A uniform all-field barrier for symmetric binary-energy shells

## Statement

For every integer `R>=1`, take `s+1` binary coordinates of weight
`2^s` for each `0<=s<=R`, and let

```math
\mathcal S_R=\{x:\sum_i w_i x_i=2^R\},\qquad
D_R=|\mathcal S_R|.
```

Suppose a crossing decomposition on this shell has pieces invariant
under independently permuting the coordinates of each equal-weight
class. Over every field, even a field depending on `R`, some single
piece has both flattening ranks at least

```math
D_R^{\,3/2-o(1)}.                                    \tag{1}
```

Thus this entire invariant shell family has `rho>=3/4` in the
width-to-factoring interface. It cannot reach the half-rate threshold
`1/sqrt(2)`. The hypothesis concerns the pieces, not just the crossing
sum: arbitrary pieces on the shell need not have this symmetry.

For the specific characteristic-two reciprocal pencil, the independent
[projected-block proof](weighted-energy-reciprocal-pencil.md) gives
the stronger restricted-family optimum `rho=1`. The present theorem
has broader scope: it allows all coefficient fields and every tensor
with the stated class symmetry, not only that pencil's coefficient
span. Neither theorem is a factoring lower bound or a new factoring
exponent.

The proof applies the [all-field middle-layer theorem](all-field-binary-word-barrier.md)
to an explicit product of balanced layers inside the shell. The
embedding preserves almost all leading label entropy.

## 1. Reserve one coordinate of every weight

Put `h=2^R` and, for sufficiently large `R`, choose

```math
L=R-2\lceil\log_2(R+2)\rceil-2.
```

In each class `0<=s<=L`, reserve one coordinate and select

```math
q_s=4\lfloor s/4\rfloor
```

of the other coordinates as free coordinates. This is possible because
`q_s<=s`, while the full class has `s+1` coordinates. Require each
label to have exactly `q_s/2` ones among those free coordinates.
All remaining nonreserved coordinates are fixed to zero.

The free coordinates contribute the same energy for every label:

```math
E=\sum_{s=0}^L2^s q_s/2 < L2^L
 \le R2^L\le \frac{R}{4(R+2)^2}h<h.                 \tag{2}
```

Use the one reserved coordinate in each class `0,...,R` to encode
the ordinary binary expansion of `h-E`. Fix these reserved bits
identically for every label. The resulting label has energy exactly
`h`. Distinct choices in the free balanced layers give distinct labels.

Consequently the shell contains a domain isomorphic to

```math
\mathcal L=\prod_{s=0}^L
 \{x\in\{0,1\}^{q_s}:|x|=q_s/2\}.                  \tag{3}
```

Its free-coordinate count is

```math
M=\sum_{s=0}^Lq_s
 =\tfrac12R^2-O(R\log(R+2)).                         \tag{4}
```

Indeed `s-3<=q_s<=s`, and `L=R-O(log R)`.
Empty free classes can be omitted.

The strict inequality in (2) is used only for sufficiently large
`L>=4`. The finitely many smaller `R` do not enter the asymptotic
statement, and the checker treats them separately.

## 2. A rank certificate supported on that domain

Restrict every tensor leg to the embedded domain (3). Both matrix
ranks can only decrease. The crossing becomes the crossing on
`mathcal L`, because the embedding is injective. Every restricted
piece remains invariant under permutations of the free coordinates
within each class: these permutations are a subgroup of the original
equal-weight class permutations, fixing the reserved and unused bits.

The product-layer theorem supplies a single linear functional on
these restricted tensors. Its value on the crossing is one. A piece
detected by this functional has both ranks at least

```math
2^{\,3M/2-O(\sum_{s=0}^L\log^2(q_s+4))}
 =2^{\,3R^2/4-O(R\log^2(R+2))}.                     \tag{5}
```

The constants are independent of the coefficient field. The same
functional detects the same summand in both cuts and across all
position classes; one does not select a different summand for each
class. Its polynomial support is balanced in every tensor leg, which
is why the restriction is valid.

## 3. Compare with the actual shell size

There are `(R+1)(R+2)/2` coordinates in the full shell description.
Thus

```math
\log_2 D_R\le (R+1)(R+2)/2=R^2/2+O(R).              \tag{6}
```

The central binomial coefficient is the largest coefficient of
`(1+z)^(q_s)`, giving

```math
\log_2 D_R\ge\log_2|\mathcal L|
 \ge M-\sum_{s=0}^L\log_2(q_s+1)
 =R^2/2-O(R\log(R+2)).                              \tag{7}
```

In particular `log2 D_R=R^2/2+o(R^2)`. Dividing the logarithm of
(5) by (6), and using (7) to control the scale, proves (1). More
explicitly, each detected rank is at least

```math
D_R^{\,3/2-O(\log^2(R+2)/R)}.
```

For the sum of the two logarithmic channel widths this is
`3 log2 D_R-o(log D_R)`, hence `rho>=3/4`.

## Verification and open scope

`tools/check_binary_energy_symmetry.py` checks the integer cutoff,
energy budget, reserved-bit embedding, and product-layer size bound.
It independently counts small shells both by their binary coordinates
and by binary partitions. The all-field rank certificate is checked
in `tools/check_all_field_binary_word_barrier.py` and proved in its
accompanying note. These finite checks do not substitute for the
asymptotic embedding and rank argument above.

A useful new shell construction must break the stated equal-weight
coordinate symmetry, change the shell, or use a different evaluation
interface. The low number of coefficient degrees alone does not
supply low matrix ranks.
