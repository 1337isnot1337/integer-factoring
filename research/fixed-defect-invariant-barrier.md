# A rank obstruction for permutation-invariant defect families

## General result

The obstruction below is uniform even when the number of defects grows.
Let `b,q` be positive integers and `t>=6b`. A label is a `b`-element
subset of `[t]`, with one of `q` colors assigned to each element. Put
`D=binom(t,b)q^b` and `s=floor(t/(2b))`.

Every exact crossing decomposition into simultaneously `S_t`-invariant
summands has a summand whose ranks in both cuts are at least

```math
(s-2)^{2b}.
\tag{U}
```

The field, coefficients, and number of summands may depend on `t,b,q`.
For fixed `q` and `t/b` tending to infinity, (U) forces a total binary
width of `(4-o(1))log2 D`. Thus no fixed contraction coefficient `rho<1`
is possible in this sparse-defect regime. Section 7 proves (U) by a
zero-sum projection that isolates one crossing-matching coefficient.

Here is a complementary finite bound that is particularly transparent
when `b` is fixed.

Fix positive integers `b,q`. A label is a `b`-element subset of `[t]`,
with one of `q` colors assigned to each of its elements. There are
`D=binom(t,b)q^b` labels. The permutation group `S_t` acts on the positions.
For `t>=4b`, every exact crossing decomposition into simultaneously
`S_t`-invariant summands has one summand whose ranks in both cuts are at
least

```math
\binom{t}{2b}-B_{4b}q^{2b}t^{2b-1},
\tag{0}
```

where `B_(4b)` is the number of set partitions of `4b` objects.
In particular both ranks are `Omega_(b,q)(t^(2b))=Omega_(b,q)(D^2)`.
No fixed contraction coefficient `rho<1` is possible for this invariant
family. The field, the coefficients, and the number of summands may
depend on `t`.

Sections 1--4 give the simpler single-defect proof with a sharper finite
constant. Section 6 proves (0) by lifting unordered labels to ordered
injective tuples. The lift only duplicates rows and columns and does not
divide by `b!`; the proof therefore includes characteristic two.

## The single-defect case

Fix an internal label dimension `q>=1`. Let the labels be pairs
`(position,color)` in `[t] x [q]`, and let `S_t` act by permuting the
position while leaving the color fixed. For every `t>=4`, consider an
exact crossing decomposition over any field:

```math
E_{qt}(x,y,u,v)=[x=v][y=u]=\sum_i A_i(x,y,u,v).
```

Suppose each summand is invariant under the same simultaneous position
permutation on all four coordinates. Then some summand satisfies

```math
\begin{aligned}
\operatorname{rank}_{(x,y)\mid(u,v)}A_i
&\ge {t\choose2}-13q^2t,\\
\operatorname{rank}_{(x,u)\mid(y,v)}A_i
&\ge {t\choose2}-13q^2t.
\end{aligned}
\tag{1}
```

For fixed `q` and sufficiently large `t`, both ranks are `Omega(t^2)`.
Their product is `Omega(t^4)`. Thus this invariant family cannot give
any contraction coefficient `rho<1` at label length
`lambda=log2(qt)+O(1)`. The conclusion allows arbitrary invariant
recombinations and arbitrary numbers of summands. The field and the
summands may depend on `t`.

This addresses a tempting endpoint of the
[three-label amplification](finite-seed-amplification.md): retain words
of length `t` with exactly one coordinate in a two-dimensional space
`U` and all other coordinates equal to `w`. There are `2t` labels.
Every type coefficient is invariant under permuting the `t` positions,
so (1) applies with `q=2`. A constant number of nonzero types would make
this an effective scalable family, but its joint rank profile cannot
give a saving.

The fixed-color conclusion does not concern label restrictions that
break this position symmetry, a growing internal color dimension, or
a defect count comparable to `t`. In particular it is not a lower bound
for arbitrary crossing
decompositions on `qt` labels.

## 1. The fifteen equality diagrams

Temporarily write only the four position coordinates as `1=x`, `2=y`,
`3=u`, `4=v`. For a set partition `pi` of `{1,2,3,4}`, define `D_pi` to
be the indicator that position coordinates in each block of `pi` agree.
Different blocks are allowed to have equal positions as well.

There are fifteen such diagrams. They form a basis for invariant scalar
functions of the four positions when `t>=4`. Indeed, the actual equality
partition of a position tuple can be any of the fifteen partitions.
Evaluating the diagrams on these fifteen orbits gives the zeta matrix
of the partition lattice: `D_pi` is one precisely when `pi` refines the
actual equality partition. This matrix is triangular with diagonal one
in a linear extension of the refinement order. It is invertible over
every field, including characteristic two.

Therefore every invariant tensor has a unique expansion

```math
A=\sum_{\pi}D_\pi\otimes C_\pi,
\tag{2}
```

where `C_pi` is a four-color tensor on `[q]^4`.
These are equality diagrams, **not indicators of exact equality
orbits**. Using exact-orbit coefficients in the argument below would be
incorrect, since an exact-orbit indicator can expand into several
diagrams of different ranks.

Name the three pairing diagrams

```math
E=(14)(23),\qquad C=(13)(24),\qquad R=(12)(34).
```

The crossing tensor is `D_E` times the color crossing tensor. Its
coefficient `C_E` in (2) is nonzero, and all its other diagram coefficients
are zero. Consequently some summand of any exact decomposition has
`C_E!=0`.

## 2. A diagram's rank is controlled by its through-blocks

For a chosen flattening, a through-block is a block of `pi` that meets
both sides of the cut. If a diagram has `b` through-blocks, its position
matrix factors through at most `t^b` coordinates: assign a position to
each through-block, and check all equalities separately on the two sides.
Blocks lying entirely on one side contribute only a factor on that side.

The color matrix has rank at most `q^2`. Thus

```math
\operatorname{rank}(D_\pi\otimes C_\pi)\le q^2t^b.
\tag{3}
```

For the row cut `(1,2)|(3,4)`, exactly two diagrams have two
through-blocks: `E` and `C`. Every other diagram has at most one.
Hence the sum of the other thirteen terms in (2) has row rank at most
`13q^2t`.

For the column cut `(1,3)|(2,4)`, the only two diagrams with two
through-blocks are `E` and `R`. The corresponding remainder again has
rank at most `13q^2t`.

## 3. A nonzero crossing coefficient is quadratic in both cuts

Suppose `C_E!=0`. In the row cut, keep just the two leading terms

```math
L_{\rm row}=D_C\otimes C_C+D_E\otimes C_E.
```

Restrict the position pairs to distinct coordinates. For each unordered
pair `{a,b}` with `a!=b`, the two row orders `(a,b),(b,a)` and the two
column orders form a separate block. Including the color coordinates,
that block is

```math
\begin{pmatrix}
 (C_C)_{\rm row}&(C_E)_{\rm row}\\
 (C_E)_{\rm row}&(C_C)_{\rm row}
\end{pmatrix}.
\tag{4}
```

Its size is `2q^2` by `2q^2`. It is nonzero because `C_E` is nonzero, so
its rank is at least one. These `binom(t,2)` blocks have disjoint row and
column coordinates. It follows that

```math
\operatorname{rank}L_{\rm row}\ge {t\choose2}.
```

Subtracting a matrix of rank at most `13q^2t` can reduce the rank by at
most that amount. This proves the first inequality in (1).

In the column cut, the same argument uses `C_R` and the column
reshuffling of `C_E`. Reshuffling cannot turn a nonzero color tensor into
zero. The corresponding block is again nonzero, proving the second
inequality in (1). No positivity or semisimplicity is used.

Since the same coefficient `C_E` controls both arguments, the large row
and column ranks occur in the same summand.

## 4. Consequence for binary channel widths

For fixed `q`, the right side of (1) is at least a positive constant
times `t^2` once `t` is large. Every exact factorization of this summand
therefore needs binary width at least `2log2 t-O_q(1)` in each cut.
At `lambda=log2(qt)+O(1)` their sum is at least

```math
4\lambda-O_q(1).
```

It cannot satisfy `4rho lambda+o(lambda)` for any fixed `rho<1`.
For the two-color single-defect construction one can take
`t=2^(lambda-1)` to obtain exactly `2^lambda` labels, so no interpolation
or padding ambiguity is present.

The argument also applies after a fixed change of the one-label basis,
provided the resulting coefficient tensors retain the simultaneous
position symmetry. Recombining invariant coefficients preserves that
symmetry. Breaking the symmetry of the label domain is a separate
possibility left open by this theorem.

## 5. Exact checks

```sh
python tools/check_fixed_defect_invariants.py --output evidence/fixed-defect-invariants.json
```

The checker verifies the fifteen-diagram basis over several fields,
counts through-blocks in both cuts, checks all fifteen small position
matrix ranks, and replays the exact block-rank formula (4) with nontrivial
color coefficients. The theorem is the uniform argument above; it does
not extrapolate those finite ranks.

The checker also verifies the ordered two-defect lift in characteristic
two, the independence of its twenty-four full-through diagrams, and
the zero-sum projection used in the uniform refinement.
The source-specific [single-defect seed audit](../tools/check_one_defect_seed.py)
is separate. The present obstruction requires only invariance and the
exact crossing sum.

## 6. Any fixed number of defects

Replace each unordered colored label by an ordered list of its `b`
distinct positions, carrying the colors along with the ordering. Forgetting
the order is a surjection onto the original label set. Pull back every
tensor on each of its four legs along this map.

In either flattening this operation merely duplicates rows and columns.
It leaves rank unchanged: the original matrix is recovered by selecting
one canonical ordering for each unordered label, and the lifted matrix
is obtained from the original by row and column replication. This uses
no averaging and no factorial denominators.

There are now `4b` labeled position slots. Their actual equality
partitions are exactly those partitions in which no block contains two
slots of the same external leg. Call these partitions admissible. Every
admissible partition is realized when `t>=4b`. Refining an admissible
partition preserves admissibility, so the zeta matrix on admissible
partitions is still unitriangular. Every invariant pulled-back tensor
therefore has a unique expansion in the admissible equality diagrams,
with color tensors on `[q]^(4b)`.

For either cut there are `2b` slots on each side. A diagram with `h`
through-blocks has rank at most `q^(2b)t^h`, by the same separate-side
factorization as (3). A diagram with exactly `2b` through-blocks must
pair each row slot with exactly one column slot. Thus the leading
diagrams are precisely the `(2b)!` bijections between those slots.
All remaining terms have total rank at most the error term in (0).

The pulled-back crossing is the sum over permutations `sigma,tau` of
the diagrams pairing the `x` slots with the `v` slots according to
`sigma`, and the `y` slots with the `u` slots according to `tau`, together
with the corresponding color equalities. These `(b!)^2` diagrams are
distinct and have nonzero color tensors. The displayed sum is exact in
every field: because each label is injective, at most one permutation
can match any given pair of its ordered position lists. There is no
multiplicity of `b!` at an entry.

Choose one such crossing-matching diagram. By uniqueness of the
expansion, some summand has a nonzero color coefficient at this diagram.
It is a full-through diagram in both cuts.

In the row cut, restrict attention to ordered row and column pairs whose
combined `2b` positions are distinct. The sum of all full-through terms
is block diagonal according to the common `2b`-element support `W`.
Fix one ordered assignment of `W` to the row slots. Different bijections
send it to different ordered assignments of `W` to the column slots.
Consequently an entry of the chosen nonzero color coefficient cannot
cancel with any other full-through diagram. The block for every `W`
is nonzero and has rank at least one.

There are `binom(t,2b)` supports `W`. Subtracting the lower-through
remainder proves (0) in the row cut. The same chosen diagram is full
through in the column cut; its color tensor is only reshuffled. The
identical argument proves (0) there as well.

Finally, fixed `b,q` give `log2 D=b log2 t+O_(b,q)(1)`. The two binary
widths are each at least `2b log2 t-O_(b,q)(1)`, so their sum is
`4log2 D-O_(b,q)(1)`. This completes the fixed-defect theorem.

## 7. A projection without a Bell-number error term

Use the ordered lift and the unique admissible-diagram expansion from
Section 6. Choose one crossing-matching diagram whose color coefficient
is nonzero in a summand. There are `2b` matched pairs of slots. Allocate
to each pair a distinct position block of size `s=floor(t/(2b))`, and
restrict each slot to the block assigned to its pair.

All within-leg injectivity constraints now hold automatically, because
different slots of the same leg use different blocks. Any diagram that
equates slots from different blocks vanishes. Consequently the only
diagrams that can survive this restriction are refinements of the chosen
matching: each of its pairs may either remain joined or split into two
singleton blocks.

On each slot apply the zero-sum difference basis

```math
e_1-e_s,\ e_2-e_s,\ \ldots,\ e_{s-1}-e_s.
```

A split pair contributes independent constant functions on its slots;
the zero-sum projection annihilates it. A joined pair contributes the
Gram matrix `I_(s-1)+J_(s-1)`. This matrix has rank at least `s-2` over
every field: it is the identity plus a rank-one matrix.

Only the chosen matching diagram remains. In either crossing cut, all
`2b` of its pairs connect the two sides. The projected matrix is, up to
row and column permutations, a tensor product of `2b` such Gram matrices
and the corresponding flattening of the chosen nonzero color tensor.
The latter has rank at least one in both cuts. Hence both projected
ranks are at least `(s-2)^(2b)`.

Restriction and multiplication by the difference-basis matrices cannot
increase rank. The ordered lift preserves rank. These facts prove (U)
for the original summand in both cuts.

For fixed `q`, put `L=t/b` and suppose `L` tends to infinity. Elementary
binomial bounds give

```math
b\log_2 L\le\log_2 D
\le b\log_2 L+O_q(b).
```

Also `log2(s-2)=log2 L-O(1)`. The sum of the two necessary binary widths
is at least `4b log2 L-O(b)=(4-o(1))log2 D`, as claimed.
There is no factorial division or assumption about semisimplicity in
this argument.
