# The exponent frontier for rectangular grids

## Result and scope

The limiting Appendix B channel profile has an exact exponent frontier at
`1/2`. The theorem permits unequal splits chosen separately at each
recursion node within the stated limiting model.

The result is a minimax theorem about a specified analysis model. It is
**not** a lower bound for factoring, matrix rank, or the existing evaluator.
Smaller factorizations, cancellations, and new ways to evaluate a term can
escape this model.

The source is `anthropics/formal-math`, revision
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Appendix H, display (77), gives the channel profile used here.
The source declarations are:

- [`AppHalfPieces.Rrow` and `Rcol`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
- [`GeneralPieces.Bits.Rrow_le`, `Rcol_le`, and `bits_le`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/Bits.lean).
- [`ClauseBound.Alg.totalG`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Alg/Data.lean),
  which counts the row and column volumes on a rectangular grid.

No Lean source changes or builds accompany this note.

## 1. The channel profile contains more information than its sum

Let `h` be the least integer with `p(h) >= 2^lambda`, where `p` counts
ordinary partitions. The Appendix B piece index satisfies `-2h <= s <= 2h`.
Its declared channel counts are

```math
R_r(h,s)=\sum_{i=0}^{2h}p(s+i),\qquad
R_c(h,s)=\sum_{i=0}^{2h}p(2h-s-i).
```

As in the source, `p(j)=0` for negative integers. The widths are the
ceiling binary logarithms of these counts, with the optional trivial
factorization at width `2lambda`.

The largest summand bounds each sum from below. The number of summands
bounds it from above:

```math
p(2h+s)\le R_r(h,s)\le(2h+1)p(2h+s),
```

and the corresponding inequality holds with `2h-s` for `R_c`.
The source's partition estimates give

```math
\log_2p(E)=\kappa\sqrt E+O(\log(E+2)),\qquad
\lambda=\kappa\sqrt h+O(\log(h+2)).
```

For any fixed `t` in `(-1,1)`, choose integer indices with `s/(2h) -> t`.
The normalized declared widths then tend to

```math
r=\sqrt{\frac{1+t}{2}},\qquad
c=\sqrt{\frac{1-t}{2}},\qquad r^2+c^2=1.
\tag{1}
```

The normalization divides each width by `2lambda`.
The endpoints follow by a limit. Every point of the positive quarter
circle occurs as a limiting declared profile. This concerns the supplied
factorization lengths. It makes no claim that those lengths equal ranks.

On a `p` by `d` grid, set `N=(pd)^2 lambda`. The exact volume formulas are

```math
U_r=p(p-1)d\lambda+2\sum_q c_r(q),\qquad
U_c=d(d-1)p\lambda+2\sum_q c_c(q).
\tag{2}
```

There are `binom(p,2)binom(d,2)` rectangles `q`.
When both dimensions and `lambda` tend to infinity, a constant profile
gives `U_r/N -> r` and `U_c/N -> c`.
For varying profiles, the average vector lies in the unit quarter disk.
This follows from the triangle inequality for the Euclidean norm.
The approximation requires both grid dimensions to grow.

For finite dimensions, let `q=(1-1/p)(1-1/d)`. A constant normalized
profile gives the exact geometry

```math
\frac{U_r}{N}=\frac{1-1/p}{d}+q r,\qquad
\frac{U_c}{N}=\frac{1-1/d}{p}+q c.
\tag{2a}
```

Replacing this affine profile by the unit circle removes both positive
straight-edge terms and the factor `q<1`. These changes have opposite
effects. Thus the unit-circle game below is not a universal relaxation
of every finite rectangular grid. In particular, a fixed side of two
vertices falls outside the limiting approximation.

## 2. One rectangular split: an exact minimax identity

For positive real costs `A,B`, define

```math
F(A,B)=\sup_{r,c\ge0,\ r^2+c^2\le1}\min(Ar,Bc).
```

Then

```math
F(A,B)=\frac{AB}{\sqrt{A^2+B^2}}
       =\frac1{\sqrt{A^{-2}+B^{-2}}}.
\tag{3}
```

**Proof.** If both products exceed `z`, then
`r > z/A` and `c > z/B`. Hence
`1 >= r^2+c^2 > z^2(A^(-2)+B^(-2))`.
This proves the upper bound in (3). Equality holds at

```math
r=\frac{B}{\sqrt{A^2+B^2}},\qquad
c=\frac{A}{\sqrt{A^2+B^2}}.
```

Thus the identity is sharp throughout the allowed profile set.

Take `a=log p`, `b=log d`, and propose a recursive coefficient
`C/(log K)^alpha`. Grant free grouping, no rounding, and destination
sizes exactly `p,d`. The costs of the two orientations are

```math
CN\frac r{a^\alpha},\qquad CN\frac c{b^\alpha}.
```

Applying (3) gives the sharp worst-case coefficient

```math
\frac{C}{\sqrt{a^{2\alpha}+b^{2\alpha}}}.
\tag{4}
```

For `alpha=1/2`, this is exactly

```math
\frac{C}{\sqrt{a+b}}
=\frac{C}{\sqrt{\log(pd)}}.
\tag{5}
```

Every aspect ratio has the same leading potential in this model.
The square grid is one choice that achieves it.

If `alpha>1/2`, then `a,b>0` imply

```math
a^{2\alpha}+b^{2\alpha}<(a+b)^{2\alpha}.
```

The coefficient in (4) therefore exceeds the proposed parent coefficient.
The extremizing profile has

```math
r^2=\frac{a^{2\alpha}}{a^{2\alpha}+b^{2\alpha}},\qquad
c^2=\frac{b^{2\alpha}}{a^{2\alpha}+b^{2\alpha}}.
```

It is an allowed limiting declared profile by (1).
A favorable profile can give a faster step. A uniform proof must also
cover this unfavorable profile.

## 3. Arbitrary split trees have the same frontier

The following game retains the limiting unit-circle profile and allows
arbitrary unequal splits. It ignores all construction costs. The conclusion
applies to this game, rather than every finite grid in (2a).

Each internal node has a parameter `K=pd`, where `p,d>=2`.
Its children have cost coefficients `A,B`, measured per unit volume.
The profile can be any point in the quarter disk from (1).
The algorithm chooses the less expensive orientation. Equation (3)
gives the node's worst coefficient.

Define the reciprocal square `H=1/F^2`. Equation (3) becomes

```math
H_{\mathrm{parent}}=H_{\mathrm{row}}+H_{\mathrm{column}}.
\tag{6}
```

Consequently, for every finite binary tree,

```math
F_{\mathrm{root}}^{-2}
=\sum_{v\text{ a leaf}}F_v^{-2}.
\tag{7}
```

This identity does not depend on the tree's balance or depth.
The proof is induction on the tree, using (6) at its root.

Suppose every leaf parameter is at least `K0>1`, and every leaf
coefficient is at least `b0>0`. The product of all leaf parameters equals
the root parameter `K`. If the number of leaves is `m`, then

```math
m\le\frac{\log K}{\log K_0},\qquad
F_{\mathrm{root}}^{-2}\le\frac m{b_0^2}.
```

Therefore

```math
F_{\mathrm{root}}
\ge b_0\sqrt{\frac{\log K_0}{\log K}}.
\tag{8}
```

Brute-force contraction on a clique of at least two vertices has
coefficient `(K-1)/(2K) >= 1/4` under the normalization `N=K^2 lambda`.
Thus the usual terminal rule meets the positive-leaf assumption.
Using fewer groups than available vertices only reduces the child
parameters. Such a restriction does not improve this game.

The same lower bound holds when the algorithm chooses each split from
previous profiles. At each internal node the available child subgames
obey (8). The profile that equalizes their two costs gives (3).
Induction applies regardless of the selected split. The argument requires
termination and forbids parameter-one splits, which contain no rectangles.

This is a barrier for proofs that replace the declared profile by its
unit-circle limit and charge only these recursive choices.
It does not cover arbitrary recursions through bounded grid dimensions.
It also makes no claim that one nonzero factoring instance realizes
all adversarial choices simultaneously.

## 4. Repeated two-row splits also retain exponent one-half

The finite geometry in (2a) needs a separate argument. This section covers
repeated splits with `p=2`. It does not assert a theorem for every finite
aspect ratio.

First let `d` tend to infinity. Let `A>0` be the fixed coefficient for
the two-vertex child, and let `0<B<=A` be the column-child coefficient.
The limiting minimax is

```math
G(A,B)=\sup_{r,c\ge0,\ r^2+c^2\le1}
\min\left(\frac{Ar}{2},\frac{B(1+c)}2\right)
=\frac{B}{1+(B/A)^2}.
\tag{9}
```

To prove this, put `z=B/(1+(B/A)^2)` and choose
`r=2z/A`, `c=2z/B-1`. These numbers are nonnegative and satisfy
`r^2+c^2=1`. Both costs equal `z`.
If both costs exceeded `z`, both coordinates would exceed these values.
Their squared norm would then exceed one. This proves (9).

For repeated limiting steps, set `x_j=A^2/B_j^2`.
Equation (9) gives the exact recurrence

```math
x_{j+1}=x_j+2+1/x_j.
```

Since `x_0>=1`, induction gives

```math
x_0+2j\le x_j\le x_0+3j.
\tag{10}
```

Thus `B_j=Theta(j^(-1/2))`. A two-row split doubles the vertex parameter,
so this remains a square-root logarithmic saving.

Now retain finite `d>=2`. The exact profile is

```math
\left(\frac1{2d}+\frac{1-1/d}{2}r,
      \frac{1-1/d}{2}(1+c)\right).
```

At the extremizer from (9), the column cost equals `(1-1/d)z`.
The row cost is `(1-1/d)z+A/(2d)`. Therefore its minimax `G_d` obeys

```math
G_d(A,B)\ge(1-1/d)\,G(A,B).
\tag{11}
```

For a chain with `d_j>=2^(j+1)`, let `B_(j+1)=G_(d_j)(A,B_j)`.
The row coordinate never exceeds `1/2`, so `B_(j+1)<=A/2`.
Consequently `x_j>=1`, and (11) gives

```math
x_{j+1}\le(1-1/d_j)^{-2}(x_j+3).
```

Unrolling this recurrence yields

```math
x_j\le\left(\prod_{i<j}(1-1/d_i)^{-2}\right)(x_0+3j)
\le16(x_0+3j).
\tag{12}
```

For the last inequality, the first factor of
`prod_(k>=1)(1-2^(-k))` is `1/2`. The remaining product is at least
`1-sum_(k>=2)2^(-k)=1/2`.
The finite-product inequality follows by induction from
`(1-u)(1-v)>=1-u-v` for nonnegative `u,v` with `u+v<=1`.
Every finite prefix product is therefore at least `1/4`.

Equation (12) implies

```math
B_j\ge\frac{A}{4\sqrt{x_0+3j}}.
\tag{13}
```

This closes the fixed-side-two loophole for this finite affine-profile
game. It still supplies no runtime lower bound for the factoring machine.

## 5. What can evade the barrier

A stronger rate needs information absent from this model. Concrete
possibilities include smaller factors for the same pieces, a narrower
joint profile, or a method that combines several terms before recursion.
Another possibility is a third evaluation method whose cost is lower
on the equalizing profiles.

The theorem also identifies the useful quantity for a new construction.
If its normalized width vectors lie in a smaller disk of radius `eta<1`,
then (3) gains the factor `eta`. Repeated square steps then target

```math
\alpha=\frac12+\log_2(1/\eta).
```

This last formula is conditional on a uniform construction with the
required size and operation bounds. No such smaller disk is proved here.

## 6. Reproduction

Run the exact arithmetic checker:

```sh
python3 tools/check_rectangular_profile.py
```

It verifies rational extremizers, the square-root potential, the tree
identity, and the two-row recurrence. Its negative controls use
exponents `3/4` and `1`.
All comparisons use integers or rational numbers. The checker supports
the ordinary proofs above and does not test an asymptotic factoring claim.

The minimax analysis and checker used substantial OpenAI assistant support.
The original pieces and counting engine retain their source attribution.
