# A fixed q8 quartic projection retains a three-quarter profile floor

## Question, result and stopped scope

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363154)
keeps the cube pairing fixed and substitutes

```math
(X,Y,Z)=(x^2,y^2,xy)                                  \tag{1}
```

in the [degree-four additive-color quartic tensor](quartic-charge-fusion.md).
This changes the three-axis transfer and its natural scaling, while avoiding
the [invalid assignment-dependent regrouping](two-axis-cube-regrouping.md).
The old generic-Z-unit order-five argument does not describe this joint
X,Z-cold specialization. All doubled degree and paired-letter costs must
be restored before comparing with the two-axis target `1/sqrt(2)`.

**For q=8, the common-scale Smith-profile envelope of (1) is at least
3/4 for every common paired physical type, including correlated types.**
Thus this fixed curve does not cross the target. No complete profile,
larger alphabet, other curve or scale search follows this certificate.
The [stop checkpoint](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363322)
was recorded before packaging. This is an obstruction to the specified
profile estimate, not an individual coefficient-rank lower bound or a
factoring lower bound. The complete proven score remains **1/2**.

## 1. Fixed wiring, endpoint and degree

Group the two old z-neighbor vertices into each superletter indexed by
the remaining x,y bits. Each superletter has alphabet `{0,...,7}^2`, of
size 64. The antipode maps it to the opposite x,y superletter and swaps
its two internal colors. Swapping the two internal positions on the two
x=1 superlegs normalizes this endpoint to the four-leg crossing. This
is one fixed permutation for the whole tensor, before choosing any piece.
It preserves occupancy and changes neither cut's Smith profile.

At x=0 in (1), the old kernel has only its Y-direction matching, wholly
inside each x half. The quartic vanishes. Hence the cold x flattening has
rank one. The y statement is symmetric. The tensor is polynomial, and
each entry has degree eight per old occupied vertex after (1). In
particular its maximum entry degree is 64, not the old 32.

The Gaussian one-particle orders are `(6,6,18,18)` in each new cut.
These are the a=2,b=1 case of the existing
[fixed-projection formula](fixed-cube-projection-profile-floor.md#1-exact-one-particle-profiles-for-every-projection);
the checker recomputes the determinantal prefixes `(6,12,30,48)`.
The Gaussian-only [block-mode floor](block-gaussian-two-axis-profile-floor.md)
already excludes that baseline. The quartic is the non-Gaussian change
tested here; no new Gaussian projection is claimed.

## 2. Complete low jets at the fixed core

Work over `F_2[[x]]` after y=1. The row/column occupancy grading recovers
the generic matrix by unit powers of y on rows/columns and x replaced
by x/y. Thus these are exact generic cut orders, with no finite-field
specialization of a remaining parameter.

Remove all pair and quartic creation gates supported wholly on one side
of the cut. They commute with the other creations, and each is invertible
with polynomial inverse. This operation preserves Smith orders. Every
remaining pair has valuation at least six; every quartic has valuation
twelve. An entry with occupied count w therefore has valuation at least
3w. Through order eighteen only the vacuum, at most three cross pairs,
or one quartic and at most one cross pair can contribute. Two quartics
start at order twenty-four.

Total color XOR is conserved. Each of the eight charge blocks has 512
rows and columns, with all four old vertices retained. Invertible linear
maps of the color space `F_2^3` preserve equality, distinctness and XOR
zero, and act transitively on its seven nonzero charges. Thus one zero
and one nonzero charge determine the entire low prefix.

The [checker](../tools/check_fixed_curve_quartic_floor.py) gives:

| Tensor | Charge zero | Each nonzero charge |
| --- | --- | --- |
| Gaussian control | `0^1, 12^7, 18^42` | `6^2, 12^12, 18^38` |
| Quartic tensor | `0^1, 12^7, 18^30` | `6^2, 12^12, 18^32` |

These list all orders through eighteen. In particular the quartic leaves
the order-twelve layer unchanged. The whole 4096-dimensional profile has
the exact prefix

```math
0^1,\qquad6^{14},\qquad12^{91},\qquad18^{254}.          \tag{2}
```

All repeated colors and four-occupied rows are included. For a matrix
with Smith orders s_j, the order-p Toeplitz truncation has rank
`sum_j max(p-s_j,0)`. Successive rank differences therefore certify (2).
Independent column construction with low pivots checks the row-streaming
ranks at truncations 7,13,18,19. Literal quartic-support/Pfaffian
enumeration checks 61,376 entry controls against the faster matching
formula, including repeated colors and the complete binary word panel.

## 3. A tail bound without computing the remaining orders

Put n=4096 and `Z(x)=sum_j x^(s_j)`. The endpoint is invertible, so all
n orders are finite. Every full-determinant term has degree

```math
8\sum_{rows,columns}\operatorname{occ}(word)=56n.
```

Indeed each of four old positions is nonzero on seven of its eight
letters. Thus the sum of the Smith orders is at most 56n. The 360
known orders in (2) sum to 5748. The other 3736 orders have mean at most

```math
\frac{56\cdot4096-5748}{3736}=\frac{55907}{934}<60.
```

Convexity and monotonicity of x^v, for 0<x<1, now give

```math
Z(x)\ge L(x):=1+14x^6+91x^{12}+254x^{18}+3736x^{60}.   \tag{3}
```

## 4. Exact all-scale certificate and paired entropy

Substitute x=u^3. The difference between L(x) and
`(1+7x^(32/3))^3`, divided by u^18, is

```math
P(u)=14-21u^{14}+91u^{18}+254u^{36}
       -147u^{46}-343u^{78}+3736u^{162}.                \tag{4}
```

All 201 coefficients of P in the degree-200 Bernstein basis are strictly
positive. The checker acquires them exactly by

```math
b_k=\sum_{j\le k}a_j\frac{\binom{k}{j}}{\binom{200}{j}},
```

where a_j are the monomial coefficients in (4), and reconstructs every
monomial coefficient, including the zero coefficients through degree 200.
The [receipt](../evidence/fixed-curve-quartic-floor.json) records the
positive minimum, its index and the complete coefficient-list hash.
Consequently

```math
Z(x)>(1+7x^{32/3})^3,\qquad0<x<1.                     \tag{5}
```

Let nu be any common type on the 64 paired letters, with entropy H(nu)>0
and mean occupied count E_nu in [0,2]. A word of length m on four such
legs has degree `32m E_nu`. The symmetric two-cut envelope is

```math
R_\nu(x)=\frac{2\log Z(x)-32E_\nu\log x}{4H(\nu)}.    \tag{6}
```

The paired-letter partition function with one-color weight x^(32/3) is
`(1+7x^(32/3))^2`. The ordinary entropy bound, without assuming independent
components of the pair, yields

```math
H(\nu)\le2\log(1+7x^{32/3})-\frac{32}3E_\nu\log x.
```

Combining this with (5) proves `R_nu(x)>3/4`. Taking the required infima
gives a floor at least 3/4, and `9/16>1/2` excludes the growth-base-two
target. No uncharged reuse of the old three-axis normalization is made.

## Verification and credit

```sh
make verify-fixed-curve-quartic
```

This reproduces the exact jets, independent entry/column controls and
positive Bernstein certificate. The [full default receipt](../evidence/fixed-curve-full-verification.json)
and [complete compressed stdout](../evidence/fixed-curve-full-verification.stdout.txt.xz)
record exit-zero `make verify` at clean frozen
`705b3e34b8564a290b6f451f623515c62a1f74db`, tree
`744418a3779e9ef2c2b79ff9a20acc3d5c78db32`, in 1981.39 seconds.
A fresh marker probe checked private /tmp, read-only host root and only
the frozen checkout writable; network was unshared, and tracked state
stayed clean. The log includes the new checker. Later main integration
and the q99 numerical replay correction have separate focused checks,
including documentation links in 1202 Markdown files; those changes are
outside this full-source pin. The receipt records the initial q99 scale
identity failure and its pressure-equation correction, without rebuilding
the q99 profile.

These are author-checked ordinary arguments and finite exact linear
algebra, without new Lean or independent mathematical review. The Gaussian,
color-class and quartic constructions retain their
linked predecessor credit; the fixed projection of the Gaussian was
already proved in the cited note. The new quartic fixed-core certificate
was developed by Sol worker 2 with substantial OpenAI Codex assistance.
Further progress must change the retained early structure, the coefficient
rank estimate or the global evaluator; expanding this stopped profile is
not that progress.
