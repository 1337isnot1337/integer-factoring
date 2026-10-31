# Joint determinant certificates for actual tensor-power ranks

## Result and scope

The [constant intertwiner](transpose-intertwiner-barrier.md) is only one
way to rule out misleading finite gains. A common determinant certificate
can force a large rank in **both** cuts of the same coefficient piece,
even when the local matrices have no constant transpose intertwiner.

The resulting ordinary theorem closes an entire class:

> A binary-alphabet quadratic with pure or zero values at two distinct
> non-endpoint projective points, or a cubic with pure or zero values at
> three such points, cannot achieve leading `rho < 4/5` through its
> full tensor-power coefficient span. This holds over every field for
> which the selected distinct points exist, and allows arbitrary
> coefficient recombinations, Taylor jets and infinity reversal.

Thus all **135,005,697** ternary pure-value cubic coefficient seeds from
the earlier catalogue are excluded as post-half improvements by an
ordinary uniform proof. The theorem concerns full binary word alphabets.
It does not exclude new physical code restrictions, mixing spaces from
different seeds, or low-rank values that are not pure. It is not a lower
bound on factoring. The complete factoring score remains **1/2**.

The key representation-theoretic dependency is a Weyl filtration, not
semisimplicity. The argument therefore retains characteristics two and
three. No Lean proof is added.

## 1. Determinant multiplicity in every characteristic

Let `V` have dimension `n`. For each integer `t>=0`, the tensor power
`V^(tensor nt)` has a Weyl filtration. The one-dimensional polynomial
representation `det^t`, of highest weight `(t,...,t)`, occurs with
multiplicity

```math
h_n(t)=\frac{(nt)!\prod_{j=0}^{n-1}j!}
                   {\prod_{j=0}^{n-1}(t+j)!}.        \tag{1}
```

For the filtration, see Stephen Doty,
[*Schur-Weyl duality in positive characteristic*, arXiv:math/0610591v3](https://arxiv.org/pdf/math/0610591v3),
Lemma 13 and the formal-character discussion before Proposition 14,
pages 13--14. The characters of the Weyl factors are the characteristic-
independent Schur characters. Expanding the character
`(x_1+...+x_n)^(nt)` in that basis gives the usual standard-tableau
multiplicity (1) for the rectangle `(t^n)`.

These are integer multiplicities of separate filtration factors. We do
not divide by `(nt)!` in the coefficient field or reduce `h_n(t)` as a
field scalar. In particular, a prime dividing a factorial in (1) does
not erase the rank bound below. For a finite base field, first extend
to its algebraic closure to use the filtration. Matrix ranks and
nonvanishing of the base-field coefficient functionals are unchanged.

For completeness, (1) counts standard tableaux of rectangular shape.
The reflection-principle expression is
`(nt)! det(1/(t-i+j)!)_(0<=i,j<n)`, with negative factorial reciprocals
zero. Factoring the factorials and evaluating the resulting Vandermonde
determinant gives (1), as an integer identity. The elementary bound

```math
h_n(t)\ge\frac{n^{nt}}
 {(nt+1)^n(t+n)^{n(n-1)/2}}                         \tag{2}
```

is sufficient here. The balanced multinomial coefficient
`(nt)!/(t!)^n` is largest among its type class and is at least
`n^(nt)/(nt+1)^n`. Its ratio to (1) is a product of at most
`n(n-1)/2` factors bounded by `t+n`, after discarding numerator
factorials at least one. Thus `h_n(t)=n^(nt)/poly(t)` for fixed `n`.

### Rank consequence for a polynomial matrix

Let `M(z,w)` be an `n` by `n` homogeneous polynomial matrix, and let
`f=det M`. Take a scalar linear functional `L` on homogeneous polynomial
coefficients. In a constant basis respecting the Weyl filtration,
`M^(tensor nt)` is block triangular and has `h_n(t)` separate diagonal
blocks equal to `f^t`. Coefficient application preserves this triangular
form. Therefore

```math
L(f^t)\ne0\quad\Longrightarrow\quad
\operatorname{rank}L(M^{\otimes nt})\ge h_n(t).      \tag{3}
```

More generally, tensor on `r` additional copies. The selected diagonal
blocks become `f^t M^(tensor r)`. If `b` is any scalar linear
combination of entries of `M^(tensor r)` and `L(f^t b)!=0`, each selected
block is nonzero, so (3) still holds for `nt+r` positions. A minor
selecting one nonzero entry from each of those triangular blocks proves
the rank bound. This avoids any assumption of a direct-sum decomposition.
The filtration is stable under generic matrices; its coefficient
identities extend polynomially to the displayed matrix `M(z,w)`.

## 2. A single polynomial witnesses both cut ranks

Let `K(z,w)` be a degree-`D` homogeneous four-leg tensor on `q` letters,
with strict crossing at a selected endpoint. Both cut matrices have
order `n=q^2`; write their determinants as `f,g`.

Suppose

```math
f=h a,\qquad g=h b,                                 \tag{4}
```

where `h` is nonzero at the endpoint and each of `a,b` is a scalar
linear combination of the entries of `K`. In particular these quotient
forms have degree `D`. Then

```math
p=f b=g a                                          \tag{5}
```

is nonzero at the endpoint. For word length `m=(n+1)t+r`, with
`0<=r<n+1`, choose an entry `ell` of `K` that is nonzero at the endpoint.
The common polynomial

```math
p^t\ell^r=f^t b^t\ell^r=g^t a^t\ell^r              \tag{6}
```

is a linear combination of entries of `K^(tensor m)` in either cut.

Consider any pieces `A_i=L_i(K^(tensor m))` summing to a nonzero scalar
multiple of the crossing on `q^m` labels. Apply an entry covector
representing (6) to that tensor identity. Since (6) is nonzero at the
endpoint, some piece satisfies `L_i(p^t ell^r)!=0`. Equation (3), using
`nt` determinant positions and `t+r` remaining positions, applies to
**both cuts of that same piece**. Consequently

```math
\boxed{\max_i r(A_i)c(A_i)\ge h_n(t)^2,\qquad
       t=\lfloor m/(n+1)\rfloor.}                  \tag{7}
```

This does not depend on the number of pieces. Equations (2) and (7)
force leading `rho >= n/(n+1)`, hence `rho>=4/5` for binary letters.
The remainder `r` shows that a selected subsequence of word lengths
does not evade the bound.

If `f=g`, the common polynomial `f^t ell^r` works already at word length
`m=nt+r`, so the stronger leading floor is **`rho>=1`**. No constant
intertwiner is required for this case.

The same argument permits quotient forms represented by `R` positions:
if `f b'=g a'`, where `a',b'` are entry combinations of `K^(tensor R)`
and the common polynomial is nonzero at the endpoint, the leading floor
is `n/(n+R)`. This provides a practical certificate to test before
trying to prove an actual-rank upper bound.

## 3. Why pure quadratic and cubic values always satisfy a certificate

Take `q=2`, so `n=4`. A rank-one value in both cuts is a pure tensor.
At a pure or zero point the determinant of either cut vanishes to order
at least three: constant row/column changes put its rank-at-most-one
constant matrix in standard form, and at least three determinant factors
must then come from positive orders.

For a cubic with three pure/zero projective values, let `H` be the
product of their three homogeneous linear point factors. Then
`h=H^3` divides both determinants, has degree nine, and is nonzero at
the crossing endpoint. The quotient forms have degree three.

Write the interpolation values as `A,B,C,E`, with `E` the crossing.
The crossing's cut rank is four, so it cannot lie in the span of the
three rank-one tensors `A,B,C`. There are two cases:

- If `A,B,C` are independent, the four value tensors are independent.
  Projective Lagrange interpolation therefore makes the local entry
  span the whole four-dimensional space of cubic binary forms. Both
  determinant quotients belong to that entry span, giving (4)--(7).
- If `A,B,C` are dependent, their span has a basis of at most two pure
  tensors. In this case the two determinant polynomials are identical,
  giving the stronger floor one.

For the second claim, write `K=sE+uA+vB` with `A,B` pure. In either
cut this is an invertible permutation matrix times a scalar, plus two
rank-one updates. The determinant lemma expresses its determinant using
the two individual traces and their two-by-two Gram minor. Both traces
agree between the cuts. If
`A=a_x b_y c_u d_v` and `B=e_x f_y g_u j_v`, the products of cross-Gram
entries in the two cuts are the same scalar:

```math
(c\cdot f)(d\cdot e)(g\cdot b)(j\cdot a)
=(b\cdot g)(d\cdot e)(f\cdot c)(j\cdot a).
```

Thus the entire determinant expressions agree, over every field.
This also covers fewer or zero pure updates, and the identity is
polynomial even at values where `s=0`.

The quadratic proof is the same with two pure/zero values. The common
factor has degree six, the quotients degree two, and either the three
interpolation values span all quadratics or the two pure values are
dependent and the determinants agree. These arguments include actual
lower-degree polynomials represented in the indicated homogeneous degree.

## 4. Exact evidence and the constructive opening left by the theorem

Run `make verify-joint-determinant`. The checker verifies (1) by an
independent tableau-path count, checks the elementary exponential bound,
and tests the rank lemma on literal tensor powers in characteristics
two, three and five, including cases where factorials vanish in the field.
The uniform filtration remains the cited mathematical dependency;
finite tests do not replace it.

For all 200 retained cubic witnesses it expands both complete determinants,
divides out the point factors exactly, and retains the actual entry
covectors for the quotient certificates. There are **114 equal-
determinant certificates** and **86 joint quotient certificates**.
All 60 witnesses with deficient local entry span have equal determinants,
as predicted. Every witness is covered.

Over `F_3`, exactly **132,599,808** of the pure-value cubic coefficient
seeds have full local entry span. The other **2,405,889** have dependent
pure values and equal determinants. The checker counts the latter from
literal finite spans and the already independently checked first-value
stabilizer orbits. The ordinary two-case proof establishes the theorem
over arbitrary fields, beyond this ternary count.

This suggests a more specific constructive opening. Quadratic values
with cut-rank products two can have a common determinant factor of only
degree four. Its quotients have degree four and need two positions to
represent them from a full quadratic entry span. The resulting certificate
only gives `rho>=4/6=2/3`, leaving room below `1/sqrt(2)`. This is an
open interval for an upper-bound construction, not an achieved rate.
It points toward deliberately asymmetric values and a uniform mechanism
for actual coefficient cancellation, rather than pure-value interpolation.

The 321 retained asymmetric quadratic witnesses give a concrete next
test set: 126 have equal determinants, 47 have common determinant degree
four, and 148 have common degree five. Every unequal-determinant witness
has full quadratic entry span. The degree-four common point factor
therefore gives a two-position quotient certificate with floor `2/3`.
The saved checker identifies example indices; these finite class counts
do not assert a complete actual-rank search within any Smith signature.

The subsequent [common-cofactor certificate](common-cofactor-barrier.md)
closes this opening for the full ternary quadratic rank-product-two class:
all 32,265 normalized coefficient seeds have a stronger projected-minor
certificate. The determinant floor above remains valid; its apparent room
does not supply a construction. Other coefficient fields and new physical
restrictions are not settled by that finite existence classification.

This proof and its exact computations used substantial OpenAI Codex assistance.
