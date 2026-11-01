# Constructive extraction from border and nonunit crossing certificates

## Result and scope

A finite polynomial border certificate can be converted into exact
crossing seeds with only a polynomial loss in rank sum. The same
construction works at explicitly bounded precision over a truncated
polynomial ring. It gives a finite, checkable route from a nonunit
leading scalar to an ordinary exact certificate.

This is an application of coefficient extraction to the two
flattenings used in the factoring reduction. It supplies no new
numerical rank certificate. In particular, it does not establish a
factoring score above `1/2`.

The useful distinction is between a formal parameter and an integer
prime. Taking one coefficient of a polynomial is linear over its
coefficient field. Taking one base-`p` digit of an integer has carries
and is not linear. Section 5 gives an exact example where the latter
operation increases a matrix's channel requirement exponentially in
the precision.

## 1. A finite certificate over a truncated polynomial ring

Let `F` be a field, let `d>=2`, and write

```math
E_d(x,y,u,v)=[x=v][y=u].
```

Fix integers `H>t>=0` and the ring `R=F[z]/(z^H)`. Suppose tensors
`A_1,...,A_M` over `R` satisfy

```math
\sum_i A_i=z^t E_d+O(z^{t+1})\pmod {z^H}.          \tag{1}
```

The notation means that all coefficients below degree `t` vanish and
the degree-`t` coefficient is exactly `E_d`. No condition on the higher
coefficients is required. Suppose that each tensor has supplied exact
factorizations over `R` through `r_i` and `c_i` channels in the two
cuts `(x,y)|(u,v)` and `(x,u)|(y,v)`, respectively. All these data are
finite polynomial tables, and every asserted equality is verified
modulo `z^H`.

Set

```math
S=\sum_i\sqrt{r_i c_i}.
```

**Finite-precision extraction.** For every positive integer `k` with
`kt<H`, one can construct an exact crossing decomposition on `d^k`
labels over `F`, with at most `(k+1)^M` pieces and rank sum at most

```math
(kt+1)S^k.                                         \tag{2}
```

These are factor-channel upper bounds over `R`; no ambiguous notion
of matrix rank over a ring is needed.

### Construction and exactness

For each type `nu=(nu_1,...,nu_M)` with `sum nu_i=k`, define

```math
\begin{aligned}
B_\nu(z)&=\sum_{w:\operatorname{type}(w)=\nu}
 A_{w_1}(z)\otimes\cdots\otimes A_{w_k}(z),\\
C_\nu&=[z^{kt}]B_\nu(z).
\end{aligned}                                    \tag{3}
```

Use any polynomial lifts of the entries to compute the coefficient.
Because `kt<H` and all exponents are nonnegative, the result is
independent of those lifts. Every word occurs in exactly one type, so

```math
\sum_\nu C_\nu
=[z^{kt}]\left(\sum_i A_i(z)\right)^{\otimes k}
=E_d^{\otimes k}=E_{d^k}.                          \tag{4}
```

Only the leading term in (1) contributes to coefficient `kt`.

### Rank accounting

Let `D_nu=k!/(nu_1!...nu_M!)`. Tensoring the given row factorizations
and taking their direct sum gives a factorization of `B_nu` modulo
`z^H` through

```math
R_\nu=D_\nu\prod_i r_i^{\nu_i}
```

channels. Write it as `U_nu(z)V_nu(z)`. Its coefficient at degree `kt`
is

```math
C_{\nu,\mathrm{row}}=
\sum_{a=0}^{kt} U_{\nu,a}V_{\nu,kt-a}.
```

Thus `rank_row C_nu <= (kt+1) R_nu`. The same argument in the other
cut gives

```math
\operatorname{rank}_{\mathrm{col}}C_\nu
\le(kt+1)D_\nu\prod_i c_i^{\nu_i}.
```

Summing the square roots of the products and using the multinomial
theorem proves (2). Discarding zero pieces can only improve it.

An entry in (3) is the coefficient of `X^nu z^(kt)` in

```math
\prod_{j=1}^{k}
\left(\sum_{i=1}^{M}
 A_i(x_j,y_j,u_j,v_j;z)X_i\right).                 \tag{5}
```

A finite dynamic program computes it. There are at most
`(k+1)^M(kt+1)` relevant states; higher powers of `z` are discarded.
This charges the coefficient computation explicitly.

For a potential strict half-score improvement, it suffices to verify

```math
(kt+1)S^k<(d^k)^{\sqrt2}.                          \tag{6}
```

As in the ordinary finite-seed theorem, one should certify a strict
inequality using rational exponents and exact algebraic comparisons,
not a floating-point estimate. The condition `kt<H` is essential to
the rank estimate: the truncated factor identities must determine the
individual coefficients used there. For larger `k`, a chosen polynomial
lift still gives the crossing sum in (4), but the supplied truncated
factorizations no longer justify the claimed ranks of its pieces.

## 2. A fixed polynomial border certificate

Suppose now that the tensors `A_i(z)` have polynomial entries over a
field `F`, and

```math
\sum_i A_i(z)=z^t E_d+O(z^{t+1}).                   \tag{7}
```

Let `r_i,c_i` be their two ranks over the rational-function field
`F(z)`. A polynomial matrix of rank `r` over `F(z)` factors through
`r` channels already over `F[z]`. Indeed, its column module is a
free rank-`r` submodule of `F[z]^N`. To see freeness directly, project
any nonzero submodule to a coordinate on which it is nonzero. Its
image is a principal ideal `(g)`. Choose a module vector with that
coordinate equal to `g`, subtract its polynomial multiples to clear
the coordinate of every other vector, and proceed on the remaining
coordinates. This gives a polynomial basis. Expressing the original
columns in that basis gives the claimed factorization. Polynomial
Smith reduction is an equivalent effective construction.

Apply this fact separately to the two flattenings of every `A_i`.
Define `C_nu` by (3) for every positive `k`. There is now no precision
cutoff: all tensor and factor identities are polynomial identities.
The coefficient argument in Section 1 proves that the exact seed has

```math
\sum_\nu
\sqrt{\operatorname{rank}_{\mathrm{row}}C_\nu\,
      \operatorname{rank}_{\mathrm{col}}C_\nu}
\le (kt+1) S^k.                                   \tag{8}
```

The polynomial factors prove the rank bounds. The algorithm can
compute the extracted tensors directly from their polynomial entries
by (5), then find their exact field factorizations by elimination.
It does not need an oracle for infinite series.

For binary label width `lambda`, choose the least `k` with
`d^k>=2^lambda` and restrict the four labels to the first `2^lambda`
words. Put

```math
\rho=\frac{\log S}{2\log d}.
```

Both ranks decrease under restriction. After rounding their two
binary channel widths, every nonzero extracted piece has width sum

```math
\begin{aligned}
\lceil\log_2R\rceil+\lceil\log_2C\rceil
&\le 2k\log_2S+2\log_2(kt+1)+2\\
&=4\rho\lambda+O(\log(\lambda+2)).                \tag{9}
\end{aligned}
```

The number of types is polynomial in `lambda`, since `M` is fixed.
For a fixed finite coefficient field, the entry dynamic program is
polynomial in `lambda`. Building all tables and computing their exact
factorizations costs `2^O(lambda)` bit operations. Thus this is an
effective all-width family at the same asymptotic rank-sum quality
as the border certificate, using the usual finite-seed factoring
transfer. It also yields arbitrarily close finite exact seed qualities.

This proof is constructive even when the border family has no useful
single evaluation at a nonzero parameter: it uses a coefficient of
its tensor powers.

Rational-function border curves can be reduced to this case. Multiply
all pieces by one common polynomial denominator, including a power of
`z` if necessary, and normalize the nonzero leading scalar. Their
generic ranks are unchanged. The leading exponent `t` may increase
by a fixed amount, which does not affect the asymptotic quality in (9).

## 3. Coefficient fields and the meaning of a finite certificate

For a polynomial border certificate over an arbitrary field, only
finitely many coefficients occur. Generate a ring over `Z` by these
coefficients and specialize it to a finite-field residue quotient as
in [finite-field-rank-specialization.md](finite-field-rank-specialization.md).
The coefficient identities (7) remain exact. Every minor larger than
the asserted generic rank is a zero polynomial, so remains zero
after specialization. Generic ranks can only decrease.

If exact generic ranks are desired, for each cut choose a nonzero
coefficient of one maximal-minor polynomial and adjoin its inverse
before specialization. The selected minor then remains a nonzero
polynomial. None of this requires evaluating or representing the
original field inside the factoring machine.

A truncated-ring certificate is distinct from a rational-function
border certificate. Two factorizations modulo `z^H` need not extend
to coherent factorizations at every higher precision. Section 1
therefore promises only the exact finite outputs allowed by `kt<H`.
Section 2 applies when a polynomial curve and its generic rank bounds
are genuinely supplied.

## 4. Weighted cut budgets also hold for nonunit chain-ring identities

The nonunit scalar does not evade the weighted cut inequalities.
Here is a direct proof that covers finite precision as well as formal
border curves.

Let `R` be a commutative local principal ideal ring with maximal ideal
`(pi)`, `pi^H=0`, and `pi^(H-1)!=0`. This includes
`F[z]/(z^H)` and `Z/p^H`. Its residue field is `K=R/(pi)`.
Suppose tensors `A,B` satisfy

```math
A+B=\pi^t(E_d+\pi D),\qquad 0\le t<H,             \tag{10}
```

where `A` factors through `r` row channels and `B` through `c`
column channels. For every integer `1<=s<=d`,

```math
r+s^2c\ge sd.                                     \tag{11}
```

### Two elementary consequences of Smith reduction

Over such a ring, row and column operations reduce a matrix to a
diagonal matrix with powers of `pi` on its nonzero diagonal. One can
obtain this by repeatedly choosing an entry with least valuation:
it divides every remaining entry, so ordinary elementary operations
clear its row and column.

First, a submodule generated by `m` vectors in `R^d` is contained in
a free direct summand of rank at most `m`. Reduce the `d`-by-`m`
generator matrix as above, then take the free span of the basis
vectors containing its nonzero diagonal entries.

Second, the matrix `pi^t I_q` cannot factor through fewer than `q`
channels when `t<H`. If `pi^t I_q=XY` with `Y` having fewer than
`q` rows, the same diagonal reduction gives a primitive vector `v`
in the kernel of `Y`: at least one coordinate of `v` is a unit.
Then `pi^t v=0`, which is impossible because a unit times `pi^t`
is nonzero. Multiplying by invertible matrices does not change this
channel lower bound.

### Applying the quotient to the crossing

In row-matrix notation the column factorization of `B` has the form

```math
B=\sum_{j=1}^{c}P_j\otimes Q_j.
```

Let `U` be the free span of the first `s` standard basis vectors of
`R^d`. The vectors `P_j u`, for its `s` basis vectors and `1<=j<=c`,
generate a submodule contained in a free direct summand `L` of rank
`ell<=cs`. Let

```math
\Pi:R^d\longrightarrow R^{d-\ell}
```

be the corresponding quotient map. It kills those vectors and its
reduction modulo `pi` has full rank. Consequently

```math
(\Pi\otimes I)(A+B)|_{U\otimes R^d}
=(\Pi\otimes I)A|_{U\otimes R^d}
```

factors through at most `r` channels.

After removing its displayed common scalar `pi^t` in (10), the
reduction of the remaining matrix is the quotient of the swap.
It has rank

```math
q=s(d-\ell)
```

over `K`: it maps `u tensor v` to `Pi(v) tensor u`. A `q`-by-`q`
minor of that remaining matrix is therefore invertible over `R`.
The corresponding submatrix of the original matrix is `pi^t W`
for an invertible `W`. The second elementary fact forces `r>=q`,
and hence `r>=s(d-cs)`, which is (11).

For any multi-piece nonunit identity, take `A` to be the sum over a
chosen subset and `B` its complement. Direct sums of the supplied
factorizations give (11) with

```math
r=\sum_{i\in I}r_i,\qquad c=\sum_{i\notin I}c_i.
```

The same proof over a discrete valuation ring uses a saturated
lattice in place of the finite-ring diagonal reduction. It applies
to the normalized formal border `E_d+O(z)` as well. Thus previously
established consequences of these weighted budgets remain necessary
for this border/nonunit route.

More explicitly, choose the subset separately at each `s` to obtain
`sum_i min(r_i,s^2 c_i)>=sd` and its transpose. The ordinary budgets
`sum r_i>=d^2` and `sum c_i>=d^2` follow from the same channel lemma
applied to the full normalized crossing matrix. Factor dimensions may
be capped at `d^2` using the trivial ambient factorization. The
[weighted interpolation bound](https://github.com/1337isnot1337/integer-factoring/blob/27339f61d41af6167dd5ebe947260907b97f604c/research/weighted-cut-rank-budget.md)
therefore gives

```math
S\ge d(2H_{2d}-H_d-1).
```

This excludes the strict target `S<d^sqrt(2)` through 41 labels even
for these border or nonunit certificates. The stronger 46-label
exact-seed refinement in that source uses an additional extremal-piece
argument; its extension to this setting is not claimed here.

## 5. Why replacing the formal parameter by an integer prime is costly

The finite-precision extraction in Section 1 uses the `F`-linear
coefficient map `[z^j]`. There is no analogous additive map from
`Z/p^H` to `F_p` that sends `p^t` to one when `t>0`: every additive
map to `F_p` annihilates `p` times every element.

The ordinary digit map can have a large rank cost. Fix `e>=1` and
put `P=p^e`. Index a matrix by `0<=a,b<P`, with entries

```math
M_{ab}=a+b\pmod {p^{e+1}}.
```

This matrix factors through two channels over `Z/p^(e+1)`. Take its
base-`p` digit at position `e` using canonical representatives. Since
`a+b<2P<=p^(e+1)`, the resulting matrix over `F_p` is

```math
D_{ab}=[a+b\ge P].
```

Its zeroth row and column vanish. Restricting to `1<=a,b<P` and
reversing the column order gives the matrix `[a>=j]`, which is unit
lower triangular. Therefore

```math
\operatorname{rank}_{F_p}D=P-1=p^e-1.              \tag{12}
```

Thus a two-channel integer residue matrix can require exponentially
many prime-field channels after extracting one high digit. The
formal-parameter transfer does not justify such an operation with
only a polynomial loss. This example does not rule out another
construction that handles those carries explicitly.

## 6. Finite verification

`check_border_nonunit_transfer.py` verifies seven explicit border
extractions over `F_2`, `F_3`, and `F_5`, including three tensor copies
and cancellation outside the crossing support. It checks both exact
flattening ranks against the coefficient bounds. Its example seeds
are only transcription checks and do not improve the rank objective.

The same checker verifies 1,800 weighted budget inequalities for
sampled nonunit identities over `Z/4`, `Z/8`, `Z/9`, and `Z/27`,
using exact Smith reduction to compute channel dimensions. It also
checks seven instances of the digit-extraction rank formula (12).
The output is recorded in `border-nonunit-transfer.json`.

These are finite checks of the written arguments. Their dense Python
storage does not implement or measure the asymptotic working-space
bounds of the factoring machine.
