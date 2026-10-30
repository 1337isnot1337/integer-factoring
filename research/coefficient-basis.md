# Changing the coefficient basis cannot improve this contraction

This note proves a restricted-method obstruction. For every even block count,
the compressed coefficient decomposition already minimizes the largest product
of row and column matrix ranks among all decompositions using linear
combinations of the same symmetric coefficients. The obstruction does not
apply to arbitrary tensor decompositions or to factoring algorithms in general.

All exact-rank statements in Sections 1--4 are over the rationals, or any field
of characteristic zero. The finite search uses exact integer arithmetic.

## 1. Simultaneous eigenspaces

Let `t` be the number of three-bit blocks, and retain

\[
 A_i=[T^i]\prod_{s=1}^{t}(R_s+TC_s).
\]

The one-block row matrices have simultaneous eigenvalues `(8,-6)` on an
eight-dimensional subspace and `(0,2)` on its 56-dimensional complement.
For `j` complement positions, the tensor sector has dimension

\[
 d_j=\binom tj8^{t-j}56^j.
\]

On that sector the row eigenvalue of `A_i` is

\[
 M_{j,i}=
 \begin{cases}
 2^j\binom{t-j}{i-j}8^{t-i}(-6)^{i-j},&j\le i,\\
 0,&j>i.
 \end{cases}
\]

The column eigenvalue is `M_(j,t-i)`. The matrix `M` is triangular with
nonzero diagonal, so `A_0,...,A_t` are linearly independent. In particular,
if a collection of pieces sums to `sum_i A_i`, the sums of their coefficients
in this basis must equal one at every coordinate.

For `F=sum_i c_i A_i`, let `L` and `U` be the smallest and largest indices
with `c_i!=0`. The row eigenvalue on sector `U` is `c_U M_(U,U)!=0`;
the column eigenvalue on sector `t-L` is `c_L M_(t-L,t-L)!=0`. Therefore

\[
 \operatorname{rank}_{\rm row}F\ge d_U,
 \qquad
 \operatorname{rank}_{\rm col}F\ge d_{t-L}.
\]

These endpoint sectors cannot be removed by cancellations involving other
coefficients. This is the part that a proposed Krawtchouk or interpolation
basis change must confront.

## 2. Exact central obstruction

Suppose `t=2m>=2` and `c_m!=0`. Write

\[
 r_m=\sum_{j=0}^{m}d_j.
\]

Then

\[
 \operatorname{rank}_{\rm row}F
 \operatorname{rank}_{\rm col}F\ge r_m^2,
\]

with equality only when `F` is a nonzero scalar multiple of `A_m`.

Here is a proof with explicit constants.

First, `d_j` is log-concave because

\[
 \frac{d_{j+1}}{d_j}=7\frac{t-j}{j+1}
\]

decreases with `j`. On the interval `j=m,...,t`, its minimum is therefore
at an endpoint. Moreover

\[
 \frac{d_t}{d_m}=
 \frac{7^m}{\binom{2m}{m}}\ge\frac72.
\]

For completeness, this ratio is `7/2` at `m=1`, and its ratio at consecutive
values of `m` is `7(m+1)/(2(2m+1))>1`.
It follows that `d_j>=d_m` for all `j>=m`.

On the shorter interval `j=m+1,...,t`, both endpoints are at least
`(7/2)d_m`: the first has ratio `7m/(m+1)>=7/2`, and the last was just
handled. Thus

\[
 d_j\ge\frac72d_m\quad(j\ge m+1).
\]

For `j<m`, consecutive terms going downward have ratio

\[
 \frac{d_j}{d_{j+1}}=
 \frac{j+1}{7(2m-j)}<\frac17.
\]

Consequently `r_m<(7/6)d_m`. If `L<m` or `U>m`, one endpoint sector
has dimension at least `(7/2)d_m` and the other at least `d_m`. Hence

\[
 \operatorname{rank}_{\rm row}F
 \operatorname{rank}_{\rm col}F
 \ge\frac72d_m^2
 >\frac{49}{36}d_m^2
 >r_m^2.
\]

The remaining case is `L=U=m`, when `F` is proportional to `A_m` and
both ranks equal `r_m` exactly.

## 3. Minimax optimality of the entire decomposition

Consider any finite decomposition

\[
 \sum_s F_s=2^t E^{\otimes t},\qquad
 F_s\in\operatorname{span}_{\mathbb Q}\{A_0,...,A_t\}.
\]

At least one piece has a nonzero coefficient of `A_m`. Section 2 shows
that its rank product is at least `r_m^2`, regardless of the number of
pieces, coefficient signs, or cancellations.

Conversely, the original decomposition into `A_i` attains maximum rank
product `r_m^2`. To see this, the exact ranks are
`rank_row(A_i)=r_i` and `rank_col(A_i)=r_(t-i)`, where
`r_i=sum_(j<=i)d_j`. The sequence `r_i` is log-concave. An elementary
proof follows by summing

\[
 d_{i+2}d_j\le d_{i+1}d_{j+1}\quad(j\le i),
\]

which holds because the ratios of `d_j` decrease. It gives
`d_(i+2)r_i<=d_(i+1)r_(i+1)`, equivalently
`r_(i+1)^2>=r_i r_(i+2)`. Thus the symmetric product
`r_i r_(t-i)` is maximized at `i=m`.

The minimax value is therefore exactly

\[
 \boxed{\min_{\{F_s\}}\max_s
  \operatorname{rank}_{\rm row}(F_s)
  \operatorname{rank}_{\rm col}(F_s)=r_m^2.}
\]

Since

\[
 r_m^2\ge d_m^2
 =\binom{t}{t/2}^{\!2}448^t
 =2^{(8+\log_2 7)t-O(\log t)},
\]

no such coefficient-basis change improves the leading contraction
`rho=(8+log2 7)/12`. The proof covers the whole linear span, not merely
separate refactorizations of individual coefficients.

In odd positive characteristic, extra binomial coefficients can vanish.
The endpoint-sector argument remains valid because its diagonal coefficients
are powers of two. It still forces the same leading-exponent obstruction,
though the exact characteristic-zero value `r_m^2` need not remain the
finite-dimensional minimax value. Characteristic two is excluded: the
crossing multiplier itself then vanishes.

## 4. Exact finite search

`python3 tools/check_coefficient_basis.py` searches all one-dimensional
intersections of `t` row/column eigenvalue-zero hyperplanes for `t=2,4,6`.
It uses signed integer minors, not floating-point rank decisions.

These intersections suffice for minimizing rank product with `c_m!=0`:
if the current zero constraints leave dimension at least two, some additional
eigenvalue functional is independent of `c_m` on that space. Setting it to
zero retains a vector with `c_m!=0` and cannot increase either rank.
Repeating reaches a one-dimensional intersection.

| Blocks `t` | Hyperplane subsets | Distinct tested lines with `c_m!=0` | Minimum rank product | Unique minimizing basis vector |
| ---: | ---: | ---: | ---: | --- |
| 2 | 15 | 13 | 921600 | `A_1` |
| 4 | 210 | 196 | 1750350168064 | `A_2` |
| 6 | 3003 | 2913 | 4009028576958480384 | `A_3` |

The search checks the finite rank formulas independently of the endpoint
inequality proof. It is not evidence for unrestricted tensor lower bounds.

## 5. Other square identities: direct substitutions tested

For a norm-composition identity with `D=2^b` inputs, cheap rank `D`, and
the analogous rank-`D(D-1)` complement, the same construction gives

\[
 \rho_b=\frac{2+2b+\log_2(D-1)}{4b}.
\]

The ordinary two-, four-, and eight-square identities give respectively
`1`, `0.9481203125901445...`, and `0.9006129101714669...`. Mixing these
block sizes gives a bit-length-weighted average, so it does not beat the
eight-square choice. A matching bilinear 16-square identity is unavailable
in characteristic different from two; the classical dimension restriction
is summarized in Section 1 of the primary paper cited below.

The published large-dimension alternative checked here is
[Pavel Hrubes, *A Subquadratic Upper Bound on Sum-Of-Squares Composition
Formulas*, CCC 2024](https://doi.org/10.4230/LIPIcs.CCC.2024.12).
Its proof gives the asymptotic square count `s=D^(c+o(1))`, with
`c=1+2H(a/2)=1.6146510368528273...`, where `H(a)=1/2` and
`a=0.11002786443835955...`.

A direct coefficient construction using cheap rank `s` and the generic
opposite rank `D^2` has limiting contraction
`(c+2)/4=0.9036627592132068...`, before any finite-block selector overhead.
That substitution does not improve the current value. Within that direct
cost calculation, one would need `c<4rho-2=1.6024516406858676...`.
This comparison does not rule out a different use of the paper's identities.

## 6. Why the nonsymmetric word expansion is insufficient

Expanding `prod_s(R_s+C_s)` into its `2^t` individual words leaves the
symmetric coefficient span. Every word has row/column rank product
`(8*64)^t=512^t`, suggesting the much smaller local contraction `3/4`.
This does not give a subexponential factoring improvement: choosing a word
independently at each grid rectangle creates

\[
 (2^t)^{\binom G2^2}
 =2^{(1/12-o(1))G^4\lambda}
\]

possibilities already at the top level. The exponent is proportional to
the top volume, so explicit path enumeration is exponential in the input
length. The local rank saving cannot pay for this cost.

An improved nonsymmetric decomposition would therefore need both a smaller
rank product and a sufficiently small number `M(lambda)` of pieces. For
this enumeration architecture, `log M(lambda)=o(lambda)` is necessary to
make this particular top-level cost sublinear. The coefficient-basis
obstruction does not rule out such decompositions; the word expansion
does not provide one.

## 7. Recurrence floor/ceiling audit

The separate recurrence note's integer choices were rechecked:

* `B=floor(log2 K)+1` exceeds `ln K`, as required by the concentration bound.
* `G<=K^8` follows from `12<=K`, `D<=(K)^2`, and `B<=K`; it is a coarse
  overhead bound and is not substituted for the sharper base-four growth.
* Maximality of the cutoff implies `n<K_(ell+1)^40<=P^640`, so
  `log n=O(log P)` holds, including its direction.
* `t=ceil(lambda/3)` contributes at most the `+3G^4` term used in the new
  volume recurrence.
* The finite-input fallback handles values below `48^40`; there is no claim
  that the admissibility inequality holds for the returned depth zero on
  those inputs.

No correction to the recurrence constants was needed in this pass.
