# Arbitrary mode sets, independent occupations, and exact type pressure

## Result and scope

The polynomial Schur kernel defines typed crossings on **arbitrary**
sets of binary indices. A Boolean down-set is not required. Restricting
tensor powers to fixed occupation counts at each mode gives an exact
crossing endpoint, including sectors with unequal plus and minus
occupations. The same restriction permits an independent degree height
at every boundary, species, and mode.

For Gaussian Smith certificates, optimization over those independent
occupations has an explicit product partition function. A strict
certificate can be checked using rational arithmetic, and its optimizing
occupation fractions are rational. We also record the balanced Schur
specialization with variable energy and the common monomial needed to
recover the complete Schur polynomial.

This supplies a wider constructive search space and its exact transfer
criterion. No parameter choice satisfying the stronger-than-half target
is established here. Numerical optimization remains separate from the
theorems and checks below.

The kernel comes from the
[complete Schur determinant formula](schur-fermion-queries.md) and
[sparse-mode construction](schur-downset-crossing.md). The latter's
down-set hypothesis remains necessary for its Pascal and holonomy
identities. This note uses their displayed polynomial entries without
asserting those identities for arbitrary sets. The arbitrary-type
transfer is developed further in
[the exact type-pressure note](type-optimized-crossing-pressure.md).
The extensions and checks below were developed with substantial OpenAI
assistant support across the Schur, holonomy, and pressure research lanes.

## 1. A polynomial kernel for any mode set

Let `M` be a power of two and let `S` be any nonempty subset of
`{0,...,M-1}`. Put `n=|S|`, `T=M-1`, `z=w^2`, `q=1+z`, and
`c=1+w`, so `c^2=q`. All matrices below are indexed by `S`:

```math
P_{a,b}=\mathbf1_{b\subseteq a},\quad
Q=\operatorname{diag}(q^a),\quad
Q_J=\operatorname{diag}(q^{T-a}),\quad
G_{a,b}=z^{a\mathbin{\&}b},\quad
H_{a,b}=z^{T-(a\mathbin{|}b)}.                              \tag{1}
```

The bitwise union is at most `T`, so all entries in the following
plus-to-minus matrix are polynomials in `w`:

```math
\mathcal A_S(w)=
\begin{pmatrix}
0&c z^M H&cP^TQ_J&I\\
c z^M G&0&z^M P&cG\\
cQP^T&z^M P&0&cQP^T\\
I&cH&cP^TQ_J&0
\end{pmatrix}.                                             \tag{2}
```

The boundary order is `B,X,V,Y`. A physical label is an arbitrary pair
`(A,B)` of subsets of `S`; their cardinalities may differ. For physical
labels `x,y,u,v`, use these occupations:

| Boundary | Plus modes | Minus modes |
| --- | --- | --- |
| Bra `B` | `B_u` | `A_u` |
| First parameter `X` | `B_x` | `A_x` |
| Second parameter `V` | `A_v` | `B_v` |
| Ket `Y` | `A_y` | `B_y` |

Let `T_S(w;x,y,u,v)` be the selected determinant of (2), with value zero
if the selected row and column counts differ. This is a Gaussian tensor
on `4^n` physical labels. It is defined independently of a Schur
partition interpretation.

At `w=1`, only the `B,Y` identity blocks and the `X,V` Pascal blocks
survive. Consequently

```math
T_S(1;x,y,u,v)
=\mathbf1_{y=u}\det P_{B_x,B_v}\det P_{A_v,A_x}.            \tag{3}
```

A rectangular minor is interpreted as zero. No involution identity for
`P` is used in (2) or (3). In particular, both formulas apply to
`S={0,3}`, although the down-set identity `G=PQP^T` fails for that set.

## 2. Fixed-mode endpoint theorem

Take `m` copies of `T_S` with the same variable `w`. For each mode
`a in S`, choose integers `k_A(a),k_B(a)` between zero and `m`.
Let `L_m` be the set of label words

```math
\mathbf x=((A_1,B_1),\ldots,(A_m,B_m))
```

in which mode `a` occurs exactly `k_A(a)` times among the `A_i`, and
exactly `k_B(a)` times among the `B_i`. Its dimension is exactly

```math
D_m=\prod_{a\in S}\binom m{k_A(a)}\binom m{k_B(a)}.         \tag{4}
```

There is no requirement that the two total occupation numbers agree.

**Theorem.** Restrict all four physical legs of `T_S(w)^(tensor m)`
to `L_m`. The resulting polynomial `Gamma_m(w)` satisfies

```math
\Gamma_m(1;\mathbf x,\mathbf y,\mathbf u,\mathbf v)
=\mathbf1_{\mathbf x=\mathbf v}\mathbf1_{\mathbf y=\mathbf u}.
                                                               \tag{5}
```

**Proof.** The identity blocks in (3) force `y_i=u_i` at each copy.
A nonzero Pascal minor has at least one matching of its column masks
into its row masks by subset inclusion. Such a matching cannot decrease
the numerical sum of the indices. Therefore the two Pascal factors
in (3) imply, copy by copy,

```math
\sum B_{v_i}\le\sum B_{x_i},\qquad
\sum A_{x_i}\le\sum A_{v_i}.                               \tag{6}
```

The fixed-mode counts give equal totals across the `m` copies, separately
for the `A` and `B` species. All increases in (6) must therefore be zero.
A proper bitwise inclusion strictly increases the numerical index, so
every matched inclusion is equality. Thus `A_(x_i)=A_(v_i)` and
`B_(x_i)=B_(v_i)` for every copy.

Conversely, when `x=v`, each remaining Pascal minor is principal. Sorting
its indices numerically makes it triangular with diagonal ones, so its
determinant is one. This proves (5).

The endpoint argument only needs fixed total weighted index sums for
each species. Fixed per-mode counts are stronger; they also permit all
the independent gauges below.

### Why species cardinalities alone are insufficient

Take `S={0,1}`, `x=({0},empty)`, `v=({1},empty)`, and
`y=u=(empty,empty)`. The two external labels have the same species
cardinalities, but (3) equals `P_(1,0)=1` although `x!=v`.
The fixed weighted-sum condition is doing necessary work in this example.

## 3. Arbitrary per-mode degree heights

At zero and infinity, independently assign integer heights

```math
\eta_{e,D,+}(a),\qquad \eta_{e,D,-}(a),
\quad e\in\{0,\infty\},\quad D\in\{B,X,V,Y\},\quad a\in S.
                                                               \tag{7}
```

At zero the local parameter is `tau=w`; at infinity it is `tau=w^-1`.
Multiply a plus row and minus column of the one-particle matrix by their
respective powers of `tau`. This multiplies a selected determinant by
the sum of its occupied heights.

For the two physical species, define

```math
\begin{aligned}
g_{e,A}(a)&=\eta_{e,B,-}(a)+\eta_{e,X,-}(a)
            +\eta_{e,V,+}(a)+\eta_{e,Y,+}(a),\\
g_{e,B}(a)&=\eta_{e,B,+}(a)+\eta_{e,X,+}(a)
            +\eta_{e,V,-}(a)+\eta_{e,Y,-}(a).
\end{aligned}                                                \tag{8}
```

On the restricted tensor in (5), the total gauge shift at end `e` is
the constant

```math
a_e=\sum_{a\in S}\bigl(g_{e,A}(a)k_A(a)+g_{e,B}(a)k_B(a)\bigr).
                                                               \tag{9}
```

No affine dependence on the numerical index is required. Put
`G_A=g_(0,A)+g_(infinity,A)` and similarly `G_B`.

## 4. Gaussian profiles and the product pressure

Choose opposite crossing cuts at zero and infinity. At each end, suppose
the internal Gaussian blocks on both sides of the cut are regular over
the corresponding power-series ring. Their Gaussian multiplication
operators are then invertible over that ring, since their constant
coefficient in fermion degree is one.

After removing those internal operators, the coefficient matrix is the
direct sum of the exterior powers of the one-particle cross matrix.
Let its `4n` Smith orders at end `e` be
`lambda_(e,1),...,lambda_(e,4n)`. Cross entries and their orders may be
Laurent; regularity of the internal blocks is the condition needed for
this argument. The full tensor profile consists of all subset sums of
the one-particle orders, so for `s>1` its profile sum is

```math
U_e(s)=\prod_{j=1}^{4n}(1+s^{-\lambda_{e,j}}),\qquad
B_{\mathrm{raw}}(s)=U_0(s)U_\infty(s).                     \tag{10}
```

If a gauge does not satisfy the internal regularity condition, formula
(10) is not asserted. One may instead compute the full tensor's Smith
profiles directly and use their profile sums.

Let `p_A(a)=k_A(a)/m` and `p_B(a)=k_B(a)/m`. Write

```math
H(p)=\sum_{a\in S}\bigl(h(p_A(a))+h(p_B(a))\bigr),
\qquad h(t)=-t\log t-(1-t)\log(1-t),                       \tag{11}
```

with natural logarithms and the usual endpoint convention. Then
`log D_m=m H(p)+O(n log(m+1))`. The coefficient-rank transfer at the
two ends gives

```math
\sum_j\sqrt{\operatorname{rank}_R\Gamma_{m,j}
                  \operatorname{rank}_C\Gamma_{m,j}}
\le m^{O(1)}\left[
B_{\mathrm{raw}}(s)\,
s^{\sum_a(G_A(a)p_A(a)+G_B(a)p_B(a))}
\right]^{m/2}.                                             \tag{12}
```

Here the mode set and heights are fixed as `m` grows. To see the degree
cancellation, a coefficient of degree `j` has zero-end index `j+a_0`
and infinity-end index `a_infinity-j`. The two geometric bounds have
opposite `j` exponents, leaving precisely the restoration factor (9).
The number of polynomial coefficients and the truncated unit-series
overhead are polynomial in `m`.

### Exact elimination of every occupation fraction

Fix a proposed rank-width rate `rho` and put `R=4rho`. Define

```math
Z_R(s)=\prod_{a\in S}
       (1+s^{-G_A(a)/R})(1+s^{-G_B(a)/R}).                   \tag{13}
```

The elementary binary entropy identity

```math
\max_{0\le p\le1}\{R h(p)-pG\log s\}
=R\log(1+s^{-G/R})                                        \tag{14}
```

shows that some independent occupation fractions satisfy the strict
rate condition in (12) if and only if

```math
\boxed{\qquad B_{\mathrm{raw}}(s)<Z_R(s)^R.\qquad}          \tag{15}
```

The maximizing fractions are explicit:

```math
p_A(a)=\frac{s^{-G_A(a)/R}}{1+s^{-G_A(a)/R}},\qquad
p_B(a)=\frac{s^{-G_B(a)/R}}{1+s^{-G_B(a)/R}}.                \tag{16}
```

There is no equal-count Lagrange multiplier in (13)--(16). The endpoint
theorem allows the two species to have different total occupations.

For an entirely rational certificate, take positive integers `u,v`,
put `R=u/v`, and set `s=t^u` for a rational `t>1`. Then (15) becomes

```math
B_{\mathrm{raw}}(t^u)^v
<\left[\prod_{a\in S}
       (1+t^{-vG_A(a)})(1+t^{-vG_B(a)})\right]^u.             \tag{17}
```

Both sides are rational, and the fractions in (16) are rational.
Taking `m` divisible by their denominators gives exact integer mode
counts. A strict inequality absorbs the polynomial factors in (12)
as `m` grows. In the repository's normalization, a certificate with
`R<2sqrt(2)` would give `rho<1/sqrt(2)` and pass the post-half target.
This note supplies no such parameter certificate.

The construction is effective: selected determinants and coefficient
products give complete tables in time exponential in the physical label
width, for any fixed finite mode set, offset, and rational certificate.
No count selection or endpoint verification requires an asymptotic oracle.

## 5. Balanced Schur sectors with variable energy

Return to labels with `|A|=|B|=r`, but allow

```math
\Delta=\sum A-\sum B,
\qquad h=Mr+\Delta.                                       \tag{18}
```

Their Schur arms are `A`, and their legs are `{M-1-b:b in B}`.
Fix `r` and `Delta` at every boundary. In the complete inhomogeneous
Schur construction, the polynomial kernel (2) has the exact relation

```math
\boxed{\qquad
P_h(w^2)=w^{4\Delta}\,
\left.\det\mathcal A_S(w)\right|_{r,\Delta}.
\qquad}                                                   \tag{19}
```

Indeed, the row and column factors used to obtain (2) contribute
`z^(2Mr)q^(h+r)`, whereas the complete Schur prefactor is
`z^(2h)q^(h+r)`. Their ratio is `z^(2Delta)=w^(4Delta)`.
Negative `Delta` is allowed; the equality can first be read as a Laurent
identity. The complete Schur polynomial on the left is polynomial.
The common monomial shifts the zero and infinity orders oppositely,
so it cancels from the product pressure.

Suppose the heights at each end have the affine form

```math
\eta_{e,D,+}(a)=c_{e,D,+}+\alpha_{e,D}a,\qquad
\eta_{e,D,-}(a)=c_{e,D,-}-\alpha_{e,D}a.                    \tag{20}
```

The gauge shift on this sector is `r C_e+Delta A_e`, where

```math
C_e=\sum_D(c_{e,D,+}+c_{e,D,-}),\qquad
A_e=\alpha_{e,V}+\alpha_{e,Y}-\alpha_{e,B}-\alpha_{e,X}.     \tag{21}
```

Let `C=C_0+C_infinity` and `A=A_0+A_infinity`.
For the full interval `S={0,...,N-1}`, define centered heights

```math
g_a=C/2+A\bigl(a-(N-1)/2\bigr).                            \tag{22}
```

Optimizing entropy subject to equal total `A` and `B` occupations has
the exact partition function

```math
Z_{R,\mathrm{balanced}}(s)
=\prod_{a=0}^{N-1}(1+s^{-g_a/R})^2.                        \tag{23}
```

One way to derive it is to introduce a count multiplier. The two
species' heights are initially `C/2+A a+mu` and `C/2-A a-mu`.
At `mu=-A(N-1)/2`, their lists are reversals. Their optimizing
occupations therefore have equal totals. This proves (23) by the
entropy identity (14).

### A finite counting bound for the balanced specialization

Let `d_(r,Delta)` count pairs of `r`-subsets with difference `Delta`.
The exact weighted equal-count sum is

```math
Z_{\mathrm{eq}}
=\sum_{r,\Delta}d_{r,\Delta}
                  s^{-(Cr+A\Delta)/R}.                    \tag{24}
```

In the product distribution from the centered heights, the two species'
cardinalities are independent and identically distributed: their
Bernoulli probability lists are reversals. If `p_r` is this common
cardinality distribution, then

```math
\frac{Z_{\mathrm{eq}}}{Z_{R,\mathrm{balanced}}}
=\sum_{r=0}^N p_r^2\ge\frac1{N+1}.                        \tag{25}
```

Thus the equal-count restriction loses at most a polynomial factor.
There are only polynomially many pairs `(r,Delta)`. Pigeonholing (24)
therefore selects one fixed occupation and Schur energy sector with
another polynomial loss. For growing mode sets, a strict exponential
pressure margin absorbs these losses whenever the coefficient and
construction overheads are subexponential in the label width. Formula
(23) alone, without such a margin and overhead accounting, is not a
stronger factoring certificate.

## 6. Exact checks

Run the standard-library checker:

```sh
python3 tools/check_fixed_mode_crossings.py \
  --output evidence/fixed-mode-crossings.json
```

It tests arbitrary and non-down-set mode sets, unequal species counts,
the exact dimension product (4), direct selected-determinant endpoints,
and constant shifts for independent per-mode heights. A negative control
exhibits the off-diagonal endpoint that survives when only species
cardinalities are fixed.

The checker compares positive and negative variable-energy sectors with
the independent complete Schur evaluator, including non-prime-power
weights. It verifies exact binomial type bounds, product partition
expansions, the finite equal-count bound (25), and its energy-pigeonhole
consequence using rational arithmetic. Finally, twelve full four-label
Smith calculations independently match the subset sums of their
one-particle profiles.

See the [checker](../tools/check_fixed_mode_crossings.py) and
[verification output](../evidence/fixed-mode-crossings.json).
