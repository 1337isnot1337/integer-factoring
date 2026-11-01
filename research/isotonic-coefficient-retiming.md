# Exact monotone coefficient retiming

## Result and scope

For polynomial crossing tensors whose two cuts have **constant diagonalizing
bases**, there is an exact optimization of the spacing between coefficient
indices. It changes no coefficient vector and preserves the endpoint. If
`a_i` counts row orders and `b_i` counts opposite-cut degrees, the optimal
pressure infimum is

\[
\boxed{\quad Q_*=
\left(\sum_{B}\sqrt{A_B B_B}\right)^2.\quad}
\tag{1}
\]

Here the consecutive blocks `B` are obtained by pooling adjacent profile
indices until their aggregate ratios `A_B/B_B` are nondecreasing. The
algorithm uses integer comparisons. The returned blocks include an exact
first-order optimality certificate.

There is also a uniform rational construction. At scale `x=1+1/M`, it
produces strictly increasing integer coefficient indices with pressure at
most `x^(D+1) Q_*`, where `D` is the original degree. Its degree is
`O(M log q+D)` for profiles of total multiplicity `q`. Taking polynomial
`M` therefore preserves polynomial degree when `D` and `log q` are
polynomial in the physical label width.

This generalizes the reciprocal retiming in
[the Walsh construction](bell-walsh-constructions.md). It does **not** permit
retiming a Smith-order list while retaining the old physical tensor. The
[32-letter coherent crossing](coherent-hot-compound-retractions.md) and
[eight-letter reflection retraction](eight-letter-reflection-retraction.md)
use polynomial Smith factors, rather than constant diagonalizing bases.
Section 5 gives a literal endpoint-preserving counterexample to that invalid
extension. No improved factoring score is claimed.

## 1. The retimed pressure bound

Write a polynomial crossing as `P(z)=sum_(i=0)^D A_i z^i` and assume that
fixed invertible bases diagonalize each of its two cuts. Thus every cut
eigenvalue is a scalar polynomial. Normalize away any common initial
monomial, so at least one row eigenpolynomial has nonzero constant term.
The opposite-cut coefficient of degree `D` is nonzero as well.

Choose integers `0=d_0<d_1<...<d_D` and put

\[
\widetilde P(z)=\sum_{i=0}^D A_i z^{d_i}.
\]

The endpoint is unchanged because `P(1)=sum_i A_i`. The constant
diagonalizing bases still apply. Strict monotonicity changes every scalar
row order `i` to `d_i` and every opposite-cut degree `j` to `d_j`.

Consequently the usual complete order/degree estimate becomes

\[
Q(d,x)=\left(\sum_i a_i x^{-d_i}\right)
        \left(\sum_i b_i x^{d_i}\right),\qquad x\ge1.
\tag{2}
\]

At tensor length `t`, every coefficient has rank product at most
`Q(d,x)^t`. Indeed, the row count at coefficient `j` is at most
`x^j (sum_i a_i x^(-d_i))^t`, while the opposite count is at most
`x^(-j) (sum_i b_i x^(d_i))^t`. Their product gives the claim. This
counts actual integer multiplicities, not field elements.

Retiming changes the polynomial degree and therefore the number of
coefficient pieces. A fixed finite schedule still has linearly many
pieces in `t`. The uniform degree bound below controls varying families.

## 2. Pooling gives the exact real optimum

Temporarily permit real nondecreasing indices `u_i=d_i log x` and write

\[
A(u)=\sum_i a_i e^{-u_i},\qquad
B(u)=\sum_i b_i e^{u_i}.
\]

Adding the same real number to all indices leaves `A(u)B(u)` unchanged.
Choosing that number to make the two sums equal shows

\[
\inf_u A(u)B(u)
=\frac14\left(\inf_u[A(u)+B(u)]\right)^2.
\tag{3}
\]

The second minimization is convex. For a block on which all indices
equal one value, its minimum is at

\[
u_B=\tfrac12\log(A_B/B_B),
\qquad
A_B=\sum_{i\in B}a_i,\quad B_B=\sum_{i\in B}b_i.
\tag{4}
\]

Begin with one block per nonzero profile position. Merge neighboring
blocks whenever the left ratio in (4) exceeds the right ratio. A zero
position can be assigned to a neighboring block. Because `a_0>0` and
`b_D>0`, every final block has positive totals in both profiles.

The resulting ratios are nondecreasing. For every prefix `J` inside
a returned block `B`, they also satisfy the exact integer inequality

\[
\left(\sum_{i\in J}a_i\right)B_B
\ge A_B\left(\sum_{i\in J}b_i\right).
\tag{5}
\]

These are the pooling algorithm's optimality conditions, and the
implementation checks them independently in its replay. They are preserved
by each merge. A prefix inside the left block has ratio at least that
block's old ratio, which is at least the merged ratio. For a prefix
containing the whole left block and part of the right, the right prefix
has ratio at least its block ratio. Combining the entire higher-ratio
left block with at most the full weight of the lower-ratio right block
again gives at least the merged ratio. This argument can be written as
the cross-multiplied inequalities (5), so zero denominators in intermediate
blocks cause no division. Single-position blocks start with the property.

For completeness,
at the values (4), the gradient of `A+B` is
`g_i=b_i e^(u_B)-a_i e^(-u_B)`. Equation (5) says that every gradient
prefix in a block is nonpositive. The whole block gradient is zero.
Taking the negative prefix sums as multipliers for `u_i<=u_(i+1)` gives
nonnegative multipliers, stationarity, and complementary slackness.
Convexity proves global minimality. Substitution in (3) gives (1).

When the individual ratios `a_i/b_i` already increase, no pooling is
necessary, and (1) is simply `(sum_i sqrt(a_i b_i))^2`. In the reciprocal
case `b_i=a_(D-i)`, it recovers the paired formula in the Walsh note.

## 3. Rational integer schedules attain the infimum

Fix a positive integer `M` and let `x=1+1/M`. Let `R_i` be the ratio of
the final block containing `i`, and let `R_0` be the first block ratio.
For every index choose the nearest nonnegative integer

\[
c_i\approx\frac{\log(R_i/R_0)}{2\log x},
\qquad d_i=c_i+i.
\tag{6}
\]

The nearest integer is obtained without numerical logarithms: it is the
least `c>=0` satisfying `x^(2c+1)>=R_i/R_0`. Doubling and binary search
use only comparisons of rational powers. Since the block ratios increase,
the `c_i` are nondecreasing. Thus (6) is strictly increasing and `d_0=0`.

The discrepancy between `d_i log x` and `u_i-u_0` lies in the interval
`[-log(x)/2,(D+1/2)log(x)]`. Applying these two endpoints separately to
the sums in (2) gives

\[
Q(d,x)\le x^{D+1}Q_*.
\tag{7}
\]

Letting `M` increase proves that the real optimum is also the infimum
over actual integer schedules and rational scales. If the two profile
totals are `q`, every positive block ratio lies between `1/q` and `q`.
Consequently the schedule degree is `O(M log q+D)`. All exponents and
pressures in the retained construction are exact rationals or integers.
The original coefficient bit lengths are unchanged.

## 4. Exact replay

```sh
python3 tools/check_isotonic_coefficient_retiming.py \
  --output evidence/isotonic-coefficient-retiming.json
```

The replay checks 192 independent profile pairs, including zero bins.
For each, it verifies the block totals, ratio ordering, and every gradient
prefix inequality. Integer square roots give rational lower bounds on
the optimum, against which the constructive approximation is checked.
Literal two-copy scalar products verify every retained diagonal rank
bound for eighteen general paths and one nonconstant degree schedule.

The nonconstant example uses the physical binary Bell coefficients
`(0,0,0,2)` and `(1,1,1,2)` over `F_3`. The two scalar profile lists are
`(1,3)` and `(3,1)`. The schedule `(0,6)` at scale `5/4` has exact pressure
`779135569/64000000<16`. This checks that the construction is being applied
to actual coefficient tensors, not merely to formal profile entries.

## 5. A Smith-order counterexample and the near-target seeds

Over `F_2`, consider

\[
M(z)=\begin{pmatrix}
z&1+z\\
z^2+z^3&z+z^2+z^3
\end{pmatrix}.
\tag{8}
\]

It has `M(1)=I`, a unit entry at zero, and determinant `z^3`. Its complete
Smith orders are therefore `(0,3)`.

Retiming the coefficient indices by `(0,1,3,4)` preserves `M(1)=I` and
produces

\[
\widetilde M(z)=\begin{pmatrix}
z&1+z\\
z^3+z^4&z+z^3+z^4
\end{pmatrix}.
\tag{9}
\]

Now the determinant is `z^2+z^3+z^4`, so the complete Smith orders are
`(0,2)`. Retiming the old order list would falsely predict `(0,4)`.
The nonlinear schedule has broken the cancellation between two determinant
terms of the original degree two.

The coherent and reflection seeds near the critical parameter use
polynomial triangular factors and normalized physical height gauges.
Their displayed Smith profiles do not meet the constant-basis hypothesis
of Section 1. A proposed retiming of those seeds must therefore reconstruct
the physical coefficient tensor, retain its endpoint, recompute both
normalized cuts, and account for its changed degree and height charges.
Applying (1) directly to their old Smith profiles would not establish a
new rank bound. This distinction is an algebraic limitation of the method,
not evidence that no physical retiming can improve those constructions.

## 6. Complete physical retimings of the eight-letter seed

The replay also applies the transformation to the actual eight-letter
reflection tensor, rather than to its profile list. Let `S(z)` be the
degree-32 nonorthogonal tensor in the preceding reflection note. Begin
with the neutral dilation `S(z^2)`, coefficient schedule `d_i=2i`, and
twice its natural physical heights. Its two complete Smith profiles
are exactly twice the original profiles, so its optimized pressure
parameter is unchanged.

For each `i=1,...,15` and each `epsilon` in `{-1,1}`, change exactly

\[
d_i=2i+\epsilon,\qquad
d_{32-i}=2(32-i)-\epsilon.
\tag{10}
\]

All remaining indices retain `d_j=2j`. These thirty schedules are strictly
increasing, keep total degree 64, and retain the exact crossing endpoint.
They are a finite declared family, not all coefficient retimings.

The implementation reconstructs all 4,096 entries of every candidate.
It then forms both physical cuts independently, clears any Laurent
monomial only for the Smith calculation, restores that monomial in the
reported orders, and computes every Smith order. Across the reference
and candidates this checks 126,976 endpoint entries and 3,968 complete
orders. Both cuts use the same doubled physical-height assignment as
the reference. The two common scalar shifts are `-32` and `+32`, and
their cancellation is retained explicitly.

**Every one of the thirty candidates has strictly larger complete
two-cut pressure than the reference at every finite nonzero scale.**
This conclusion is exact and does not use a numerical minimizer.

To prove it, form the integer counting functions for the two fully
recomputed order lists. Subtract the reference product from the candidate
product. Each retained difference has the form

\[
r^a(1-r)Q(r),\qquad 0<r<1,
\tag{11}
\]

where `a` may be negative because the reported orders include restored
scalar shifts. Every `Q` has a nonzero expansion

\[
Q(r)=\sum_{j=0}^{N}c_j r^j(1-r)^{N-j},
\qquad c_j\ge0.
\tag{12}
\]

The checker computes each integer coefficient `c_j` and reconstructs
the original power-basis polynomial from (12). Thus (11) is strictly
positive throughout the stated interval. The full order lists, degree
schedules, and quotient polynomials are saved in the evidence.

Because the physical heights agree, their type-restoration terms agree
for every common physical type. The thirty strict comparisons therefore
remain valid after type optimization at the retained heights. Other
height gauges, simultaneous changes to several coefficient pairs, and
physical retiming of the 32-letter coherent seed are not covered by
this finite result.

These results use ordinary algebra and exact finite controls with
substantial OpenAI assistant support. No new Lean formalization or
independent proof review is claimed. The established factoring score
remains `1/2`.
