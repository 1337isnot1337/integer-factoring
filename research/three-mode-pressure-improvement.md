# A certified three-mode pressure improvement

The raw polynomial crossing on the mode set \(S=\{0,1,2\}\), with ambient
size \(M=4\), admits

\[
\boxed{\rho<\frac{741}{1000}=0.741.}
\]

This improves the previous Gaussian numerical control
\(\rho\approx0.7602278114\). It remains above \(1/\sqrt2\), so it does not
improve the factoring exponent or the speedrun score. The result is an exact
pressure certificate for the crossing construction, not a new factoring
algorithm.

The improvement comes from allowing different occupation densities for
every mode and for the two fermion species. The three occupations of one
species are substantial; those of the other species are positive but small.
The fixed-mode endpoint theorem permits precisely this asymmetric choice.

## 1. Crossing and allowed restrictions

Use the raw polynomial kernel \(A_S(w)\) from
[the down-set construction](schur-downset-crossing.md), with \(z=w^2\) and
\(S=\{0,1,2\}\subset[0,4)\). Its rows and columns are ordered by the four
groups \(B,X,V,Y\), followed by the mode order \(0,1,2\).

In a tensor power, restrict each physical word to prescribed occupation
counts for each \(A\)-mode and each \(B\)-mode, independently. The
[fixed-mode theorem](fixed-mode-crossings.md) proves that the projected
endpoint is still a crossing. In particular, the two species need not have
equal cardinality or equal total energy. If the six mode probabilities are
\(p_i\), the logarithmic dimension per copy is

\[
h(p)=\sum_{i=1}^{6}H(p_i),
\qquad H(p)=-p\log p-(1-p)\log(1-p).
\]

The polynomial loss from choosing exact occupation counts does not affect
the exponential rate.

At zero use the row cut, with groups \((X,Y)\mid(V,B)\). At infinity use
the column cut, with groups \((X,B)\mid(V,Y)\). The gauges in the certificate
make every within-cut Gaussian block regular at its endpoint. The remaining
cross blocks may be Laurent matrices. Their exterior-compound Smith
profiles give the pressure product below.

## 2. Exact gauge certificate

The complete input is
[three-mode-pressure-refinement.json](../evidence/three-mode-pressure-refinement.json).
It uses ramification \(w=t^{4096}\). Each of its two height vectors lists
the row heights for \(B,X,V,Y\), followed by the column heights in the same
order. All entries are integers.

The four cross blocks have the following sorted one-particle orders:

| Endpoint and direction | Orders |
| --- | --- |
| Zero, first direction | \(4084,12276,20468,45044,53236,54442\) |
| Zero, reverse direction | \(4110,12302,20494,34356,50740,58932\) |
| Infinity, first direction | \(36510,44702,52894,69278,77470,85662\) |
| Infinity, reverse direction | \(-20124,-11932,-3740,4452,12644,20836\) |

The corresponding physical mean weights, in \(A\)-mode order followed by
\(B\)-mode order, are

\[
G=(-32054,\ 714,\ 33482,\ 161818,\ 136036,\ 103268).
\]

These are the sums of all relevant physical heights at the two endpoints.
The checker reconstructs them from the input; it does not trust this
displayed list. For each mode \(a\),

\[
\begin{aligned}
G_A(a)&=\sum_{\text{ends}}
  \bigl(h_{V,+}(a)+h_{Y,+}(a)+h_{B,-}(a)+h_{X,-}(a)\bigr),\\
G_B(a)&=\sum_{\text{ends}}
  \bigl(h_{B,+}(a)+h_{X,+}(a)+h_{V,-}(a)+h_{Y,-}(a)\bigr).
\end{aligned}
\]

## 3. Rational Gibbs inequality

Let

\[
R=\frac{741}{250},\qquad
r=\frac{281475046555937}{281474976710656}>1,\qquad s=r^{741}.
\]

Let \(\kappa_1,\ldots,\kappa_{24}\) denote the four displayed profiles.
The exact certificate verifies

\[
\left[\prod_{j=1}^{24}(1+r^{-741\kappa_j})\right]^{250}
<
\left[\prod_{i=1}^{6}(1+r^{-250G_i})\right]^{741}. \tag{1}
\]

Every number in this inequality is rational. Choose

\[
p_i=\frac{r^{-250G_i}}{1+r^{-250G_i}}.
\]

These probabilities are also rational and strictly between zero and one.
Numerically they are approximately

\[
(0.8796,\ 0.4890,\ 0.1114,\ 0.0000437,\ 0.000216,\ 0.00165).
\]

Write

\[
L=\sum_j\log(1+s^{-\kappa_j}),\qquad
Z_R=\prod_i(1+s^{-G_i/R}).
\]

Inequality (1) is exactly \(L<R\log Z_R\). The Gibbs identity gives

\[
R h(p)=R\log Z_R+(\log s)\sum_i p_iG_i.
\]

Consequently the two-cut pressure numerator is strictly below \(Rh(p)\).
The fixed-mode transfer therefore yields

\[
\rho
\le
\frac{L+(\log s)\sum_i p_iG_i}{4h(p)}
<
\frac{R}{4}
=\frac{741}{1000}.
\]

Ramification clears gauge denominators and changes the polynomial degree by
a fixed factor. It does not change the asymptotic rate in the tensor-power
parameter.

## 4. Exact verification

The checker is
[check_three_mode_pressure_certificate.py](../tools/check_three_mode_pressure_certificate.py).
It performs four separate checks.

1. It reconstructs the polynomial kernel over \(\mathbb F_2\), substitutes
   the integer heights, and checks every nonzero within-cut entry for
   regularity at the appropriate endpoint.
2. It computes the four weighted Smith profiles with the exact
   [separation oracle](weighted-smith-separation-oracle.md).
3. Independently, it enumerates every nonzero minor of each \(6\times6\)
   cross matrix and checks every determinantal prefix minimum against the
   returned profile.
4. It verifies (1) using outward-rounded dyadic intervals with 256-bit
   precision. Every multiplication, inversion, exponentiation, and final
   comparison uses integers. The upper endpoint on the left is strictly
   less than the lower endpoint on the right.

The interval routines also pass independent comparisons against exact
fractions for positive and negative powers and a product of powers.
Floating-point search output is not used to certify the inequality.

Reproduce the stronger certificate with:

~~~sh
python tools/check_three_mode_pressure_certificate.py \
  --certificate evidence/three-mode-pressure-refinement.json \
  --output evidence/three-mode-pressure-refinement-verification.json
~~~

The earlier, weaker exact checkpoint \(\rho<19/25=0.76\) is retained in
[three-mode-pressure-certificate.json](../evidence/three-mode-pressure-certificate.json).
It can be replayed with:

~~~sh
python tools/check_three_mode_pressure_certificate.py \
  --output evidence/three-mode-pressure-verification.json
~~~

Both verification reports record the full profiles, the exact interval
endpoints, the strict integer gap, and hashes of the relevant source files.

## 5. Search scope and next question

The exploratory
[mode-pressure driver](../tools/search_schur_mode_pressure.py) alternates a
fixed-type convex gauge search with an independent Gibbs type update. Its
minor valuations are exact, but its optimization uses floating-point
arithmetic. A temporary finite search box prevents incomplete minor
relaxations from diverging; feasibility is checked before accepting an
iterate. The search does not prove global optimality.

The best saved numerical point is approximately \(0.7409183468\). The
displayed rational certificate deliberately leaves a margin. It proves only
the stated strict bound \(0.741\).

The small positive \(B\)-occupations suggest comparing the boundary where
that species is absent and embedding the asymmetric pattern in larger mode
sets. Those are constructive extensions of this certificate. Neither this
note nor the earlier scoped sparse-bulk floor excludes them.
