# Symmetric multilevel filters and exact natural-charge gains

## Results and scope

An occupation-symmetric filter on `r` selected binary identity modes can
be synthesized in `r+1` polynomial coordinates. The full physical
construction has `2^r` occupation states. This reduction permits
multilevel filters and arbitrary nonzero endpoint colors; it does not
assume a homogeneous Boolean pattern.

The reduced transform also proves that filters on disjoint groups of
selected modes combine by multiplying their exact cost ratios. Thus
small correlated blocks can make a change proportional to the interval
length, while keeping every polynomial degree and both physical costs.

A literal four-mode filter gives a strict cost reduction relative to the
Gaussian at its natural charge `J=0`, as well as at an interior charge.
This disproves an extension of homogeneous Gaussian domination to all
polynomial filters. The full critical pressure of this witness is
positive. Its positive Gibbs root lies between `0.966` and `0.967`, far
above `1/sqrt(2)`. There is no improved factoring exponent.

## 1. Compress the complete transform by occupation weight

Let `n` be a power of two, let `P` be positive and odd, set `m=nP`, and
select any `r<=n` distinct identity modes. For `k in {0,1}^r`, put

\[
f_k(x)=\alpha_{\mathbf1-k}(x)=f_{|k|}(x),
\qquad C=1+x^m.
\]

The complete [joint transform](joint-binary-filter-synthesis.md) becomes

\[
\begin{aligned}
\beta_k(x)
 &=\sum_{\ell\le k}C^{|k|-|\ell|}f_{|\ell|}(x)\\
 &=\sum_{h=0}^{w}\binom wh C^{w-h}f_h(x)
 =:\beta_w(x),\qquad w=|k|.
\end{aligned}
\tag{1}
\]

All coefficients are reduced modulo two. Lucas's binary identity says
that `binomial(w,h)` is odd exactly when every set bit of `h` is also a
set bit of `w`. Consequently

\[
\boxed{\beta_w(x)=\sum_{h\subseteq_{\rm bit}w}
                 (1+x^m)^{w-h}f_h(x).}
\tag{2}
\]

Equation (1) is obtained by grouping the full tensor sum. It therefore
preserves every occupation state's actual polynomial, rather than only
an averaged spectrum.

For weights `0<=w<2^a`, the matrix in (2) is the tensor product of the
`a` matrices

\[
\begin{pmatrix}1&0\\1+x^{m2^j}&1\end{pmatrix},
\qquad 0\le j<a.
\tag{3}
\]

Indeed, `C^{2^j}=1+x^{m2^j}` in characteristic two. A non-power-of-two
weight range is the corresponding initial principal block. Each
matrix in (3) is its own inverse, and the full weight transform is
also its own inverse.

## 2. Exact polynomial synthesis and endpoints

Prescribe one degree bound `D_w` and one beta-order target `B_w` for
each weight `0<=w<=r`. Because (2) has coefficients in `F_2[x^m]`,
decompose every `f_w` by residue modulo `m`. For a residue `eta`, the
local bounds are

\[
d_w(\eta)=\left\lfloor\frac{D_w-\eta}{m}\right\rfloor,
\qquad
b_w(\eta)=\max\left(0,\left\lceil\frac{B_w-\eta}{m}\right\rceil\right).
\tag{4}
\]

There are at most `2(r+1)+1` distinct bound patterns, since the only
changes occur at `D_w mod m + 1` or `B_w mod m`. Each local system
uses `r+1` polynomial coordinates and degrees measured in units of `m`.
Negative local degree bounds mean that the corresponding polynomial
vanishes in that residue.

The sum of the local endpoint images is a binary vector space
`E subset F_2^(r+1)`. Choose a fixed binary extension field of order
`Q>r+1`, for example the smallest such power of two. Then an endpoint
vector with every coordinate nonzero exists over that field exactly
when every coordinate projection of `E` is nonzero. Necessity is
immediate. For sufficiency, the `r+1` coordinate-zero hyperplanes cannot
cover `E tensor F_Q`, since `r+1<Q`.

The [implementation](../tools/symmetric_binary_filter_synthesis.py)
retains coefficient lifts during exact binary elimination and constructs
literal polynomials. At `x=1`, equation (1) gives

\[
\beta_w(1)=f_w(1)\ne0.
\tag{5}
\]

Every returned beta polynomial therefore has finite order. A constant
physical endpoint repair uses the inverse of `f_w(1)` and has no degree
or Smith-order cost, as in the joint identity theorem. In the explicit
counterexample below, every endpoint is already one.

## 3. Cost evaluation uses the actual selected modes

Let the actual positive odd mode energies be
`lambda_i=2n-1-2j_i`, and put `t=P log s`. Define the coefficients

\[
\prod_{i=1}^r(1+e^{-t\lambda_i}y)=\sum_{w=0}^r c_w(t)y^w.
\tag{6}
\]

Then the normalized occupation-weight distributions at the two cuts are

\[
p_w=\frac{c_w(t)}{\sum_h c_h(t)},\qquad
q_w=\frac{c_w(t)e^{-tJw}}{\sum_h c_h(t)e^{-tJh}}.
\tag{7}
\]

For the full degrees `D_w=deg f_w` and exact orders
`B_w=ord_x beta_w`, the complete relative cost is

\[
\boxed{\mathcal R=
 \left(\sum_wq_w e^{2t(nw-B_w/P)}\right)
 \left(\sum_wp_w e^{2tD_w/P}\right).}
\tag{8}
\]

The coefficients in (6) keep each selected mode's real physical energy.
Replacing them by a binomial distribution is generally incorrect unless
all those energies are equal. The compressed profile is symmetric; the
physical occupation measure need not be.

## 4. Disjoint blocks tensorize exactly

Partition selected modes into disjoint blocks `A` and `B`. They may
use different filters, different block sizes, and different endpoint
colors. Define

\[
f_{k_A,k_B}(x)=f^A_{k_A}(x)f^B_{k_B}(x).
\tag{9}
\]

The full transform is a tensor product, so its finite sum factors:

\[
\beta_{k_A,k_B}(x)=\beta^A_{k_A}(x)\beta^B_{k_B}(x).
\tag{10}
\]

Over a field, the nonzero leading and lowest coefficients of a product
do not cancel. Hence

\[
B_{k_A,k_B}=B^A_{k_A}+B^B_{k_B},\qquad
D_{k_A,k_B}=D^A_{k_A}+D^B_{k_B}.
\tag{11}
\]

Both occupation measures are independent across the blocks. Equations
(8)--(11) therefore give the exact identity

\[
\boxed{\mathcal R_{A\cup B}=\mathcal R_A\mathcal R_B.}
\tag{12}
\]

Every infinity degree is included in (11). Each factor in (12) uses its
own block's actual mode energies.

A source block constructed at `m_0=rP_0` can be transported to an
ambient interval with `m=nP` whenever `m_0` divides `m`: substitute
`x -> x^(m/m_0)` in every polynomial. Both orders and full degrees
are multiplied by `m/m_0`, and the coefficient `1+x^m_0` becomes
`1+x^m`. For power-of-two block sizes dividing `n`, taking a common
odd multiple of the source ramifications supplies this divisibility.
Different blocks' additive coefficient models can be realized in one
sufficiently large fixed binary field.

The independent checker expands a literal four-mode filter made from
two different two-mode polynomial blocks. It verifies all sixteen
transformed polynomials, all orders, all full infinity degrees, all
endpoints, and fifteen exact rational cost factorizations. The blocks
use physical energies `(7,1)` and `(5,3)`, so the control does not rely
on equal energies or on adjacent mode indices.

## 5. An exact natural-charge relative gain

Take `n=4`, `P=1`, all four modes, and

\[
f_0(x)=1+x^8+x^{16},\qquad f_1=f_2=f_3=f_4=1.
\tag{13}
\]

The complete transformed polynomials are

| Occupation weight | `f_w(x)` | `beta_w(x)` | Full degree | Exact beta order |
|---|---|---|---:|---:|
| 0 | `1+x^8+x^16` | `1+x^8+x^16` | 16 | 0 |
| 1 | `1` | `x^4+x^8+x^12+x^16+x^20` | 0 | 4 |
| 2 | `1` | `x^24` | 0 | 24 |
| 3 | `1` | `x^8+x^24+x^28` | 0 | 8 |
| 4 | `1` | `x^8+x^24+x^32` | 0 | 8 |

All endpoints are one. The identity at weight two is particularly
useful: a polynomial beta, which becomes a single high-order monomial,
is obtained while all positive-weight `f_w` have degree zero.

Set `exp(-t)=49/50`. Every expectation in (8) is then rational at
integer charge. The [exact certificate](../evidence/symmetric-binary-filter-synthesis.json)
proves `mathcal R<1` at both `J=0` and `J=8` by integer comparison of
the numerator and denominator of that rational number. Outward-integer
logarithm and Gaussian-pressure bounds give:

| Charge | Approximate `log mathcal R` | Approximate full critical pressure | Certified positive Gibbs root |
|---:|---:|---:|---|
| 0 | -0.06063424885 | 2.87577318071 | `0.966 < rho < 0.967` |
| 8 | -0.07006014850 | 2.87394296510 | `0.966 < rho < 0.967` |

The displayed decimals summarize exact rational and interval certificates.
The true root is bracketed by evaluating the complete pressure, including
the physical partition, at both rational endpoints. The physical entropy
is positive, so pressure decreases strictly with `rho`.

The support indicator of (13) is the zero-occupation indicator
`product_i(1+k_i)`, which has terms of several degrees and a constant
term. Therefore this example does not satisfy the hypotheses of the
[global homogeneous exclusion](global-homogeneous-boolean-filter-floor.md).
It shows why that theorem's homogeneity condition cannot be dropped
when asserting relative Gaussian domination. It does not refute an
absolute pressure floor for a larger family.

## 6. The relative gain scales to every larger interval

For any ambient power of two `n>=4`, partition the modes into groups
of four. In each group use

\[
f_0=1+x^{2m}+x^{4m},\qquad f_1=f_2=f_3=f_4=1.
\tag{14}
\]

Its normalized order and degree vectors are

\[
B/m=(0,1,6,2,2),\qquad D/m=(4,0,0,0,0).
\tag{15}
\]

Put `tau=nt`, and fix a scaled charge `j=J/n`. At `tau=0`, every
four-mode occupation distribution is uniform. Differentiating (8) at
zero scale gives

\[
\left.\frac{d}{d\tau}\log\mathcal R_{\rm block}\right|_{\tau=0}
=2\left(2-\frac{25}{8}+\frac14\right)=-\frac74.
\tag{16}
\]

The derivatives of the normalized occupation probabilities sum to zero,
so neither the actual physical energies nor `j` enter this first
derivative. Combining the `n/4` blocks gives

\[
\boxed{\frac1n\log\mathcal R_{\rm tiled}
=-\frac7{16}\tau+O_j(\tau^2).}
\tag{17}
\]

For bounded `j`, the remainder is uniform in `n` and the partition:
every normalized energy lies in `[0,2]`, and the finite four-mode
expressions have uniformly bounded second derivatives on a compact
scale interval. Thus these are strict natural-charge relative gains
at all admissible lengths, and their logarithms are proportional to
`n`. This remains a relative statement. The positive Gaussian
pressure at small scale is also proportional to `n` and is not
overcome by (17).

## 7. Replay and limits

```sh
python3 tools/check_symmetric_binary_filter_synthesis.py --output /tmp/symmetric-filter-replay.json
```

The checker compares 144 residue-compressed feasibility decisions with
independent full-state coefficient systems, replays 17 feasible literal
witnesses in 78 full occupation states, checks baseline constructions
through 64 modes, and verifies the tensor and pressure certificates
described above. These finite controls support the exact formulas; the
grouping and tensor-product arguments prove them for arbitrary size.

The synthesis routine decides each stated finite bound problem exactly.
It does not enumerate every possible degree/order profile. Finite
unsuccessful pressure searches cannot exclude this multilevel family,
and a relative gain alone does not establish an improved factoring bound.
