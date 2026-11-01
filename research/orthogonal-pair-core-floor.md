# The first orthogonal paired-code candidate has an uncompetitive natural pressure

## Result and scope

Use the 92-letter F4 orthogonal paired code from
[orthogonal-differential-pair-codes.md](orthogonal-differential-pair-codes.md),
with its displayed coefficient matrix and natural physical heights.
Its balanced two-cut Smith-profile pressure satisfies

\[
\rho_{\rm prof}\ge\frac{32}{45}=0.711111\ldots
>0.710005861>\frac1{\sqrt2}.                         \tag{1}
\]

Thus this candidate cannot beat the certified hyperoval bound by the
balanced natural-profile criterion. This does not exclude other letter
gauges, asymmetric supporting weights, different coefficient matrices,
or changes to the low-height physical maps.

The proof does not require completing the expensive 92-letter profile.
It recomputes both complete profiles of the 43-letter height-at-most-20
core. The full profile run was stopped once this exact all-scale exclusion
settled the natural balanced test. It is not recorded as a completed run.

## 1. An exact shared low core

Take pairs `(0,7),(1,6),(2,5),(3,4)` and

\[
C=\begin{pmatrix}1&1&\omega&\omega\\0&1&1+\omega&1\end{pmatrix},
\qquad \omega^2=\omega+1.
\]

Let V be the full common differential kernel, of dimension 92, and let
V_low retain heights at most 20, of dimension 43. The reduced basis P and
orthogonal retract `L=(P transposed P) inverse P transposed` are formed
separately in height and occupation cells. Therefore the low rows and
columns of both full maps agree exactly with the maps constructed directly
on V_low. The checker verifies this equality.

For each natural cut, the 43-squared channel matrix is consequently a
polynomial submatrix of the 92-squared channel matrix. Let Z_low(t) be its
complete channel partition function, and let Q(t) count all 92 physical
heights. Both low cut profiles agree by the orthogonal Gram identity.

The truncated Toeplitz-rank argument from the
[common-core pressure proof](hot-boundary-natural-pressure-floor.md)
gives `Z_full(t) >= Z_low(t)` for every `0<t<1`, independently at both cuts.
This monotonicity does not require a termwise comparison of individual
Smith orders.

## 2. A finite all-scale certificate

Set `a=32`, `b=45`. The checker proves

\[
Z_{\rm low}(s^{32})^{45}>Q(s^{45})^{64}
\qquad(0<s<1).                                     \tag{2}
\]

Write `Z_low(t)=sum_i z_i t^i` and `Q(t)=sum_j q_j t^j`. For

\[
F(s)=45\log Z_{\rm low}(s^{32})-64\log Q(s^{45}),
\]

the derivative has the sign of the integer polynomial

\[
G(s)=\sum_{i,j}(i-2j)z_iq_j s^{32i+45j-32},          \tag{3}
\]

with zero summands omitted. All surviving exponents are nonnegative.
Its degree is 4816. The coefficients of

\[
(1+u)^{4816}G\left(\frac{u}{1+u}\right)
\]

have exactly two sign variations. Also `G(0)>0`, and exact integer
Horner evaluation gives opposite signs across the box

\[
s_- = \frac{990294031}{10^9},\qquad
s_+ = \frac{990294032}{10^9},\qquad
G(s_-)<0<G(s_+).
\]

Descartes' rule therefore leaves exactly one local minimum of F in this
box, preceded by one local maximum. Since `F(0)=0`, it suffices to prove
that the boxed minimum is positive.

Both partition polynomials have nonnegative coefficients. Outward 192-bit
arithmetic proves

\[
45\log Z_{\rm low}(s_-^{32})-64\log Q(s_+^{45})>0.
\]

The lower dyadic numerator, with denominator `2^192`, is
`528991607314108202477555528305982990079905401095304915735`.
A separate integer-power comparison of the same outward rational bounds
also proves the sign. This establishes (2) on the whole open interval.

## 3. Every common physical type is covered

Substitute `t=s^32` in (2). Then

\[
\log Z_{\rm low}(t)>
2\rho\log Q(t^{1/\rho}),\qquad \rho=32/45.
\]

For a physical type p with mean height E and entropy H, the Gibbs
variational inequality gives

\[
\log Q(t^{1/\rho})\ge H+(E/\rho)\log t.
\]

The two full cut partition functions each dominate Z_low. Hence their
balanced natural-pressure estimate is at least 32/45 for every scale and
every positive-entropy common physical type. Taking the infimum proves
(1). The argument applies over extensions of F4 without changing its ranks.

## 4. Reproduction

The [checker](../tools/check_orthogonal_pair_core_floor.py) reconstructs
both core profiles, checks agreement of the physical maps, and repeats
all polynomial, root-box, logarithm and integer-power comparisons. The
[receipt](../evidence/orthogonal-pair-core-floor.json) retains the exact
profiles, height inventory, polynomial coefficients, transformed signs,
source hashes and the explicit incomplete-full-profile status.

This is a limitation of one tested construction and its stated pressure
method, not a lower bound for factoring algorithms.
