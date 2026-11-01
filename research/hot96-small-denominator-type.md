# A small exact common type for the q96 certificate

The strict bound rho<710005874881191/10^15 also holds with explicit
positive common physical type counts of denominator
2^36=68,719,476,736. The counts are in
[`hot96-small-denominator-type.json`](../evidence/hot96-small-denominator-type.json).
This retains the certified bound while avoiding the much larger
denominator of its exact Gibbs probabilities.

Write A/B=2rho, z for the rational base, P=sum z^(A sigma), and
Q=sum z^(B E), as in [the interval certificate](exact-natural-pressure-intervals.md).
The exact interval computation gives a positive lower bound G for
A log Q-B log P. Its ideal Gibbs type is g_i=z^(B E_i)/Q.
For any positive rational type p, put t=A log(1/z). The entropy identity is

\[
H(p)=\log Q+\frac{t}{A/B}\mathbb E_p E-D(p\Vert g).
\]

The full natural pressure at p is
(2 log P+2t E_p E)/(4H(p)). Thus it stays strictly below rho if

\[
A D(p\Vert g)<G. \tag{1}
\]

Both physical means are included. The inequality log x<=x-1 gives

\[
D(p\Vert g)\le\sum_i\frac{(p_i-g_i)^2}{g_i}. \tag{2}
\]

Use outward power bounds and opposite denominator endpoints to compute
exact intervals 0<L_i<=g_i<=U_i. Set K=2^36. For each nonvacuum letter,
take c_i=floor(K L_i). Give the remaining mass to the vacuum:
c_0=K-sum_(i>0)c_i. Every saved count is positive and the sum is exactly K.
For p_i=c_i/K, the right side of (2) is at most

\[
\Delta=\sum_i\frac{\max(|p_i-L_i|,|p_i-U_i|)^2}{L_i}.
\]

The checker proves G-A Delta>0 using rational arithmetic. The remaining
margin is approximately 1.738226575; the stored rational numerator is
strictly positive. The minimum count is 2,813,121, so this construction
does not approach a zero-type face.

The [checker](../tools/certify_hot_rational_type_rounding.py) first repeats
the complete interval certificate and checks its source and input hashes.
It then reconstructs every count and the rational bound above. No
floating-point inequality is used. The finite-type and all-width proof
now uses the explicit common denominator K. The original exact Gibbs
certificate remains unchanged.

```sh
python tools/certify_hot_rational_type_rounding.py
```
