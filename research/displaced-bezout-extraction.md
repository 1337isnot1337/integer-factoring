# Extracting factor polynomials from displaced Bezout data

## Constructive result and remaining obligation

The [fixed regularization theorem](regularized-bezout-factor-moments.md)
turns a correctly supplied global Bezout jet into moments supported
exactly on factor rows. This note gives two further results.

First, the **optimal primitive exponent** itself distinguishes a
rectangle with factor incidence from a rectangle without incidence,
once the displacement exceeds an explicit bound. The large baseline
exponent is stored in binary.

Second, an implemented polynomial-time decoder recovers the exact
integer factor polynomial from `2R` supplied moments, where `R` bounds
the number of factor rows. It uses a declared precision bound and
handles vanishing early Hankel minors. No global jet provider is
implemented here. The finite checker explicitly constructs its small
input polynomials, and the large decoder control starts with supplied
moments. The complete factoring score remains `1/2`.

## 1. Domains, displacement, and the fixed residue polynomial

Let `A` and `B` be positive odd arithmetic progressions with common
step two and respective lengths `m=2^u` and `s=2^v`. Put

\[
 P(X)=\prod_{a\in A}(X-a),\qquad
 Q_t(X)=\prod_{b\in B}(N+t-bX),\qquad t=2^K w,
 \quad w^2+w+1=0.
 \tag{1}
\]

Assume `N` is odd. Choose `H` so that every nonzero `N-ab` has
absolute value less than `2^H`, and enlarge `H` to at least
`u+1,v+1`. Let `V_t` be the reduced inverse of `Q_t` modulo `P`.
The two derivative orders are constant on their domains:

\[
 E_A=2m-2-u,\qquad E_B=2s-2-v,\qquad E=E_A+E_B.
 \tag{2}
\]

For instance, the derivative at row index `i` is
`2^(m-1)(-1)^(m-1-i)i!(m-1-i)!`. Complementary binary digits of
`i` and `m-1-i` give (2).

Let `F={a in A:N/a in B}`, let `D=|F|`, and put

\[
 G(X)=\prod_{a\in F}(X-a),\qquad P_0=P/G.
 \tag{3}
\]

For `a in F`, write `b=N/a` and define

\[
 \alpha_a=2^{-E_A}P'(a),\qquad
 \beta_b=2^{-E_B}P_B'(b),\qquad
 c_a=\frac1{a^{s-1}\alpha_a\beta_b}.
 \tag{4}
\]

Every `c_a` is a rational dyadic unit. Define the fixed polynomial

\[
 W(X)=\sum_{a\in F}c_a\frac{P(X)}{X-a}
     =P_0(X)\sum_{a\in F}c_a\frac{G(X)}{X-a}.
 \tag{5}
\]

For `D=0`, this is the zero polynomial. The global fixed regularization
is

\[
 w\,2^{K+E}V_t(X)\equiv W(X)\pmod {2^p}
 \qquad(K\ge H+p).
 \tag{6}
\]

All polynomial congruences are coefficientwise. In `Y=X-1`, the
coefficient `d` places below top degree in `W(1+Y)` is divisible by
`2^d`: every summand in (5) is a unit times a monic product with even
roots. Thus the complete fixed-regularized polynomial modulo `2^p`
has only its first `p` reverse coefficients possibly nonzero.

## 2. The primitive exponent stabilizes after a factor-count bound

Let `v_c` denote the minimum valuation of the coefficients of a nonzero
polynomial. If `D>0`, set

\[
 \mu=v_c(W).
 \tag{7}
\]

The polynomial `W` is nonzero, since its value at a factor row is
`c_a P'(a)`. Its coefficients are integral, so `mu>=0`.
The monic polynomial `P_0` in (5) has content zero. Gauss content
multiplicativity therefore gives

\[
 \mu=v_c\left(\sum_{a\in F}c_a\frac{G(X)}{X-a}\right).
 \tag{8}
\]

Evaluate the polynomial in (8) at any `a in F`. The value is
`c_a G'(a)`. If `ell` bounds every `v_2(a-a')` for distinct row nodes,
then

\[
 \boxed{0\le\mu\le\min_{a\in F}v_2(G'(a))
                    \le(D-1)\ell.}
 \tag{9}
\]

One can take `ell=max(1,bits(max(A)-min(A)))`; when there is only one
row, any positive bound suffices.

### A coefficientwise remainder bound

There is an exact expansion

\[
 V_t=\frac{W}{t\,2^E}+R_t,
 \qquad v_c(R_t)\ge-E-H+1.
 \tag{10}
\]

For each row `a`, the column domain contains exactly one residue
`b_0` with `ab_0=N mod 2^(v+1)`. For `b!=b_0`, the valuation of
`N-ab` equals that of `b-b_0`. Their total is `E_B`. If the matched
difference is nonzero, its valuation is at most `H-1`; the entire
Lagrange summand for `1/Q_t(a)` consequently has coefficient valuations
at least `-E-H+1`.

At a factor row, separate the factor `t` in `Q_t(a)`. The remaining
product has constant term `a^(s-1)P_B'(N/a)`. Every other factor changes
by a relative quantity of valuation at least `K-H+1`. Inverting their
product and subtracting its constant term leaves the same bound
`-E-H+1` after division by `t P'(a)`. The numerator `P(X)/(X-a)` is
integral. Summing the Lagrange terms proves (10), without a factor
depending on the number of rows.

Let `r_K` be the least nonnegative exponent for a polynomial certificate
`U P+V Q_t=2^(r_K)` over the valuation ring `Z_2[w]`. Reduction modulo
the monic `P` shows that `r_K=-v_c(V_t)`, just as in the
[rectangular primitive theorem](rectangular-bezout-denominator.md).
Equations (9)--(10) imply

\[
 \boxed{\begin{array}{ll}
 D=0:&r_K\le E+H-1\quad(K\ge H),\\[2mm]
 D>0:&r_K=K+E-\mu\quad(K\ge H+\mu).
 \end{array}}
 \tag{11}
\]

In the second line the pole term in (10) has strictly smaller
coefficient valuation than the remainder. Thus, if a supplied bound
`R>=D` is known, then at

\[
 K\ge H+(R-1)\ell,
 \qquad
 \boxed{D>0\ \Longleftrightarrow\ r_K>E+H-1.}
 \tag{12}
\]

For nonempty incidence, the primitive exponent increases by exactly
one when `K` increases by one in this range. The magnitude test (12)
also covers the empty case; it does not assume an equally early
stabilization theorem for the empty rectangle.

The baseline `E` can be exponential in the interval description, but
its binary representation is short. Equations (11)--(12) are not an
algorithm for finding the primitive exponent. They identify precisely
what a fast primitive-scale computation would decide.

## 3. Supplied jets give unit-weight moments

The fixed-scale interface supplies leading coefficients of
`w 2^(K+E) V_t(1+Y)` modulo `2^p`, with `K>=H+p`.
The [moment postprocessor](../tools/regularized_bezout_factor_moments.py)
divides its reversed series by the compressed reversed row product
and changes from `(a-1)^j` to `a^j` moments. By (5)--(6), the output is

\[
 S_j=\sum_{a\in F}c_a a^j\pmod {2^p},\qquad v_2(c_a)=0.
 \tag{13}
\]

This conversion visits no original domain node. It requires only as
many supplied reverse coefficients as the requested moment prefix,
capped by the polynomial degree and the proved precision support.

For completeness, `R` moments at precision `p>(R-1)ell` already decide
whether the support is empty. If all first `D` moments vanished, apply
their linear combination with coefficients from
`prod_(b in F, b!=a)(X-b)`. It equals `c_a G'(a)`, whose valuation is
at most `(D-1)ell`. This contradicts the precision bound. The moment
weights can have signs and odd denominators; no positivity argument
is needed.

## 4. Exact polynomial recovery with an explicit precision budget

Suppose `R>=D`, every factor row is at most `A_max`, and pairwise
row-difference valuations are at most `ell`. Put

\[
 B=R\,\operatorname{bits}(A_{\max}+1)+2,
 \qquad J=R(R-1)\ell,
 \qquad p=B+J+1.
 \tag{14}
\]

Supply `S_0,...,S_(2R-1)` modulo `2^p`. For `d>=1`, let

\[
 H_d=(S_{i+j})_{0\le i,j<d}.
 \tag{15}
\]

At the true degree, the Vandermonde factorization gives

\[
 \det H_D=\left(\prod_{a\in F}c_a\right)
           \prod_{a<a'}(a'-a)^2,
 \qquad
 0\le v_2(\det H_D)\le D(D-1)\ell\le J.
 \tag{16}
\]

For every `d>D`, the exact Hankel determinant is zero. It follows that

\[
 \boxed{D=\max\{d\le R:\det H_d\not\equiv0\pmod {2^p}\},}
 \tag{17}
\]

with maximum zero for an empty set. Earlier principal minors are
allowed to vanish. In particular, stopping at the first zero determinant
would be incorrect: unit weights `1,-1` on two distinct nodes make
`H_1=0` and `H_2` nonsingular.

Write `G(X)=X^D+sum_(i<D)g_i X^i`. The moments satisfy

\[
 H_D\begin{pmatrix}g_0\\ \vdots\\g_{D-1}\end{pmatrix}
 =-\begin{pmatrix}S_D\\ \vdots\\S_{2D-1}\end{pmatrix}.
 \tag{18}
\]

Let `j=v_2(det H_D)`. Cramer's rule loses at most `j<=J` bits:
divide each determinant numerator and denominator by `2^j`, invert
the remaining odd denominator, and recover every coefficient modulo
`2^B`. This is exact modular division, with its precision charged.

The integer coefficients satisfy

\[
 |g_i|\le\prod_{a\in F}(1+a)
          \le(1+A_{\max})^R<2^{B-2}.
 \tag{19}
\]

Their centered residues modulo `2^B` therefore determine them uniquely.
This reconstructs the exact polynomial in (3), rather than a rational
or modular surrogate. The usual interval-root recovery in the
[orbit-gcd reduction](interval-orbit-gcd.md) then extracts factor rows.

## 5. Implementation and cost

The [decoder](../tools/displaced_bezout_extraction.py) provides
`extraction_precision` for (14), `decode_factor_moments` for exact
polynomial recovery, and a fraction-free determinant implementation.
The determinant computation and Cramer solves have polynomial bit cost
in `R`, `p`, and the encoded node bound. Intermediate determinant lifts
have `O(R(p+log R))` bits. The implementation examines all leading
Hankel sizes through `R`; it does not assume that nonzero principal
minors form an initial segment. It checks the recovered recurrence
against every supplied moment and checks the integer coefficient bound.

The large clearing exponent `K+E` is never expanded by the postprocessors.
Its binary value belongs to the provider contract. A provider must supply
the correctly normalized coefficients before reduction; reducing an
unscaled numerator and then dividing by a vanished power of two is invalid.

## 6. Where the global arithmetic cost remains

The matched column for one row is the representative of `N/a` modulo
`2s` in `B`. Its individual construction has polynomial bit cost.
The whole interval can nevertheless contain exponentially many rows.
The high displacement removes every nonfactor row only after the
global coefficients have been normalized correctly.

Equations (5) and (13) identify the resulting scalar sums exactly.
They do not evaluate those sums without locating or aggregating their
factor support. A direct Lagrange construction still visits every row.
A residue-tree construction must retain the correlated reciprocal
matches until the required depth. Neither the short final jet nor
the binary representation of `E` bounds that tree's width.

Computing the global regularized jet is the missing engine.
Normalization, the nonempty-support test, and exact factor-polynomial
recovery now have explicit finite-precision algorithms and bit bounds.
The [one-leading-coefficient continuation](leading-bezout-factoring.md)
uses extra interval padding to remove the factor-count bound from the
decision interface.

## 7. Exact evidence

The [checker](../tools/check_displaced_bezout_extraction.py) supplies a
separate finite provider by full rational interpolation in `Q(w)`.
It verifies (9)--(12), compares the fixed-regularized coefficients
against the exact factor-supported polynomial (5), and runs the
independent moment postprocessor and decoder.

The [receipt](../evidence/displaced-bezout-extraction.json) contains
12 complete small-interval providers, 35 generic unit-weight moment
decodings, a cancellation case with zero first Hankel minor, and a
large-integer-node decoding control. The latter begins with supplied
moments and is explicitly not a large implicit-interval provider.

```sh
python3 tools/check_displaced_bezout_extraction.py --output /tmp/displaced-bezout.json
cmp /tmp/displaced-bezout.json evidence/displaced-bezout-extraction.json
```

These checks exercise the complete conditional extraction path. They
do not close its global coefficient-provider obligation.
