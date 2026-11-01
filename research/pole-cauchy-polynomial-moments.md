# Polynomial inverse moments of implicit pole Cauchy matrices

## Constructive result

The [individual inverse-entry oracle](pole-cauchy-conditioning.md)
evaluates one entry of a large pole Cauchy inverse. This note sums one
or both complete matrix axes when the input weights are polynomials of
bounded degree. The domains remain unions of binary-encoded odd
progressions; their numerical sizes are never enumerated.

Let `A,B` be equally sized distinct positive odd node sets, each given
as a union of `c` disjoint odd progressions. Let `m=|A|=|B|`, let `N` be
odd, and set

\[
 z=N+2^k w,\qquad t=2^k w,\qquad w^2+w+1=0,\qquad
 M(a,b)=\frac{t}{z-ab}.
 \tag{1}
\]

As before, choose `k>=h>=1` with `|N-ab|<2^h` on the rectangle. Given
integer polynomials `f` and `g` of degrees `D,E`, the following have
bit cost polynomial in `c`, the binary input lengths, `k`, the requested
precision, `D`, and `E`:

- one entry of `M^(-1) f(A)`;
- the full scalar `g(B)^T M^(-1) f(A)`.

The second operation sums both matrix axes. For example,

\[
 \boxed{\displaystyle
 \boldsymbol 1^T M^{-1}\boldsymbol 1
 =\frac zt\left(1-\frac{\prod_{a\in A}a\prod_{b\in B}b}{z^m}\right).
 }
 \tag{2}
\]

No balance or conditioning hypothesis is needed for these polynomial
moments. For every pair of integer polynomials,

\[
 2^k g(B)^T M^{-1}f(A)
 \quad\text{is two-adically integral}.
 \tag{3}
\]

This does not yet evaluate the required Hadamard-square trace
`tr(M^(-1)(M o M))`. Its weights depend jointly on both nodes through
`(z-ab)^(-2)`. A polynomial representation of unbounded degree is not
covered by a polynomial-degree cost bound. The factoring score is
unchanged.

The implementation is
[pole_cauchy_polynomial_solve.py](../tools/pole_cauchy_polynomial_solve.py).

## 1. A rational interpolation formula

Put `x_a=z/a`, and introduce monic polynomials

\[
 P_X(Y)=\prod_{a\in A}(Y-x_a),\qquad
 P_B(Y)=\prod_{b\in B}(Y-b).
 \tag{4}
\]

All `x_a,b` are two-adic units. Define the expansion at zero

\[
 R(Y)=\frac{P_B(Y)}{P_X(Y)}=\sum_{j\ge0}r_jY^j,
 \qquad R_{\le d}(Y)=\sum_{j=0}^d r_jY^j.
 \tag{5}
\]

Let `u_d=M^(-1)(a^d)_(a in A)`, and encode its entries by

\[
 F_d(Y)=\sum_{b\in B}\frac{u_d(b)}{Y-b}.
 \tag{6}
\]

The defining matrix equation is

\[
 \frac t a F_d(z/a)=a^d.
 \tag{7}
\]

The unique rational function with denominator `P_B`, numerator degree
at most `m-1`, and the values (7) is

\[
 \boxed{\displaystyle
 F_d(Y)=\frac{z^{d+1}}{tY^{d+1}}
       \left(1-\frac{P_X(Y)}{P_B(Y)}R_{\le d}(Y)\right).
 }
 \tag{8}
\]

To check this, the numerator `P_B-P_X R_(<=d)` vanishes to order at
least `d+1` at zero by (5). After dividing by `Y^(d+1)`, it is a
polynomial of degree at most `m-1`. At each `Y=x_a`, the second term in
the parentheses vanishes and (7) follows. The Cauchy interpolation
matrix is invertible: its nodes are distinct, and `z-ab` is nonzero in
`Q(w)` for every integer pair.

Taking the residue of (8) at `Y=b` gives

\[
 u_d(b)=-\frac{z^{d+1}}t
       \frac{P_X(b)R_{\le d}(b)}{b^{d+1}P_B'(b)}.
 \tag{9}
\]

In particular,

\[
 \boxed{\displaystyle
 u_0(b)=\frac{(-1)^{m+1}}{t z^{m-1}}
    \frac{\prod_{v\in B}v}{b}
    \frac{\prod_{a\in A}(z-ab)}{P_B'(b)}.
 }
 \tag{10}
\]

If `S(Y)=R(Y)/r_0` and `S_(<=d)` is its degree-`d` truncation, then

\[
 u_d(b)=u_0(b)\left(\frac zb\right)^d S_{\le d}(b).
 \tag{11}
\]

Every factor multiplying `u_0(b)` in (11) is integral. Thus the known
valuation of (10) is also a valid lower bound for every polynomial
combination with integer coefficients. It need not be the exact
valuation after cancellation.

## 2. The whole bilinear form uses only two short jets

Define the expansion at infinity, written in `U=1/Y`, by

\[
 T(U)=\frac{\prod_{a\in A}(1-(z/a)U)}
            {\prod_{b\in B}(1-bU)}
      =\sum_{q\ge0}t_qU^q.
 \tag{12}
\]

Then `P_X(Y)/P_B(Y)=T(1/Y)`. On the other hand, expanding (6) at
infinity gives

\[
 F_d(Y)=\sum_{e\ge0}\frac{\mu_{e,d}}{Y^{e+1}},\qquad
 \mu_{e,d}=\sum_{b\in B}b^e u_d(b).
 \tag{13}
\]

Taking the coefficient of `Y^(-e-1)` in (8) proves

\[
 \boxed{\displaystyle
 \mu_{e,d}=\frac{z^{d+1}}t
   \left(\mathbf 1_{e=d}
     -\sum_{j=\max(0,d-e)}^d r_j t_{e-d+j}\right).
 }
 \tag{14}
\]

For all `d<=D` and `e<=E`, only `r_0,...,r_D` and `t_0,...,t_E`
are required. Their coefficients are integral because all the
denominator constant terms are units. The prefactor has valuation
`-k`, proving (3). Arbitrary `f,g` follow by the corresponding finite
linear combination of (14).

At `d=e=0`, one has `t_0=1` and

\[
 r_0=\frac{P_B(0)}{P_X(0)}
     =\frac{\prod_{a\in A}a\prod_{b\in B}b}{z^m}.
 \tag{15}
\]

This yields (2), including its sign.

## 3. Computing the short jets without listing the nodes

The normalized origin series is

\[
 S(Y)=\frac{\prod_{b\in B}(1-Y/b)}
            {\prod_{a\in A}(1-aY/z)}.
 \tag{16}
\]

Its column numerator can be computed as `prod(b-Y)/prod(b)`. For a
progression `b=b_0+2j`, its factors form an affine progression of jets
with unit base `b_0-Y` and even step `2`. Its row denominator has unit
base `1-a_0Y/z` and even step `-2Y/z`.

Likewise, the numerator of (12) is `prod(a-zU)/prod(a)`, with base
`a_0-zU` and even step `2`. Its denominator has base `1-b_0U` and
even step `-2U`. Each product therefore uses the existing
[unit-progression jet algorithm](normalized-interval-jets.md).

For clarity, the finite-precision identity behind that algorithm is

\[
 \prod_{j=0}^{L-1}(u+jv)
 =u^L\sum_{q=0}^{p-1}e_q(0,1,\ldots,L-1)(v/u)^q
 \pmod{2^p},
 \tag{17}
\]

when `u` is a unit jet and every coefficient of `v` is even. Terms
with `q>=p` vanish. The elementary symmetric functions in (17) are
computed from exact integer power sums and Newton identities before
reduction modulo `2^p`. Their bit lengths and computation costs are
polynomial in `p` and `log(L+1)`. Jet multiplication and inversion have
polynomial cost in the truncated degree and precision. The power
`u^L` uses binary exponentiation.

The implementation now uses
[precision-bounded Newton arithmetic](modular-interval-jet-arithmetic.md)
for these products. Its scalar coefficients have `O(p)` bits rather
than the `O(p log(L+1))` bits of the exact integers. Divisions by powers
of two retain an explicit precision budget; odd divisors are inverted.
Formal-variable nilpotence also shortens (17) whenever the even step
has positive formal degree. This changes the arithmetic cost, not the
identity or its scope.

The constant ratio (15) has the same unit-progression construction.
Thus (14) computes the whole bilinear form with polynomial cost in
the binary description, `p,D,E`, and the number of progression pieces.
No determinant, matrix, or node permutation is expanded.

For a single solution entry (10), the only nonunit products are
`prod_a(z-ab)` and `P_B'(b)`. Split each affine progression into its
exact binary valuation classes and normalize every nonzero factor.
There are only polynomially many classes in the endpoint bit length.
The remaining unit products again use (17). This is the same
normalization used by the individual inverse-entry oracle.

## 4. Precision and interface

The module exposes three operations:

| Operation | Returned representation |
|---|---|
| `inverse_constant_vector_entry` | Exact valuation `E` and unit of `(M^(-1)1)[b]` modulo `2^p` |
| `inverse_polynomial_vector_entry` | `2^E` times an integral residue modulo `2^p`, with error valuation at least `E+p` |
| `inverse_polynomial_bilinear_form` | `2^(-k)` times an integral residue modulo `2^p`, with error valuation at least `p-k` |

The second and third operations allow a zero or even residue. They do
not infer an exact valuation from a truncated cancellation. To recover
a bilinear form to absolute error valuation at least `q`, request
`p=q+k`. This remains a polynomial precision budget when `k` is
polynomial in the input bit length.

All polynomials in the interface use integer coefficients. Rational
coefficients can also be handled by clearing their explicit common
denominators and charging the resulting precision. The implementation
does not expand a polynomial given by a succinct arithmetic circuit;
its degree and coefficient list are explicit costs.

## 5. What remains for factoring

The existing pole identity uses

\[
 \mathrm{tr}\bigl(M^{-1}(M\circ M)\bigr)
 =\sum_{a,b}M(a,b)-\binom m2\frac{t}{z}.
 \tag{18}
\]

A bounded-degree separable polynomial matrix can now replace
`M o M` in the left side of (18) and its trace can be evaluated without
listing either axis. The unresolved step is a representation of the
actual Hadamard-square kernel, at the required absolute precision,
whose polynomial degrees, coefficient lengths, number of terms, and
construction cost meet an improved uniform budget.

The short displacement relation of a Cauchy matrix does not itself
give such a representation. Nor may one count an exponential-degree
interpolating polynomial as a bounded-degree input. Equations (2) and
(14) establish genuine compressed global sums for a restricted class
of weights and identify precisely where the factor-count trace differs.

## 6. Exact checks

[check_pole_cauchy_polynomial_solve.py](../tools/check_pole_cauchy_polynomial_solve.py)
forms independent exact inverses over `Q(w)` for 150 small matrices.
The first 100 use balanced padding; the remaining 50 use arbitrary
distinct odd nodes supplied as singleton progressions. It compares
every constant-vector output and four polynomial-vector degrees per
column, and checks four bilinear degree pairs per matrix.

A 402-bit control has domain size `2^200+13`, represented by three
progression pieces in total. The program computes degree-eight
solution queries and a degree-eight-by-degree-eight bilinear form.
It checks consistency at 12 and 16 bits of normalized precision
without enumerating a node. This verifies execution on the compressed
representation; the proof of polynomial cost is (17) and the fixed
number of short-jet operations.

The same complete bilinear form is also computed at 627 bits of
normalized precision with `k=615`. That gives absolute error valuation
at least 12, including all bits needed to compensate the factor
`2^(-k)`. The control therefore exercises an absolute precision budget,
as well as the smaller normalized checks.

The exact counts and residues are recorded in
[pole-cauchy-polynomial-solve.json](../evidence/pole-cauchy-polynomial-solve.json).
