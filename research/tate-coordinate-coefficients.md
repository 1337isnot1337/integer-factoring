# Tate coordinates reduce the coefficient engine to a quadratic composition

For every odd prime phase order q, this note constructs an algebraic series
X_q(v) over an explicit degree-(q-1)/2 finite algebra. At every odd positive
integer N, the signed divisor residues are the coefficient of t^N in

```math
X_q(\Delta(t)),\qquad \Delta(t)=\sum_{r\ge0}t^{(2r+1)^2}
\quad\text{over characteristic two}.                    \tag{1}
```

The q-division polynomial, coefficient algebra, branch and rational-diagonal
representation are uniformly constructible in polynomial time in q. A single
coefficient [v^n]X_q(v) is computable in poly(q,log(n+1)) bit operations.
The implementation performs O(log(n+1)) bounded-size diagonal steps.

The remaining target is coefficient extraction **after** the sparse quadratic
substitution in (1). The fast internal coefficient algorithm does not supply
a fast coefficient of that composition. No improved factoring exponent is
proved; the strongest complete exponent remains 1/2.

The Tate curve, division polynomials, Hensel lifting and the
Furstenberg/Cartier compiler are standard constructions. The contribution
here is their explicit composition into the existing signed divisor interface,
with a factorization-free coefficient algebra, uniform cost accounting and
exact controls. This work uses OpenAI Codex assistance and builds on the
[coefficient reduction](weight-one-divisor-coefficients.md),
[signed sparse recovery](sparse-divisor-residues.md) and
[theta-section engine](lambert-theta-sections.md).

## 1. The odd-index divisor coefficient is a Tate x-coordinate coefficient

Retain the additive phase algebra R_q=F_2[Z]/Phi_q(Z) of the theta note.
The integral Tate curve has the equation

```math
y^2+xy=x^3+a_4(t)x+a_6(t),
```

with integral coefficient series, discriminant
t product_(m>=1)(1-t^m)^24, and x-coordinate

```math
X(Z,t)=\frac{Z}{(1-Z)^2}
 +\sum_{m\ge1}\left(
 \frac{Zt^m}{(1-Zt^m)^2}
 +\frac{Z^{-1}t^m}{(1-Z^{-1}t^m)^2}
 -\frac{2t^m}{(1-t^m)^2}\right).                         \tag{2}
```

These standard formulas and the torsion embedding are stated in
[Meier--Ozornova, Appendix A.3](https://arxiv.org/html/1812.04425#A3),
Theorem A.16 and the following construction of mu_q -> Tate.
They are first integral identities, and are then reduced modulo two. In
particular, we do not divide the displayed integral formula for a_6 by twelve
inside a characteristic-two algebra.

Since z/(1-z)^2=sum_(d>=1) d z^d, reduction of (2) gives, for N>=1,

```math
[t^N]X(Z,t)=\sum_{\substack{d\mid N\\d\text{ odd}}}
             (Z^d+Z^{-d}).                              \tag{3}
```

For odd N this is exactly S_N(Z), the signed divisor polynomial. For even
N, (3) is the signed coefficient for its odd part; the factoring endpoint
already returns factor two on even inputs. The constant in (2) is

```math
x_0=Z/(1+Z)^2.
```

Evaluating at a nontrivial q-th root gives a point of order q. The integral
torsion construction is defined after inverting q; q is odd, so its
reduction modulo two is valid. All constant denominators 1-Z are units in
R_q. Positive t-degree denominators in (2) have unit constant term.

## 2. A normalization that keeps x unchanged

In characteristic two, translating the y-coordinate by a_4(t) gives

```math
y^2+xy=x^3+v,\qquad v=a_6(t)+a_4(t)^2=\Delta(t).          \tag{4}
```

The discriminant formula for a_1=1,a_2=a_3=0 is Delta=a_6+a_4^2 modulo two.
The x-coordinate is unchanged by this translation. No square root or field
extension is used.

To identify Delta, differentiate Jacobi's triple product in its independent
phase variable and then set that phase to one. Pairing its two equal
triangular exponents gives

```math
\prod_{m\ge1}(1-t^m)^3
 =\sum_{r\ge0}(-1)^r(2r+1)t^{r(r+1)/2}.
```

Modulo two, Frobenius therefore yields

```math
\Delta(t)=t\prod_{m\ge1}(1-t^m)^{24}
 =t\sum_{r\ge0}t^{8r(r+1)/2}
 =\sum_{r\ge0}t^{(2r+1)^2}.                             \tag{5}
```

The product identity is the derivative specialization of
[Jacobi's triple product](https://dlmf.nist.gov/17.8.E1); the argument above
also fixes its normalization.

## 3. Construct the division polynomial without factoring a large integer

For E_v:y^2+xy=x^3+v, use the standard division polynomials, with

```math
\psi_0=0,\quad\psi_1=1,\quad\psi_2=x,\quad
\psi_3=x^4+x^3+v,\quad\psi_4=x^6+vx^2.
```

The recurrences specialize to

```math
\psi_{2m+1}=\psi_{m+2}\psi_m^3+
             \psi_{m-1}\psi_{m+1}^3,
\qquad
x\psi_{2m}=\psi_m(\psi_{m+2}\psi_{m-1}^2+
                         \psi_{m-2}\psi_{m+1}^2).        \tag{6}
```

Division by x is exact. If m is even, psi_m already contains x; if m is
odd, the two even-index factors in the bracket contain x, and are squared.
The code multiplies the entire numerator before exact division.

The usual division-polynomial interpretation, also documented for
[general Weierstrass equations in Sage](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_generic.html#sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic.division_polynomial),
implies psi_q(X(Z,t),Delta(t))=0. The equation is a polynomial identity over
the integral torsion construction and hence remains valid after reduction.

For odd q, deg_x psi_q=(q^2-1)/2. Standard Weierstrass weights, or (6) before
setting a_1=1, give 2i+6j<=q^2-1 for its monomials x^i v^j. Thus its
degrees, coefficient count and construction cost are polynomial in q.
The implementation stores binary bivariate coefficients, not a supplied
factorization of the original coefficient index N.

## 4. The intrinsic coefficient algebra and a simple branch

The recurrences imply the nodal rational identity

```math
\psi_n\left(\frac{z}{(1+z)^2},0\right)
 =\frac{z^{n(n-1)/2}(1+z^n)}{(1+z)^{n^2}}.              \tag{7}
```

It holds for n=1,2,3,4 directly. Substitution into (6) proves the induction.
For example the even step reduces to

```math
(1+z^{m+2})(1+z^{m-1})^2+
(1+z^{m-2})(1+z^{m+1})^2
 =z^{m-2}(1+z)^4(1+z^m).
```

The odd step follows by expanding the two cubes; its bracket equals
z^(m-1)(1+z)^3(1+z^(2m+1)). No finite-field point count is used in this
all-index argument.

For the odd prime q put s=q(q-1)/2. Equation (7) shows

```math
\psi_q(x,0)=x^s C_q(x),\qquad
\deg C_q=(q-1)/2,\quad C_q(0)=1.                         \tag{8}
```

The nonsingular roots are x=z/(1+z)^2 for z^q=1,z!=1. Opposite roots z and
z^-1 give the same x and are the only duplication. Their derivative in z
is nonzero: in characteristic two, dx/dz=(1+z)^-2. Equation (7) then shows
that all these roots of C_q are simple. Consequently C_q is squarefree.

Construct the finite etale algebra

```math
A_q=\mathbb F_2[c]/C_q(c).                               \tag{9}
```

The known x^s factor in (8) is removed by an exact coefficient shift;
constructing (9) does not factor a polynomial. In A_q, c is a unit and

```math
\partial_x\psi_q(c,0)=c^s C_q'(c)
```

is a unit. Its inverse is found by the polynomial Euclidean algorithm.
Hensel's recursion therefore produces a unique series

```math
X_q(v)\in A_q[[v]],\qquad X_q(0)=c,\quad\psi_q(X_q(v),v)=0.
                                                               \tag{10}
```

The embedding A_q -> R_q sends c to x_0=Z/(1+Z)^2. It is injective because
the phase roots cover all the distinct roots of C_q. Both the embedding and
the phase-to-signed-vector conversion have polynomial cost in q.

After v=Delta(t), (10) and the Tate x-coordinate have the same constant and
the same simple equation. Hensel uniqueness proves their equality, giving
(1) and (3) at every index. This is the uniform identity; the prefix checks
below are controls of its implementation.

## 5. A polynomial-size coefficient compiler in v

Set W_q(v)=X_q(v)-c and F(v,w)=psi_q(c+w,v). Then F(0,0)=0 and
F_w(0,0) is a unit. Furstenberg's formula represents W_q as a rational
diagonal with

```math
W_q(v)=\operatorname{Diag}\frac{P(s,w)}{Q(s,w)},\qquad
P=wF_w(sw,w),\quad Q=F(sw,w)/w.                          \tag{11}
```

Here Q(0,0)=F_w(0,0) is a unit. Formula (11) and its degree bounds are
given in [Bostan--Christol--Dumas, Section 3.1](https://arxiv.org/html/1602.00545#S3.SS1).
Although A_q can be a product of fields, the identity applies componentwise;
the implementation never needs to construct that product decomposition.

Let sigma square coefficients in A_q. For the current diagonal target n,
with b=n mod 2, the characteristic-two update is

```math
(P,Q,n)\longmapsto
(\Lambda_{b,b}(PQ),\ \sigma(Q),\ \lfloor n/2\rfloor).     \tag{12}
```

Indeed Q(s,w)^2=sigma(Q)(s^2,w^2), so multiplication by Q and simultaneous
parity selection in both variables proves (12). At target zero the result
is P(0,0)/Q(0,0). The constant coefficient of X_q is added separately.

Unlike the original nome sections, all polynomials here have bounded degree.
For d=(q^2-1)/2, the weight bound in Section 3 gives

```math
\deg_s(P),\deg_s(Q)\le\lfloor d/3\rfloor,\qquad
\deg_w(P),\deg_w(Q)\le d.
```

The same rectangle is preserved by (12). It contains O(q^4) monomials,
each with O(q) coefficient bits. A naive step uses polynomially many
coefficient-ring operations, and there are exactly bitlength(n) steps.
Constructing the equation, inverses and diagonal data is also polynomial
in q. Thus [v^n]X_q has complete bit cost poly(q,log(n+1)), without an
initial expansion through n. Integer-index handling in the literal Python
loop is polynomial in its bit length; no unit-cost large-integer assumption
is needed.

This is a new-variable coefficient algorithm. Replacing v by Delta(t)
does not preserve the meaning of a coefficient index. Applying the
algebraic-series algorithm directly to the original t-series would require
a new hypothesis, and the
[common-linear-space obstruction](lambert-theta-sections.md#6-exponential-dimension-for-a-common-linear-section-space)
already excludes the proposed small common linear closure at q=3.

## 6. Independent high-index control at q=3

Here A_3=F_2, c=1, and psi_3=x^4+x^3+v. Writing x=1+w gives

```math
v=w(1+w)^3.
```

Lift this equation to the unique integral formal branch. Lagrange inversion,
in the normalization of [DLMF 1.10](https://dlmf.nist.gov/1.10.vii), gives
for n>=1

```math
[v^n]w(v)=\frac{(-1)^{n-1}}{n}\binom{4n-2}{n-1}.
```

Its parity is one exactly when

```math
s_2(n-1)+s_2(3n-1)-s_2(4n-2)=v_2(n),                    \tag{13}
```

using v_2(m!)=m-s_2(m). This is an independent O(poly(log(n+1))) control
of the compiled diagonal values, without constructing a huge binomial.
It concerns the internal algebraic index n. The polynomial-time square test
for the original q=3 coefficient is a different formula.

## 7. A cheap nonlinear trace after the Delta substitution

There is one uniform family for which the original composition can be
evaluated cheaply. Let Tr:A_q -> F_2 be the algebra trace, and put
T_q(t)=Tr(X_q(Delta(t))). Then

```math
T_q(t)^2+T_q(t)=\Delta(t)+\Delta(t^q).                    \tag{14}
```

To prove this, take the trace in (3). Its conjugates correspond to all the
inverse pairs of nontrivial q-th roots. For each odd divisor d, the trace
of Z^d+Z^-d is the sum of Z^(ad) over 1<=a<q. It is one if q does not divide
d and zero otherwise. Thus, for N>=1,

```math
[t^N]T_q=[N/(2^{v_2(N)}q^{v_q(N)})\text{ is a square}].   \tag{15}
```

The positive-degree part is therefore

```math
\sum_{a,b\ge0}\sum_{\substack{m>0\text{ odd}\\q\nmid m}}
 t^{2^a q^b m^2}
 =\sum_{a\ge0}\bigl(\Delta(t^{2^a})+\Delta(t^{q2^a})\bigr).
```

The equality follows by removing the odd squares divisible by q and
telescoping in b. Squaring shifts the latter sum from a>=0 to a>=1,
proving (14); the constant trace is a bit and cancels. Equations (14)--(15)
give a polynomial-time original-index evaluator for this trace, for every q.
For odd original indices the removal of the power of two is unnecessary.

The trace can be constructed by taking the trace of the multiplication
matrix in the basis 1,c,...,c^(h-1), h=(q-1)/2. This needs no factorization
of C_q. Because A_q is etale, its trace pairing is nondegenerate. Search
this basis for eta with Tr(eta)=1; this works even when h is even and
Tr(1)=0. The decomposition

```math
X_q(\Delta(t))=\eta T_q(t)+
 \bigl(X_q(\Delta(t))-\eta T_q(t)\bigr)                   \tag{16}
```

has an explicitly cheap first term and a second term with trace zero.
The latter has h-1 coefficient coordinates. It is the entire unresolved
part when q>3. At q=3, h=1 and (16) recovers the full square-test engine.
This trace result does not supply the other coordinates or a new factoring
bound.

## 8. Evidence and the remaining construction

Run `make verify-tate-coordinate`, or

```sh
python3 tools/check_tate_coordinate_compiler.py
```

The [compiler](../tools/tate_coordinate_compiler.py),
[checker](../tools/check_tate_coordinate_compiler.py) and
[source-pinned receipt](../evidence/tate-coordinate-compiler.json) retain:

- 65 exact nodal identities (7), with rational denominators cleared;
- intrinsic nodal algebras, simple branches and embeddings for nine prime
  levels through 31;
- 221 internal coefficients checked against Newton/Hensel prefixes;
- additional dense bivariate division controls of the diagonal formula;
- 221 finite Delta-composition coefficients checked against literal Tate
  Lambert coefficients, including 108 odd original indices;
- 128 small exact binomial controls and 30 internal-index controls at q=3
  through 1,024 bits, checked independently by (13).
- the trace identity (14), trace-one projections including even algebra
  degrees, and 54 large original-index trace controls.

These are ordinary proofs and finite source checks, not new Lean verification
or independent acceptance of an improved factoring bound.

The explicit remaining engine is: extract [t^N]X_q(Delta(t)) at odd
binary-described N and prime levels q in the earlier polynomial level
budget. The inner algebraic coefficients are cheap; enumerating all of
them through N, or expanding the full composition through N, is still
exponential in bitlength(N). A nonlinear quadratic-composition method with
its full cost charged is the next target.
