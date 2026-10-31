# Known class roots give an integral conic parametrization

The prescribed-square targets in the arithmetic program already have known
class square roots. Consequently rational solvability of their conics is
available in polynomial bit time, without factoring the discriminant. The
remaining problem is finding a **smooth denominator**, with a precise
discriminant bound on cancellation. This note proves the bridge in both
directions and implements a direct gcd certificate for ramified splitting.
It does not supply the missing search algorithm or improve the score `1/2`.

Parametrizing a conic from a known rational point is classical; see the
primary [PARI documentation for `qfparam`](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Vectors__matrices__linear_algebra_and_sets.html)
and [Magma's conic documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/1538).
The contribution here is an explicit integral matrix from ideal squaring,
its determinant exactly equal to the absolute order discriminant, the
resulting cancellation/divisibility contract, and checked conversion code.
No novelty claim is made for conic parametrization or ambiguous-form
factor extraction.

## 1. An integral matrix of determinant |Delta|

Let F=(a,b,c) be a primitive positive definite form of discriminant
Delta<0. The form need not be reduced. Let Q be the canonical reduced
representative of [F]^2. There is an explicitly computable integer
three-by-three matrix C such that

```math
C\begin{pmatrix}u^2\\uv\\v^2\end{pmatrix}
 =\begin{pmatrix}P_0(u,v)\\P_1(u,v)\\F(u,v)\end{pmatrix},
\qquad \det C=-\Delta,                                  \tag{1}
```

and, for three independent variables r_0,r_1,r_2,

```math
Q\bigl((Cr)_0,(Cr)_1\bigr)-(Cr)_2^2
 =-\Delta(r_1^2-r_0r_2).                                \tag{2}
```

Thus Q(P_0,P_1)=F^2. All coefficients and the construction cost are
polynomial in the bit length of F. There is no search over represented
integers and no factorization of Delta.

Here is the integer construction. Put s=Delta mod 2,
k=(s^2-Delta)/4, and t=(b+s)/2. In O=Z[omega], with
omega^2=s omega-k, use the invertible ideal I=(a,omega-t). An element

```math
\alpha=au+(t-\omega)v
```

has norm aF(u,v). The polynomial composition routine in
[factored-norm search](factored-norm-search.md) computes d=gcd(a,b)
and the primitive ideal J=I^2/d. Write its norm form before reduction
as Q_0=(A,B,C_0), where A=a^2/d^2, and put h=(-B-s)/2. The ideal J has
basis (A,h+omega). The horizontal and vertical coefficients of alpha^2
in the basis (1,omega) are

```math
H=(au+tv)^2-kv^2,\qquad V=-2auv-bv^2.
```

Both are divisible by d. Set

```math
Y=-V/d,\qquad X=\frac{H/d+hY}{A}.                       \tag{3}
```

These are integral quadratic polynomials because alpha^2/d belongs to
J. The sign in Y makes its norm form Q_0(X,Y), rather than Q_0(X,-Y).
The norm identity gives Q_0(X,Y)=F(u,v)^2. Applying the inverse integral
reduction matrix to (X,Y) gives P_0,P_1 for the reduced Q.

Before reduction, the coefficient rows of X,Y,F are

```math
\left(d,\frac{2d(t+h)}a,\frac{d(t^2-k+bh)}{a^2}\right),
\quad (0,2a/d,b/d),\quad(a,b,c).
```

Their determinant is

```math
2ac-b^2+2b(t+h)-2(t^2-k+bh)
 =4ac-b^2=-\Delta.
```

The reduction matrix has determinant one, so it preserves this
determinant. Expanding Q_0(X,Y)-F^2 on independent monomials yields
(2), also preserved by reduction. These arguments cover conductor
primes and even discriminants as well as maximal orders.

## 2. Every rational conic point, and its exact cancellation

For coprime integers u,v, let

```math
g(u,v)=\gcd(P_0(u,v),P_1(u,v),F(u,v))>0.
```

The adjugate identity and (1) show that g divides each of
Delta*u^2, Delta*uv, and Delta*v^2. Their gcd is |Delta| because
gcd(u,v)=1. Therefore

```math
g(u,v)\mid|\Delta|.                                    \tag{4}
```

In particular all common cancellation is supported on ramified primes;
it is not an uncharged arbitrary smooth factor.

Conversely, suppose p=(x,y,z) is a primitive integer point with z>0
and Q(x,y)=z^2. Compute the integer vector

```math
r=\operatorname{adj}(C)p.
```

Equation (2) gives r_1^2=r_0r_2. If r_0 is nonzero, reduce the pair
(r_0,r_1) by its gcd to get (u,v); if r_0=0 take (u,v)=(0,1).
The vector r is projectively (u^2,uv,v^2). Since F is positive definite
and z>0, the projective scale has positive sign. It follows that

```math
p=\frac{(P_0(u,v),P_1(u,v),F(u,v))}{g(u,v)},\qquad
g(u,v)=\frac{|\Delta|}{\gcd(r_0,r_1,r_2)}.               \tag{5}
```

This proves surjectivity and gives an integer inverse. Nonprimitive input
points are first divided by gcd(x,y,z). Construction, inversion, and
content extraction all have polynomial bit cost in the supplied data.

## 3. The fixed-square problem is a denominator-divisor problem

Fix M>0, supplied in factored form if desired. Equations (1)--(5) prove
the exact equivalence

```math
Q(x,y)=M^2\text{ for some integers }x,y
\quad\Longleftrightarrow\quad
\exists\,\gcd(u,v)=1:
       \frac{F(u,v)}{g(u,v)}\mid M.                     \tag{6}
```

For the forward direction normalize (x,y,M) to a primitive point. Its
last coordinate is a positive divisor of M. Invert that point by (5).
For the reverse direction scale the primitive point from (1) by
M/(F/g). This produces the required coordinates exactly.

If M is Y-smooth and coprime to Delta, (6) requires a form value

```math
F(u,v)=g\,d,\qquad g\mid|\Delta|,\quad d\mid M,
\quad\gcd(d,\Delta)=1.                                 \tag{7}
```

The exact content condition g=g(u,v) remains part of the contract;
arbitrary divisors g of Delta cannot simply be substituted. In
particular it is unnecessary to factor Delta to compute g from supplied
parameters or a supplied conic point.

For the [universal-square targets](universal-smooth-square.md),
Q=red(g_i^(2K)) has the explicit root F=red(g_i^K). A rational point and
the entire parametrization are therefore immediately accessible. For
every query in the [decision reduction](smooth-square-decision-reduction.md)
one can likewise maintain a known root by updating the generator powers
used to shift that query. Neither a rational point nor the promise that
some denominator divides M finds the needed parameters. General conic
solvability should not be counted as the missing oracle.

The [fixed-strip conditional theorem](fixed-strip-smooth-square.md)
can make M polynomial in bit length under its explicit analytic inputs.
It does not choose the parameters in (6). Conversely, the large
unconditional universal M retains its charged encoding length in this
model. Rewriting its conic does not replace that length by the original
discriminant bit length.

## 4. Direct splitting certificates

For Delta=-dN, a supplied point gives the integer content g in (5).
If 1<gcd(g,N)<N, that gcd is already a certified proper factor.
There is no need to factor g first. If the gcd is trivial or all of N,
the point supplies no proper split by this test.

For Delta=-4N, N odd, and an ambiguous root F, Q=(1,0,N).
Both unit points (1,0,1) and (-1,0,1) are known. Applying (5) to them
recovers the usual ambiguous-form factors. The three reduced boundary
cases make this explicit:

| Root F | The two contents, in either order | Product |
| --- | --- | --- |
| (A,0,C) | A, C | N |
| (2a,2a,C) | 2a, 2(2C-a) | 4N |
| (A,2v,A) | 2(A-v), 2(A+v) | 4N |

The first case includes the identity. The exceptional nonsplitting form
(2,2,(N+1)/2), when primitive, gives contents 2 and 2N. These are exactly
the exceptional classes retained by
[early-stop extraction](early-stop-class-factoring.md). Thus this is a
uniform matrix-and-gcd implementation of the same splitting map, not a
new source of ambiguous classes. General supplied-square witnesses can
also yield useful content; uniform discovery of such witnesses remains
open.

The concurrent [conic denominator audit](conic-denominator-descent.md),
integrated from `014ba9f`, gives a sharper specialization for Delta=-4N.
After preparing an odd coprime leading coefficient, its content is
exactly a parity factor times one linear gcd with N. Unless that gcd
already factors N, the remaining denominator belongs to the root norm
class or its product with U. It also recovers a square-root class from
a supplied square representation without factoring its denominator and
gives explicit failures of monotone denominator descent. That audit and
the general-discriminant matrix and adjugate proof here are complementary
ordinary results; neither supplies smooth-denominator search.

## 5. Reproduction and evidence

Run `make verify-class-square-conic`, or

```sh
python3 tools/check_class_square_conic.py > evidence/class-square-conic.json
```

The [implementation](../tools/class_square_conic.py) uses the polynomial
ideal composer. The [checker](../tools/check_class_square_conic.py) compares
its squared classes with the independent finite composer, verifies (1)
and the full three-variable identity (2), tests primitive parameters and
inverse maps on independently enumerated conic points, and compares unit
point factors with the independent ambiguous boundary formulas. It also
replays the two-part tops in the existing 60/80-bit certificates, includes
large unreduced roots, and records a known-root point whose denominator
is not 1000-smooth. Exact counts are in the
[receipt](../evidence/class-square-conic.json).

The ordinary proof and finite checks have not received independent
mathematical review. This work used substantial OpenAI Codex assistance
and adds no Lean code.
