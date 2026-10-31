# Rational conic points expose the same smooth-value search

The fixed-square relation targets are supplied as squares g^2, so a
rational point on the corresponding conic is easy to construct. This
note audits the proposed next step: parameterize other rational points
and cancel their denominators until only the prescribed split primes
remain.

For discriminant -4N with N odd, the cancellation has an exact simple
form. With notation below, its odd part is just

```math
\gcd(ar-bs,N).
```

A proper value already factors N by that linear gcd. If it does not,
the resulting allowed denominator is a represented value of the original
square-root class, or its product with the single unhelpful ambiguous
class U. Thus rational parameterization reduces the remaining task to
finding smooth values of one or two binary quadratic forms. It supplies
no automatic decrease of the denominator or of its largest prime factor.

The result is an exact algebraic description of this proposed search
route, with a checked obstruction to monotone numerical descent. It does
not rule out an effective nonmonotone sieve or another arithmetic method.
It also does not supply the unresolved smooth-square decision algorithm.

## 1. Preparing a rational point from the supplied square root

Let g be a given proper class of discriminant Delta=-4N. The input is its
reduced form, not a short word in the factor base. We can either obtain a
proper factor of N or choose an equivalent representative of g with odd
leading coefficient a coprime to N.

Indeed, if the reduced leading coefficient is odd, take its gcd with N.
A nontrivial gcd is proper because the leading coefficient is below N.
If the leading coefficient is even, the middle coefficient is even and
primitivity makes the final coefficient odd. Swap the outer coefficients
and negate the middle coefficient by a determinant-one substitution,
then take the gcd of the new leading coefficient with N. For a reduced
form with even leading coefficient A>=2, its final coefficient satisfies

```math
C\le N/A+A/4\le 4N/(3A)\le 2N/3<N,
```

so this gcd also either equals one or gives a proper factor.

Squaring the resulting ideal class by the polynomial composition routine
produces the primitive equivalent target

```math
H=(a^2,2b,C),\qquad a^2C=b^2+N,\qquad \gcd(a,2N)=1.       \tag{1}
```

No factorization of a is needed. The point (x,y,z)=(1,0,a) satisfies
H(x,y)=z^2 and gives the rational point (1/a,0) on H=1. Its denominator
is generally not smooth. A recorded determinant-one reduction matrix
transfers every point below to the canonical reduced representative of
the target without changing its denominator or coordinate content.

## 2. Parameterization and exact content

For coprime integers r,s, put

```math
w=ar-bs,
\qquad
x=Cs^2-r^2,
\qquad
y=2sw,
\qquad
z=ar^2-2brs+aCs^2.                                      \tag{2}
```

Then H(x,y)=z^2. These are the usual line-through-a-point formulas:
substituting z=ax+(r/s)y into H(x,y)=z^2 gives (2). The projective value
s=0 and the finite parameter with y=0 include the endpoint cases, so the
parameterization covers the rational conic.

The denominator before cancellation is the positive quadratic-form value

```math
z=F(r,s),\qquad F=(a,-2b,aC).                             \tag{3}
```

The form F has discriminant -4N and represents the inverse of the
original square-root class g. Inverting a form does not change its set
of integer represented values: negate one coordinate. Thus even before
examining cancellation, the denominator polynomial is a norm form of
the known square root.

Let c=gcd(x,y,z)>0. Because gcd(r,s)=1, equation (2) gives gcd(s,c)=1.
The identity

```math
4a^2x=-\Delta s^2-(2ar-2bs)(2ar+2bs)
```

then shows c divides |Delta|=4N. There is at most one factor of two in c.
To see this, if x and z are even then s is odd, and
C= b+1 mod 2 follows from (1). The equation x=0 mod 2 forces
r=C mod 2. Therefore w=ar-bs is odd, and y=2sw has two-adic valuation
exactly one.

Write d=gcd(w,N). The odd part of c divides d because c divides 2N and
c divides y=2sw, with gcd(s,c)=1. Conversely, d divides all three
coordinates by

```math
az=w^2+Ns^2,\qquad
a^2x=Ns^2-w(w+2bs),\qquad y=2sw.
```

Here a is invertible modulo d. Consequently the content is exactly

```math
c=\varepsilon\gcd(ar-bs,N),\qquad
\varepsilon=
 \begin{cases}
 2,&F(r,s)\text{ is even},\\
 1,&F(r,s)\text{ is odd}.
 \end{cases}                                             \tag{4}
```

The normalized integer point has coordinates (x/c,y/c) and denominator
z/c. Its two coordinates are coprime: a common divisor would also divide
z/c by the norm-square identity, contradicting removal of the full
projective content.

If 1<gcd(ar-bs,N)<N, the parameter already gives a proper factor by a
single linear gcd. Computing the larger quadratic coordinates and
cancelling them does not create an additional source of factors.

## 3. The remaining cancellation cases preserve the root norm class

Assume that the linear gcd does not factor N. Equation (4) leaves exactly
four possible contents: 1, 2, N, and 2N.

### Cancellation by all of N

If N divides the content, then ar-bs=0 mod N. Set
r_N=b a^-1 mod N and write r=Nt+r_Ns. Substituting this parameter lattice
into F and dividing by N gives the integral primitive form

```math
F_N=
\left(aN,
      2ar_N-2b,
      \frac{ar_N^2-2br_N+aC}{N}\right).                  \tag{5}
```

It has the same discriminant as F. Its middle coefficient is congruent
to that of F modulo 2a and to zero modulo 2N. Coprime composition
therefore identifies its class with

```math
[F_N]=[F]\,[(N,0,1)]=[F],                               \tag{6}
```

because (N,0,1) is equivalent to the principal form (1,0,N).
Primitivity can also be checked directly. A common prime divisor over a
prime dividing a is excluded by gcd(a,b)=1; at a prime dividing N, the
last coefficient of (5) is congruent to a^-1 modulo that prime.

Thus cancellation by the whole N is a determinant-N change of parameter
lattice followed by division of its quadratic value by N. It preserves
the proper root norm class. It does not provide a new small class or a
smoothness guarantee.

### Cancellation by two

This step may be applied either to F or to F_N. Write the current
primitive norm form as J=(A,B,C_0), with A odd and B even, and set
r_2=C_0 mod 2. Its even values lie on the parameter lattice
r=2t+r_2s. Dividing the value by two gives

```math
J_2=
\left(2A,
      2Ar_2+B,
      \frac{Ar_2^2+Br_2+C_0}{2}\right).                 \tag{7}
```

The middle coefficient is 2 mod 4. The discriminant identity shows that
the last coefficient is odd when N=1 mod 4 and even when N=3 mod 4.
In the first case J_2 is primitive, and coprime composition gives

```math
[J_2]=[J]\,U,\qquad U=(2,2,(N+1)/2).                    \tag{8}
```

This U is precisely the unique nonidentity ambiguous class that fails
to factor N in [early-stop extraction](early-stop-class-factoring.md).
In the second case every coefficient of J_2 is even, so the normalized
denominator is still even. It cannot be supported on the permitted
unramified split primes, because two divides Delta and is excluded from
that support.

It follows that, unless the parameter already exposes a proper gcd,
every admissible denominator is an integer represented by g or, when
N=1 mod 4, by gU. Conversely, this is an arithmetic search in the same
root-form classes. It is not a generic denominator-reduction operation.

## 4. Recovering and comparing square roots without factoring the denominator

There is a useful direct check on any returned representation. Suppose
a target Q represents d^2 primitively, with gcd(d,Delta)=1. Extend its
primitive representing vector to a determinant-one integer basis. The
transformed target is

```math
(d^2,B',C').
```

Then

```math
q=[(d,B',dC')]
```

is an explicitly computable square root of [Q]. Its discriminant is
Delta, and composing it with itself gives Q because the leading norm
and the compatible root lift are d^2 and B'. This construction uses
extended gcd and reduction; it does not factor d.

If g was the supplied square root of Q, the quotient qg^-1 is ambiguous.
The early-stop formulas either find a proper factor or identify it as
the identity or U. A smooth-square oracle can therefore sometimes factor
immediately by returning a representation with a different helpful root.
It may also return the original root, so this observation does not
replace the relation construction in a worst-case algorithm.

For the free parameterization (2), the preceding calculation is sharper:
a helpful root difference requires the proper linear gcd already
exhibited in (4). Thus merely generating rational conic points does not
supply new ambiguous directions beyond those gcd events.

One can track the orientation, rather than only the represented norm,
as follows. Put omega=sqrt(-N), I=(a,omega-b), and
beta=w+s omega in I. The compatible lift in (1) gives
I^2=(a^2,omega-b). The element encoding the target representation is

```math
\alpha=a^2x+(b+\omega)y=-\overline\beta^{\,2}.             \tag{9}
```

The integral ideal J=(beta)I^-1 has norm z. The ideal encoded by the
target point is its conjugate square, so its recovered root class is g
when no content is removed. If N is removed, beta/omega is integral and
still lies in I. Integrality follows from N dividing w; membership in I
follows because multiplication by omega is invertible modulo I, as
gcd(a,N)=1. Equation (9) divided by N is again a square of that conjugate
element, up to a unit. The root class is unchanged.

When two is removed and N=1 mod 4, the ideal
R_2=(2,omega-1) is the invertible ideal of class U and satisfies
R_2^2=(2). The odd values of w and s put beta in I R_2. Dividing the
conjugate square ideal by (2) therefore multiplies its root class by U.
The same reasoning applies after division by N, because N is odd.
These arguments use only principal division at N and the proper
norm-two ideal when it exists; they do not assume that N is squarefree.
For an admissible denominator all remaining prime ideals are away from
Delta, so their square-root ideal is unambiguous. This proves the stated
root-class conclusion for the actual point, as well as for its norm.

## 5. Monotone numerical descent can stop before a smooth denominator

A reduced positive form (A,B,C) has every nonzero integer represented
value at least A. This follows from

```math
Ar^2+Brs+Cs^2
 \ge A(r^2-|rs|+s^2)\ge A.
```

If the initial denominator is already the smallest leading coefficient
among the permitted root classes g and gU, a strict decrease cannot find
a different admissible denominator. Smoothness can require a larger
value even when the small prime support generates the full class group.

For a particularly small control, take Delta=-188, N=47, and
g=(7,6,8). The prime form at 3 generates all five classes. There is no
U alternative because N=3 mod 4. The minimum root-form value is seven,
while g(1,-1)=9 is the first denominator supported on {3}. The target
square is represented by H=(49,20,3), and parameter (r,s)=(2,1) gives

```math
(x,y,z)=(-1,8,9),\qquad H(-1,8)=81.
```

After reducing H to (3,-2,16), the representing vector is (5,1).
A successful supported-denominator search increases seven to nine.

There is also a composite control with both prime factors outside the
chosen support. Take

```math
N=299=13\cdot23,\quad \Delta=-1196,\quad
g=(7,-6,44),\quad P=\{3,5\}.
```

The prime forms at 3 and 5 generate all 24 classes. The known-square
target is H=(49,-34,12), reduced to (12,10,27). Parameter (-3,1) gives

```math
(x,y,z)=(3,-8,45),\qquad c=1,
```

so 45=3^2*5 is an allowed denominator in the original root class.
It is its first allowed denominator. For a nonzero second coordinate,
completion of the square bounds g(r,s) below by (299/7)s^2, and the
first values with |s|=1 are 44 and 45. With s=0, every value has a
factor seven. There is no smaller supported value.

By contrast, parameter (-3,2) has w=13 and produces

```math
(x,y,z)=(39,52,195),\qquad c=13,
\qquad (x/c,y/c,z/c)=(3,4,15).
```

Its recovered square root is (15,-2,20), whose quotient with g is
(13,0,23). This factors N, exactly as the earlier linear gcd
gcd(13,299) already did. The two controls distinguish a nonmonotone
smooth-value search from cancellation that immediately solves factoring.

These examples obstruct a claim of guaranteed strictly decreasing
denominators. They do not establish a lower bound for arbitrary conic
algorithms, nor do they exclude progress under a different potential
function. A certified strategy for finding the required smooth values
remains open.

## 6. Exact implementation and replay

The implementation is
[conic_denominator_descent.py](../tools/conic_denominator_descent.py).
The independent finite controls are in
[check_conic_denominator_descent.py](../tools/check_conic_denominator_descent.py).
Run

```console
python3 tools/check_conic_denominator_descent.py > evidence/conic-denominator-descent.json
```

The [evidence](../evidence/conic-denominator-descent.json) checks all 127
odd N from 3 through 255. Among 991 reduced root forms, 119 immediately
expose a factor while preparing an odd coprime leading coefficient.
For the remaining 872 forms it checks 35,752 primitive parameter points
from a fixed 41-point stencil. The exact content formula finds 4,092
proper linear-gcd factors. All four remaining content cases occur.

The checker verifies 26,558 unramified root recoveries and 5,102 even
imprimitive residual forms in the N=3 mod 4 parity case. It also verifies
the full class-group generation and exact coordinates in both descent
controls above, including the helpful root difference for N=299.
The computations validate the arithmetic identities and implementation;
they are not a smooth-value density theorem or a search-time bound.
