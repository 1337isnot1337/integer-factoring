# Why the small moment quotient cannot directly invert the pole

## Result and scope

The [polynomial-circuit moment construction](pole-cauchy-circuit-moments.md)
represents every integral polynomial contraction in a quotient whose
dimension is polynomial in the working precision. Directly inverting
`z-XY` in that quotient does not extend the construction to the
factor-count kernel. On balanced power-of-two intervals, the quotient
has an explicit artificial evaluation point on `z-XY=0` at every
precision up to half the matrix dimension.

More precisely, let both domains consist of `m=2^s` consecutive positive
odd integers, with `s>=2`. Let `N` be positive and odd, set
`z=N+2^k w`, where `w^2+w+1=0`, and assume `k>=s+1`. At every precision

```math
2\le P\le m/2,
```

the moment quotient admits the evaluation `X -> z, Y -> 1`. Hence no
polynomial class `R` in that quotient can satisfy

```math
(z-XY)R=2^r u,\qquad 0\le r<P,
```

where `u` is a scalar unit. The same statement holds with any positive
power of `z-XY` on the left.

This includes the separated factoring regime, where every actual node
denominator has valuation at most `k`. Taking `k<P<=m/2` therefore
separates two facts: the physical denominators need at most `k` powers
of two, whereas inversion in this particular polynomial quotient needs
at least `P` coefficient-denominator bits. Raising a polynomial working
precision does not close that gap when `m` is exponential in its binary
description.

This is an obstruction to direct quotient inversion. Common-factor
cancellation changes the quotient, and a scalar-specific extension
after such cancellation remains outside this theorem. Section 5 gives
an exact example in which cancellation removes the artificial point.
No factoring bound is improved or ruled out in general.

## 1. A complete-residue coefficient bound

Work in the unramified quadratic valuation ring
`O=Z_2[w]/(w^2+w+1)`. Suppose `u_1,...,u_m` reduce to a complete residue
system modulo `2^s`, with `m=2^s`. For every integer `ell>=1`,

```math
v_2\!\left(\sum_i u_i^\ell\right)\ge s-1.                 \tag{1}
```

For the representatives `0,...,2^s-1`, this follows by induction on `s`.
Pair `i` with `i+2^(s-1)`. Their power sum is congruent to `2*i^ell`
modulo `2^(s-1)`, and the induction hypothesis applies to the half-size
sum. Replacing the representatives by the given `u_i` changes the sum
by a multiple of `2^s`.

Let `e_j` be the elementary symmetric functions of `2u_1,...,2u_m`,
with `e_0=1`. Then

```math
v_2(e_j)\ge j+s-1-v_2(j),\qquad 1\le j\le m.              \tag{2}
```

Indeed the power sums of these even elements have valuation at least
`ell+s-1`. Newton's identity is

```math
j e_j=\sum_{\ell=1}^j(-1)^{\ell-1}e_{j-\ell}p_\ell.
```

The elementary bound `v_2(e_(j-ell))>=j-ell` makes every summand
divisible by `2^(j+s-1)`. Dividing by the integer `j` proves (2).

For a domain of `m` consecutive odd integers, the numbers `(b-1)/2`
form such a complete residue system. The reciprocal roots `z/a` have
the same property: because `k>=s+1`, they reduce modulo `2^(s+1)` to
`N/a`, and inversion followed by multiplication by odd `N` permutes all
odd residue classes. Thus (2) applies to both root-correction lists
`b-1` and `z/a-1`.

## 2. Both packets vanish at one

Set `J=P-1`; the assumed range of `P` ensures `J<m`. The moment packets
are the degree-`J` polynomials

```math
F(Y)=\sum_{j=0}^J(-1)^j e_j(b-1)(Y-1)^{J-j},
\qquad
G(Y)=\sum_{j=0}^J(-1)^j e_j(z/a-1)(Y-1)^{J-j}
                                                               \tag{3}
```

over `O/2^P`. These are exactly the monic packets used by the existing
unit-jet constructor; expanding the original product around `Y-1`
gives its defining packet congruence.

Since `1<=J<2^(s-1)`, one has `v_2(J)<=s-2`. Equation (2) therefore
gives

```math
v_2(e_J)\ge J+s-1-v_2(J)\ge J+1=P.
```

Evaluating (3) at one leaves only that last coefficient. Consequently,

```math
F(1)=G(1)=0\quad\text{in }O/2^P.                          \tag{4}
```

There is often a larger common power. A uniform lower bound is

```math
(Y-1)^r\mid F,G,\qquad
r=\min\{P-1,\max\{0,s-1-\lfloor\log_2(P-1)\rfloor\}\}. \tag{5}
```

The cap ensures `r<=J`. For every index `j>=P-r`, (2) makes its
coefficient zero modulo `2^P`.
The remaining terms in (3) have exponent at least `r`. The one-point
obstruction needs only (4).

## 3. The evaluation certificate rules out scaled inversion

The moment quotient is

```math
\mathcal A_P=(O/2^P)[X,Y]/(Q_A(X),Q_B(Y)),
\quad Q_B(Y)=F(Y),\quad
Q_A(X)=\frac{X^J G(z/X)}{G(0)}.                            \tag{6}
```

Both quotient polynomials are monic and `G(0)` is a unit. Equation (4)
gives

```math
Q_B(1)=0,\qquad Q_A(z)=z^J G(1)/G(0)=0.
```

Thus substitution `X=z,Y=1` defines a ring homomorphism from
`A_P` to `O/2^P`. It sends `z-XY` to zero but sends `2^r u` to a
nonzero element for every scalar unit `u` and every `r<P`. This proves
the claimed impossibility, without a degree restriction on `R` or a
restriction on how its polynomial circuit is written.

Arbitrary integral lifts of the quotient coefficients do not supply
a low-denominator inverse either. If an inverse exists over the
quadratic fraction field for such a lift, and its reduced polynomial
representative has coefficient denominators cleared by `2^r`, monic
polynomial reduction gives the corresponding integral identity modulo
`2^P`. The evaluation certificate then forces `r>=P`. This statement
also applies to positive powers of the denominator.

For `m=2^s` and input lengths `O(s)`, every polynomial precision lies
below `m/2` for sufficiently large `s`. The first precision not covered
by this obstruction is `P=2^(s-1)+1`. This is a precision lower bound
for the specified inversion strategy, not a general lower bound for
factor-count extraction.

## 4. A separated control with one genuine factor pair

Take `A=B={1,3,...,127}`, `N=10201=101^2`, and `k=14`. Here `m=64`
and every integer difference satisfies `|N-ab|<2^14`. There is exactly
one factor pair in the rectangle: `(101,101)`.

At precision `P=16`, the degree-15 packet quotient has the evaluation
certificate above. It cannot solve `(z-XY)R=2^14`, even though the
physical reciprocal at every node has denominator valuation at most
14. Using the exact trace identity in
[pole-cauchy-content.md](pole-cauchy-content.md), a literal node calculation gives

```math
\operatorname{tr}(M^{-1}(M\circ M))\equiv1\pmod4,
\qquad M_{a,b}=2^{14}w/(z-ab).                            \tag{7}
```

There is also a misleading residual check to avoid. In this quotient
`t^2=2^28 w^2` is already zero modulo `2^16`, so `R=0` satisfies
`(z-XY)^2 R=t^2`. Nevertheless its contraction is zero, while the true
unscaled contraction is one modulo four. Multiplying a rational
identity by a nonunit and checking its residual does not identify the
required rational value.

The same failure occurs for a proper-factor rectangle with different
domains. Take `A={3,5,...,129}`, `B={131,133,...,257}`,
`N=13837=101*137`, `k=15`, and `P=17`. Its unique factor pair is
`(101,137)`, and its trace is again one modulo four. The artificial
evaluation still rules out the physically sufficient scaled inverse.

## 5. Cancellation can remove the artificial point

If `F=C F_0` and `G=C G_0`, where the common monic factor `C` has unit
constant term, then both moment-generating ratios are unchanged:

```math
F/G=F_0/G_0,\qquad G^*/F^*=G_0^*/F_0^*.
```

The polynomial functional therefore factors through the smaller
quotient constructed from `F_0,G_0`. This applies to the common powers
of `Y-1` in (5), which can be divided out by monic synthetic division.

The effect can be substantial. Over the quadratic field, consider

```math
F=(Y-1)(Y-3),\qquad G=(Y-1)(Y-5).
```

The original quotient contains the point `(X,Y)=(z,1)`. After canceling
`Y-1`, the quotient has `X=z/5,Y=3`, so

```math
z-XY=2z/5,\qquad (z-XY)(5/z)=2.                           \tag{8}
```

The original evaluation obstruction has disappeared. This example is
an algebraic cancellation control, not a claim that these particular
packets arise from a factoring interval.

### A canonical lift still fails on the proper-factor interval

For the actual `N=13837` interval instance in Section 4, canceling all
four common powers of `Y-1` leaves degree-12 packets modulo `2^17`.
Their quadratic coordinates are zero. Lift each coefficient to the
balanced integer range `[-2^16,2^16)`. In ascending coefficient order,
the resulting monic polynomials `f,g` are

```text
f: [-4863,21684,-22462,44004,-22545,-56728,-5220,
    -34712,16111,27876,-49086,-12364,1]
g: [17409,-57932,-55998,-284,56815,-21400,-19044,
    52840,25071,-5660,18754,22196,1]
```

These polynomials are coprime and squarefree over the rationals. The
reduced characteristic-zero quotient therefore permits exact rational
functional calculus. Let its virtual row nodes be `z/x`, where `g(x)=0`,
and its virtual column nodes be the roots of `f`. Distinctness,
nonzero constant terms, and coprimality make its 12-by-12 Cauchy matrix
invertible.

Its normalized polynomial moment functional agrees with the original
64-by-64 functional modulo `2^17` for **every** integral polynomial.
Indeed `(Y-1)^4 f` and `(Y-1)^4 g` recover the original packets modulo
`2^17`, so both origin and infinity ratios agree identically.

Nevertheless its rational pole trace differs. The exact trace identity
reduces its value to

```math
S_{\rm virtual}=\frac{2^{15}w}{z}K,
\quad
K=\operatorname{Tr}_{\mathbb Q[Y]/g}
       \left(Y f'(Y) f(Y)^{-1}\right)-\binom{12}{2}
 =\frac{2274178272512455570345056001700682549607049910138}
        {210931526610747294728550895089720540989500093909}. \tag{9}
```

The numerator has valuation one and the denominator is odd. Thus
`v_2(S_virtual)=16`, whereas the original trace is one modulo four.
The checker verifies a rational Bezout identity for the inverse of
`f modulo g`, verifies squarefreeness, and computes the displayed trace
both from Newton sums and from the diagonal of the multiplication
matrix. No algebraic roots are approximated numerically.

This disproves the proposed rule of canceling the common packet,
choosing these canonical coefficient lifts, and then applying exact
rational functional calculus. It does not rule out a different
reconstruction that uses additional information about the original
integer domains. A larger precision can change the lifted quotient;
a valid method would need a proved precision and reconstruction bound.

What cancellation preserves is the polynomial functional. Extending
it to the singular pole-square weight still needs a separate proof.
An exact reason is supplied by the following elementary observation.
For `H=2^s` consecutive odd nodes `a_i=a_0+2i`, define

```math
\nu(h)=\sum_{i=0}^{H-1}(-1)^i\binom{H-1}{i}h(a_i,b).
```

Every integral polynomial `h` satisfies

```math
\nu(h)\equiv0\pmod {2^L},\qquad L=2H-2-s,              \tag{10}
```

by the step-two finite-difference identity. If `N=a_(i_0)b` and the
displacement separates all other nodes, the pole-square weight is the
indicator of `i=i_0` modulo two. All binomial coefficients in row
`H-1=2^s-1` are odd, so its value under `nu` is odd.

Thus even *all* integral polynomial moments at the precision in (10)
do not determine that rational contraction for general integral
atomic functionals. The zero functional and `nu` provide an explicit
pair with identical polynomial moments and different rational values.
These auxiliary functionals are not asserted to be Cauchy inverse
functionals. Additional Cauchy structure might support a scalar-only
extension; polynomial-moment agreement by itself does not prove one.

## 6. Exact checks

Run:

```sh
python3 tools/check_pole_cauchy_localization.py --output /tmp/pole-cauchy-localization.json
cmp /tmp/pole-cauchy-localization.json evidence/pole-cauchy-localization.json
```

The checker compares compressed packets with independently multiplied
finite root lists, checks the evaluation certificates in separated
rectangles, and evaluates (7) from the literal pole sum. It also checks
a large binary-described interval without enumerating its nodes,
the exact cancellation control (8), the rational lift failure (9), and
the finite-difference controls behind (10). The receipt records which controls are literal and which
are executions of the proved compressed construction.

This note analyzes the direct global extraction attempt suggested by
the existing moment quotient. It supplies no new local evaluation
oracle and no new Lean formalization.
