# Retaining the arithmetic correlations as an interval orbit gcd

## Exact interface, unresolved algorithm

The [inverse-section framework](dyadic-inverse-sections.md) preserves
the arithmetic relation that arbitrary clique tables discard. Its
[complete-orbit aggregation](dyadic-orbit-aggregation.md) is effective
for uniform families, but the actual sections have correlated parameters.
Here is a second exact representation of the same factor query, without
replacing those correlations by an average.

For positive integer intervals `I=[L,U)`, `J=[V,W)`, define

```math
P_I(X)=\prod_{a\in I}(X-a),\qquad
Q_{N,J}(X)=\prod_{b\in J}(N-bX).
```

Their monic gcd over `Q` is exactly

```math
G_{N,I,J}(X)=
\prod_{\substack{d\mid N\\d\in I,\ N/d\in J}}(X-d).  \tag{1}
```

Both input polynomials have distinct rational roots, and a root `d`
of `P_I` is a root of `Q` precisely when `N=bd` for some `b in J`.
This proves (1), including empty intervals. The gcd is a monic integer
polynomial. If its degree is `D`, every coefficient has absolute value
at most `(N+1)^D`, since the sum of their absolute values is at most
the product of `1+d` over its roots. Its dense output therefore has
`O((D+1)^2*(n+1))` bits when `N<2^n`.

Equation (1) is an elementary exact reformulation. The implemented
reference expands the inputs and performs primitive pseudo-remainders;
it has no compressed-input or improved factoring bound. A binary
description of an interval product is not an inexpensive arithmetic
circuit for that product.

## 1. A precise algorithmic target and its consequence

Consider the following **unproved primitive**: given `N,I,J` by binary
endpoints of `O(n)` bits, compute the monic gcd (1) deterministically
in `(n+D+1)^C` bit operations, for one absolute constant `C`. The cost
must include all generation, zero testing, normalization, scalars,
and coefficient output. It must work for degree zero as well.

If this primitive existed, it would imply unconditional deterministic
factoring in `exp(O(sqrt(n)))` bit operations. This is a conditional
transfer theorem, not an implementation or an improved bound.

To prove the implication, trial-divide up to `2^b`, stopping on a
divisor. If none exists, every prime factor exceeds `2^b`, so the
number of prime factors with multiplicity is at most `n/b`. Hence

```math
\tau(N)=\prod_p(e_p+1)\le2^{\sum_p e_p}\le2^{n/b}.
```

Every gcd query has `D<=tau(N)`. The dyadic rectangles and binary
search in the inverse-section note recover the least divisor using
polynomially many calls to the primitive; over `Q`, their exactness
also follows directly from (1). Thus the total cost would be

```math
2^b\operatorname{poly}(n)+
 \operatorname{poly}(n,2^{n/b})=
\exp(O(b+n/b+\log n)).                            \tag{2}
```

Choose the integer `b=ceil(sqrt(n))`, using trial division alone for
the finitely many small boundary cases. Equation (2) gives the claim.
For every fixed scoring exponent `alpha`, `sqrt(n)` is eventually
smaller than `n/(log n)^alpha`. Consequently a proof of this primitive
would go far beyond a small improvement on `1/2`. There is currently
no proof that the primitive can be constructed with the stated cost.

Even a weaker coprimality decision procedure is relevant: a uniform
cost `F(n)` for deciding `G=1` would give factoring in
`poly(n)*(F(n)+1)` by interval bisection. This version avoids paying
to write all the divisor roots, but its existence is also open here.

## 2. The inputs have short boundary recurrences

These products are not arbitrary high-degree polynomials. Cancelling
the overlapping factors gives the exact identity

```math
(X-U+1)P_I(X+1)=(X-L+1)P_I(X).                    \tag{3}
```

For `h=W-V` and `S_N(X)=NX/(N+X)`, the reciprocal product satisfies

```math
(N+X)^h\bigl(N-(W-1)X\bigr)Q_{N,J}(S_N(X))
=N^h\bigl(N-(V-1)X\bigr)Q_{N,J}(X).              \tag{4}
```

Indeed `N/S_N(X)=N/X+1`, so (4) is (3) in the reciprocal coordinate,
with the homogeneous degree factors restored. The two relevant
actions are an additive translation and its conjugate by inversion.
The open constructive question is whether their common finite orbit
factor can be eliminated from these boundary recurrences while paying
for the small output rather than the long input products.
The [unit-jet algorithm](interval-unit-jets.md) now provides all local
Taylor data at a suitable quadratic 2-adic point with polynomial bit
cost in the endpoint length, precision, and jet degree. Extracting the
global common factor from this data and the recurrences remains open.
The [global packet circuits](interval-packet-circuits.md) further give
small ordinary circuits for both input polynomials modulo `2^p`, with
polynomial generation cost. Large artificial modular common factors
persist in an explicit family, so a rational-factor certificate or
effective multiplicity removal is still required.
The [normalized integer-center evaluator](normalized-interval-jets.md)
now separates an enormous local content exponent from a short unit
jet with polynomial bit cost. It does not select a common root or
provide the missing global factor certificate.

Equations (3)--(4) specify a new elimination problem. They do not
authorize treating a long product, a shifted evaluation, or a large
integer power as a unit-cost operation. In particular, a short
description of the recurrence is not a proof of a fast solution.

## 3. A global rational-part projection and its specialization trap

Let `G` be the gcd in (1) and put `R=(P'/P)(Q'/Q)`. Both input
polynomials have simple rational roots. Thus

```math
R(X)=\sum_{a\in I,\ b\in J}
 \frac1{(X-a)(X-N/b)}.
```

If `a!=N/b`, a summand has only simple poles, and its antiderivative
is a difference of logarithms. If `a=N/b`, its antiderivative is
`-1/(X-a)`. Consequently the unique proper rational part of an
antiderivative of `R` is exactly

```math
\mathcal H(R)=-\sum_{\substack{a\in I,b\in J\\ab=N}}\frac1{X-a}
=-G'(X)/G(X).                                    \tag{5}
```

This describes the desired output as a linear rational-part
projection, with denominator degree equal to the number of genuine
factor witnesses. Dense Hermite reduction would compute it, but
charges for the long input degrees. A fast algorithm for these
implicitly described inputs is the new question in
[Discussion #53](https://github.com/1337isnot1337/integer-factoring/discussions/53).
Equation (5) is an identity, not an algorithm for that projection.
The [global pole-trace continuation](interval-pole-trace.md) moves the pole
to a separate product parameter, proves a polynomial-precision normalized
filter, and closes both boundary translations in a bilinear state family.
Its exact dyadic reference still has an enumerative state count; a compact
state description is not an improved elimination bound.

Treating `N` as an indeterminate does not solve the problem. Over
`Q(N)`, none of the two root families coincide, so the rational part
is identically zero. Specialization can create common roots exactly
when

```math
\prod_{a\in I,b\in J}(N-ab)=0.                   \tag{6}
```

For `I=J=[2,4)`, this exceptional polynomial is
`(N-4)(N-6)^2(N-9)`. At `N=5` the gcd is one; at `N=6` it is
`(X-2)(X-3)`, so the rational part becomes
`-(2X-5)/((X-2)(X-3))`. Generic integration followed by specialization
therefore loses precisely the resonant terms sought by factoring.
A constructive method must handle these specializations, rather than
assume that the generic rational-part projection commutes with them.

## 4. What existing algebraic results actually supply

Kaltofen's [Single-Factor Hensel Lifting, 1987](https://userpages.cs.umbc.edu/phatak/691a/supl/hensel-kaltofen-1987.pdf)
proves a bound on the straight-line complexity of factors and gcds in
terms of input program size and output degree. Its introduction
explicitly says that polynomial-time construction is lost when input
degrees are unbounded, including because required scalars may have
exponential size. It therefore does not establish the uniform bit-cost
primitive above, even if short input programs were supplied.

Bostan, Gaudry, and Schost's
[Linear Recurrences with Polynomial Coefficients, 2007](https://specfun.inria.fr/bostan/publications/BoGaSc07.pdf)
give effective recurrence and integer-factorization algorithms.
The recurrence bound still depends essentially on the square root of
the numerical index, rather than polynomially on its bit length.
Applying that algorithm to (3) does not give the proposed primitive.

Bhattacharjee, Kumar, Ramanathan, Saptharishi, and Saraf's
[2025 deterministic circuit-factorization theorem](https://arxiv.org/abs/2504.08063v2),
Theorem 1.2, assumes constant depth and polynomial degree and size in
the number of variables, over the rationals (with an extension to
sufficiently large characteristic). Our exponentially high-degree
circuits over `Z/2^p` do not meet those hypotheses. This theorem is
not a ready-made extraction step for the packet circuits.

There is another obstruction to a shortcut through small primes:
reducing the two polynomials modulo a small prime can create common
roots that are absent over `Q`. For example, `N=35`, `I=[2,5)`, and
`J=[2,9)` have gcd one over `Q`, but `x=2,y=4` give
`xy congruent N mod 3`. Any modular implementation needs a proved
separation and reconstruction cost; a single small-prime gcd is not
an exact factor decision.

## 5. Replay and current scope

```sh
python3 tools/check_interval_orbit_gcd.py --output /tmp/interval-orbit-gcd.json
cmp /tmp/interval-orbit-gcd.json evidence/interval-orbit-gcd.json
```

The [dense reference](../tools/interval_orbit_gcd.py) never calls a
factoring oracle. The [independent checker](../tools/check_interval_orbit_gcd.py)
compares its integer polynomial coefficients with trial-divisor
products on 1,024 selected rectangles and all 688 dyadic factor boxes
for inputs 2 through 256. It also checks both boundary recurrences
at exact rational arguments and the small-prime false-root example.
It independently compares 192 literal pole decompositions and their
rational parts with logarithmic derivatives of the dense gcd, and
retains the `N=5,6` specialization control above.
The [receipt](../evidence/interval-orbit-gcd.json) records those finite
controls. This is an exact interface and a source audit; the required
compressed elimination algorithm remains unresolved. The strongest
complete factoring score remains `1/2`.
