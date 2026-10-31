# Arithmetic inverse sections and exact interval aggregation

## Status and reason for changing the representation

The strongest complete factoring score remains `1/2`. The owner requested
a framework shift after the crossing searches reached repeated structural
obstructions. [PR #41](https://github.com/1337isnot1337/integer-factoring/pull/41)
checkpoints and releases that Bell lane. An obstruction for its coefficient
or pressure method is not a lower bound on factoring.

This note starts directly with the arithmetic relation

```math
xy\equiv N\pmod{2^m},\qquad x,y\text{ odd}.
```

It retains the inverse function and the integer interval restrictions.
It does not construct arbitrary clique tables. The positive results are
exact polynomial and linear-fractional descriptions of all inverse
fibres, a closed interval recursion, and an evaluator for an entire
affine fibre by Euclidean descent. These are elementary identities and
an implemented reference representation, not a claimed new number
theoretic breakthrough. The unresolved task is aggregation across fibres.

The current evaluator costs `2^(ceil(m/2))*poly(m)` in the worst case.
For factoring this is weaker than the existing half-rate algorithm;
its purpose is to expose a different arithmetic operation to improve.
No new Lean formalization is used.

## 1. Every fibre is one short polynomial

Let `M=2^m`, `B=2^s`, `Q=2^(m-s)`, with `1<=s<=m`. Fix an odd
`0<a<B`, and let `r=a^(-1) mod M`. Put

```math
b=Nr\bmod B,\qquad c=(Nr\bmod M-b)/B.
```

For every `0<=t<Q`, the inverse partner of `x=a+Bt` is

```math
y=b+B\bigl(P_a(t)\bmod Q\bigr),\qquad
P_a(t)=c+\sum_{j=1}^{\lfloor(m-1)/s\rfloor}
 (-1)^jNr^{j+1}B^{j-1}t^j\pmod Q.                 \tag{1}
```

Indeed, the geometric series for `(a+Bt)^(-1)` truncates as soon as
`sj>=m`. Subtract its constant low residue `b` and divide by `B`.
The nonconstant coefficients in (1) have exactly the divisibilities
`2^((j-1)s)` until they vanish modulo `Q`.

At `s=1`, the entire odd inverse graph is represented by one polynomial
of degree at most `m-1` with `O(m^2)` coefficient bits. This description
does not make an interval count inexpensive: the modular reduction and
the interval indicator remain part of the problem.

At `2s>=m`, equation (1) is affine:

```math
P_a(t)=c+h t\pmod Q,\qquad h=-Na^{-2}\pmod Q.       \tag{2}
```

The implementation constructs (1) and (2) without visiting any of the
`Q` high-bit completions.

## 2. Linear-fractional sections close exactly

A more economical description before the affine stage follows by
writing `x=a+Bt`, `y=b+Bu`. With `c_0=(N-ab)/B`, which is an integer,

```math
u\equiv\frac{c_0-bt}{a+Bt}\pmod Q.                 \tag{3}
```

The denominator is odd. Its integer matrix
`[[-b,c_0],[B,a]]` has determinant exactly `-N`.

More generally, consider `F(t)=(At+H)/(Ct+D) mod 2^r`, with `C` even,
`D` odd, and `AD-HC` odd. It is a permutation because

```math
F(t)-F(u)=\frac{(AD-HC)(t-u)}{(Ct+D)(Cu+D)}.
```

The denominators and determinant are units, so two values have the same
2-adic distance as their arguments. For an input digit `e in {0,1}`,
let `f=F(e) mod 2`. Then

```math
F(e+2t)=f+2F_e(t)\pmod{2^r},
\qquad
F_e(t)=
\frac{(A-fC)t+(Ae+H-f(Ce+D))/2}{2Ct+Ce+D}
\pmod{2^{r-1}}.                                  \tag{4}
```

The constant numerator is integral by the choice of `f`. The
unnormalized matrix in (4) has the **same integer determinant**
`AD-HC`. A projective normalization to denominator constant one, as in
the implementation, multiplies this determinant by an odd square.
The invariant after normalization is its unit square class.

Thus the recursion stays within four integer parameters. For the odd
inverse graph, after extracting the first input and output bits it
starts with

```math
F(t)=\frac{(N-1)/2-t}{1+2t}\pmod{2^{m-1}}.
```

Successive sections multiply the denominator's variable coefficient
by two. After a total of `s` low bits of `x` have been fixed, this
coefficient is divisible by `2^s`, while the remaining precision is
`m-s`. This recovers the affine stopping rule `2s>=m`.

## 3. Ordinary interval restrictions also close

For an interval `I=[L,U)` of residues, conditioning on a low digit `e`
leaves the interval

```math
I_e=\left[\left\lceil\frac{L-e}{2}\right\rceil,
           \left\lceil\frac{U-e}{2}\right\rceil\right)
       \cap[0,2^{r-1}).
```

The exact count therefore satisfies

```math
\mathcal C(F,r,I,J)=
\sum_{e=0}^1\mathcal C(F_e,r-1,I_e,J_{F(e)\bmod2}). \tag{5}
```

An empty interval gives zero. If either interval is the full residue
range, permutation of `F` gives the other interval's length. These
terminal rules use exact arithmetic, not equidistribution estimates.

At an affine leaf, let `t in [L,U)` and require
`V<=(ht+c mod Q)<W`. Define

```math
S_Q(n,h,c)=\sum_{i=0}^{n-1}\left\lfloor\frac{hi+c}{Q}\right\rfloor.
```

For `n=U-L` and `c'=c+hL`, the number with residue below `T` is

```math
n+S_Q(n,h,c')-S_Q(n,h,c'+Q-T).                    \tag{6}
```

Take the difference of (6) at `T=W,V`. The floor sum is computed by
Euclidean descent: remove the integral parts of `h/Q,c/Q`, then swap
the residual slope and modulus. The modulus strictly decreases, so
there are `O(m)` divisions on `O(m)`-bit integers. School arithmetic
already makes the whole leaf cost polynomial in `m`.

The recursion uses depth `O(m)` and polynomial bit space without a
memoization table. Stopping at (2) gives at most `2^(ceil(m/2)-1)`
low roots and polynomial work per root. Directly enumerating the
distinct low roots in `I` gives a second implementation of the same
count, without recursion. Neither implementation removes this root
enumeration cost.

## 4. Exact transfer back to integer factors

Let `N>=3` be odd and `m` its bit length, so `N<M=2^m`. For each
power of two `A>=2`, consider

```math
I_A=[A,\min(2A,\lfloor\sqrt N\rfloor+1)),\qquad
J_A=[\lfloor N/(2A)\rfloor+1,\lfloor N/A\rfloor+1).
```

Every divisor `2<=d<=sqrt(N)` lies in exactly one `I_A`, and its
partner `N/d` is in `J_A`. Conversely, for any positive integers in
this rectangle, `xy<2N<N+M` and `xy>0>N-M`. Consequently
`xy congruent N mod M` implies **integer equality** `xy=N`.
There are no wraparound false factors.

Summing the rectangle counts therefore counts the nontrivial divisors
up to `sqrt(N)`. Restricting `I_A` to an additional query interval
preserves this statement. Binary search recovers the least divisor,
which is prime. Even inputs return two immediately; a zero count on
an odd input returns `N`. With the current evaluator the full uniform
bit cost is `2^(ceil(m/2))*poly(m)`. Polynomially many queries and the
sizes of every integer used above are included in this bound.

## 5. The next mathematical obligation

Equation (5) is a closed arithmetic recursion. A small formula for
each state does not bound the number of distinct states that must be
evaluated. The useful new target is to sum families of its sections
while retaining the interval selectors, with complete generation and
bit cost below `exp(O(m/sqrt(log m)))`. The determinant `-N`, its
unchanged square class under normalization, and the increasing
denominator divisibility are structure available for that aggregation.

In particular, evaluating one high-degree polynomial fibre quickly,
or replacing its count by an asymptotic density with error exceeding
one, would not by itself give the required all-input algorithm.
The next experiment is exact merging and algebraic aggregation of
linear-fractional sections, including all interval boundary data.

Related earlier repository work concerns
[Fourier multiplication-gate queries on PR #23](https://github.com/1337isnot1337/integer-factoring/blob/research/beyond-half-20261030/research/fourier-multiplication-gates.md)
and [p-adic digit/carry queries on PR #26](https://github.com/1337isnot1337/integer-factoring/blob/research/further-bounds-20261030/research/padic-digit-carry-transfer.md).
Those do not aggregate the present inverse-section recursion.
For analytic background, Shparlinski's
[Modular Hyperbolas, arXiv:1103.2879v4](https://arxiv.org/abs/1103.2879v4)
surveys the distribution problem and its relation to Kloosterman sums.
No estimate from that survey is assumed in the elementary identities
or cost proofs above.

## 6. Reproduction

```sh
python3 tools/check_dyadic_hyperbola.py --output /tmp/dyadic-hyperbola.json
cmp /tmp/dyadic-hyperbola.json evidence/dyadic-hyperbola.json
make verify
```

The [reference evaluator](../tools/dyadic_hyperbola.py) is independent
of the clique evaluator. The [checker](../tools/check_dyadic_hyperbola.py)
compares floor sums with literal integer sums, sections with direct
modular inversion, both rectangle evaluators with complete scanning,
and divisor-prefix and factor recovery with trial division. The
[saved evidence](../evidence/dyadic-hyperbola.json) gives the exact
finite ranges and counts. It also counts twelve full affine fibres
of size `2^64` using permutation, without expanding them. These are
controls for the displayed proofs, not evidence of a new exponent.
