# Fourier queries for every original multiplication gate

## Result and scope

The original factoring multiplication network admits an exact Fourier
representation in which **one entry of every original gate tensor is
computable in polynomially many ring operations**. This includes the
prefix restrictions, a monomial weight on the first factor, and all carry
gates. A deterministic `exp(O(sqrt(n)))` preprocessing step supplies a
characteristic-two extension field. The transformed contraction has
exactly the original weighted value.

This removes the local-setting enumeration cost for an individual
original gate. It does not remove the cost of contracting a cluster of
gates, writing a dense table, or summing the transformed edge frequencies.
There is **no improved factoring exponent in this note**. The positive
primitive is a different exact representation of the arithmetic network,
which a more structured evaluator could use.

The source is the factoring reference at commit
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1):

- [`FactoringAux/Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean)
  defines the two-stage transforms and the arithmetic modulus.
- [`FactoringAux/MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
  defines all local settings, messages, prefix tests, and carry tests.
- [`FactoringAux/EngineCostRate.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostRate.lean)
  proves the relation between the source parameter and input length.

No new Lean proof is used. The identities below have ordinary proofs and
an exact finite-field checker,
[`check_fourier_multiplication_gates.py`](../tools/check_fourier_multiplication_gates.py).

## 1. Parameters and the character field

Write the source parameters as

```math
S=2^s,\qquad h=S,\qquad A=2S,\qquad d=2S^2=Ah,
\qquad B=2^d,\qquad M=B^4+1.
```

Thus `A²=2d`. An input gate holds `h` base-`B` digits. An arithmetic
wire holds a residue modulo `M`, represented by its least nonnegative
integer. A carry input ranges over `0<=c<2B⁴`. For the source's minimal
parameter, `n<=d²<=16n+16`; hence `d=O(sqrt(n))`.

Use the extension field

```math
\mathbb F=\mathbb F_{2^{8d}},\qquad
L=2^{8d}-1=B^8-1.
```

Choose a generator `zeta` of its multiplicative group, and put
`Omega=zeta^(L/M)=zeta^(B⁴-1)`. The element `Omega` has order `M`.
The source arithmetic transform root, which is `16 modulo M`, is a
different object; it is not this character root.

All assertions hold after embedding `F` into any nonzero commutative
`F`-algebra `R`. In particular, one may take
`R=F[X]/(X^r-1)` for the weighted nonzeroness reduction in
[`sparse-divisor-nonzeroness.md`](sparse-divisor-nonzeroness.md).

For any odd `m` dividing `L` and any root `chi` of order `m`,

```math
\sum_{\alpha=0}^{m-1}\chi^{\alpha(a-b)}
=[a\equiv b\pmod m].
\tag{1}
```

Indeed, the usual geometric sum is zero for a nontrivial field root
and is the scalar `m=1` in characteristic two for the trivial root.
This identity remains true in `R`. It uses no inverse of an arbitrary
element of `R`.

Replace every arithmetic-wire equality by (1) with `m=M`. Replace
every carry-to-carry equality by (1) with `m=L`. Both carry endpoints
are in `[0,L)`: an outgoing carry is less than `3B³`, which is less
than `L`. Thus equality modulo `L` is literal equality here.

Give the two endpoints of an edge opposite frequency signs. The
transformed gate tensor is

```math
\widehat T_v((\alpha_e)_e)=
\sum_{\text{local settings }u}
 W_v(u)\prod_{e\ni v}\chi_e^{\pm\alpha_e\,m_{v,e}(u)}.
\tag{2}
```

Expanding (1), interchanging finite sums, and grouping by vertices
proves that the contraction of (2) equals the original network value.
There is no normalization factor because all character-group orders
are odd and hence equal to one in the coefficient field.

### Deterministic preprocessing

An elementary `2^O(d)` algorithm suffices to construct the field and
roots. Enumerate monic degree-`8d` binary polynomials and test
irreducibility, even by trial polynomial factors. Factor `2^(8d)-1`
by trial division. Enumerate field elements and test whether their
orders are `L` by exponentiation for each prime divisor of `L`.
Finite fields exist and their multiplicative groups are cyclic, so
these searches terminate. Their crude costs are all `2^O(d)`.

The resulting irreducible polynomial and generator have `O(d)` bits.
Polynomial arithmetic modulo that polynomial implements a field
operation in `poly(d)` bit operations. The preprocessing cost is
`exp(O(sqrt(n)))` and the extension increases the representation of a
cyclic-ring element only by an `O(sqrt(n))` factor.

## 2. Prefix and weighted input gates

The messages of an `FA` gate are linear forms in its input digits:

```math
y_b=\sum_{i=0}^{h-1}C_{bi}u_i\pmod M.
```

For a fixed tuple of signed character frequencies, set
`alpha_i=sum_b alpha_b C_(bi) modulo M`.
The source prefix condition fixes a final interval of bit positions
inside each digit. Consequently its allowed digit values are exactly

```math
u_i\in[a_i,a_i+2^{t_i}),\qquad
0\le t_i\le d,\qquad 2^{t_i}\mid a_i.
\tag{3}
```

The unweighted second-factor gate has `a_i=0,t_i=d`. An invalid
global prefix makes the gate zero, as in the original source test.

Allow a multiplicative digit weight `W_i^(u_i)`. For the divisor
polynomial one takes `W_i=X^(2^(d position(i)))`, reducing the
exponent modulo `r` before constructing the cyclic-ring monomial.
Equation (2) then factors as

```math
\widehat T_{\rm FA}
=\prod_{i=0}^{h-1}z_i^{a_i}S_{2^{t_i}}(z_i),
\qquad
z_i=\Omega^{\alpha_i}W_i,\qquad
S_m(z)=\sum_{j=0}^{m-1}z^j.
\tag{4}
```

Both `z^m` and `S_m(z)` can be computed together by binary recursion:

```math
S_{2m}(z)=S_m(z)(1+z^m),\qquad
S_{2m+1}(z)=S_{2m}(z)+z^{2m}.
\tag{5}
```

This uses `O(log(m+1))` ring operations and no division. It remains
valid when `1-z` is a zero divisor, as can occur in the cyclic
coefficient rings. Thus (4) is computable in `poly(n)` operations;
the direct sum over `B^h` local settings is unnecessary.

## 3. Inverse-transform and middle gates

For an inverse-transform gate, let `T` be its arithmetic linear
map over `Z/MZ`, with input setting `u in (Z/MZ)^A`.
Equation (2) gives

```math
\widehat T_{\rm IB}(a,b)
=\sum_u\Omega^{a\cdot u+b\cdot Tu}
=[a+T^\top b=0].
\tag{6}
```

The equality is coordinatewise character orthogonality.

For a middle gate, the source map has the form

```math
(u,v)\longmapsto F_k\bigl((Du)\odot(Dv)\bigr).
```

Here `D` is the length-`A` arithmetic DFT and `F_k` is an
inverse-DFT/twiddle map. Both are invertible over `Z/MZ`. The source
root-sum identities give explicit inverses; `A` is a power of two
and is a unit modulo the odd integer `M`.

Change summation variables to `Du,Dv` and transform the incident
frequencies to

```math
a=D^{-T}\alpha,\qquad b=D^{-T}\beta,\qquad
c=F_k^\top\gamma.
```

The transformed tensor is a product of `A` scalar sums

```math
G_M(a,b,c)=
\sum_{u,v\bmod M}\Omega^{au+bv+cuv}.
\tag{7}
```

Let `g=gcd(c,M)`. Then

```math
G_M(a,b,c)=
\begin{cases}
0,&g\nmid a\text{ or }g\nmid b,\\
\Omega^{a u_0},&g\mid a,\ g\mid b,\quad cu_0=-b\pmod M.
\end{cases}
\tag{8}
```

To prove (8), sum over `v` first. Orthogonality leaves the solutions
to `cu=-b`. If there are solutions, they are
`u=u_0+(M/g)t` for `0<=t<g`. The sum over `t` is zero unless
`g|a`, and otherwise is the odd scalar `g=1` in characteristic two.
The value is independent of the chosen solution `u_0`.

A solution is obtained by the Euclidean algorithm after dividing
`b,c,M` by `g`. Hence (8) works for composite `M` and nonunit
`c`. It does not assume the arithmetic residue ring is a field.
All matrix products, gcds, modular inversions, and character
exponentiations take `poly(n)` operations. A middle-gate entry is
therefore zero or an `M`th root of unity, computed without enumerating
the `M^(2A)` local settings.

## 4. Carry gates also have compact exact queries

Fix one required base-`B` output digit `m in [0,B)`. An interior
carry gate has a coefficient `0<=v<B⁴+1` and incoming carry
`0<=c<2B⁴`. It requires `v+c=m modulo B` and emits
`t=(v+c-m)/B`. Its transformed tensor is the sum of

```math
\Omega^{\alpha v}\zeta^{\beta c+\gamma t}
```

over these allowed settings.

Put `H=B³` and define

```math
U=\Omega^{\alpha B}\zeta^\gamma,\qquad
Q=\zeta^{\beta B+\gamma},\qquad
Z=\Omega^{-\alpha}\zeta^\beta,
```

```math
J=\Omega^{\alpha m}S_{m+1}(Z)
+\Omega^{\alpha(m+B)}\zeta^\gamma
 Z^{m+1}S_{B-m-1}(Z).
```

Then the interior tensor is exactly

```math
\widehat T_{\rm carry}
=S_{2H}(Q)\left(S_H(U)J+U^H\zeta^{\beta m}\right).
\tag{9}
```

For a proof, write `c=Bq+s` with `0<=q<2B³` and `0<=s<B`.
The allowed coefficient is `v=Ba+(m-s mod B)`, normally with
`0<=a<B³`. The additional value `a=B³` is allowed only for
`s=m`, when `v=B⁴`. The outgoing carry is
`q+a+[s>m]`. Summing over `a,q` gives the geometric factors
in (9); splitting `s<=m` from `s>m` gives `J`. The extra
value `v=B⁴` gives the last term in parentheses.

At the first carry gate the extra condition is `c=0`. Its value is

```math
\Omega^{\alpha m}S_H(U)+[m=0]U^H.
\tag{10}
```

At the final carry gate the condition is `t=0`, and its value is

```math
\Omega^{\alpha m}S_{m+1}(Z).
\tag{11}
```

If both endpoint conditions are imposed, the value is
`Omega^(alpha m)`. All these expressions use a constant number of
geometric sums of lengths at most `2B³`, so (5) gives `O(d)` ring
operations apart from computing the character phases.

## 5. A cancellation primitive for one reciprocal phase

The simple class consisting of zero and one root phase is not closed
under contractions. For example, if `M=5`,

```math
\sum_{c=1}^4\Omega^{c-c^{-1}}=\Omega+\Omega^4.
```

The right side has multiplicative order three in characteristic two,
so it is nonzero and is not an `M`th root of unity.

There is nevertheless a more useful exact formula. Suppose `a,b`
are units modulo an odd `M`. For the unit Kloosterman sum,

```math
K_M(a,b)=
\sum_{c\in(\mathbb Z/M\mathbb Z)^\times}
 \Omega^{ac+bc^{-1}},
```

one has the characteristic-two identity

```math
K_M(a,b)=
\sum_{\substack{c\in(\mathbb Z/M\mathbb Z)^\times\\ac^2=b}}
 \Omega^{2ac}.
\tag{12}
```

The involution `c -> b/(ac)` preserves the phase. Every nonfixed
pair cancels in characteristic two. Its fixed points are exactly
the roots `ac²=b`, and at such a point the phase is `2ac`.
This proves (12) without analytic estimates.

For an odd prime power, a unit quadratic congruence has at most two
roots. Indeed, if `x²=y²` modulo `p^e`, the two factors
`x-y,x+y` cannot both be divisible by `p` because `2x` is a unit;
one factor is therefore divisible by `p^e`. The Chinese remainder
theorem gives at most `2^omega(M)` roots for general odd `M`.

If `M` has `t` distinct prime divisors, their product is at least
`t!`, so `t=O(log(M)/log log(M))`. Here this yields

```math
2^{\omega(M)}
=\exp\!\left(O\!\left(\frac{\sqrt n}{\log n}\right)\right).
\tag{13}
```

All root lists may be precomputed by squaring every unit modulo
`M` in `2^O(d)` time and storing lists keyed by their squares.
Thus, after the allowed preprocessing, (12) can be evaluated in
`2^omega(M) poly(n)` bit operations plus ring additions.

The unit hypotheses matter. The involution used in (12) need not be
defined on units when `a` or `b` is a nonunit. Formula (8) still
handles those cases at an individual middle gate, but (12) as stated
does not supply their subsequent contraction.

## 6. The remaining closure question

A point query is not a contraction algorithm. The Fourier network still
contains sums over many coupled frequency variables. Linear inverse
gates can be eliminated with (6), but that changes which frequencies
occur together in the other factors.

The most immediate obstruction to using (12) repeatedly is the input
and carry weights. A weighted reciprocal sum has terms such as

```math
\sum_{c\in(\mathbb Z/M\mathbb Z)^\times}
 S_B(\Omega^{uc})\,\Omega^{ac+bc^{-1}}
=\sum_{t=0}^{B-1}K_M(a+ut,b).
\tag{14}
```

The new factor need not be invariant under `c -> b/(ac)`.
The right side has `B` summands, and some values `a+ut` may be
nonunits. Evaluating it by expanding the interval is not the
polynomial-time gate query from Sections 2–4, and repeated use can
restore a large global search.

A successful evaluator would need a uniform representation closed under
the required weighted contractions, with the number and cost of those
contractions proved. Equations (4), (6), (8), and (9) give explicit initial
data for that question. Equations (12)–(13) show that at least one
otherwise exponential-looking frequency contraction has a much smaller
exact representation. No closure claim, width recurrence, or factoring
rate beyond the established half exponent follows yet.

## 7. Exact checks

The companion checker constructs binary fields from explicit irreducible
polynomials, verifies their irreducibility and full primitive-root cycles,
and compares formulas with direct sums. It covers composite moduli and
nonunit bilinear coefficients, digit-prefix intervals, all endpoint
choices for carry gates, geometric sums in rings with zero divisors,
and complete small product-network contractions. It also checks the
failure of the one-root phase class under a reciprocal contraction.

These finite computations test the formulas and conventions. The
identities and their uniform operation bounds are justified by the
ordinary proofs above, not inferred from the finite examples.
