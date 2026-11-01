# Compact arithmetic-gate Fourier queries in odd characteristic

## Result

The compact Fourier gate representation can be used in every fixed odd
characteristic without constructing an exponentially long field element.
The key is to change the **integer arithmetic modulus** of the FFT
network together with the character field.

Fix an odd prime `p`. Keep the source's binary digit parameters

```math
S=2^s,\quad h=S,\quad A=2S,\quad d=2S^2,\quad B=2^d,
\qquad A^2=2d.
```

Replace its arithmetic modulus and root by

```math
M=\frac{p^{4d}+1}{2},
\qquad \omega=p^4\pmod M.
\tag{1}
```

Then the same gate topology computes the same integer multiplication
relation. Its exact Fourier tensors can be queried in `poly_p(n)`
ring operations after `exp(O_p(sqrt(n)))` preprocessing. The
character field has degree only `8d=O(sqrt(n))` over `F_p`.
The preprocessing and queries over this field or its polynomial-dimension
tag algebras use polynomial space. For a general coefficient algebra,
the representation of its live elements must additionally be charged.

This also supports the Walsh signs and bit tags in
[bit-tagged-witness-recovery.md](bit-tagged-witness-recovery.md).
It is an exact initial representation and a local query bound.
It supplies no faster algorithm for contracting arbitrary blocks of
these tensors, and does not improve the factoring exponent.

The reference topology and transform identities are those in
[`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean)
and
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
at `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
The modulus substitution and the formulas below are ordinary arguments;
no Lean formalization is claimed.

## 1. The changed arithmetic modulus remains exact

The integer `M` is odd. Indeed `p^(4d)=1 modulo 8`, so
`M=1 modulo 4`. Also `gcd(p,M)=1` and

```math
\omega^d=p^{4d}=-1\pmod M.
\tag{2}
```

The transform length is `2d=A²`, a power of two. It is a unit
modulo `M`. The source's root-sum proof uses (2) and the binary
geometric-product identity, so the forward and inverse DFT remain
inverse over `Z/MZ`. The inverse normalization is obtained by the
integer Euclidean algorithm.

These arithmetic objects must not be confused with the coefficient
field. In particular `p^4` in (1) is an integer residue modulo
`M`; it is not an element used as a root in characteristic `p`.

An input has at most `d` base-`B` digits. Its ordinary convolution
coefficients are less than

```math
dB^2<2^{4d}<M,
\tag{3}
```

where the last inequality holds for every fixed odd prime `p>=3`
and `d>=1`. Thus the inverse arithmetic transform returns the exact
integer convolution coefficient when its residue is read in `[0,M)`.
There is no coefficient wraparound.

Let

```math
f=\lceil\log_2 M\rceil,\qquad C=2^f.
\tag{4}
```

Use `0<=c<C` for the incoming carry domain and `f` bits for
every arithmetic or carry message. We have `M<C<2M`. An outgoing
carry satisfies
`floor((v+c)/B)<C` whenever `v<M,c<C,B>=2`, so every local
message fits this width.

### The witness bijection is retained

Given the two input digit strings, all forward, middle, and inverse
transform settings are determined by their adjacent equalities.
Equation (3) identifies the resulting residues with the ordinary
convolution coefficients. Starting from carry zero, each next carry
is then uniquely determined. If a carry is less than `M`, the next
one is also less than `M`, because `B>=2` and the coefficient
is less than `M`. Thus every true product fits the carry domain.

The output-digit tests and final zero carry are equivalent, by
telescoping the carry equations, to equality of the integer product
with `N`. Conversely, any two integers whose product is `N`
give exactly this setting. For `N>=2` neither factor is zero.
If `N<2^(d²)`, every positive divisor and its cofactor fit the
input digit domains. This proves the same bijection between valid
settings and divisors that the original network uses.

Prefix conditions and bit-multivariate weights can be imposed on the
same input digits. They preserve this one-setting-per-divisor property.

All network size estimates retain their asymptotic form:
`f=O_p(d)`, the largest original setting has `O_p(Ad)` bits,
and the gate topology is unchanged. Repeating the explicit grouping
counts in [half-top-constructor.md](half-top-constructor.md) merely
changes constants depending on the fixed prime. In particular the
top information volume remains `O_p(n)` and top-table queries cost
`2^O_p(n/P)` when `P^40<=n`. No extra asymptotic parameter is
introduced by (1).

## 2. A polynomial-degree character field

Use

```math
\mathbb F=\mathbb F_{p^{8d}},\qquad
L=p^{8d}-1,
```

with a primitive element `zeta` of order `L`. Since

```math
L=2M(p^{4d}-1),
```

the element `Omega=zeta^(L/M)` has order `M`.
Use `Omega` on arithmetic wires and `zeta` on carry-to-carry
wires. All actual message values are smaller than `L`, so the
latter equality is literal integer equality.

For either order `m=M` or `m=L`, the scalar `m` is a unit
in `F_p`. Character orthogonality is

```math
[a=b\bmod m]
=m^{-1}\sum_{\alpha=0}^{m-1}\chi^{\alpha(a-b)}.
\tag{5}
```

As before, form each gate tensor by summing its weight times the
incident character phases, with opposite signs at the two ends of
an edge. The complete transformed contraction must additionally be
multiplied by the product of the edge normalizations `m^(-1)`.
This is a computable, nonzero prime-field scalar. It can equivalently
be incorporated into the tensor at one endpoint of each edge.

Every identity holds over an arbitrary commutative `F`-algebra,
including square-zero tag algebras and cyclic polynomial rings.

Enumerating degree-`8d` polynomials to find an irreducible one,
factoring `p^(8d)-1` by trial division, and searching for a
primitive field element all take `p^O(d)=exp(O_p(sqrt(n)))` time.
These searches can be sequential. Polynomial-basis field arithmetic
and each search counter use `poly_p(d)` bits, so the preprocessing
does not require exponential space.

## 3. The input, inverse, and middle gates

For an ordinary prefix-weighted input gate, the division-free formula
from [fourier-multiplication-gates.md](fourier-multiplication-gates.md)
is unchanged:

```math
\widehat T_{\rm FA}
=\prod_i z_i^{a_i}S_{2^{t_i}}(z_i).
\tag{6}
```

Each `z_i` combines the arithmetic-wire character phase with the
monomial input weight. Geometric sums use the binary recursion and
do not divide by `1-z_i`.

For an inverse gate with arithmetic linear map `T`, its raw Fourier
tensor is

```math
\widehat T_{\rm IB}(a,b)
=M^A[a+T^\top b=0].
\tag{7}
```

The scalar in (7) was one in characteristic two; it is retained here.

After the invertible arithmetic DFT changes of variables, a middle
gate again factors into scalar bilinear sums. If `g=gcd(c,M)`, then

```math
\sum_{u,v\bmod M}\Omega^{au+bv+cuv}
=
\begin{cases}
0,&g\nmid a\text{ or }g\nmid b,\\
Mg\,\Omega^{au_0},&
g\mid a,\ g\mid b,\quad cu_0=-b\pmod M.
\end{cases}
\tag{8}
```

Sum over `v` first, then over the `g` solutions for `u`.
The first orthogonality factor is `M` and the second is `g`.
Since `g|M` and `p` does not divide `M`, both are units in
the coefficient field. Numeric gcds and a modular inverse after
dividing by `g` compute (8) in polynomial time.

Equations (6)–(8) therefore have the same polynomial point-query
complexity as their characteristic-two versions.

## 4. A carry formula for arbitrary rectangular ranges

Changing `M` means that its remainder modulo the binary digit base
is no longer always one. The carry sum is still a constant number of
geometric products.

The following identity works for any `B>=2`, `M,C>=B`. Write

```math
M=BH+r,\quad 0\le r<B,\qquad
C=BQ+c_0,\quad0\le c_0<B.
```

Fix an output digit `0<=m<B`. For `0<=s<B`, define

```math
\kappa_s=[s>m],\qquad
t_s=m-s+B\kappa_s,\qquad
\delta_v(s)=[t_s<r],\qquad
\delta_c(s)=[s<c_0].
```

Writing `c=Bq+s`, the allowed coefficient values are
`v=Ba+t_s` with
`0<=a<H+delta_v(s)`, while
`0<=q<Q+delta_c(s)`. The outgoing carry is
`q+a+kappa_s`.

Set

```math
U=\Omega^{\alpha B}\zeta^\gamma,\quad
V=\zeta^{\beta B+\gamma},\quad
Z=\Omega^{-\alpha}\zeta^\beta.
```

The interior carry tensor is therefore

```math
\sum_{s=0}^{B-1}
\Omega^{\alpha(m+B\kappa_s)}\zeta^{\gamma\kappa_s}
Z^s\,S_{H+\delta_v(s)}(U)\,
S_{Q+\delta_c(s)}(V).
\tag{9}
```

This is a compact query despite the displayed `B`-term sum.
Clamp the following integers to `[0,B]` and sort the distinct results:

```math
0,\quad B,\quad m+1,\quad c_0,\quad
m-r+1,\quad m+B-r+1.
\tag{10}
```

They partition the integer values of `s` into at most five half-open
intervals. On each interval all three functions
`kappa_s,delta_v(s),delta_c(s)` are constant.
Replace `sum_(s=a)^(b-1) Z^s` by `Z^a S_(b-a)(Z)` on that
interval. Thus (9) uses at most five geometric products with exponents
of `O_p(d)` bits, and is computable in `O_p(d)` ring operations
apart from the initial character exponentiations.

At the initial carry gate, `c=0` gives

```math
\Omega^{\alpha m}S_{H+[m<r]}(U).
\tag{11}
```

At the final gate, outgoing carry zero gives

```math
\Omega^{\alpha m}S_{m+1}(Z).
\tag{12}
```

If both endpoint conditions hold, the value is
`Omega^(alpha m)`. These endpoint identities use `M,C>=B`.
For the original modulus `M=B⁴+1` and range `C=2B⁴`,
equations (9)–(12) reduce to the earlier carry formulas.

## 5. Walsh signs and bit tags still give compact input queries

In the odd-characteristic witness decoder, each input bit receives an
independent sign `s_j in {1,-1}` and a tag `1+epsilon_j`.
For a digit interval `u=a+v` with `0<=v<2^t` and `2^t|a`,
the high bits are fixed and the low `t` bits are independent.
Let `z` be the character phase times any cyclic monomial weight.

Up to the known product of the high-bit signs, the scalar digit sum is

```math
z^a\prod_{j=0}^{t-1}(1+s_jz^{2^j}).
\tag{13}
```

For a free-bit tag `epsilon_j`, replace its one factor in (13) by
`s_j z^(2^j)`. A fixed high-bit tag contributes either zero or
the full scalar sum. Prefix and suffix products compute all these
tagged components with polynomially many ring operations.
No division by any factor of (13) is used.

This supplies the exact original-gate queries for the signs
`s_j=(-1)^(v dot H_t(e_j))` in the polynomial-space Walsh recovery.
Every bit is tagged separately; integer powers of one tag are not
substituted for a digit's binary tags.

## 6. Scope of the construction

The changed FFT modulus, character-field construction, generic carry
sum, and Walsh-tagged input sum give a fully explicit initial tensor
family. They can be used with the same scale-one modular crossings,
because all coefficient operations occur in commutative algebras over
the fixed prime field.

They do not prove that a useful class of tensors is closed under the
required block contractions. In particular the characteristic-two
involution cancellation for a Kloosterman sum is absent in odd
characteristic. Polynomial-time original-gate queries and polynomial
space do not by themselves imply a better global contraction exponent.

The companion checker tests the new integer DFT at actual source
parameter sizes, the nontrivial normalization factors in (5), (7),
and (8), the general carry formula including both remainders, and the
Walsh/tagged digit product. Its finite comparisons supplement the
uniform proofs above.
