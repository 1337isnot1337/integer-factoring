# An explicit polynomial inverse for binary character coefficients

This note gives a Gaussian-free inverse for the character map already proved
injective in [the weight-one reduction](weight-one-divisor-coefficients.md#41-one-bit-per-coordinate-suffices-at-every-level).
It uses the same supplied coefficient bits, preserves the full nilpotent
algebra, and avoids storing a linear system. It also identifies exactly which
character orders are necessary for decoding arbitrary binary residue vectors.

The new result concerns the decoder. The earlier note supplies the divisor
reduction, coefficient normalization, and conditional factoring transfer; the
[paired sparse reconstruction](sparse-divisor-residues.md) supplies recovery
of the divisors. None of these results supplies a faster composite-index
coefficient engine. The complete factoring saving exponent remains **1/2**.
This construction and its checks were developed with OpenAI Codex assistance.

## 1. The input map

Let q be an odd prime and write

```math
q-1=2eu,\qquad e=2^a,\qquad u\text{ odd},\qquad h=eu.
```

Set phi(T)=Phi_u(T) modulo two and use the full algebra

```math
R=\mathbb F_2[T]/(\phi(T)^e)
 =\mathbb F_2[T]/(\Phi_{q-1}(T)).                         \tag{1}
```

The second equality follows from the cyclotomic identity in characteristic
two. In R, T^h=1. Given a binary polynomial P(X) of degree below h, the
observations are

```math
M_j=P(T^j)\in R,\qquad j\mid u.                         \tag{2}
```

These are exactly the inputs of the existing decoder. For the primitive
generator g modulo q and divisor counts c_k at residue g^k, take

```math
P(X)=\sum_{k=0}^{h-1}(c_k+c_{k+h}\bmod2)X^k.
```

Every j in (2) is odd. Since T^h=1, the full weight-one coefficient
sum_(d|N) chi_j(d), reduced modulo two in its integral power basis, equals
P(T^j). Divisors divisible by q contribute zero. Decoding P therefore recovers
the signed residues modulo two by copying its coefficients to the two
opposite residue positions. No divisors are enumerated by the decoder.

The algebraic inverse below works for every even character order 2eu, whether
or not 2eu+1 is prime. Primality enters only when interpreting the coordinates
as multiplicative residues modulo q.

## 2. Obtain all odd-order evaluations by substitutions

Choose integers a,b modulo h with

```math
a=0\pmod e,\quad a=1\pmod u,\qquad b=1-a\pmod h.
```

The case u=1 uses a=0. Put omega=T^a and eta=T^b. Then

```math
T=\omega\eta,\quad \omega^u=1,\quad \eta^e=1,
\quad \eta^u=T^u.                                      \tag{3}
```

For 0<=r<u let k_r=(ar+b) mod h and j_r=gcd(k_r,u). We can choose
an integer ell_r, coprime to 2eu, such that

```math
\ell_r=k_r/j_r\pmod{eu/j_r}.                            \tag{4}
```

Indeed k_r/j_r is a unit modulo eu/j_r. This includes e=1 and k_r=0,
when that modulus can be one. A unit modulo a divisor of 2eu lifts to a
unit modulo 2eu: retain its nonzero residues at the primes in that divisor
and choose nonzero residues at the remaining primes. Trying at most 2j_r
successive lifts in (4) also gives a direct polynomial-time construction.

Substitution sigma_r:T -> T^ell_r is an automorphism of R. It permutes the
primitive u-th roots and preserves the multiplicity e; an inverse exponent
modulo h gives its inverse. Apply it to the supplied observation M_(j_r):

```math
B_r=\sigma_r(M_{j_r})=P(T^{k_r})=P(\omega^r\eta).
                                                               \tag{5}
```

Thus the tau(u) observations supply all u evaluations in (5), using only
modular polynomial composition. This step uses their canonical power-basis
coordinates, including nilpotent coordinates.

## 3. Fourier inversion separates the u residue classes

For an integer k, the geometric sum in R satisfies

```math
\sum_{r=0}^{u-1}\omega^{rk}=[u\mid k].                   \tag{6}
```

If u divides k, the sum is the odd integer u, hence one. Otherwise
1-omega^k is a unit: its value at every primitive u-th root is nonzero.
Multiplying the sum by this unit gives 1-omega^(uk)=0, proving (6).
This argument applies to the full algebra (1), even when it has nilpotents.

Write P(X)=sum_(i=0)^(h-1) p_i X^i. Equations (3), (5), and (6) give,
for 0<=s<u,

```math
S_s=\eta^{-s}\sum_{r=0}^{u-1}\omega^{-rs}B_r
   =\sum_{l=0}^{e-1}p_{s+ul}T^{ul}.                      \tag{7}
```

There is no division by an even integer. In particular, one-bit input
precision suffices for this explicit inverse, just as it does for the
previous Gaussian decoder.

## 4. Recover all nilpotent coefficients by exact division

Let V=T^u and H(T)=(T^u+1)/phi(T). Because u is odd, T^u+1 is squarefree,
so H is a unit modulo phi^e. The element V+1=phi H has nilpotent order
exactly e. Consequently (7) has a unique expansion

```math
S_s=C_s(V+1)=\sum_{l=0}^{e-1} c_l(V+1)^l,\qquad c_l\in\mathbb F_2.
                                                               \tag{8}
```

Here is an explicit coefficient extraction. At precision phi^k, start with
the canonical representative A of the current remainder. Reducing A modulo
phi gives its constant digit c in {0,1}. Replace A by

```math
\left((A-c)/\phi\right)H^{-1}\pmod{\phi^{k-1}}.           \tag{9}
```

Division in (9) is exact polynomial division. If A is changed by a multiple
of phi^k, the quotient changes by a multiple of phi^(k-1), so the update is
well defined. Equation (9) removes the first digit in powers of V+1.
Repeating e times recovers every c_l; H^{-1} is obtained by the polynomial
Euclidean algorithm. The inverse computed modulo phi^e works at every lower
precision too.

Finally expand C_s(V+1) in powers of V. From
(V+1)^(2^t)=V^(2^t)+1, the coefficient at V^k is the XOR of c_l over
indices l containing every binary digit of k. The usual successive subset
summations compute this transform in O(e log(e+1)) binary additions. These
coefficients are exactly p_(s+uk) in (7).

Equations (5), (7), and (9) are a complete inverse of (2), at every length e
and every odd u. They are the uniform argument; the finite replays below
check the implementation.

## 5. Every sampled character order is necessary for the full binary map

For j|u put d=u/j. The kernel, as an ideal of F_2[X], of substitution

```math
P(X)\longmapsto P(T^j)\pmod{\phi(T)^e}
```

is exactly (Phi_d(X)^e). To see this, work temporarily over an algebraic
closure. The j-th powers of primitive u-th roots run through all primitive
d-th roots. This follows by lifting a unit modulo d to a unit modulo u.
Since j is odd, substitution preserves each root multiplicity. Vanishing
to order e at every primitive u-th root after substitution is therefore
equivalent to divisibility of P by Phi_d^e before substitution.

The same statement holds for any odd index k, with d=u/gcd(k,u). Thus all
indices with the same d carry the same information, up to the substitutions
of Section 2. Distinct d are precisely the factors in

```math
X^h+1=\prod_{d\mid u}\Phi_d(X)^e.                        \tag{10}
```

If the index j_0=u/d_0 is omitted from (2), the nonzero polynomial

```math
P_*(X)=\frac{X^h+1}{\Phi_{d_0}(X)^e},\qquad \deg P_*<h,
                                                               \tag{11}
```

has zero observations at every retained index. Its observation at j_0
is nonzero. This proves necessity of every divisor-index observation.
More generally, sampling complete odd-character coefficient vectors requires
at least one index for each d|u to decode **arbitrary binary paired vectors**.
The tau(u) indices in (2) meet this bound.

This is an information statement for the stated full binary map. It is not
a lower bound for actual divisor instances, compressed linear observations,
other coefficient engines, or factoring algorithms.

## 6. Algorithm, cost, and implementation

The [implementation](../tools/cyclotomic_character_decoder.py) performs only
binary polynomial arithmetic and small integer gcd computations. It
constructs Phi_u modulo two by dividing X^d+1 by previously constructed
cyclotomic factors at the divisors d|u. The sum of their degrees is u, so
integer cyclotomic coefficients are never constructed. It does not factor
polynomials or replace R by its reduced quotient.

There are tau(u) calls to the supplied coefficient provider, each returning
e phi(u) bits. A naive implementation of (7) uses u^2 modular compositions
of polynomials of degree O(q), plus polynomially many modular products,
powers, gcd computations, and exact divisions. Schoolbook arithmetic makes
the entire decoder polynomial-time in q. The implementation recomputes
B_r for each s to keep workspace at

```math
O(q\tau(u)+q\log(q+1))\text{ bits},                     \tag{12}
```

including stored observations, output, polynomial factors and loop data.
This avoids the stored Gaussian system. No practical runtime advantage over
the previous bit-packed Gaussian implementation is asserted or measured.

The original decoder remains available as an independent comparison. The
new wrapper also reconstructs every supplied observation after decoding and
rejects inconsistent input, rather than assuming a genuine coefficient
provider. Prime-level validation, a primitive generator, and conversion to
residue positions have polynomial cost in q.

Substituting this decoder into the earlier reduction retains its conditional
factoring time: a coefficient provider costing F(b)poly(b,q) gives the same
F(b)poly(b,q) bound after decoding. The actual computation of those
coefficients at the original composite index remains the missing step.

## 7. Exact evidence and reproduction

Run `make verify-cyclotomic-decoder`, or

```sh
python3 tools/check_cyclotomic_character_decoder.py
```

The [checker](../tools/check_cyclotomic_character_decoder.py) generates
observations with the older exact **integer** cyclotomic root-power routine,
independently of the new binary polynomial arithmetic. It compares the new
decoder with the existing Gaussian decoder and the literal input vectors.
The [receipt](../evidence/cyclotomic-character-decoder.json) records:

- 1,036 binary paired vectors across 17 prime levels, with every vector
  enumerated at q=3,5,7,11,13,17,19;
- 120 additional vectors at general even character orders through 80;
- 504 actual divisor-residue conversions for N=2,...,64, including q|N;
- 114 explicit missing-character witnesses (11), for orders through 96;
- 11 rejected invalid inputs, plus the special q=2 zero-output control.

The large nilpotent control q=257 retains all 128 layers; the mixed control
q=337 has u=21 and e=8. Both are finite controls of the same all-order proof.
The receipt pins the decoder, checker and prior comparison source by SHA-256.
It does not certify a faster coefficient provider, an independently reviewed
mathematical result, a Lean theorem, or an improved factoring exponent.
