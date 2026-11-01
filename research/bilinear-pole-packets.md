# Global modular packets for the complete bilinear pole product

## Constructive theorem and remaining normalization

The [pole trace](interval-pole-trace.md) uses
`C_(I,J)(Z)=sum_(u in I,v in J) 1/(Z-uv)` and notes that writing it as
`R'(Z)/R(Z)` does not supply a small product circuit. This note supplies
such a circuit **modulo a power of two**, including the whole bilinear
state family used by that trace.

Let `I,J` be nonnegative integer half-open intervals, and let

```math
q(u,v)=Auv+Bu+Cv+D,\qquad
R_q(Z)=\prod_{u\in I,v\in J}(Z-q(u,v)).             \tag{1}
```

The four coefficients may be signed. Suppose their absolute values and
all endpoints have at most `n` bits, with `n>=1`. Given `p>=1`, a uniform
deterministic algorithm constructs an ordinary addition/multiplication
circuit for `R_q` in `(Z/2^p)[Z]` with `O(n+p)` nodes. Every constant in
the circuit has at most `p` bits, binary powering is expanded into
ordinary nodes, and full generation has polynomial bit cost in `n,p`.
It enumerates neither interval nor their product.

This is a polynomial identity, valid under substitution into every
commutative algebra over `Z/2^p`, not an interpolation from a list of
values. It applies to `q(u,v)=uv`, to both boundary translations, and to
every dyadic child in the pole-trace note.

There is also a polynomial-cost evaluator for the complete trace at
unit arguments, including all truncated jets at the quadratic point
`w^2+w+1=0`. The argument `N+2^k*w` used to extract a factor count is
outside this unit region. Computing its normalized trace remains open;
the factoring exponent remains `1/2`.

## 1. Two root-residue packets suffice

For `c in {0,1}`, let

```math
S_c=\{(u,v)\in I\times J:q(u,v)\equiv c\pmod2\},\quad
m_c=|S_c|,\quad r_{u,v}=(q(u,v)-c)/2\quad((u,v)\in S_c).
```

Let `e_(c,j)` be the elementary symmetric function of degree `j` of the
multiset of these `r` values. Repeated values are allowed. Set
`J_c=min(m_c,p-1)`. Ordinary expansion gives

```math
R_q(Z)\equiv
\prod_{c=0}^1
 (Z-c)^{m_c-J_c}
 \sum_{j=0}^{J_c}(-2)^j e_{c,j}(Z-c)^{J_c-j}
 \pmod {2^p}.                                    \tag{2}
```

Every omitted term has an explicit factor `2^j` with `j>=p`.
An empty packet contributes one. This truncation concerns polynomial
coefficients and requires no division in the output ring.

Each inner sum is evaluated by Horner's rule in `O(p)` nodes. Each
remaining large power has exponent at most `|I||J|`, whose binary length
is `O(n)`, so binary powering needs `O(n)` nodes. Thus only the generation
of the first `p` elementary functions needs justification.

## 2. Compute the moments without visiting a pair

Split the input box into its four parity classes. Write
`u=e+2x`, `v=f+2y`, where `e,f` are zero or one. Ceiling division converts
`I,J` to ordinary nonnegative intervals for `x,y`. On that box, put

```math
c\equiv Aef+Be+Cf+D\pmod2,\qquad 0\le c<2.
```

Then the integer polynomial whose moments are needed is

```math
\frac{q(e+2x,f+2y)-c}{2}
 =2Axy+(B+Af)x+(C+Ae)y+
   \frac{Aef+Be+Cf+D-c}{2}.                       \tag{3}
```

Its `j`th power has degree at most `j` in each variable and at most
`(j+1)^2` monomials. The sum of every monomial separates:

```math
\sum_{x\in I_e,y\in J_f}x^a y^b
 =\left(\sum_{x\in I_e}x^a\right)
  \left(\sum_{y\in J_f}y^b\right).               \tag{4}
```

Exact prefix power sums through degree `p-1` follow from the triangular
identity

```math
t^{j+1}=\sum_{a=0}^{j}\binom{j+1}{a}
             \sum_{x=0}^{t-1}x^a .              \tag{5}
```

The division by `j+1` in solving (5) is exact over the integers. Taking
differences gives arbitrary interval sums. Summing the parity boxes
with the same value of `c` yields
`s_(c,j)=sum_(S_c) r_(u,v)^j`, including `s_(c,0)=m_c`.

Newton's identities now compute the required elementary functions:

```math
e_{c,0}=1,\qquad
e_{c,k}=\frac1k\sum_{j=1}^k(-1)^{j-1}e_{c,k-j}s_{c,j}.      \tag{6}
```

These divisions are again exact integer divisions. Only afterward are
the circuit constants `(-2)^j e_(c,j)` reduced modulo `2^p`.

All values of `q` in the box have `O(n)` bits. There are at most
`2^(2n)` roots, counted with multiplicity. A moment or retained
elementary function has `O(np)` bits; the polynomial-expansion
intermediates satisfy the same polynomial bound. There are polynomially
many integer operations: the bivariate powers and their sums use
`O(p^3)` coefficient operations, and (5)--(6) use `O(p^2)` more, apart
from ordinary integer powering. Standard binary arithmetic therefore
gives polynomial generation cost, including exact division, all
constants, and writing the `O(n+p)` output nodes.

For `q=uv`, the moments can alternatively be obtained from separate
ordinary and odd-only power sums of the two original intervals. The
checker uses this specialized formula as an independent large-input
control for the general four-box implementation.

## 3. The complete trace is cheap where both residue denominators are units

In any commutative `Z/2^p`-algebra where `Z-c` is a unit for every
nonempty packet, geometric expansion gives

```math
\mathcal C_q(Z)
 =\sum_{u\in I,v\in J}\frac1{Z-q(u,v)}
 \equiv
 \sum_{c=0}^1\sum_{j=0}^{p-1}
       \frac{2^j s_{c,j}}{(Z-c)^{j+1}}
 \pmod {2^p}.                                    \tag{7}
```

Terms from an empty packet are omitted. The error in each individual
geometric series is divisible by `2^p`, and every denominator used in
the calculation is a unit. No bound on the number of pairs is needed
to preserve the congruence after summation.

Take the unramified quadratic extension with `w^2+w+1=0` and substitute
`Z=w+T` in its jet ring modulo `T^(d+1)`. Both `w` and `w-1` have nonzero
residue in `F_4`, so (7) computes every jet coefficient through degree
`d` in polynomial bit time in `n,p,d+1`. The circuit (2) can also be
differentiated and evaluated there; its value is a unit and `R'_q/R_q`
agrees with (7).

The same argument extends to any fixed number of variables and fixed
polynomial degree: finitely many parity boxes reduce the moments to
separable power sums. This statement does not give a polynomial bound
uniform in a growing number of variables or growing polynomial degree.

## 4. Why this does not yet evaluate the factor-count filter

For factoring, the required normalized argument is
`Z=N+2^k*w`, with `N` an integer. Its residue lies in `F_2`, so the
denominator `Z-c` for the matching parity packet is not a unit. Neither
(7) nor division by the reduced product can be used there.

The lost content can be exponentially large even in a rectangle with
no factor pair. Let `H=2^r`, `r>=1`, take `I=J=[1,H+1)`, `q=uv`, and
`N=(H+1)^2`. Every product is strictly below `N`, so the exact factor
count is zero. Nevertheless, all `(H/2)^2` pairs of odd integers have
`N-uv` even. For `k` larger than the valuations of all these nonzero
differences,

```math
v_2\bigl(R_q(N+2^k w)\bigr)
 =\sum_{u,v}v_2(N-uv)\ \ge H^2/4.                \tag{8}
```

The lower bound is exponential in the binary input length. In the small
control `H=16`, both `N=289` (zero pairs) and `N=225` (one pair) give a
zero raw product modulo `2^32` at sufficiently large `k`. A normalized
pole trace distinguishes them, but raw zero testing does not.

Thus the new circuits close a modular product-generation obligation.
They do not compute the content in (8), divide it away at low precision,
or extract a rational common factor. A successful next step must handle
that singular normalization while charging the total bit cost.

## 5. Exact checks

Run:

```sh
python3 tools/check_bilinear_pole_packets.py --output /tmp/bilinear-pole-packets.json
cmp /tmp/bilinear-pole-packets.json evidence/bilinear-pole-packets.json
```

The [implementation](../tools/bilinear_pole_packets.py) writes explicit
addition/multiplication node lists. The
[checker](../tools/check_bilinear_pole_packets.py) compares full dense
polynomial coefficients on small signed bilinear inputs, exact parity
moments with literal sums, and unit jets with literal resolvent sums.
Large cases compare the general moment generator with the independent
separable formula for `q=uv`; no enormous box is expanded. The
[receipt](../evidence/bilinear-pole-packets.json) records exact counts,
circuit sizes, and hashes.

The construction extends the truncation and Newton method in
[interval-packet-circuits.md](interval-packet-circuits.md) from one
interval product to the complete bilinear image. The proof of its
global polynomial identity is independent of the finite checks. No new
Lean formalization or stronger factoring bound is claimed.
