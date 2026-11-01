# Factoring from parity counts of one quadratic congruence

## Result and limits

The characteristic-two specialization of
[unitary-orbit-factor-recovery.md](unitary-orbit-factor-recovery.md)
admits a one-scalar verifier. After removing powers of two and testing
perfect powers, use the congruence

```math
z^2\equiv1\pmod N,\qquad 2\le z\le\lfloor N/2\rfloor,              \tag{1}
```

for odd `N>1`. Its number of solutions is exactly

```math
2^{\omega(N)-1}-1.                                                \tag{2}
```

It is odd whenever `N` has at least two distinct prime divisors.
Binary prefix parity queries recover one solution, and
`gcd(N,z-1)` is a proper divisor.

The unique arithmetic witness is

```math
z^2=Ny+1,\qquad 0\le y<N.                                         \tag{3}
```

It uses one squaring gadget, one multiplication gadget whose first
operand is the fixed input `N`, and a deterministic addition-by-one
chain. The total setting-bit and ordered-wire volume is `O(n)`, with
maximum gate load `O(n^(3/4))` under the source gadgets.

This preserves the conditional evaluator-to-factoring transfer with
no sparse-detector overhead. It does not prove a faster evaluator,
a larger logarithmic saving exponent, or a sublinear workspace bound.
The structural target is now precise: improve prefix counting for
a square and a fixed-operand linear multiplication network.

## 1. Exact count of square roots of one

Let `N` be odd and write `N=prod r^(a_r)` for its distinct odd prime
divisors `r`. If `z^2=1 (mod r^a)`, then
`r^a | (z-1)(z+1)`. The two factors have greatest common divisor
dividing two. Since `r` is odd, the entire prime power divides exactly
one factor. Thus

```math
z\equiv1\pmod{r^a}\quad\hbox{or}\quad z\equiv-1\pmod{r^a}.
```

The two residues are distinct. The Chinese remainder theorem therefore
gives exactly `2^omega(N)` roots modulo `N`.

Negation modulo `N` pairs these roots. It has no fixed point:
`2z=0 (mod N)` would imply `z=0 (mod N)` because `N` is odd, whereas
a root of one is a unit. Each pair has exactly one representative
in `[1,floor(N/2)]`. The trivial pair `{1,N-1}` contributes the
representative one. Deleting it proves (2).

If `z` satisfies (1), then `gcd(N,z-1)>1`: otherwise `z-1` would be
a unit modulo `N` and `N|(z+1)`, contradicting `2<=z<=N/2`.
The gcd is smaller than `N` because `0<z-1<N`. Thus Euclid's ordinary
integer algorithm produces a proper factor from any recovered root.

## 2. Perfect powers and complete recovery

First divide out all factors of two. For every remaining current odd
integer `M>1`, perform the exact integer-root loop from the unitary note:
test all exponents `2<=e<=floor(log_2 M)` by bounded integer powers and
binary search. If `M=b^e`, factor `b` recursively and multiply the output
multiplicities by `e`. This has polynomial overhead.

If `M` is not a perfect power, it is prime exactly when `omega(M)=1`.
Equation (2) therefore vanishes in characteristic two exactly for primes.
For a composite, its initial value is one.

For any binary prefix `w` of the `n`-bit word `z`, let `C(w)` count
solutions of (1) with that prefix. The equality
`C(w)=C(w0)+C(w1)` holds in every characteristic-two coefficient ring.
From a nonzero parent, select the left child when its count is nonzero;
otherwise select the right child. After at most `n` queries, a single
valid `z` remains. The witness `y` in (3) is unique, so the leaf count
is one.

Recover `gcd(M,z-1)` and recurse on it and the complementary factor.
Both factors must be processed: the root need not yield the smallest
prime factor. The complete recursion uses `O(n)` nodes and `O(n^2)`
counting queries, with polynomial ordinary arithmetic overhead.
For input `N=2^a`, dividing out the powers of two already finishes.

## 3. Exact verifier and resource budget

Give `z,y` `n` bits each. Use a source multiplication gadget with
requested product length `2n` for `z*z`, with the two operand inputs
connected to the same bit string by explicit equality wires.
Use a second such gadget for `N*y`, pinning the first operand to the
known binary expansion of `N`.

Connect the product strings by the deterministic addition-by-one chain
from the unitary note, with initial carry one and final carry zero.
Require every output bit of the chain to equal the corresponding bit
of `z*z`. Integer comparisons enforce `2<=z<=floor(N/2)` and `y<N`.
The nonnegative `n`-bit convention supplies `y>=0`.

Every root of (1) has precisely one quotient
`y=(z^2-1)/N` in the stated range. Conversely, (3) and the comparisons
imply (1). The source multiplication internals, equality copies,
comparison states, and carry states are unique for this `z,y`.
Therefore the verifier counts each root once.

Each multiplication contributes `O(n)` setting bits and wire bits.
The shared squaring operand is copied twice, the fixed operand is
pinned once, and each comparison or carry chain has `O(n)` constant-size
states and wires. These costs yield

```math
W=O(n),\qquad D=O(n^{3/4}),                                        \tag{4}
```

for the total weighted volume and maximum gate load defined in
[weighted-network-top-compiler.md](weighted-network-top-compiler.md).
That compiler supplies the usual `P`-vertex clique interface with
volume `O(n)` and table cost `2^O(n/P)` when `P^40<=n`.

An exact prefix-count evaluator with charged cost `A(Cn)` therefore
gives deterministic complete factoring in `poly(n,A(Cn))` time.
In particular, constant rescaling does not impose a new floor on a
hypothetical evaluator bound `exp(O(n/(log n)^alpha))`.
The source table-construction costs remain part of this hypothesis.

## 4. Fuse the verifier into one source-shaped network

The two products in (3) can share their transform and carry layers.
This gives a more specific evaluation target with exactly the original
source's local setting types and graph shape.

Use source parameters `L=2^(2s+1)`, `A=2^(s+1)`, `A^2=2L`,
`B=2^L`, and arithmetic modulus `M=B^4+1`. Let `n0=L^2` be the
input-word capacity, choose it at least `bit_length(N)`, and put
`Q=B^L-1=2^n0-1`. Introduce an unrestricted `n0`-bit word `u`.
Then the desired equation is equivalent to

```math
z^2+Nu=NQ+1,\qquad 2\le z\le\lfloor N/2\rfloor,
\qquad 0\le u\le Q.                                               \tag{5}
```

Indeed, set `y=Q-u`. Equation (5) forces `Ny=z^2-1`.
Its root interval makes `0<=y<N`, so every valid root has exactly
one such `u`. There is no need to add a separate lower bound on `u`.

Write `z_i,u_i,N_i` for the base-`B` digits. The target convolution
coefficients are

```math
c_j=\sum_{a+b=j}(z_a z_b+N_a u_b).                                 \tag{6}
```

All coefficients are nonnegative, and
`c_j<=2L(B-1)^2<B^4=M-1` for every source scale `L>=2`.
Thus taking a residue modulo `M` does not lose any coefficient.
The natural incoming carry bound `2L(B-1)` is preserved, since
`(2L(B-1)^2+2L(B-1))/B=2L(B-1)`. This fits inside the
source's `4L+1` carry bits.

The source FA layer is unchanged, with the two variable families now
holding digits of `z` and `u`. Precompute the transform of the fixed
digit string `N`. At each frequency replace the pointwise product by

```math
\widehat z(t)^2+\widehat N(t)\widehat u(t).                          \tag{7}
```

Apply the unchanged inverse transform and carry chain, now comparing
its output digits with the fixed integer `NQ+1`. The convolution
theorem and linearity give exactly (6). Its strict coefficient bound
justifies using the natural residue representative as the integer
coefficient. Keep zero initial and final carries. Arbitrary local input
pairs can make the sum exceed the `2n0`-bit output capacity, although
the desired fixed target fits; the final condition is part of the exact
output interface.

There are still two FA families of `A` gates, one MID family of `A`
gates receiving two `A`-vectors, one IB family of `A` gates, and
`A^2` carry gates. The setting domains and message widths are unchanged.
Each MID map is now quadratic in one vector and linear in the other,
with precomputed coefficients. It remains polynomial-time computable.
This is an exact construction with the same `O(n0)` volume and
`O(n0^(3/4))` largest-gate bounds as one source multiplication network.

Even the interval can use the existing prefix interface. Decompose
`[2,floor(N/2)+1)` into `O(n0)` disjoint dyadic prefix blocks.
Their parity counts sum to (2), so one block has odd count when the
input requires splitting. Scan the blocks to find one, then refine its
prefix. Every query pins only the `z` family through the ordinary source
prefix condition; `u` remains completely unrestricted. Thus this form
needs neither an extra comparator network nor a new unrestricted-input
oracle.

The primary formulas are `faOut`, `midOut`, `ibOut`, and `ib_mid_fa` in
[`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean),
and `mS`, `mV`, `mW` in
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean).
The fused verifier is proved here by linearity and the explicit no-wrap
bound, rather than being asserted as an existing source theorem.

This removes a constant number of complete arithmetic gadgets. It does
not itself improve the asymptotic contraction rate. A stronger evaluator
would need to exploit the quadratic-plus-fixed-linear MID map in (7).

## 5. Idempotents work without the odd-input restriction

A closely related statement holds for every `N>1`, including even
integers. The roots of

```math
x^2\equiv x\pmod N
```

are the idempotents of `Z/NZ`. For each prime power dividing `N`,
the equality `x(x-1)=0` forces `x=0` or `x=1`, because consecutive
integers are coprime. Hence there are again exactly `2^omega(N)` roots.

The complement map `x -> 1-x (mod N)` is free on these roots. A fixed
point would give `2x=1` and `x^2=x` modulo `N`. Multiplying the latter
equality by four and substituting `2x=1` gives `1=2` modulo `N`,
contradicting `N>1`.

The trivial orbit is `{0,1}`. For every other root `2<=x<=N-1`, its
partner is `N+1-x`. Exactly one of the pair belongs to
`[2,floor(N/2)]`; the midpoint is never an idempotent. Thus

```math
\#\{x:2\le x\le\lfloor N/2\rfloor,\ x(x-1)=0\pmod N\}
=2^{\omega(N)-1}-1.                                               \tag{8}
```

The unique quotient system is `x(x-1)=Ny`. A valid root gives the proper
divisor `gcd(N,x)`: an invertible idempotent is one, and a zero idempotent
would have representative zero.

This gives an alternative verifier with one variable-variable product
and one fixed-operand product. Its first operands differ by one, whereas
the odd-input verifier (3) uses an actual square. The latter is the
cleaner target when investigating a symmetric arithmetic contraction.

## Verification

`check_quadratic_congruence_recovery.py` enumerates every scalar residue
for small inputs. It checks the idempotent and square-root counts,
free involutions, the canonical intervals, unique quotients, and proper
gcd recovery. It runs complete parity-based factoring against independent
trial division and includes prime powers and even inputs.

The checker also replays the literal fused FA/MID/IB transforms,
including the composite auxiliary modulus at source scale one, and
compares the resulting coefficients with ordinary integer convolution.
The finite checker validates the transcription. The arguments above
and the source-compatible compiler provide the all-size statements.
No Lean formalization is added.
