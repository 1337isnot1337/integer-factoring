# Factoring from unitary orbit counts in every fixed characteristic

## Result and scope

There is a shorter alternative to the power-divisor systems in
[fixed-characteristic-factor-recovery.md](fixed-characteristic-factor-recovery.md).
Count pairwise coprime factorizations, and use ordinary integer roots to
remove perfect powers before calling the counting evaluator.

For every fixed integer `k>=2`, the number of nontrivial cyclic orbits of
pairwise coprime ordered `k`-factorizations of an integer `N>1` is exactly

```math
H_k(N)=k^{\omega(N)-1}-1,                              \tag{1}
```

where `omega(N)` is the number of distinct prime divisors. All cyclic
orbits have size `k`, including when `k` is composite or `N` is a perfect
power. A canonical representative has a unique bounded inverse witness
for every required coprimality test.

In characteristic `p`, take `k=p`. The count is zero for prime powers
and is `-1` for every integer having at least two distinct prime factors.
After deterministic perfect-power preprocessing, it therefore separates
primes from composites. Binary prefix queries recover a proper divisor.

More generally, take `k=m` over any fixed nonzero commutative ring of
positive characteristic `m`. Equation (1) has the same values `0,-1`.
Require a nonzero unital ring, with effective arithmetic and zero testing
whose bit costs are charged. A fixed finite ring is the simplest case.
Zero divisors do not affect the proof or prefix recovery. In
characteristic zero, `k=2` gives a nonzero integer for every input with
at least two distinct prime factors; this does not make arithmetic on
large characteristic-zero coefficients a constant-cost operation.

The full verifier uses `3(k-1)` ordinary multiplication gadgets, integer
comparisons, and deterministic addition-by-one chains. For fixed `k`,
its total setting-bit and ordered-wire volume is `O_k(n)`, where
`n=bit_length(N)`. The largest gate load is `O_k(n^(3/4))` under the
source multiplication construction. The
[weighted network compiler](weighted-network-top-compiler.md) applies.

Thus a suitable exact counting evaluator gives deterministic complete
factoring in `poly(n,A(C_k n))` time, without a separate sparse-detector
cost or a higher-exponent-digit loop. This is a reduction theorem, not
a faster crossing construction. It proves no new unconditional
factoring exponent and makes no sublinear workspace claim. Constants
depend on fixed `k`; growing characteristic requires new accounting.

## 1. Unitary factorizations and the exact orbit count

Write the prime factorization only for the proof:

```math
N=\prod_{r\mid N}r^{a_r},\qquad \omega=\#\{r:r\mid N\}.
```

An ordered unitary `k`-factorization is a tuple of positive integers

```math
(v_1,\ldots,v_k),\qquad
\prod_i v_i=N,\qquad \gcd(v_i,v_j)=1\quad(i\ne j).     \tag{2}
```

Every prime power `r^(a_r)` must lie wholly in one coordinate. Conversely,
assigning each prime power to one of the `k` coordinates produces exactly
one tuple satisfying (2). There are therefore `k^omega` ordered tuples.

The cyclic group rotates the coordinates. Suppose a tuple is fixed by a
nonidentity rotation. Its coordinates repeat with some proper period;
every coordinate value appears at least twice. Pairwise coprimality then
forces every coordinate to equal one. This contradicts `N>1`.
Consequently the action is free, with `k` tuples per orbit.

There is one trivial orbit: one coordinate is `N` and all others are one.
Every other orbit has at least two coordinates greater than one. Hence
the number of nontrivial orbits is `k^omega/k-1`, proving (1).
The division by `k` is an integer counting argument. The algorithm counts
canonical representatives; it never divides by `k` in its coefficient
ring.

For `k=2`, a canonical nontrivial representative is simply

```math
1<a<b,\qquad ab=N,\qquad \gcd(a,b)=1.
```

These are the unordered proper unitary divisors. Their number is
`2^(omega-1)-1`, which is odd whenever `omega>=2`. A square input causes
no fixed point: a nonunit pair `(a,a)` is not coprime.

## 2. Prime powers are handled without a primality oracle

Before querying the orbit count for a current integer `M>1`, test whether
`M=b^e` for integers `b>=2` and `2<=e<=floor(log_2 M)`.
Use exact integer binary search for each `e`, and compare bounded integer
powers with `M`. Repeated multiplication may stop once the product exceeds
`M`. All temporary integers can then be bounded by a polynomial number
of bits in `bit_length(M)`, and the full loop takes polynomial time.

If a representation is found, factor the smaller base `b` recursively
and multiply all output multiplicities by `e`. Any valid representation
suffices; neither a maximal exponent nor a prime exponent is required.
The base has at most half the logarithmic size of `M`.

If no representation is found, `M` is not a proper prime power. It then
has exactly one distinct prime divisor if and only if it is prime.
By (1), the modular orbit count is zero exactly in that prime case.
A nonzero count certifies `omega(M)>=2` and admits divisor recovery.

For the full recursion, a perfect-power node produces one smaller base;
a recovered proper factor produces two smaller factors. The total number
of nodes is `O(n)`, for example by charging binary splits and power
reductions to the at most `floor(log_2 N)` prime factors counted with
multiplicity. The sum of the polynomial root-testing overheads is still
polynomial in `n`. No primality oracle or unproved distribution of primes
is used.

## 3. Coprimality has a unique bounded arithmetic witness

Let `a,b>=1` and put `delta=[b>1]`. Consider integers `u,v` satisfying

```math
0\le u<b,\qquad 0\le v<a,\qquad au=bv+\delta.         \tag{3}
```

There is exactly one such pair if `gcd(a,b)=1`, and no pair otherwise.

If `b=1`, the first bound forces `u=0`; (3) then forces `v=0`. This
is the unique witness, and `gcd(a,1)=1` as required.

Suppose `b>1`. Any witness gives `au-bv=1`, so `a,b` are coprime.
Conversely, when they are coprime, let `u` be the unique inverse of `a`
modulo `b` in `[0,b)`. This inverse is nonzero. Then

```math
v=(au-1)/b
```

is an integer in `[0,a)`. Uniqueness of the inverse and the equality
give uniqueness of both witnesses. This argument includes `a=1`, when
`u=1,v=0` for `b>1`.

The strict range bound is essential: without `u<b`, inverse witnesses
can repeat as `u+t b,v+t a`, changing the modular count.

For a tuple (2), introduce prefix products

```math
R_1=v_1,\qquad R_i=R_{i-1}v_i\quad(2\le i\le k),
\qquad R_k=N.                                        \tag{4}
```

Pairwise coprimality is equivalent to
`gcd(R_(i-1),v_i)=1` for every `2<=i<=k`.
For one direction, a prime dividing the product `R_(i-1)` divides one
of its factors; for the other, any common factor of two coordinates
appears in the corresponding prefix-product test.

Attach one copy of (3) with `a=R_(i-1),b=v_i` at each step. All prefix
products and all inverse witnesses are uniquely determined by the tuple.
Moreover,

```math
au<ab=R_i\le N,\qquad bv<ab=R_i\le N.                 \tag{5}
```

Thus the actual products in the inverse checks need at most `n` bits,
not a growing precision along the chain. Using a conservative `2n`-bit
product gadget also preserves every asymptotic claim.

## 4. Canonical representatives and nonzero-prefix recovery

Use the lexicographically least rotation of each tuple as its unique
representative. Since the cyclic action is free, a tuple is canonical
exactly when it is strictly smaller than each of its `k-1` nonidentity
rotations. Equal coordinates, including repeated ones, are permitted;
the complete rotated tuples are distinct.

Represent each coordinate by exactly `n` binary bits, most significant
first. Comparing the concatenated `kn` bits agrees with tuple
lexicographic order. A deterministic three-state comparison chain
starts in `EQ`; at the first unequal bit it enters `LT` or `GT` and
then retains that state. Requiring final state `LT` enforces one strict
comparison. Given the input bits, every state has a unique value.
Hence comparison states introduce no multiplicity.

The trivial orbit's representative is

```math
t_N=(1,\ldots,1,N).
```

For a binary prefix `w` of the tuple encoding, query the count of all
canonical unitary tuples extending `w`, then subtract one if `t_N`
extends `w`. This gives exactly the nontrivial count `C(w)`.
One may instead compile exclusion of `t_N` as a deterministic Boolean
comparison; subtracting its known contribution is simpler.

For every prefix,

```math
C(w)=C(w0)+C(w1).                                    \tag{6}
```

Start from nonzero `C(empty)`. Query the left child. If its value is
nonzero, choose it; otherwise choose the right child, which is nonzero
by (6). This argument uses only addition and equality with zero. It
holds over a ring with zero divisors as well as over a field.

After `kn` such queries the encoding is a single nontrivial tuple.
Its arithmetic and comparator extensions are unique, so the leaf count
is one. At least two components exceed one, and any component greater
than one is a proper divisor of `N`. Recursively factor that component
and its complementary quotient. This obtains the complete prime
factorization; it does not assume that the recovered component is the
least prime divisor.

There are `O_k(n^2)` counting queries across the complete recursion.
All root tests, known-trivial-tuple prefix checks, divisions, and output
handling have polynomial cost.

## 5. Linear-volume verifier compilation

Use the multiplication output ports proved in
[weighted-network-top-compiler.md](weighted-network-top-compiler.md).
Every requested multiplication scale is at most a constant multiple
of `n`. A gadget has unique internal settings given its operands and
product, `O(n)` total setting bits and ordered wire bits, maximum gate
load `O(n^(3/4))`, and polynomial-time local predicates.

There are `k-1` products in (4), and two products for each of the `k-1`
inverse equations (3): exactly `3(k-1)` multiplication gadgets.
Their combined volume is `O(kn)` and their largest gate load remains
`O(n^(3/4))` up to constants depending on the source gadget.

The checks `0<=u<b` and `0<=v<a`, the positivity of every tuple
coordinate, and `delta=[b>1]` use deterministic comparison chains of
`O(n)` states per inverse test. To check `au=bv+delta`, expose the
two product digit strings and run a binary addition-by-one chain.
Start its carry at `delta`, use

```math
z_j=(y_j+c_j)\bmod2,\qquad
c_{j+1}=\lfloor(y_j+c_j)/2\rfloor,
```

and compare every `z_j` to the corresponding bit of `au`. Require
zero final carry. A fixed-length chain and these initial/final
conditions give exactly one extension when the integer equality is
true and none otherwise. It has `O(n)` constant-size gates and wires.
Allowing `2n` product bits is harmless if the arithmetic gadget uses
that convention; unused high bits are fixed to zero.

The `k-1` lexicographic comparisons use `k n` positions each, for
`O(k^2 n)` gates and copied operand bits. Every repeated use of a
coordinate is connected by an explicit equality wire. Copies to the
prefix-product multipliers, inverse multipliers, bounds, and rotation
comparators likewise cost `O(k^2 n)` bits in total. Wire refinement
between different multiplication digit layouts is charged by bit length,
as in the compiler note. No full tuple is introduced as one giant gate.

Combining these costs gives

```math
W=O(k^2 n),\qquad D=O_k(n^{3/4}).                     \tag{7}
```

Here `W` includes setting bits and ordered wire widths; `D` is the
largest per-vertex setting and incident-wire load. All auxiliary values
are uniquely determined by a valid canonical tuple. Applying the
weighted compiler with fixed `k` and `P^40<=n` produces the same
`P`-vertex clique interface with volume `O_k(n)` and table cost
`2^O_k(n/P)`.

## 6. Transfer statement and verification

Fix the coefficient ring and the constant tuple length `k` chosen above.
Suppose an exact evaluator handles every prefix-pinned verifier from
Sections 3--5 with time `A(C_k n)`, including the required construction
and representation costs. Then deterministic complete factoring takes

```math
\operatorname{poly}\bigl(n,A(C_k n)\bigr).            \tag{8}
```

The hypothesis concerns an actual uniform evaluator. In particular,
the compiler's `2^O_k(n/P)` table cost cannot be discarded when stating
a faster target. For regular exponential scales, constant rescaling
and polynomial overhead preserve a bound such as
`exp(O(n/(log n)^alpha))` for any fixed `alpha>0`, or
`exp(O(n^beta))` for fixed `0<beta<1`, if the corresponding evaluator
bound is supplied. No such improved evaluator is constructed here.

`check_unitary_orbit_recovery.py` uses independent trial division and
literal ordinary factor tuples to verify the orbit counts, free cyclic
actions, bounded-inverse uniqueness, and complete recovery. It also
checks the deterministic comparator and addition-by-one traces, including
unit coordinates and composite tuple lengths. Negative controls expose
the need for perfect-power preprocessing, coprimality, bounded inverse
witnesses, and removal of the trivial orbit. These are finite checks of
the transcription; the arguments above prove the all-size statements.

This note introduces no Lean formalization.
