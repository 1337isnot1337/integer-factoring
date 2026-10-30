# Rough inputs remove the sparse-transfer exponent ceiling

## Result and scope

The fixed-characteristic divisor-polynomial transfer can use a preliminary
trial-division step to reduce its additional cost from
`exp(O(n/log n))` to `exp(O(sqrt(n)))`.
This removes the transfer's previous ceiling at logarithmic exponent one.
It does not construct a faster arithmetic evaluator or improve the
repository's current unconditional factoring exponent.

**Conditional transfer theorem.** Fix a prime `p`. For every `n>=2`,
`2<=N<2^n`, and divisor prefix `(k,x)`, suppose the polynomial

```math
C_{N,k,x}(z_0,\ldots,z_{n-1})
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
  \prod_{i=0}^{n-1}z_i^{d_i},
\qquad d=\sum_i d_i2^i,
\tag{1}
```

has a uniformly generated arithmetic circuit over `F_p` of size and
generation cost at most `A(n)`. The bound includes the circuit encoding
and its generation on a finite-alphabet machine. The circuit uses only
addition, multiplication, and constants from the fixed field, and is
correct after substitution into every commutative unital `F_p`-algebra.
Inputs of fewer than `n` bits are allowed with leading-zero padding.

For every integer `1<=b<=n`, a deterministic finite-alphabet machine
computes the smallest prime factor of `N` in time

```math
T(n,b)\le
\exp\!\left(O\!\left(
\log(A(n)+2)+b+\frac nb+\log n
\right)\right).
\tag{2}
```

The same bound holds for complete prime factorization. The construction
can use at most

```math
S(n,b)\le
\exp\!\left(O\!\left(
\log(A(n)+2)+\frac nb+\log n
\right)\right)
\tag{3}
```

working bits. These are ordinary circuit-storage and coefficient-vector
space bounds, not a logarithmic-space factoring claim.

In particular, the integer choice `b=ceil(sqrt(n))` gives

```math
T(n),S(n)\le
\exp\!\left(O\!\left(\log(A(n)+2)+\sqrt n\right)\right).
\tag{4}
```

Consequently, for every fixed real `alpha>0`, the hypothesis

```math
A(n)\le\exp\!\left(O\!\left(
\frac n{(\log n)^\alpha}\right)\right)
```

gives deterministic factoring in

```math
\exp\!\left(O\!\left(
\frac n{(\log n)^\alpha}\right)\right).
\tag{5}
```

There is no replacement of `alpha` by `min(alpha,1)` in (5).
For each fixed `alpha`, the limit
`(log n)^alpha/sqrt(n) -> 0` absorbs the extra term in (4).
The implicit constant may depend on the fixed evaluator and on `alpha`.

## 1. Trial division gives a smaller support bound

Test the integers

```math
2,3,\ldots,\min(2^b,N)
```

in increasing order. If one divides `N`, return the first such integer.
It is prime: a proper prime divisor of a composite first divisor would
have occurred earlier in the same loop. This also handles a prime `N`
that is no larger than the threshold.

The loop takes `2^b poly(n)` bit operations. It stores one counter and
ordinary division work registers, requiring `poly(n)` bits. It does not
store a table indexed by all trial divisors.

If the loop finds no divisor, every prime factor of `N` is strictly
greater than `2^b`. Write `Omega(N)` for the number of prime factors
with multiplicity. Then

```math
b\,\Omega(N)<\log_2N<n.
```

Since `b Omega(N)` is an integer,

```math
\Omega(N)\le h:=\left\lfloor\frac{n-1}{b}\right\rfloor.
\tag{6}
```

For each positive exponent `e`, the inequality `e+1<=2^e` gives

```math
\tau(N)=\prod_{q^e\parallel N}(e+1)
\le2^{\Omega(N)}\le K:=2^h.
\tag{7}
```

The bound is used only after the trial loop fails. It is false in
general without that condition; for example, `N=36`, `n=6`, `b=3`
has nine divisors while the expression in (7) is two.
That input is correctly handled by the trial divisor two.
If `h=0`, the rough branch cannot occur for `N>=2`, since the trial
range then already includes `N`.

## 2. A small cyclic-ring family detects each remaining prefix

For the rough branch, define the divisor-prefix polynomial

```math
G_{N,k,x}(X)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}X^d
 -[\lfloor1/2^k\rfloor=x]X.
\tag{8}
```

It removes the divisor one and preserves the divisor `N`.
Its nonzero coefficients are one, its degree is below `2^n`, and
its support has at most `tau(N)-1<=K` elements.

Here is the relevant sparse-polynomial argument, included to make
the dependence on `K` explicit. For a nonzero polynomial

```math
F(X)=\sum_{d\in S}a_dX^d,\qquad
|S|\le K,\quad 0\le d<2^n,\quad a_d\in\mathbb F_p\setminus\{0\},
```

fix any `d0 in S`. If the reduction of `F` modulo `X^r-1` is zero,
another support exponent must have residue `d0 mod r`.
Thus a prime `r` with zero reduction divides

```math
P=\prod_{d\in S\setminus\{d_0\}}|d-d_0|.
```

The empty product is one. Otherwise `1<=P<2^(n(|S|-1))`, so it
has fewer than `n(|S|-1)` distinct prime divisors. Therefore any

```math
M=n(K-1)+1
\tag{9}
```

distinct odd primes form a nonzeroness family for every such polynomial.
There are at least `M` odd primes at most

```math
Q=4(M+1)^2.
\tag{10}
```

An elementary central-binomial-coefficient proof of (10) and an explicit
trial-division prime generator appear in
[Sparse divisor nonzeroness, Sections 2--3](sparse-divisor-nonzeroness.md).
That generator takes a fixed polynomial in `Q` bit operations and can
stream the primes, retaining only a counter and one candidate.

For each prime `r` in the family, use the commutative algebra

```math
R_r=\mathbb F_p[X]/(X^r-1)
```

and substitute `z_i=X^(2^i mod r)` in (1). The circuit then produces
the first sum in (8), represented by its `r` coefficients. Subtract
the known divisor-one monomial when its prefix matches. The resulting
coefficient vector is nonzero for at least one family member exactly
when that prefix contains a divisor greater than one.

This argument works if `r=p`, too. Nilpotents in that quotient do not
affect its coefficient basis or the support-collision argument.
It uses no reconstruction of an integer divisor count and no Chinese
remainder theorem.

## 3. Recover the least prime factor and the full factorization

Start with `(k,x)=(n,0)`. The allowed prefix contains `N`, so it
contains a divisor greater than one. If `k>0`, test the lower child
`(k-1,2x)` using Section 2. Choose that child if it is nonempty;
otherwise choose `(k-1,2x+1)`. Decrease `k` by one.

The selected interval remains nonempty. The two children are consecutive
intervals, so choosing the first nonempty child preserves the smallest
eligible divisor. After `n` steps the prefix is a single integer: the
least divisor of `N` greater than one. It is prime by the argument in
Section 1. For prime `N`, that integer is `N` itself.

For complete factorization, repeat the smallest-factor procedure on the
current quotient and divide by the returned prime until the quotient is
one. Pad each quotient to the original `n` bits and use the same `b`.
The hypothesis therefore supplies the same `A(n)` bound at each call.
There are at most `n-1` calls, since the product of `n` prime factors is
at least `2^n`. The extra factor is polynomial in `n` and is absorbed
by (2). Emitting the prime factors needs only `O(n)` total output bits
apart from delimiters.

## 4. Bit cost and uniformity

The sizes in the rough branch satisfy

```math
K\le2^{n/b},\qquad
M=O(nK),\qquad Q=O(n^2K^2).
```

Represent a ring element by its `r` coefficients in the fixed field.
Addition costs `O(r)` elementary field operations and cyclic convolution
costs `O(r^2)` using the direct algorithm. Generating the monomial inputs
requires repeated doubling modulo `r` and a polynomial amount of work
in `n+r`. The field is fixed, so coefficient bit costs are bounded by
a constant.

Generate one prefix circuit within its assumed cost `A(n)`. Evaluate
it in one ring at a time, storing its description and coefficient-vector
gate values. A finite-alphabet machine can access those records by
sequential scans; even this simple implementation has a fixed polynomial
cost in `A(n)+n+Q`. The `n` prefix steps, the `M` moduli, and their
generation retain a fixed polynomial total cost in those quantities.
Combining this with the trial loop gives

```math
T(n,b)\le 2^b\operatorname{poly}(n)
        +\operatorname{poly}(A(n)+n+2^{n/b}),
```

which proves (2). Retained gate values, the generator state, counters,
and arithmetic registers have a fixed polynomial total bit length in
`A(n)+n+2^(n/b)`, proving (3). No data-dependent circuit division or
comparison is used.

The schedule `b=ceil(sqrt(n))` is computed with integer arithmetic.
It satisfies `b<=sqrt(n)+1` and `n/b<=sqrt(n)`, giving (4).
All loops, encodings, and prime choices are specified uniformly; no
asymptotic divisor bound or ineffective number-theoretic constant enters
the algorithm.

## 5. Relation to the current evaluator interface

The required circuit is the same weighted, algebra-valid circuit in
[the earlier sparse transfer](sparse-divisor-nonzeroness.md), not an
arbitrary program that happens to count divisors over the integers.
The exact divisor-to-network-setting bijection is supplied by the
reference multiplication network's
[`MulCount.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean)
and network construction. The existing note explains where to place
the digit weights and why the coefficients are one.

[Scale-one modular transfer](scale-one-modular-transfer.md) constructs
such circuits conditionally from effective scale-one crossings over a
fixed finite field. Its stated range `1/2<=rho<1` gives only
`alpha=log2(1/rho)<=1`, so (5) does not enlarge that source's parameter
range. It does remove the additional obstruction that the previous
sparse-polynomial step would impose on any stronger evaluator obtained
by another route.

The precise improvement here is the conditional transfer cost (4).
Neither the exact divisor-maximum oracle nor this rough-input step
reduces the internally added interpolation precision in the current
OpenAI Cauchy-deformation evaluator. A larger achieved logarithmic
exponent still requires a new evaluator theorem.

## 6. Reproducible finite controls

Run

```sh
python3 tools/check_rough_sparse_transfer.py \
  --output evidence/rough-sparse-transfer.json
```

The checker uses explicit divisor enumeration as an independent reference
evaluator. It is not an implementation of the hypothesized fast circuit.
The recorded run verifies:

- 6,141 smallest-factor cases: every integer from 2 through 2,048 in
  each of characteristics 2, 3, and 5;
- 1,801 cases covering every threshold `1<=b<=bit_length(N)` for
  every integer from 2 through 256;
- 511 complete prime factorizations, from 2 through 512;
- 2,104 uses of the rough-input support bound;
- 18,099 binary-prefix queries, with 523,072 cyclic reductions;
- 99,999 integer schedule controls, for `2<=n<=100000`.

The controls include a negative test showing why (7) needs prior trial
division and prime inputs for which the prefix recovery must return `N`.
The evidence file records the checker hash and representative rough
semiprime and prime cases. These finite controls support the algorithmic
details; the general complexity and transfer statements follow from the
proof above.
