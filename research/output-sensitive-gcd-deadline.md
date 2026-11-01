# An output-sensitive orbit gcd would give polynomial-time factoring

## Conditional theorem

The [interval orbit-gcd target](interval-orbit-gcd.md) asks for a uniform
deterministic algorithm that computes

```math
G_{N,I,J}(X)=
\prod_{\substack{d\mid N\\d\in I,\ N/d\in J}}(X-d)
```

in time polynomial in the binary input length and the output degree
`D=deg G`. Under exactly that proposed bound, deterministic factoring
would be **polynomial-time**. The earlier `exp(O(sqrt(n)))` implication is
valid but weaker than necessary: a bounded simulation removes the need
to compute a large nonempty output.

This is a strengthened conditional implication. No such orbit-gcd
algorithm is supplied, and the achieved factoring exponent remains `1/2`.

## 1. A deadline decides whether the output is empty

Suppose `A` is one fixed deterministic machine, correct on every valid
orbit-gcd input. Let `s` be that input's binary encoding length. Assume
there is a known computable polynomial `B(s)` such that

```math
D=0\quad\Longrightarrow\quad
 A\text{ halts within }B(s)\text{ steps and returns }1.       \tag{1}
```

Whenever `A` halts, assume its output is the correct monic gcd. The
output-sensitive hypothesis

```math
T_A(s,D)\le K(s+D+1)^C
```

implies (1) with `B(s)=ceil(K)(s+1)^ceil(C)`, after enlarging fixed
constants if needed. The algorithm, constants, and polynomial are
independent of the integer being factored. In a complexity-class
existence statement they can be hardwired; an executable reduction must
be supplied with an actual proved deadline.

Run `A` for at most `B(s)` simulated steps.

1. If it halts, answer whether its output differs from the constant
   polynomial one.
2. If it has not halted when the deadline expires, answer that the gcd
   has positive degree.

The second answer follows by the contrapositive of (1). A degree-zero
input cannot time out. The first answer follows from correctness of
`A`. Thus this is a deterministic, always-correct existence test.

The comparison with one is also inexpensive: a machine that halts within
`B(s)` steps can have written only `O(B(s))` output symbols. Inspect its
complete output only after halting; a partial output is not a result.
Adding a clock and simulating a fixed finite-alphabet machine has
polynomial overhead, so the existence test takes polynomial time in
`s`. It never needs to write a large positive-degree gcd.

Total termination of `A` on positive-degree instances is part of the
original target, but is not needed for this argument. It is enough that
any completed answer is correct and every degree-zero instance halts
within the proved deadline. A heuristic deadline, expected running-time
bound, or unproved empirical cutoff does not satisfy (1).

## 2. Recover the smallest prime factor with ordinary intervals

Let `N>=2`, set `b=floor(sqrt(N))`, and fix the half-open interval
`J=[2,N+1)`. Consider queries with `I` contained in `[2,b+1)`. All
endpoints have `O(n)` bits for `n=bit_length(N)`, and the complete query
encoding has length `O(n)`.

For these intervals,

```math
\deg G_{N,I,J}>0
\quad\Longleftrightarrow\quad
\text{some integer }d\in I\text{ divides }N.                \tag{2}
```

Indeed, if `2<=d<=b` divides `N`, then `2<=N/d<=N`, so its cofactor is
in `J`. Conversely, every root of this orbit gcd is such a divisor.
The equality is over the rationals, so there is no modular wraparound.

If `[2,b+1)` is empty, return `N`. Otherwise query that whole interval
with the deadline test. If it contains no divisor, `N` is prime and the
smallest prime factor is `N`.

If it contains a divisor, maintain a half-open interval known to contain
the least divisor in `[2,b+1)`. Bisect the interval. Query the lower
half. Retain that half if the answer is positive; otherwise retain the
upper half. At most `ceil(log_2(b))=O(n)` queries leave a singleton.
Its member is the least nontrivial divisor of `N`, and is therefore
prime. Integer square root, endpoint arithmetic, and the final exact
division check have polynomial bit cost.

Every call has input length `O(n)` and a polynomial deadline. Including
simulation, output inspection, and all bisection queries, the complete
finite-alphabet factoring machine therefore has polynomial running time.
No trial-division prepass or bound on the total divisor count is needed.

## 3. The useful weaker research target

The argument applies more generally when `B(s)` is any explicit
time-constructible deadline for the zero-degree case. A fixed polynomial
simulation overhead and `O(n)` calls give

```math
T_{\rm factor}(n)
\le \mathrm{poly}(n,B(cn))                              \tag{3}
```

for a fixed encoding constant `c`. In particular,

```math
B(s)=\exp\!\left(O\!\left(\frac{s}{(\log s)^\beta}\right)\right)
```

would preserve the same logarithmic saving exponent `beta`. Small input
lengths can be handled separately. The positive-degree behavior of `A`
need not have that time bound.

This isolates a weaker sufficient target than full output-sensitive gcd
construction: an exact procedure with a proved fast termination bound on
coprime interval products and no incorrect completed answers. Establishing
that one-sided bound remains open. The deadline argument supplies no way
to obtain it from the existing modular packet circuits or unit jets.

## Scope and source

This proof uses the exact root characterization and the degree-zero
runtime requirement in [interval-orbit-gcd.md](interval-orbit-gcd.md).
It sharpens the consequence of the proposed primitive; it does not
construct the primitive or change an unconditional factoring bound.
The reasoning is a finite-machine reduction, so a mock implementation
with an artificially delayed positive output would add no mathematical
evidence. No new Lean formalization is included.
