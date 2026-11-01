# All-input minimum-factor adapters for Sanftenberg's network

## Scope and source

These are ordinary, conditional reductions to the benchmark's
smallest-prime-factor output: one uses exact prefix counts, and one places a
min-plus objective in the network. Neither claims that the current D1 counter
computes the required values for arbitrary inputs or that the resulting
evaluator improves the factoring score.

Sources were read at
[`dbsanfte/RiemannGaussian` commit `f5fd3ddb407972f9ca65298fbe9be3225fe06635`](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635).
The public graph's accepted assignments correspond to ordered nontrivial
factor pairs for every positive `N`; see
[`SemiprimePublicAcquiredExecution.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimePublicAcquiredExecution.lean).
The closed D1 caller instead emits a two-bit root count, and its completeness
theorem is semiprime-scoped; see
[`SemiprimeD1PublicCounter.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PublicCounter.lean)
and
[`SemiprimeD1PublicCountConsequences.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PublicCountConsequences.lean).

## Exact prefix counts recover the smallest factor on every input

Let `ell = bitLength(N)`, so `N < 2^ell`. For a fixed list `pi` of wire-bit
constraints, let `C_pi` be the number of accepted assignments extending
`pi`. The source theorem identifies the number of accepted assignments with
the number of ordered pairs `(a,b)` satisfying `a,b >= 2` and `a*b=N`. A
first operand `a` determines at most one `b`, so the total number of pairs,
and therefore every restricted prefix count, satisfies

```math
0 <= C_pi <= N < 2^ell.
```

Consequently, a counter returning `C_pi mod 2^ell` returns the exact count;
zero is equivalent to there being no completion. An exact Boolean-semiring
contraction would provide the same zero/nonzero information without carrying
an integer count.

Assume one of those exact zero/nonzero prefix-query procedures is available.
Run zero-first search over the public wire bits, but order the operand-0 input
bits first in descending numerical significance, followed by every remaining
wire. At each bit, query whether the current prefix extends with that bit set
to zero; keep zero when the count is nonzero, and otherwise set the bit to
one. The exact count invariant makes this ordinary lexicographic minimization
of the satisfying assignments.

The source decoder represents operand 0 by radix digits, with binary weights
`2^(radixBits*i + bit)` for digit position `i` and bit position `bit`. Ordering
these input wires by descending weight makes lexicographic order agree with
numeric order of the first operand. The selected witness therefore has the
least nontrivial divisor `a` of `N`. If `N` is composite, that least divisor
is prime: a composite `a` would have a smaller prime divisor that also divides
`N`. It is thus the smallest prime factor. If there is no accepted assignment,
`N` is prime and the algorithm returns `N`.

This gives an all-input correctness adapter using at most `M+1` prefix
queries for `M` wires, plus polynomial-time operand decoding. It still needs
an efficient exact prefix counter. Raising the current D1 root count from
modulo four to modulo `2^ell` changes the coefficient precision; replacing
the counter with Boolean-semiring contraction changes the evaluation
semantics. Neither change has an established whole-recursion cost here.

## A one-contraction min-plus encoding

There is a more direct way to put the output objective into the network. Set
`B=2^ell` and use the finite tropical semiring

```math
S_ell={0,1,...,B},
\qquad x \oplus y=\min(x,y),
\qquad x \otimes y=\min(B,x+y).
```

Here `B` is the infinity sentinel. The additive identity is `B`, the
multiplicative identity is `0`, and both operations are associative,
commutative, and distributive as required of a commutative semiring. A
coefficient occupies `ell+1` bits.

Replace each original local Boolean table by a table over `S_ell`. A rejected
local pattern receives weight `B`. An accepted pattern receives weight zero,
except at the input vertices owning operand-0 bits: there add the bit's
numeric weight `2^(radixBits*i+bit)` whenever that bit is one. Each operand
bit is charged at exactly one input vertex. This only changes table values;
it adds no vertex or wire. Thus a complete rejected assignment has product
weight `B`; a complete accepted assignment has product weight equal to its
first operand `a`, which is less than `B`.
Contracting the network with `oplus` over assignments and `otimes` over local
tables returns

```math
\min\{a : a,b\ge2,\ ab=N\},
```

or `B` when this set is empty. For composite `N` the least divisor above one
is prime, so the finite result is the smallest prime factor; sentinel `B`
means `N` is prime and the output is `N`. This avoids count wraparound,
semiprime assumptions, and a sequence of prefix queries.

The deterministic boundary contractor is proved for commutative semirings,
so its algebraic contraction identity accepts this coefficient system. A
concrete port still has to construct every local table, implement saturated
addition and minimum on `ell+1`-bit words, and charge those operations. More
critically, the public root graph's width estimate below remains linear, and
the current D1 recursion uses jet/ring operations rather than this tropical
semiring. A recursive min-plus decomposition with sublinear acquired width
has not been established.

## Direct contraction of the public root graph is too wide

The acquired root graph has the right incidence shape for the deterministic
boundary contractor: its equality subdivision is bipartite and loopless,
and each wire occurs in exactly two local scopes. But its full-root parameter
size leaves the contractor's bound linear in the input length.

Write `T = transformLength(ell)` and `D = digitCount(ell)`. The exact edge
count in
[`SemiprimeMacroGraph.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroGraph.lean)
is

```math
M = 2 ((13T+1)(T+1)+4D-2).
```

The source proves `ell <= M <= 6400 ell` and bounds the maximum local degree
by `Delta <= 512 ell^(3/4)` in
[`SemiprimeMacroParameters.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroParameters.lean).
Applying the acquired-width theorem
[`acquired_width_cuberoot`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1BoundaryPartitionWidthRate.lean)
gives

```math
W <= M/2 + 3(M^2 Delta)^(1/3)
  <= 3200 ell + 8400 ell^(11/12).
```

Thus a direct whole-root contraction certifies only `2^(O(ell))` work. That
does not meet any positive Big O Speedrun saving exponent. The contractor
remains promising inside a smaller recursive leaf, where the relevant `M`
must first be shown sublinear and the recursive arithmetic must preserve exact
prefix existence or counts.

## Status

The prefix adapter proves a correctness bridge once exact prefix queries are
available; the min-plus adapter proves the desired output of one exact
tropical contraction. Neither extends the current D1 runtime theorem from
semiprimes to all inputs, and neither improves the repository's complete
score. The upstream whole primitive-machine bit-time theorem also remains
open. No Lean rebuild was run; this note applies the pinned source statements
by ordinary argument.
