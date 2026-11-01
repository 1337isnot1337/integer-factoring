# Factor recovery in every fixed prime characteristic

## Result and scope

The [characteristic-two power-divisor reduction](parity-power-divisor-recovery.md)
extends to every fixed prime characteristic. It gives a deterministic
reduction from complete integer factorization to exact counts modulo a
fixed prime `p`. Every query is a system of multiplication equations and
binary comparisons with a total bit budget `O_p(n)`. The reduction uses
`O_p(n^2)` queries, and each arithmetic witness has a unique extension to
the compiled network.

Consequently, a uniform evaluator for these networks with bit complexity
`A(n)` gives

```math
T_{\rm factor}(n)\le\operatorname{poly}_p(n)\,(1+A(C_p n)),       \tag{1}
```

where `C_p` is a fixed constant. If `A` instead bounds circuit size and
generation time, ordinary finite-tape simulation gives
`poly_p(n,A(C_p n))`. There is no sparse-polynomial detection overhead
such as `exp(O(sqrt(n)))` in this recovery reduction.

Here `p` is fixed once for the machine. The constants are not asserted
to be uniform when the characteristic grows with the input. An extension
of the fixed prime field is permitted: all witness counts lie in its
prime subfield, and the same zero tests apply.

This is a conditional evaluator-to-factoring result. It constructs no
improved evaluator or crossing family and does not increase the proved
unconditional saving exponent `1/2`.

## 1. Detect higher base-p digits of prime exponents

Write the prime factorization only for the proof:

```math
N=\prod_r r^{a_r},\qquad n=\operatorname{bitlength}(N),\quad N\ge2.
```

For `j>=1`, let `k=p^j` and count ordered `(p-1)`-tuples of positive
integers satisfying

```math
(u_1u_2\cdots u_{p-1})^k e=N,
\qquad (u_1,\ldots,u_{p-1})\ne(1,\ldots,1).                    \tag{2}
```

The cofactor `e` is uniquely determined by the tuple. At a prime `r`,
the exponents of the `p-1` coordinates are arbitrary nonnegative integers
whose sum is at most `B_r=floor(a_r/k)`. Stars and bars gives the exact
number of permitted exponent tuples as

```math
\binom{B_r+p-1}{p-1}.
```

Therefore the count in (2) is

```math
Q_j(N)=\prod_r\binom{\lfloor a_r/p^j\rfloor+p-1}{p-1}-1.       \tag{3}
```

For every nonnegative integer `B`,

```math
\binom{B+p-1}{p-1}
=\frac{(B+1)(B+2)\cdots(B+p-1)}{(p-1)!}
\equiv
\begin{cases}
1& B\equiv0\pmod p,\\
0& B\not\equiv0\pmod p.
\end{cases}                                                   \tag{4}
```

The denominator in (4) is a unit modulo `p`. If `B` is not zero modulo
`p`, exactly one of the `p-1` consecutive numerator factors is zero
modulo `p`; if it is zero, numerator and denominator have the same
nonzero residues. This proves (4) without a division by `p`.

It follows that `Q_j(N)` is either zero or `-1` modulo `p`. It is nonzero
exactly when at least one prime exponent has a nonzero `j`th base-`p`
digit. Since every `a_r<=n-1`, it suffices to test

```math
1\le j\le\lfloor\log_p(n-1)\rfloor.                           \tag{5}
```

An empty range is allowed. If all tests vanish, every prime exponent is
less than `p`. Call this property *p-free*; it does not assert that `N`
is prime or squarefree.

If a test is nonzero, its recovered tuple in (2) has some `u_i>1`.
Put `D=product_i u_i`. Then `D>=2` and `D^k<=N`, with `k>=p>=2`.
Every nonunit coordinate is consequently a proper divisor of `N`.

## 2. Nonzero counts support binary self-reduction over any field

The recovery step does not require a parity invariant. Given a finite
set of witnesses whose total count is nonzero in a field, split it by
the next bit of its witness encoding. If the zero-child count is
nonzero, retain that child. If it is zero, retain the one child: its
count equals the nonzero parent count. Continue until the witness is
fully specified.

At a leaf there is at most one allowed witness, so a nonzero leaf count
means that the witness is valid. This last assertion depends on unique
extension of a fixed tuple to the rest of the arithmetic network.
It would not be justified by a system with uncontrolled witness
multiplicities.

For (2), fix a width `b=ceil(n/k)` for each `u_i`, concatenate their
binary encodings, and restrict an initial segment of that string.
Every permitted coordinate fits in this width because
`u_i<=D<2^(n/k)`. The all-ones tuple has one known encoding. For each
prefix, subtract its contribution if that encoding has the prefix.
The remaining cofactor and all multiplication intermediates are unique.

The same procedure will be applied to the orbit representatives in
Section 3. It returns a proper factor after reading `O_p(n)` witness
bits. It does not claim that the first factor is the least prime factor.

## 3. A p-free composite has a nonzero cyclic-orbit count

Assume the tests in (5) all vanish, so `0<a_r<p` at every prime dividing
`N`. Let

```math
\mathcal F_p(N)=\{(v_1,\ldots,v_p)\in\mathbb Z_{>0}^p:
                         v_1\cdots v_p=N\}.
```

The number of ordered factorizations is

```math
F_p(N)=\prod_{r\mid N}\binom{a_r+p-1}{p-1}.                    \tag{6}
```

Act on this set by cyclic rotation of its `p` coordinates. Since `p`
is prime, a tuple fixed by a nonidentity rotation is constant. Such a
tuple would make `N` a `p`th power. A positive p-free integer greater
than one is not a `p`th power. Every orbit therefore has exactly `p`
members.

Select the lexicographically least rotation of each tuple as its unique
representative. Exclude the one orbit with a single nonunit coordinate;
its representative is `(1,...,1,N)`. The number of remaining representatives
is the integer

```math
H_p(N)=F_p(N)/p-1.                                           \tag{7}
```

The quotient in (7) belongs to the proof. The counting algorithm does
not divide a field element by its characteristic. It directly counts
one uniquely specified representative from each orbit.

For `1<=a<p`,

```math
\binom{a+p-1}{p-1}
=\frac{p(p+1)\cdots(p+a-1)}{a!}
```

has exactly one factor of `p`, and

```math
\frac1p\binom{a+p-1}{p-1}\equiv a^{-1}\pmod p.                \tag{8}
```

Indeed, after removing the initial `p`, the numerator reduces to
`(a-1)!` and the denominator is `a!`, both units modulo `p`.

If `N` has at least two distinct prime factors, (6) has at least two
factors of `p`; hence `H_p(N)=-1` modulo `p`. If `N=r^a`, with
`1<=a<p`, then

```math
H_p(N)\equiv a^{-1}-1\pmod p.                                 \tag{9}
```

This is zero exactly when `a=1`. We have proved

```math
\boxed{\quad N\text{ p-free and }N\ge2:
       \quad H_p(N)\ne0\pmod p
       \iff N\text{ is composite}.\quad}                     \tag{10}
```

For a nonzero count, Section 2 recovers an admitted representative.
It lies outside the excluded orbit and therefore has at least two
nonunit coordinates. Any nonunit coordinate is a proper divisor.

For `p=2`, a least cyclic rotation is the smaller-first ordered pair,
and (7) counts proper divisors below `sqrt(N)`. Thus Section 3 specializes
to the squarefree interval step of the characteristic-two reduction.

## 4. Counting least rotations without multiplicity

Give every coordinate `v_i` a fixed `n`-bit encoding. Lexicographic
comparison of tuples is then ordinary lexicographic comparison of their
concatenated `pn`-bit strings. A tuple is the least rotation exactly when
it is strictly less than each of its `p-1` nonidentity rotations.

Each comparison uses a three-state deterministic chain, with states
`equal`, `less`, and `greater`. Start in `equal`. Scan the two bit strings
from most to least significant. A non-equal state remains unchanged;
from `equal`, the two current bits either leave the state equal or set
the indicated comparison result. Accept only a final `less` state.

A comparison gate stores its incoming state and the two current bits.
Its outgoing state is a deterministic function of those data. Equality
wires connect adjacent states and copy the compared bits from their
operand-input ports. The initial and final states are local weight
conditions. For fixed operand tuples, all copied bits and all states are
forced, so the comparison chain has either zero or one accepted setting.
It cannot multiply the arithmetic witness count.

There are `(p-1)pn` comparison positions, each of constant setting and
wire width. Every operand bit participates in at most `2(p-1)` comparisons.
The total added comparison setting length and wire width are `O(p^2 n)`.
Since `p` is fixed, these are `O_p(n)`. Directing each copied bit from
its existing input gate, and concatenating bits between the same pair
of gates, preserves the maximum gate-load bound `O_p(n^(3/4))`.

Strict comparison is appropriate because the p-free hypothesis rules
out fixed rotations. When a tuple has repeated coordinates but is not
constant, its `p` rotations are still distinct, and exactly one passes
every strict comparison. The finite checker includes such tuples.

The exceptional representative `(1,...,1,N)` is subtracted after the
network count. Its compatibility with any requested witness prefix is
decidable from `N` and that prefix in polynomial time.

## 5. The multiplication systems retain linear total length

Use the deterministic-output port interface from
[the weighted-network compiler](weighted-network-top-compiler.md).
Each multiplication gadget has unique internal settings for its fixed
operands and output, total setting and ordered-wire weight linear in its
requested length, maximum gate load of order length to the `3/4`, and
polynomial-time local queries. These properties come from the pinned
Rudich--Rousseau multiplication network; the compiler explains the port
modification and its source dependencies.

### 5.1 Higher-digit queries

For (2), first form `D=u_1...u_(p-1)`, then compute

```math
x_0=D,\qquad x_i=x_{i-1}^{p}\ (1\le i\le j),\qquad x_j e=N.
```

Let `b_i=ceil(n/p^(j-i))`. Every partial product forming `D` is at
most `D`, and every intermediate product forming `x_i=x_(i-1)^p` is at
most `x_i`, because all factors are positive. Bound their inputs and
outputs by `b_0` and `b_i`, respectively. Construct a `p`th power using
`p-1` ordinary multiplications by the same base. Each reused base is
linked by explicit equality wires to its operand occurrences.

The requested multiplication scales sum to at most

```math
(p-2)b_0+(p-1)\sum_{i=1}^j b_i+n=O_p(n).                       \tag{11}
```

The geometric sum of the nonrounded terms is at most `pn/(p-1)`,
and the rounding contributes only `j=O(log n)` terms. There are no
negative terms when `p=2`: forming `D=u_1` needs no multiplication.
Input-prefix restrictions and zero padding are local conditions.
All repeated-base and inter-gadget copies have total width `O_p(n)`.

The inputs `u_i`, their product and powers, and `e` determine a unique
accepted network setting exactly when (2) holds. The only removed setting
is the known all-ones tuple and its unique extension.

### 5.2 Orbit queries

Form `v_1...v_p=N` with a chain of `p-1` multiplication gadgets.
Every positive intermediate product is at most `N`, so each gadget
uses requested length `n`. Their total scale is `(p-1)n=O_p(n)`.
The operand input ports supply the bit copies used by Section 4.

Product intermediates, carries, transform values, comparison bits, and
comparison states are all uniquely determined by the tuple. A tuple
passes the network exactly when it is a lexicographically least cyclic
representative of a factorization of `N`. The tuple prefix is imposed
at the operand input ports.

### 5.3 Apply the general top compiler

Both constructions satisfy

```math
W=\sum_v\left(\sigma_v+\sum_w c_{vw}\right)=O_p(n),\qquad
\max_v\left(\sigma_v+\sum_w c_{vw}\right)=O_p(n^{3/4}).          \tag{12}
```

The weighted-network compiler therefore writes a value-preserving
`P`-vertex clique with `P^2 lambda=O_p(n)` and local setting budgets
`O_p(n/P)`, in time `2^O_p(n/P)`, for `P^40<=n`. It preserves the count
over the chosen characteristic-p field. The proof does not substitute
field-element bit length for an unrelated integer-height parameter.

Equation (12) and the top-clique compiler are the concrete interface for
using a generic crossing evaluator. A faster hypothetical evaluator must
pay for this constructor or provide its own faster accepted-input
constructor. Formula (1) does not omit that cost.

## 6. Complete recovery and complexity

At a current input `N>=2`, run the higher-digit tests (5). A nonzero
test yields a proper factor by Sections 1--2. If they all vanish, count
the least-rotation representatives of Section 3. A nonzero answer yields
a proper factor; a zero answer proves that `N` is prime by (10).

Whenever a proper factor is recovered, recurse on it and its cofactor.
There are at most `n` prime leaves with multiplicity, hence at most
`2n-1` nodes. At a node, the root tests number `O_p(log n)`, and the
selected witness encoding has `O_p(n)` bits. The entire factorization
therefore needs `O_p(n^2)` exact field-count queries and polynomial
additional integer arithmetic. Taking the smallest resulting prime
gives the repository's required output.

For a correctly compiled arithmetic circuit over a fixed field, every
operation has constant scalar bit cost. Circuit construction and ordinary
finite-tape evaluation have polynomial overhead in their charged size.
This proves the stated forms of (1). For example, each of the following
hypothetical evaluator bounds retains its exponent under the reduction:

```math
\exp(O(n/(\log(n+2))^\alpha))\quad(\alpha>0\text{ fixed}),
\qquad
\exp(O(n^\beta))\quad(\beta>0\text{ fixed}).
```

The input replacement `n -> C_p n` changes only constants in these
forms. This statement is about the recovery stage and its admitted
networks. It does not claim that the existing top schedule or crossing
construction achieves every hypothetical displayed rate.

Ordinary binary witnesses, factors, and prefixes can use linear or larger
working storage. The reduction does not itself preserve a simultaneous
sublinear-space guarantee. A field whose prime characteristic grows with
`n` would require a separate uniform bound for the `p`-dependent networks
and comparisons; it is outside this theorem.

## 7. Exact checks and attribution

Run

```sh
python3 tools/check_fixed_characteristic_recovery.py
```

The checker enumerates small integer power-divisor tuples and cyclic
factorization orbits independently of their binomial formulas. It checks
the residue predicates, free-orbit counts, unique least rotations,
deterministic comparator states, binary nonzero self-reduction, and
complete prime factorizations in characteristics `2,3,5,7`. Negative
controls cover division by the characteristic without free orbits,
composite rotation length, missing comparator initialization, and failure
to remove the trivial orbit.

The checks are finite exact references. They do not run the asymptotic
clique evaluator. The ordinary arguments above establish the all-input
reduction, subject to the stated evaluator interface. No new Lean theorem
is claimed. The underlying multiplication construction retains the
Rudich--Rousseau source attribution and Apache-2.0 license. The new
recovery proof and checker used substantial OpenAI assistant support.
