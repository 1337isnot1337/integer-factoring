# Deterministic factoring with a square-root logarithmic saving

## Theorem and evidence level

There is one deterministic finite-alphabet Turing machine which computes
the smallest prime factor of every binary integer `N >= 2` and satisfies

```math
T(n)\le C\exp\!\left(\frac{Cn}{\sqrt{\log n}}\right)
\qquad(n\ge2),
\tag{1}
```

where `n` is the input bit length and `C` is an absolute constant. In the
site's scoring convention this gives exponent `alpha=1/2`.

This is a complete ordinary mathematical composition argument. It uses
the existing Appendix B clique-evaluation theorem of Claude (Anthropic)
and the existing multiplication-network reduction of Isaac Rudich and
Louis-Martin Rousseau. The new work supplies and pays for their written-table
connection. It does not claim a new crossing decomposition or a new clique
evaluator. The complete composition is now formalized by
[`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](../formal-appendix/HalfTheorem.lean).
Independent review and leaderboard acceptance remain separate from that
formal theorem; neither is asserted by this note.

The two source revisions are fixed:

- Factoring reference:
  [`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
- Appendix B evaluator:
  [`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).

The top constructor is proved in [half-top-constructor.md](half-top-constructor.md).
The evaluator, encoding, and explicit pad are audited in
[half-evaluator-audit.md](half-evaluator-audit.md). The proof below states
the resulting single algorithm and accounts for all its costs.

## 1. The prior evaluator, with its actual input charged

The assembled declarations `AppRunningTime.thmHalfPieces` and
`AppRunningTime.thmRunningTime` provide absolute constants `c,C0` and,
after fixing an integer `x0>=18` with `2^x0>=c+5`, one deterministic
finite-symbol machine fragment `run`. The public statement is `PartB`.
Its input consists of an explicit unary pad, a unary level, a unary label
width, and all signed integer tables. It returns their exact integer
clique contraction.

Use the source's schedule

```math
K_0=2^{x_0},\qquad
g_d=2^{3d+x_0}K_d,\qquad K_{d+1}=g_d^2.
\tag{2}
```

For `P=K_ell`, positive widths `lambda,B`, and an explicitly written
`P`-vertex clique whose rows have `2^((P-1)lambda)` entries bounded in
absolute value by `2^B`, the sufficient pad exponent obeys

```math
\operatorname{need}(\ell,\lambda,B)
\le e_\ell P^2\lambda+C_0P^2\log_2(P(\lambda+B)),
\qquad e_\ell\le0.85\,2^{-\ell/2}.
\tag{3}
```

One fixed polynomial in the complete encoded input length bounds the
fragment's actual machine steps. In particular, (3) is not a unit-cost
oracle for unwritten tables. The source's recursion proof pays for its
signed intermediate values, denominators, groupings, and every recursive
level. The exact input encoding and the assembled dependency chain are
given in the evaluator audit. This use requires the assembled proof,
not merely the statement file.

All constants and the fragment are fixed once. Their existence does not
permit input-dependent advice. Their potentially enormous values do not
give a practical factoring implementation.

## 2. Select one level from the input length

For `n>=K_0^40`, choose the largest `ell` with `K_ell^40<=n`, and put

```math
P=K_\ell,\qquad r=2^{\lfloor\ell/2\rfloor}.
\tag{4}
```

For smaller lengths, use trial division and stop. This is a fixed finite
initial range. Handle `N=0,1` explicitly as `2,1` if retaining the
reference's all-natural-number convention.

The source proves

```math
\begin{aligned}
\log_2 K_d&=(3x_0+6)2^d-6d-2x_0-6,\\
K_d&<K_{d+1}\le K_d^{10}.
\end{aligned}
\tag{5}
```

Consequently

```math
P^{40}\le n<K_{\ell+1}^{40}\le P^{400},\qquad r\le P,
\tag{6}
```

and

```math
\log_2n<80(3x_0+6)2^\ell
\le160(3x_0+6)r^2.
\tag{7}
```

Successive integer powers and comparisons find `ell`. The next candidate
has polynomially many bits in `n`; (5)--(6) bound its bit length by
`O(log n)` above the fixed initial range. No real logarithms or powers
are computed by this algorithm. Equation (7) implies
`n/r=O(n/sqrt(log n))`.

## 3. A prefix count becomes a written clique

For a query `0<=k<=n`, `0<=x<2^(n-k)`, define

```math
D(N,k,x)=\#\{d\mid N:\lfloor d/2^k\rfloor=x\}.
```

The top-constructor lemma provides a deterministic procedure
`Construct(N,P,k,x)` that writes an instance `I` and positive integers
`lambda,B` with

```math
\operatorname{WFe}(P,\lambda,B,I),\quad
\operatorname{zval}(P,I)=D(N,k,x),\quad
P^2\lambda\le4717n,\quad B\le6769n/P.
\tag{8}
```

Its complete bit cost and written input size are `2^O(n/P)`, with constants
independent of the prefix, `P`, and the level.

The construction retains the factoring reference's top grouping
`G=3P^8`; `G` is distinct from the recursive `g_d` in (2). It searches
the `P^G` colorings deterministically. A valid grouping exists by the
reference's original theorem. Its local tables are built by enumerating
the original gate settings in each final group and testing equality of
their deterministic messages. An internal edge contributes a 0/1 test,
so this procedure does not perform an uncharged sum over internal channels.

The grouping and materialization cost is bounded explicitly by

```math
n^a\left(P^G+P\,2^{(P-1)\lambda+H}+1\right),
\qquad H=B-1\le6768n/P,
\tag{9}
```

for a fixed integer `a`. Because `P^40<=n`,
`log2(P^G)=3P^8 log2P<=3n/P`; (8)--(9) give `2^O(n/P)`.
Binary arithmetic and all table output are included.

At each vertex the table contains only the `P-1` off-diagonal edge blocks.
The source's zero self block is removed; Appendix B inserts an empty
self block. Each pair uses the same original crossing-edge order at
both endpoints, and added channel bits are forced to zero. The codec
is a bijection between the two summation domains. It preserves the exact
integer value, including the unique extension of a globally valid gate
setting to its boundary labels. This proves the value equality in (8),
not just a size comparison.

## 4. Write a sufficient computable pad

The existential public theorem does not assert that its function `need`
is computable. Use the following integer instead:

```math
q=\left\lfloor\frac{P^2\lambda}{r}\right\rfloor+1
  +C_0P^2\operatorname{size}\bigl(P(\lambda+B)\bigr),
\tag{10}
```

where `size` is binary bit length. Since `e_ell<=1/r`,
`log2(t)<=size(t)` for positive integers, and the floor plus one dominates
the unrounded quotient, (3) implies `need(ell,lambda,B)<=q`.

This pad is affordable. Equation (8) and `P>=2` imply
`P(lambda+B)<=9128n<2^14 n`, so

```math
\operatorname{size}(P(\lambda+B))\le\operatorname{size}(n)+14.
```

From (6), `size(n)<=400 size(P)<=400P`. Thus

```math
rP^2\operatorname{size}(P(\lambda+B))
\le414P^4\le414n.
```

Multiplying (10) by `r`, using (8) and `r<=P<=n`, gives the explicit bound

```math
qr\le(4718+414C_0)n.
\tag{11}
```

The algorithm writes exactly `2^q` pad units, using a binary counter or
repeated doubling. Its time is polynomial in the written pad length.
Neither `need` nor a real-number approximation is an instruction.

## 5. Run the evaluator and account for its full input

Let `J=2^((P-1)lambda)` and `F=2^q`. The actual four-symbol input encoding
has length at most

```math
14+4(F+\ell+\lambda+P)+(3B+10)PJ.
\tag{12}
```

This follows from the source's explicit list, product, unary, and signed
binary integer encodings. There is no arbitrary encoding choice at this
input boundary. Equations (8), (11), and `r<=P` make (12)
`2^O(n/r)`. Building the encoding, running the fixed polynomial-time
fragment, and reading its integer output therefore cost `2^O(n/r)`.

The result is exactly `D(N,k,x)` by (3) and (8). It is a nonnegative
integer at most the number of divisors, hence at most `N`. The signed
binary output format does not change the value or introduce an oracle.

The source fragment has finite control, finitely many stacks, and a
four-symbol stack alphabet. It can be attached to the explicit
constructor and arithmetic routines and simulated by one ordinary
finite-alphabet multitape Turing machine with a fixed polynomial
overhead. A fixed polynomial changes only the constant multiplying
`n/r`. This supplies actual bit-machine cost for the composition.

## 6. Recover the smallest prime factor

Start with the empty prefix, numerically `p=0`. For `k=n-1,n-2,...,0`:

```text
y = 2*p
count = CountPrefix(N, k, y) - [floor(1/2^k) = y]
if count > 0:
    p = y
else:
    p = y+1
return p
```

The invariant is that some divisor greater than one has the current
prefix, and this is the smallest such prefix. Initially `N` is one such
divisor. The subtraction removes divisor one exactly. If the left
extension has no remaining divisor, the right extension must have one.
Choosing the left extension whenever possible preserves minimality.
After `n` queries, `p` is the smallest divisor greater than one.
It is prime: a composite value would have a smaller divisor greater
than one. Thus the output is `Nat.minFac N`.

There are exactly `n` prefix counts, each with the same parameter schedule
chosen from the original bit length. All queries satisfy the range in
Section 3. The factor `n`, binary-prefix bookkeeping, and any fixed
polynomial machine simulation are absorbed into `2^O(n/r)` because
`r<=P<=n^(1/40)`. Equation (7) proves (1).

This is one algorithm for every input. Neither the level nor the chosen
input prefix selects a new machine or a new polynomial degree.

## 7. Relation to the reference conjecture and the score

The statement of
[`SqrtLogLogConjecture`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/Challenge.lean)
uses the integer denominator
`D(n)=floor(sqrt(floor(log2 n)))+1` and the time bound
`C*2^floor(C*n/D(n))`. This denominator is within fixed positive factors
of `sqrt(log n)` for large `n`. Choose an integer `C` large enough to
absorb the fixed runtime constants and a factor two for rounding the
exponent down. Increase `C` again to cover the finite fallback and every
small input length. Thus the ordinary argument also proves that exact
integer-bound form. It has not produced a Lean term of that proposition.

For the leaderboard, for every fixed `epsilon>0`, (1) is eventually at
most

```math
2^{n/(\log n)^{1/2-\epsilon}}.
```

After taking logarithms, the ratio of the leading term from (1) to the
allowed exponent is a constant times `(log n)^(-epsilon)`, which tends
to zero. The same machine works for every `epsilon`; only the threshold
changes. This establishes score `1/2` under the
[published rule](https://github.com/bigospeedrun/site/blob/main/data/integer-factoring.json).

It improves the prior certified `1/30` and the sharper handwritten k6
exponent `log(12/11)/log12`. The assertion compares proved upper-bound
forms. It supplies neither a practical timing claim nor a lower bound
on the runtime of an older implementation.

## 8. Checks, limitations, and reproduction

The proof obligations are discharged as follows:

| Obligation | Evidence |
| --- | --- |
| Uniform evaluator, signed arithmetic, all recursion levels | Pinned assembled Appendix B theorem; [source audit](half-evaluator-audit.md) |
| Exact divisor-prefix network | Pinned factoring reference's `mulNet_value` |
| Deterministic grouping and bounded explicit tables | [Top-constructor proof](half-top-constructor.md), including (8)--(9) |
| Codec, orientation, padding, integer value | Bijection proof and finite adapter checks with negative control |
| Computable sufficient pad and complete encoded input | (10)--(12), with the source's fixed encodings |
| One machine, every input, full factor recovery | Sections 2, 5, and 6; finite fallback paid by a fixed constant |
| Exact logarithmic rate and score | (5)--(7) and Section 7 |

Run the finite implementation and arithmetic checks with

```sh
make verify
```

The checkers do not materialize the enormous evaluator or benchmark a
practical factoring implementation. They test finite adapters and arithmetic;
the infinite-input result is proved by the pinned source theorems and the
Lean composition. [HalfTheorem.lean](../formal-appendix/HalfTheorem.lean)
proves the reference `Factoring.SqrtLogLogConjecture`, including the
polynomial-time padded algorithm, its binary input and output, the
finite-alphabet machine, and the discrete square-root-log bound. The final
axiom report contains `propext`, `Classical.choice`, and `Quot.sound`, with
no `sorryAx`. Run `make verify-half-appendix` to rebuild the pinned sources
and audit the formal composition. Independent review is still outstanding.

The existing half-pieces channel counts have a proved asymptotic barrier
at `1/2` within that particular interface; see the evaluator audit.
That is not a lower bound on the ranks of all possible factorizations
or on deterministic factoring. A stronger bound needs a further idea.
