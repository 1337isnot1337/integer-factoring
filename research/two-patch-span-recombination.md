# Exact recombination controls for two small trace spans

The [two-patch trace crossing](two-patch-trace-crossing.md) has degree
`2t-1` and an exact crossing endpoint. The raw coefficient rank sums fail
the finite target. The following complete finite calculation proves a
stronger statement at two widths: **every binary recombination and every
number of pieces in the same coefficient span has rank cost at least
`d^2`**. The single crossing tensor attains that cost.

This is a finite statement for the specified field models and trace bases.
It is not an all-size obstruction, and it does not cover recombinations
with nonbinary coefficients in a larger coefficient field.

## 1. The exact finite optimization

Write the polynomial tensor as

```math
T_t(z)=\sum_{j=0}^{2t-1}z^j P_j,\qquad
E=T_t(1)=\sum_j P_j.
```

For a binary coefficient mask `a`, let `Q_a=sum_j a_j P_j`. The tested
coefficient tensors are independent, so the masks identify the span
without collisions. The target mask is the all-one mask.

The optimization includes every finite decomposition

```math
E=\sum_i Q_{a_i},\qquad
S=\sum_i\sqrt{r(Q_{a_i})c(Q_{a_i})}.
```

Zero pieces can be removed. Every such decomposition is a path from the
zero mask to the target mask in the binary additive group: adding a piece
changes the current mask by XOR with its coefficient mask. Conversely,
every such path specifies a decomposition. Repetition, grouping, and
arbitrary binary basis changes are therefore included.

For each nonzero mask define the integer lower cost

```math
L(a)=\left\lfloor\sqrt{r(Q_a)c(Q_a)}\right\rfloor.
```

The checker computes shortest-path distances using these integer costs.
It then independently checks the complete finite potential certificate

```math
\ell(0)=0,\qquad
\ell(x\mathbin{\oplus}a)\le\ell(x)+L(a)
\quad\text{for every state }x\text{ and nonzero mask }a.
```

Telescoping along any decomposition proves `S>=ell(target)`. This lower
bound requires no comparison or rounding of sums of radicals.

## 2. Exact outcome

| `t` | Field polynomial | Span dimension | Tensors checked | Certified minimum `S` |
| ---: | ---: | ---: | ---: | ---: |
| 3 | `X^3+X+1` | 6 | 64 | 64 |
| 4 | `X^4+X+1` | 8 | 256 | 256 |

The trace coordinates are `Tr(alpha^j z)`, with `alpha=X` and `0<=j<t`.
In both cases the integer lower distance is exactly `d^2`. The target
tensor has both ranks `d^2`, so its one-piece decomposition supplies the
matching upper bound. Thus the displayed minima are exact, not intervals.

The proof includes all 69,312 potential inequalities and all 320 span
tensors, including zero. Every tensor's row and column ranks are computed
by binary convolution blocks and checked independently on the full
four-leg flattenings: 640 dense cut comparisons in total. The evidence
contains every rank pair, every integer lower cost, and every potential
value, along with source hashes.

The raw costs were approximately 165.774 at `t=3` and 865.793 at `t=4`.
Recombination can reduce those costs to the trivial one-piece values, but
cannot improve on those values within these two binary spans. No larger
width or other trace basis is excluded by this calculation.

## 3. A concrete next synthesis target

The next target should change the coefficient space while retaining an
exact endpoint and only linearly many coefficients. One inexpensive candidate is to
assign weight two to one trace coordinate. Choose a patch `ell in {0,1}`
and a coordinate `j_0`. Replace its factor

```math
1+z\tau_{j_0}(B_\ell)
\quad\text{by}\quad
1+z^2\tau_{j_0}(B_\ell),
```

leaving every other factor unchanged. Evaluation at `z=1` is unchanged,
so the resulting tensor still equals `E`. Its degree is exactly `2t`:
the impossible all-one trace pattern removes degree `2t+1`, while a
pattern missing only one untagged coordinate realizes degree `2t`.
Such a pattern is realizable because any distinct pair `(B_0,B_1)`
comes from `e=B_0+B_1`, `S=B_0/e`, and `w=S+e`.

Seek an exact binary recombination of these at most `2t+1` coefficients
whose total rank cost satisfies `S<d^sqrt(2)`. This is a proposed target,
not a claimed improvement. It uses the higher Boolean degrees necessary
for a crossing, and it uses no ordinary-kernel multiplication, so neither
the low-degree span barrier nor the full-cube composition cancellation
excludes it.

A bounded next experiment at `t=4` has at most 512 span tensors for each
choice of tagged coordinate. The same complete XOR-potential calculation
can test whether the optimum drops strictly below `d^2=256`. That is a
diagnostic for useful new cancellation, not a claim that the factoring
threshold has been reached. Any larger-width or asymptotic conclusion
would still need its own identity, exact ranks, and resource accounting.

## Reproduction

[Exact checker](../tools/check_two_patch_span.py) and
[complete finite certificate](../evidence/two-patch-span-verification.json).

```sh
python3 tools/check_two_patch_span.py --output /tmp/two-patch-span-verification.json
```

The checker uses the Python standard library and the existing binary-field
implementation. It regenerates all tensors, both rank calculations, and
the complete lower certificates. The optional rational refinement in the
code is unnecessary for these two instances: integer floors already
certify both exact optima.
