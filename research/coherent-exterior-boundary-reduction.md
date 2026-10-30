# Coherent exterior boundaries and a common marginal projection

## Results and limits

The exterior-monodromy family admits arbitrary normalized auxiliary
boundary functionals. This note proves a common projection that reduces
their parameter space without worsening either cut rank by more than a
polynomial factor. For binary letters and the two-particle auxiliary
space, the reduction leaves eight parameters instead of thirty-six.

It also strengthens the fixed-occupation rank-sum certificate by removing
its factor `1/(t+1)`. Consequently every binary fixed occupation is
excluded from the improving finite-seed criterion at **every** positive
width, not merely asymptotically.

The remaining coherent family is not excluded asymptotically. An exact
finite screen covers all 128 normalized reduced boundaries over `F_2`
at four positions, including arbitrary binary recombinations of their
coefficients. Each has minimum rank sum exactly `256`, the cost of the
crossing itself. This finite result does not extend to larger widths or
extension-field coefficients without further proof.

The construction and the fixed-boundary triangular restrictions come
from [exterior monodromy](https://github.com/1337isnot1337/integer-factoring/blob/bb3bf763ed054f6fdaf882649c037ad810e16b97/research/exterior-monodromy-and-occupation-graphs.md).
The common projection and the union-of-levels refinement below are new
arguments. The code extends the repository's exact exclusion walk and
packed elimination, with an independent exterior-operator replay.
Substantial OpenAI assistant support was used. No new factoring exponent
or Lean result is asserted.

## 1. The full boundary interface

Fix `d>=2` and `0<=r<=d^2`, and work in characteristic two. Let
`D=d^2`. The auxiliary basis consists of `r`-element occupations of
the pair alphabet `[d] times [d]`. For `p=(x,y)`, `q=(v,u)`, put

```math
K(z;p,q)=[p=q]I+(1+z)c_q a_p.
```

For an arbitrary linear functional `ell` on the auxiliary matrix algebra,
define

```math
P_{t,\ell}(z;\mathbf x,\mathbf y,\mathbf u,\mathbf v)
 =\ell\!\left(K(z;p_1,q_1)\cdots K(z;p_t,q_t)\right).
```

The order of the positions is retained. If `ell(I)=1`, then

```math
P_{t,\ell}(1)=E_{d^t}.
```

Thus every normalized boundary gives a polynomial crossing with at most
`t+1` coefficients. Fixed `d,r,ell` permit polynomial entry construction
and `d^{O(t)}` dense factor construction, as in the preceding note.

Write `w_(S,T)=ell(E_(S,T))`. An entry of the row-vector auxiliary walk
begins at `S`, ends at `T`, and is multiplied by `w_(S,T)`. Diagonal
steps retain the occupation. An off-diagonal step replaces occupied
`q` by empty `p`. Therefore every surviving path satisfies

```math
\mathbf 1_T-\mathbf 1_S=\sum_{i=1}^t(e_{p_i}-e_{q_i}).\tag{1}
```

Equation (1) is an equality of integer occupation counts. It remains
valid for the support of the matrices over characteristic two.

For an occupation `S`, let `mu(S)` be its two external-letter marginals:
the `d` row counts and the `d` column counts of its occupied pairs.
Define the reduced boundary

```math
w^0_{S,T}=[\mu(S)=\mu(T)]w_{S,T}.\tag{2}
```

All diagonal entries are retained, so `ell^0(I)=ell(I)`.

## 2. One projection controls both cuts

Let `type(x)` be the histogram of letters in an external word `x`.
Define the linear tensor projection

```math
(\Pi F)(x,y,u,v)
 =[\operatorname{type}(x)=\operatorname{type}(v)]
  [\operatorname{type}(y)=\operatorname{type}(u)]F(x,y,u,v).\tag{3}
```

Taking the two marginals in (1) proves the exact polynomial identity

```math
\boxed{\Pi P_{t,\ell}=P_{t,\ell^0}.}\tag{4}
```

Indeed the conditions in (3) hold exactly when the two auxiliary endpoint
occupations have equal marginals. Every path with such endpoints already
satisfies (3); every path with other endpoints is removed.

There are

```math
h_t=\binom{t+d-1}{d-1}
```

possible histograms of a length-`t` word. In the row cut, (3) extracts
the matching blocks indexed by `(type(x),type(y))` on one side and
`(type(v),type(u))` on the other. In the column cut it extracts the
matching blocks indexed by `(type(x),type(u))` and
`(type(v),type(y))`. In either arrangement it is a sum of `h_t^2`
two-sided coordinate projections. Rank subadditivity gives

```math
r(\Pi F)\le h_t^2 r(F),\qquad c(\Pi F)\le h_t^2 c(F).\tag{5}
```

This is not a claim that block-diagonal extraction never increases rank.
The polynomial factor in (5) pays for that possible increase.

The crossing is fixed by `Pi`. If `E=sum_j L_j(P_(t,ell))` is any
decomposition in the coefficient span, then the same coefficient
functionals give a decomposition in the reduced boundary's span. Equations
(4)--(5) imply

```math
\sum_j\sqrt{r(\Pi F_j)c(\Pi F_j)}
 \le h_t^2\sum_j\sqrt{r(F_j)c(F_j)}.\tag{6}
```

For fixed `d,r`, `log h_t=O(log t)`. Therefore the infimum leading
rank-sum exponent among **all** normalized fixed boundaries is the same
as the infimum among reduced boundaries (2). One inequality follows
from (6), and the other follows because reduced boundaries are a
subfamily. The reduction also applies input by input if the boundary
coefficients change with `t`, provided the original effectiveness and
field-cost hypotheses are supplied.

It does not preserve a strict finite-seed inequality without the factor
in (6), and it does not average over an auxiliary group or divide by an
auxiliary dimension.

## 3. Binary two-particle boundaries have eight parameters

For `d=r=2`, encode an occupied pair `(x,y)` by `x+2y`. In increasing
occupation-mask order, the six states are

```text
3, 5, 6, 9, 10, 12.
```

The four rectangular occupations `3,5,10,12` have distinct marginals.
The two perfect-matching occupations `6,9` have the same marginals.
Consequently (2) retains six diagonal entries and the two directed
off-diagonal entries `6->9` and `9->6`. Every other boundary entry is
removed.

The reduced boundary algebra is therefore the vector-space direct sum
of four scalar blocks and one `2 by 2` block. Over `F_2`, the condition
`ell(I)=1` is one nonzero linear equation in eight coordinates, leaving
exactly `2^7=128` boundaries. Coherence in the matching block survives
the projection; deleting it would be an additional unsupported step.

## 4. A stronger fixed-occupation rank-sum inequality

Keep a fixed occupation `S` of size `r`, with `0<r<D`. Define its
column occupation graph on the pair alphabet as follows: for
`a=(x,u)` and `b=(v,y)`, there is an edge `b->a` precisely when
`a!=b`, `(v,u)` is occupied in `S`, and `(x,y)` is unoccupied.
Let `T` be a vertex subset on which this graph is acyclic, with
`n=|T|` and `s=|T intersect S|`. For a coefficient functional `L`, put

```math
I(L)=\{k:0\le k\le t,\ L(z^k)\ne0\}.
```

The row order puts every unoccupied letter before every occupied letter;
the column order is one fixed topological order of the occupation graph
on `T`. These orders are chosen before `k` or `L` and do not depend on
the diagonal level. The first-differing-position argument therefore makes
the entire polynomial restriction triangular in each cut. More explicitly,
all earlier equal pair letters have left the initial occupation unchanged.
A first unequal row-cut pair can act only from an occupied to an
unoccupied letter, which gives the chosen triangular order. A first
unequal column-cut pair can act only along an edge of the displayed
graph, which gives its topological order. On a diagonal word with `k`
occupied letters the polynomial entry is `z^k`. Retaining
**all** nonzero diagonal levels at once takes one principal submatrix of
that same triangular matrix, so it remains triangular with nonzero
diagonal. Hence

```math
r(L(P))\ge\sum_{k\in I(L)}a_k,\qquad
 c(L(P))\ge\sum_{k\in I(L)}b_k.\tag{7}
```

where

```math
a_k=\binom tk r^k(D-r)^{t-k},\qquad
 b_k=\binom tk s^k(n-s)^{t-k}.
```

Cauchy--Schwarz, applied only to the selected levels, gives

```math
\sqrt{r(L(P))c(L(P))}
 \ge\sum_{k\in I(L)}\sqrt{a_k b_k}.\tag{8}
```

Suppose arbitrary coefficient-span tensors `L_j(P)` sum to the crossing.
For each `k`, a row diagonal word with exactly `k` occupied letters
forces `sum_j L_j(z^k)=1`. At least one `I(L_j)` therefore contains
that `k`. Summing (8) over pieces and then summing the binomial expansion
proves the exact finite certificate

```math
\boxed{\mathcal S_t\ge\left(\sqrt{rs}+\sqrt{(D-r)(n-s)}\right)^t.}\tag{9}
```

There is no `1/(t+1)` loss. The argument permits arbitrarily many
pieces, arbitrary coefficient recombinations, and every
characteristic-two extension field.

For rectangular occupations or their complements the full occupation
graph is acyclic, so (9) gives `S_t>=D^t`. For the two diagonal binary
occupations, the three-letter acyclic restrictions give

```math
\mathcal S_t\ge(2+\sqrt2)^t.
```

Thus every binary fixed occupation fails `S_t<2^(sqrt(2)t)` at every
positive width. The 126 ternary half-occupation witnesses already in
the repository similarly strengthen to

```math
\mathcal S_t\ge(\sqrt{12}+\sqrt{20})^t.
```

These remain restrictions on fixed occupation boundaries. Equation (9)
does not assert the same inequality for a coherent boundary.

## 5. Complete finite binary screen

At `t=4`, every normalized reduced boundary was tested against every
nonzero binary coefficient functional. There are `128*31=3968`
rank pairs. Native elimination uses the exact histogram blocks; the
independent checker forms ordinary exterior creation/annihilation matrices
and full small cut matrices.

The optimum is certified without floating-point comparisons. For each
boundary, construct a graph on the `32` coefficient-functional masks.
An edge from `a` to `a xor b` has weight

```math
\left\lfloor2^{16}\sqrt{r_b c_b}\right\rfloor.
```

These integer weights are lower bounds on the true rank cost multiplied
by `2^16`. Exact Dijkstra distances from zero to the all-ones functional
are `256*2^16` for all 128 boundaries. Any binary coefficient
decomposition traces such a path. If two coefficient masks represent
the same tensor, their difference has zero rank and zero edge cost, so
coefficient dependencies are included rather than assumed absent.

The crossing itself has ranks `(256,256)` and cost `256`. It supplies
the matching upper bound. Thus the minimum binary coefficient-span
rank sum is exactly `256` for every screened boundary.

```sh
python3 tools/check_coherent_exterior_boundaries.py \
  --output evidence/coherent-boundary-verification.json \
  --screen-output evidence/coherent-t4-all-span.json
```

The native source is `tools/coherent_exterior.cpp`. The output retains
every boundary, functional, and exact rank pair. The checker also verifies
the common projection on arbitrary full boundaries, exercises off-diagonal
terms removed by the projection and matching coherence that survives it,
and checks the simultaneous triangular-level inequalities (7).

The unresolved problem is a useful large-width rank bound or a competitive
finite certificate in the remaining coherent family. The reduction makes
that search smaller; the four-position screen does not settle it.
