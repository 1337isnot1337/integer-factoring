# Component-path refinement of channel-walk bounds

## Constructive result

The [channel-walk transfer](twisted-channel-walk-transfer.md) bounds
the rank sum using one set of reciprocal piece weights for the two
channel graphs. A single weighting can lose information about which
strongly connected components one walk can actually visit.

This note removes part of that loss. Cover each graph by subgraphs on
maximal paths in its condensation directed acyclic graph. Use a separate
reciprocal weighting and positive-vector certificate for **each pair**
of a row path and a column path. The largest of their certified bases
is a valid amplification base. The number of path pairs is fixed by
the finite seed, so it contributes only a constant factor.

The actual tensor pieces remain the original piece-type sums. This is
a refinement of their rank-count bound; it does not require a new
tensor identity or a simultaneous row/column decomposition by paths.
It also applies to the [terminal-rank automata](scalarized-channel-automata.md)
obtained by scalarizing ordered product spaces.

An exact abstract graph example improves the certified base from two
to one. No improving crossing tensor is supplied by that example.

## 1. A finite cover of all walks

Let `G` be a finite directed graph or multigraph. Contract each
strongly connected component to one vertex. The resulting condensation
graph is acyclic: a cycle of components would itself be strongly
connected.

Every walk visits a directed path of components, with no return to a
component after leaving it. Extend this path to a maximal path in the
finite condensation graph. For each maximal path `pi`, let `G_pi`
be the subgraph induced by the union of its components. Edges that
skip intermediate components may be retained. Then every original
walk lies in at least one `G_pi`.

The cover may count one walk more than once. Only an upper bound is
needed. If `R_nu(k)` counts row channel words of length `k` and piece
type `nu`, and `R_nu^pi(k)` counts them in `G_pi`, then

```math
R_\nu(k)\le\sum_{\pi\in\mathcal P_R}R_\nu^\pi(k).
\tag{1}
```

The same holds in the column graph. Initial and terminal weights are
retained on each subgraph. For a deterministic product-space automaton,
use its original initial state and terminal bilinear ranks. Subgraphs
with zero contribution can be omitted.

## 2. Independent weights for each path pair

The unrefined factor-count bound has the form

```math
S_k\le c_0\sum_{|\nu|=k}\sqrt{R_\nu(k)C_\nu(k)},
\tag{2}
```

where `c_0=b` for the original matrix-entry scalarization, and `c_0=1`
when terminal bilinear ranks have already been included in the counts.
From (1) and `sqrt(sum a_i)<=sum sqrt(a_i)` for nonnegative terms,

```math
S_k\le c_0\sum_{\pi\in\mathcal P_R,\ \tau\in\mathcal P_C}
              \sum_{|\nu|=k}
                \sqrt{R_\nu^\pi(k)C_\nu^\tau(k)}.
\tag{3}
```

For each path pair `(pi,tau)`, choose its own positive piece weights
`z_i^(pi,tau)`. Weighted Cauchy--Schwarz gives

```math
\sum_\nu\sqrt{R_\nu^\pi(k)C_\nu^\tau(k)}
\le\sqrt{P_R^\pi(k,z^{\pi,\tau})
          P_C^\tau(k,(z^{\pi,\tau})^{-1})}.
\tag{4}
```

The counts in different terms of (3) need not use the same weights.
They are separate nonnegative numerical bounds, not different
definitions of the original tensor pieces.

Supply rational positive-vector inequalities for each of these two
restricted transition matrices, as in the channel-walk theorem.
They give a constant `C_(pi,tau)` and a base

```math
L_{\pi,\tau}=\sqrt{\beta_R^{\pi,\tau}
                            \beta_C^{\pi,\tau}}
```

such that the right side of (4) is at most
`C_(pi,tau) L_(pi,tau)^k`. In the vertex-word convention the power
is initially `k-1`; positivity of the bases absorbs that one-step
shift into the constant. Therefore

```math
\boxed{S_k\le C L_*^k,\qquad
       L_*=\max_{\pi,\tau}L_{\pi,\tau},\qquad
       C=c_0\sum_{\pi,\tau}C_{\pi,\tau}.}
\tag{5}
```

All path families, rational data, and constants are finite parts of
the seed. The original polynomial number of type pieces, effective
entry queries, and `2^O(lambda)` table construction remain unchanged.
Thus `L_*<d^sqrt(2)` is another sufficient certificate for the strict
finite-seed target after sufficiently many copies.

## 3. Comparison with the single-weight bound

At the optimization level the refined sufficient base is

```math
\max_{\pi,\tau}\inf_{z_i>0}
 \sqrt{\operatorname{spr}T_R^\pi(z)
       \operatorname{spr}T_C^\tau(z^{-1})}.
\tag{6}
```

As before, a strict upper bound is implemented by positive rational
vectors, not by an unverified eigenvalue computation. Since the
number of pairs is finite, strict rational certificates can be chosen
for all pairs at once.

The refined infimum is never worse than the original one. Every
restricted transition matrix is entrywise bounded by the full matrix
after zero padding. A positive-vector certificate for the full graph
therefore restricts to each subgraph. Using the same piece weights
and bounds in every pair recovers the old base; independent choices
can only improve the infimum.

The construction can be expensive as a function of the graph size:
a directed acyclic graph can have exponentially many maximal paths.
Here the graphs are fixed finite certificate data. This is not a claim
that the refinement is polynomial time in arbitrary input graph size.

## 4. A strict and exactly verified improvement

Use piece labels `A,B,C`. The row graph has three isolated vertices:
two parallel `A` loops at the first, two parallel `B` loops at the
second, and one `C` loop at the third. The column graph has a
two-cycle with two parallel `A` edges in one direction and two
parallel `B` edges in the other, plus an isolated `C` loop.

### The original optimized base is two

Give count one to each of the four `A/B` edges in both graphs,
and zero to the `C` loops. The row and column flows are stationary,
have common total `D=4`, and have piece totals `(2,2,0)`. Each
flow has the integer ratio `K=2^4` from the
[flow dual certificate](channel-flow-dual-certificates.md).
Consequently every original positive-vector certificate satisfies
`L^8>=2^8`, so `L>=2`.

Conversely all piece weights and all positive-vector entries equal
to one give `beta_R=beta_C=2`. Hence two is the exact optimized
base of the original method.

### The component-path base is one

There are three row paths and two column paths. Their six pairings
admit the following rational certificates. Entries in the `z` column
are ordered `(z_A,z_B,z_C)`. The column vector is restricted to the
displayed column component; unused entries can be set to one.

| Row component | Column component | `z` | Row beta | Column beta | Column positive vector |
| --- | --- | --- | ---: | ---: | --- |
| `A` | `AB` cycle | `(1/4,4,1)` | `1/2` | `2` | `(4,1)` |
| `B` | `AB` cycle | `(4,1/4,1)` | `1/2` | `2` | `(1,4)` |
| `C` | `AB` cycle | `(4,4,1)` | `1` | `1/2` | `(1,1)` |
| `A` | `C` loop | `(1/4,1,1)` | `1/2` | `1` | `(1)` |
| `B` | `C` loop | `(1,1/4,1)` | `1/2` | `1` | `(1)` |
| `C` | `C` loop | `(1,1,1)` | `1` | `1` | `(1)` |

The row positive vector is one in every case. Every product of
the two betas is at most one, proving `L_*<=1`. The `C/C` pair
has a matching stationary flow with ratio one, so its optimum is
at least one. Thus the refined optimum is exactly one.

The exact common-type count also equals one for every walk length
at least two edges. The row `A/B` components give monochromatic
words, while the column cycle gives words containing both labels.
Only the all-`C` type occurs in both graphs. This verifies the
mechanism of the improvement independently of the spectral bounds.

These are abstract counting graphs, not a crossing certificate for
label size at least two. A genuine factoring application must still
supply the matrix-valued crossing and all factor/zero identities
from which its graphs arise.

## Verification

```sh
python3 tools/check_component_channel_walks.py
```

The checker derives the strongly connected components and maximal
path covers, checks all six rational certificates, verifies the two
integer dual bounds, and enumerates the matching types. A connected
two-component control prevents an invalid refinement that splits a
walk into separate components even when it can travel between them.
Evidence is in `evidence/component-channel-walks.json`.
