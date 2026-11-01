# Complete recombination screens for the weight-four Schur pieces

## Result

Over `F_2`, neither complete weight-four Schur coefficient span gives
an improving finite crossing seed, even after arbitrary linear
recombination and allowing any number of pieces. This is a finite
restriction result for the two supplied spans. It is not an obstruction
to larger weights, auxiliary matrix constructions, polynomial Smith
compression, or factoring algorithms in general.

The integral construction is proved in
[integral-schur-crossing-pr30.md](integral-schur-crossing-pr30.md). Here `h=4`
and the external dimension is `q=p(4)=5`. The two families are the
inhomogeneous source duals `T(s_lambda)` and the full homogeneous
dual family at the prime-power weight four.

| Family | Dimension of coefficient span | Tensors checked | Minimum total rank cost interval |
| --- | ---: | ---: | ---: |
| Complete inhomogeneous Schur | 14 | 16,384 | `[20.491933,20.491934]` |
| Complete homogeneous Schur | 8 | 256 | `[20.491933,20.491934]` |

For a scalar tensor `A`, the cost is
`w(A)=sqrt(rank_row(A) rank_column(A))`. The optimization asks for the
smallest `sum_j w(A_j)` subject to `sum_j A_j=E` and every `A_j` in
the specified span. The number of summands is unrestricted.

An explicit upper-bound witness in both spans is the usual trace
split. Up to exchanging the two cuts, its rank pairs are `(1,25)`
and `(24,10)`, so its cost is

```math
5+4\sqrt{15}=20.491933384829\ldots.
```

The rigorous lower bound is already larger than `q^(3/2)`, and hence
larger than `q^sqrt(2)`, the rank-sum threshold for beating the current
half exponent through finite-seed amplification. The substantial
reduction from each family's raw shift sum does not supply a new seed.

## Why the optimization covers every decomposition

Choose a binary basis of the complete coefficient span. A tensor has
a unique mask in `F_2^m`, and addition is XOR. Build the complete
weighted graph on these masks, with edge cost `w(u XOR v)` between
vertices `u` and `v`. Every decomposition of the target tensor is a
path from zero to the target mask. Conversely, the edge differences
of such a path give a decomposition.

All nonzero tensors have positive cost. Removing repeated vertices
from a path cannot increase its cost, so the finite shortest-path
problem covers decompositions of arbitrary length.

The native checker computes both matrix ranks exactly over `F_2`
for every tensor in the span. It replaces each edge cost by the
integer

```math
\left\lfloor 10^6\sqrt{r c}\right\rfloor.
```

The square root is rounded using integer comparisons. Complete
Dijkstra search then gives an integer shortest-path distance.
Dividing this distance by `10^6` is a rigorous lower bound for the
original radical-cost optimization. The recovered path is checked
separately with each edge rounded upwards, giving the stated upper
bound. The search does not rely on floating-point comparisons for
its exclusion.

## Replay

The source files are `schur_span_recombination.cpp` and
`check_schur_span_recombination.py`. The Python replay independently
builds the complete spans, verifies the endpoint, checks the rank
pairs and sum of the returned decomposition, and records the hashes
of the input matrices and both checkers. Its source matrices are the
integrally normalized weight-four crossings, before any recombination.

```sh
g++ -O3 -std=c++20 tools/schur_span_recombination.cpp \
  -o /tmp/schur_span_recombination
python3 tools/check_schur_span_recombination.py \
  --native /tmp/schur_span_recombination \
  --pieces evidence/homogeneous-h4-pieces.json evidence/full-h4-pieces.json \
  --output evidence/full-schur-span-recombination.json
```

Every assertion in this note concerns these finite spans. Polynomial
coefficient amplification can behave differently from a finite scalar
decomposition, so the result does not replace the Smith-profile or
auxiliary-path searches.
