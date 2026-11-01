# Coordinate restrictions of the four-permutation seed

## Result and scope

Let `S` be any nonempty subset of the 24 permutations in `S_4`,
and write `d=|S|`. Restrict all four label legs of the identity-kernel
[split-permutation coefficient family](split-permutation-seed.md)
to these coordinates. Let `V_S` be the resulting coefficient span
over any field of characteristic two.

For every exact decomposition

\[
E_S=\sum_j B_j,\qquad B_j\in V_S,
\]

the strict finite rank-sum condition fails:

\[
\boxed{\qquad
\sum_j\sqrt{r(B_j)c(B_j)}\ \ge\ d^{\sqrt2}.
\qquad}                                                       \tag{1}
\]

For `d>=2` the inequality obtained below is strict. Thus deleting
coordinate labels cannot make this entire four-permutation
coefficient span pass the
[finite rank-sum amplification criterion](finite-seed-amplification.md).
This allows arbitrary recombination inside the restricted span.

The statement concerns coordinate restriction, with the same selected
permutations on all four legs. It does not apply to general dual
linear embeddings, a changed mode kernel, arbitrary crossing seeds,
or every collective recombination of tensor powers.

## 1. Single-cycle blocks survive exactly when both labels survive

Use the full-cycle graph from
[permutation-span-cycle-obstruction.md](permutation-span-cycle-obstruction.md).
Its vertices are the 24 permutation labels. Join two vertices when
their union is one alternating cycle through all eight bipartite
vertices. Each permutation has six neighbors, corresponding to the
six relative four-cycles. Hence the graph has 72 edges.

An edge uses only its two permutation labels in each cut. Its union
has exactly two alternating perfect matchings, so distinct edges
give distinct matrix blocks. The polynomial blocks are

\[
\begin{pmatrix}P_*&P_*Q_+\\P_*Q_+&P_*\end{pmatrix},
\qquad
\begin{pmatrix}P_*&Q_+\\Q_+&P_*\end{pmatrix},                    \tag{2}
\]

where `P_*=P_1P_2P_3P_4` and `Q_+=prod_j(1+Q_j)`.
The first and second matrices are the row and column cuts,
respectively. These formulas depend only on being a full cycle.

Deleting `24-d` vertices removes at most `6(24-d)` edges.
For `15<=d<=24` there are consequently at least

\[
M=72-6(24-d)=6d-72>0                                         \tag{3}
\]

surviving two-by-two blocks in each cut.

The common diagonal-entry functional `a(B)` from (2) is defined
on the restricted tensor span: choose any surviving block and read
its diagonal entry. It is independent of the coefficient-functional
representation of `B`, and it has the same value on every surviving
block. The crossing satisfies `a(E_S)=1`. In a decomposition of
the crossing, some summand therefore has `a(B_j)!=0`. Every
surviving block of that summand has rank at least one in both cuts.
The parity support still separates these blocks after restriction,
so the same distinguished summand obeys

\[
r(B_j)\ge M,\qquad c(B_j)\ge M.                               \tag{4}
\]

No assumption is made about which subset of a given cardinality
was retained.

## 2. A general rank budget with one distinguished piece

Here is a useful scalar consequence of (4). Suppose positive rank
pairs have both total ranks at least `N`, and one distinguished
pair has both ranks at least `M`, with `1<=M<=N-1`.
Then their cost is at least

\[
\Phi(N,M)=\min\left\{
N,\quad
\sqrt{M(N-1)}+\sqrt{N-M},\
M+2\sqrt{N-M-1}
\right\}.                                                     \tag{5}
\]

To prove this, let `k` be the number of pairs. If
`M+k-1>N`, the distinguished cost is at least `M` and every
other cost is at least one, giving a total greater than `N`.
Otherwise decrease the row and column totals to exactly `N`,
preserving their respective lower bounds `M,1,...,1`.

Minimize by concavity first in all row ranks and then in all column
ranks. Each surplus can be placed at a single coordinate. If the
surpluses land on the same coordinate, the cost is `N`.
If one lands on the distinguished coordinate and the other elsewhere,
the cost is

\[
\sqrt{M(N-k+1)}+\sqrt{N-M-k+2}+k-2,\qquad k\ge2.               \tag{6}
\]

If they land on two different nondistinguished coordinates, it is

\[
M+2\sqrt{N-M-k+2}+k-3,\qquad k\ge3.                           \tag{7}
\]

Both expressions are nondecreasing on their allowed ranges
`k<=N-M+1`. For (6), the derivative is

\[
1-\frac{\sqrt M}{2\sqrt{N-k+1}}
 -\frac1{2\sqrt{N-M-k+2}}\ge0,
\]

because the denominators are at least `sqrt M` and one.
For (7), the derivative is
`1-1/sqrt(N-M-k+2)>=0`.
Taking `k=2` and `k=3`, respectively, gives (5).
This is a real-rank relaxation, so it applies to integer matrix
ranks without further assumptions on the number of pieces.
At `M=N-1`, the `k>=3` case is empty; including its formal
candidate `N-1` in (5) merely weakens the valid bound.

## 3. The exact finite inequalities

The crossing restricted to `d` labels has rank `N=d^2` in both
cuts. Equations (3)--(5) therefore give

\[
\sum_j\sqrt{r(B_j)c(B_j)}
\ge \Phi(d^2,6d-72)>d^{\sqrt2},
\qquad 15\le d\le24.                                         \tag{8}
\]

The following entries are downward-rounded rational lower bounds
from the exact checker:

| Retained labels `d` | Forced blocks `M` | Lower bound for the rank sum |
|---:|---:|---:|
| 15 | 18 | 46.7054000 |
| 16 | 24 | 54.3973682 |
| 17 | 30 | 62.1247568 |
| 18 | 36 | 69.8821486 |
| 19 | 42 | 77.6651090 |
| 20 | 48 | 85.4699878 |
| 21 | 54 | 93.2937654 |
| 22 | 60 | 101.1339276 |
| 23 | 66 | 108.9883704 |
| 24 | 72 | 116.8553228 |

For each integer `d`, the checker replaces every square root
`sqrt z` in (5) by `isqrt(10^14 z)/10^7` and calls the resulting
rational lower bound `L_d`. It verifies

\[
L_d^{408}>d^{577},\qquad 577^2-2\cdot408^2=1.
\]

These integer inequalities prove the strict comparison in (8);
no floating-point decision is used.

For `2<=d<=14`, the universal
[minimum-cut rank-budget exclusion](minimum-cut-rank-budget.md)
already excludes the finite criterion for every crossing
decomposition over every field. At `d=1`, the scalar crossing has
at least one nonzero summand, so the cost is at least one.
This proves (1) for all nonempty coordinate subsets.

At the full 24-label domain, (5) also strengthens the earlier
rank-sum lower bound 96 to more than 116.8553. It does not determine
the optimal rank sum in the full coefficient span.

## Reproduction

The checker constructs the complete 24-vertex cycle graph, verifies
all 576 scalar entries of its two polynomial blocks, checks its
degree/deletion accounting, and verifies the ten exact comparisons
in (8). The universal small-label statement has its own proof and
checker in the linked minimum-cut note.

```sh
python tools/check_permutation_coordinate_restrictions.py \
  --output evidence/permutation-coordinate-restrictions.json
```
