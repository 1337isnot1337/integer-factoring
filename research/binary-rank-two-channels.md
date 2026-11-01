# Complete binary four-label normalized rank-two search

## Result and exact scope

Consider every channel over `F_2` of the normalized form

\[
 U(X)=P\left[X^{\mathsf T}(I+AB)+AX^{\mathsf T}B\right]P^{-1},
 \qquad A,B\in\mathrm{Mat}_4(F_2),\quad AB=BA,
 \quad P\in\mathrm{GL}_4(F_2).                    \tag{1}
\]

For every such channel for which `S=I+U` is regular nilpotent, the
existing resolvent Smith-profile method fails to improve the factoring
bound. An exhaustive exact search establishes a stronger finite fact:
every rank-two case has at least eleven row-Smith orders equal to one,
in addition to its two orders equal to zero. Consequently its profile
cost is at least 197 at every scale, and its normalized contraction
proxy is greater than `19/20`.

This covers all commuting matrix modules in (1), including noncyclic
ones. It does not cover a normalization requiring a coefficient-field
extension, nonunital or non-trace-preserving channels, larger channel
ranks, or other coefficient decompositions of the same tensors. The
[all-field normal-form theorem](rank-two-channel-normal-form.md)
can require a scalar extension; that possibility is deliberately not
identified with the finite binary space searched here.

The result is an exact finite computer-assisted classification with
ordinary proofs of its reductions. It is not a new factoring exponent
or a Lean formalization.

The separate [determinant-budget theorem](regular-resolvent-determinant-budget.md)
now excludes every regular rank-two source over every field and at every
alphabet size. The finite census here remains stronger about the actual
first-order pivots in its declared binary normalized space: it proves
eleven order-one pivots and the cost 197, rather than only failure of the
half-score threshold.

## 1. Reducing A by similarity preserves the whole search

The channel fixes the identity and preserves trace. Its crossing
reshuffling has rank at most two. To reduce its parameters, change
the physical matrix basis by `X -> CXC^-1`. Conjugating `U` by this
operator replaces

\[
 A\longmapsto C^{\mathsf T}AC^{-\mathsf T},
 \quad B\longmapsto C^{\mathsf T}BC^{-\mathsf T},
 \quad P\longmapsto C^{-1}PC^{-\mathsf T}.         \tag{2}
\]

Substitution into (1) proves (2) directly. The transformed `P` ranges
over the whole general linear group, and transformed `B` ranges over
the centralizer of transformed `A`. This physical change preserves
nilpotent Jordan type and both crossing ranks and Smith profiles.
Thus it suffices to choose one representative of each similarity
class of `A`, then every commuting `B` and every invertible `P`.

The native search constructs these representatives without assuming
a classification table. It scans all 65,536 binary matrices in order.
For each first unseen matrix, it marks its complete orbit under all
20,160 invertible matrices. The result has 34 orbits. The independent
checker reconstructs every orbit, checks that they are disjoint, and
verifies that their union is all 65,536 matrices.

For each representative, the native search tests every binary `B`
for `AB=BA`. The checker instead constructs the kernel of the linear
commutator map `B -> AB+BA`, enumerates that kernel from a basis, and
checks its size. These are different enumeration methods for the same
centralizer.

## 2. Exact finite census

There are 136,960 pairs `(A,B)` after choosing the 34 representatives
of `A`. The two scalar representatives `A=0,I` each contribute 65,536
pairs, but (1) then reduces to a single transpose channel. The
remaining 5,888 pairs contain 58 further rank-one channels and 5,830
distinct rank-two base channels. Each is combined with all 20,160
invertible `P`, for

\[
 5830\cdot20160=117532800                         \tag{3}
\]

exact rank-two candidate tests.

Regular nilpotence is tested by

\[
 S^{16}=0,\qquad S^{15}\ne0.                     \tag{4}
\]

In dimension sixteen, (4) is equivalent to a single nilpotent Jordan
block: a vector on which `S^15` is nonzero supplies a chain of length
sixteen. The native program uses repeated squaring; the independent
checker recomputes the successive powers of each retained case.
All arithmetic is binary integer arithmetic.

Exactly 5,888 parameter cases pass (4). These are counted as parameter
cases, without claiming that all their resulting operators are distinct.
Their numbers of order-one row pivots are:

| Number of order-one pivots | Regular parameter cases |
|---:|---:|
| 11 | 480 |
| 12 | 3,792 |
| 13 | 784 |
| 14 | 832 |

Every retained case has row reshuffling rank exactly two.

## 3. The first two coefficients certify the failure

For the [nilpotent resolvent crossing](nilpotent-flag-crossings.md),
the reversed row polynomial begins

\[
 H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr)
     =H_0+wH_1+O(w^2),
 \quad H_0=\Phi(U),\quad H_1=\Phi(U+U^2).         \tag{5}
\]

If its Smith orders are `mu_i`, then

\[
 \mathrm{rank}\begin{pmatrix}H_0&0\\H_1&H_0\end{pmatrix}
 =\sum_i\max(0,2-\mu_i)
 =2\#\{i:\mu_i=0\}+\#\{i:\mu_i=1\}.            \tag{6}
\]

This is multiplication by `H` on truncated length-two power series;
invertible Smith row and column operations remain invertible after
that truncation. Thus an exact binary rank calculation on a 32 by 32
matrix determines the number of order-one pivots. The native search
and the independent checker use separate rank implementations.

The column orders for (4) are `(0,16,...,16)`, with fifteen copies
of sixteen. The sufficient profile cost at `s>=1` is therefore

\[
 B(s)=\left(\sum_i s^{15-\mu_i}\right)(1+15s^{-16}).
 \tag{7}
\]

Each zero order contributes `s^15+15/s>=16`. Each order one contributes
`s^14+15/s^2>=15`, by the elementary scalar inequality proved in
the [partial-transpose rank bound](partial-transpose-power-rank.md).
The census and (6) now imply

\[
 \boxed{B(s)\ge2\cdot16+11\cdot15=197.}          \tag{8}
\]

In particular, `197>64>16^sqrt(2)`, so the desired sufficient threshold
cannot hold. More quantitatively, `197^5>4^19` gives the claimed
proxy bound `log(B)/(4 log 4)>19/20`.

The omitted rank-one cases also fail: the existing general power-rank
bound gives `B>=16+15*14=226` for a regular source. Therefore omitting
them from the 117,532,800 tests does not leave a gap in the stated
failure result for (1).

## 4. A retained minimum witness and full profile

One case attaining the minimum eleven order-one pivots is encoded
in row-major binary matrix masks by

\[
 A=72,\qquad B=2068,\qquad P=4788.               \tag{9}
\]

The independent checker constructs each entry of (1) from the two
transpose-channel factors, verifies the endpoint normalizations and
(4), and constructs the entire degree-fifteen crossing polynomial.
Its complete row orders at infinity are

\[
 (0,0,\underbrace{1,\ldots,1}_{11},2,2,9).        \tag{10}
\]

The full series computation uses precision 257. This is sufficient
to resolve all orders because a full-rank 16 by 16 polynomial matrix
of degree at most fifteen has determinant degree at most 240. Its
crossing endpoint proves full rank. The retained profile contains
no unresolved cutoff values.

## 5. Reproduction and certificate scope

The [native enumerator](../tools/search_cyclic_rank_two_channels.cpp)
has a smaller default cyclic-matrix mode. Use `--all-commuting` for
the complete finite space in this note. The
[independent checker](../tools/check_binary_rank_two_channels.py)
verifies the orbit partition, all centralizers, every retained regular
case, its first two coefficient ranks, and the complete minimum witness.
The [receipt](../evidence/binary-rank-two-channels.json) retains all
5,888 regular cases, the 34 class sizes, both census summaries, and
source hashes.

```sh
python3 tools/check_binary_rank_two_channels.py
python3 tools/check_binary_rank_two_channels.py --replay-search
```

The second command additionally compiles and reruns all 117,532,800
native candidate tests, then compares the complete ordered census
with the retained receipt. A witness-only replay does not reclassify
rejected candidates; the command and output distinguish these scopes.
The complete native run was executed when producing this checkpoint.
