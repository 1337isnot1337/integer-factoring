# Disjointness minors from paths and cycles

## Result and scope

Let `A,B` be partial matchings with the same number `k` of edges. Their
supports need not agree. Over `F_2`, the determinant

```math
\det D(A,B),\qquad D(A,B)_{e,f}=[e\cap f=\varnothing],
```

is computed with `O(k)` graph operations. The algorithm reads the
components of an intersection graph. It does not build a dense matrix.
This applies to the disjointness minors in the homogeneous matching
construction, including minors whose two remaining vertex sets differ.

This is a query algorithm. It does not bound the rank of a crossing
flattening. The distinction matters: the full even-particle matching
Gram matrix has factorial dimension and is invertible, despite each
entry having this short computation.

## 1. The intersection graph

Make a bipartite graph whose left vertices are the edges of `A` and
whose right vertices are the edges of `B`. Join two vertices when the
corresponding physical edges meet. This is a simple graph: a common
physical edge gives one adjacency, although its endpoints meet twice.
Each graph vertex has degree at most two. Every component is a path,
an isolated vertex, or an even cycle.

Let `C` be its bipartite adjacency matrix, so `D=J+C` over `F_2`, where
`J` is the all-ones matrix. Classify the components as follows.

| Component | Left and right sizes | Contribution |
|---|---|---|
| Balanced path | `r,r` | Invertible block, with inverse-entry sum `ceil(r/2)` modulo two |
| Left-major path | `r+1,r` | One left null vector, all ones on its left vertices |
| Right-major path | `r,r+1` | One right null vector, all ones on its right vertices |
| Cycle | `r,r`, `r>=2` | Nullity one, with left and right null vectors all ones |

Isolated vertices are the unbalanced paths with `r=0`.
A common physical edge is a balanced path with `r=1`.

There are equally many left-major and right-major paths. Write that
number as `u`, and write `v` for the number of cycle components. Then

```math
\operatorname{nullity} C=u+v.
\tag{1}
```

This follows by summing the ranks of the component blocks.

## 2. The complete determinant formula

If `u+v>=2`, the answer is zero: adding the rank-one matrix `J` cannot
repair a nullity of at least two.

If `v=1,u=0`, let `r` be the number of vertices on either side of the
unique cycle. Then

```math
\det(J+C)=r\pmod2.
\tag{2}
```

If `u=1,v=0`, let `a` be the left size of the left-major path and let
`b` be the right size of the right-major path. Then

```math
\det(J+C)=ab\pmod2.
\tag{3}
```

If `u=v=0`, every component is a balanced path. If their sizes on one
side are `r_1,...,r_t`, then

```math
\det(J+C)=1+\sum_{j=1}^t\left\lceil\frac{r_j}{2}\right\rceil
\pmod2.
\tag{4}
```

The empty determinant has value one, as included in (4).

To verify these formulas, order a balanced path from one end. Its block
is a unit upper or lower bidiagonal matrix. Its inverse has
`r(r+1)/2` ones, whose parity is `ceil(r/2)`. The matrix determinant
lemma gives (4).

A cycle block is `I+S`, where `S` is a cyclic permutation. Its rank is
`r-1` and every entry of its adjugate is one. Other balanced blocks
have determinant one, so the rank-one update contributes `r^2=r`
modulo two. This proves (2).

An unbalanced path has a maximal minor equal to one. Its only left or
right null vector is the all-ones vector on its larger side. With one
path of each orientation, the full matrix has corank one and its
adjugate is the outer product of these null vectors. The all-ones
update therefore contributes `ab`, giving (3).

## 3. Construction cost and a useful limitation

Store the owner of each endpoint in `B`. Each edge of `A` then finds at
most two neighbors. Deduplicate them when the edge occurs in both
matchings. A component traversal examines `2k` graph vertices and at
most `2k` graph edges. Equations (1)--(4) then use only component
sizes and parities. This proves the stated linear graph-operation
bound. With binary endpoint identifiers, the bit cost includes their
length and dictionary operations.

For two perfect matchings on the same support there are no path
components except common edges. The formula reduces to Section 5 of
[the matching-code lift](matching-code-lift.md). The companion
[Jack-integrality theorem](matching-gram-unipotence.md) proves that its
Gram matrix is invertible at every even particle count. Thus a fast
entry oracle alone cannot justify a small matrix-rank claim.

A different matrix, which tests whether the union is one Hamiltonian
cycle, has an explicit `2^(k-1)`-dimensional factorization over `F_2`;
see [Cygan, Kratsch, and Nederlof, Theorem 3.4 and Corollary 3.5](https://arxiv.org/pdf/1211.1506).
Applying that result to complete crossing pieces would require an
additional factorization. The component formula above does not supply
that missing argument.

## 4. Exact verification

Run:

```sh
python3 tools/check_partial_matching_minors.py --output evidence/partial-matching-minors.json
```

The checker compares the component algorithm with independent Gaussian
elimination on every equal-size pair of partial matchings through
eight physical vertices: 257,664 pairs. It also checks 1,000 fixed-seed
samples through forty vertices and rejects malformed matchings.
Every component case with a possible nonzero value is tested at both
parities. These checks validate the implementation; the proof above
establishes the formula for all sizes.
