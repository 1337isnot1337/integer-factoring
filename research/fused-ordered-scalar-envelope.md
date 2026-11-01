# A complete rank-sum obstruction in the ordered binary envelope

Every binary decomposition of the nine-letter crossing inside the full
ordered two-copy scalar-projection envelope of the literal q3, theta4
fused seed has cost at least `3*sqrt(56)=sqrt(504)`. This exceeds the
critical finite-block cost `9^sqrt(2)`. The result covers any number of
pieces in this 33-dimensional space. It does not cover higher tensor
powers, other coefficient fields or physical type restrictions.

The proof combines two exhaustive threshold audits with an integer
bound on the number of affordable pieces. The original audits covered
all `2^33=8,589,934,592` coordinate vectors, using different thresholds
on the two halves. They are local author computations with frozen source
and data fingerprints, not an independent replay or mathematical review.
The complete factoring score remains `1/2`.

## 1. The complete ordered space

Use `T_0,T_1,T_2`, `chi` and the literal fused matrices from
[the symmetric-envelope certificate](fused-scalar-projection-envelope.md).
The full ordered space is

```math
V=\operatorname{span}_{\mathbb F_2}
\{\chi(a T_i\otimes T_j):i,j\in\{0,1,2\},\
a\in\{1,2,4,8\}\}.
```

Its 36 displayed generators have dimension 33. Projection is applied
after forming the whole two-copy tensor. The physical alphabet is nine
and each flattening is 81 by 81. The crossing `E_9` has coordinate
integer `1118481` in the retained basis.

Keep the functional `f(x)=x_1+x_19+x_20` from the symmetric audit. It
detects `E_9`, but the larger ordered space has low-product vectors with
`f=1`. Thus the earlier simple low-product-span argument does not by
itself prove the present result.

## 2. Two complete threshold audits

Let `p(X)=r_R(X)r_C(X)` be the product of the two binary ranks.

| Affine half | Vectors visited | Threshold | Eligible nonzero vectors |
|---|---:|---:|---:|
| `f=1` | `2^32` | `p<=501` | 21 |
| `f=0` | `2^32` | `p<=55` | 0 |

Both halves are divided into four disjoint shards. If `p` is the least
nonzero coordinate of `f`, the affine-half generators are
`e_i+f_i e_p`, for `i!=p`; offsets are `e_p` and zero. Fixing the last
two free coordinates gives the four shards. Gray-code traversal visits
every remaining free vector exactly once. The zero vector is visited
and excluded from the nonzero eligible list.

Every basis matrix has the exact 25-block occupation/charge splitting
used in the symmetric audit. The code verifies that no generator has
off-block entries. Exact block ranks add, and the product of partial
rank sums can only increase. A threshold rejection is therefore sound
even when the full ranks of that rejected vector are not computed.

The smallest rank product in the complete affine-one eligible catalogue
is 162. All other affine-one vectors have product at least 502. Together
with the affine-zero audit, this proves

```math
p(X)\ge56\quad\text{for every nonzero }X\in V.
```

The checker reconstructs all 21 eligible matrices independently as full
81 by 81 bit matrices, and checks both ranks by a different elimination
algorithm. It also computes both ranks of all crossing complements
`E_9+X`. Their minimum product is 3969. These are 84 independent full
matrix rank calculations, without using the block decomposition.

## 3. Every decomposition costs at least sqrt(504)

Consider `E_9=sum_j X_j`, with all pieces in `V`, and discard zeros.
Each remaining piece has cost at least `sqrt(56)`.

If there are at least three pieces, their cost sum is at least
`3*sqrt(56)`. With one piece, that piece is the crossing, whose two
ranks are both 81, so its cost is larger still.

For two pieces, exactly one has `f=1`. If that piece is absent from the
21-vector catalogue, its rank product is at least 502, and

```math
\sqrt{502}+\sqrt{56}>3\sqrt{56},\qquad502>4\cdot56.
```

If it is listed, its other piece is the uniquely determined crossing
complement, of rank product at least 3969. The same inequality follows
from `3969>4*56`. This covers every decomposition, including repeated
pieces or an unrestricted number of summands.

Finally the comparison with the critical cost uses integers only:

```math
99^2>2\cdot70^2,\qquad504^{35}>9^{99}.
```

Hence `sqrt(504)>9^(99/70)>9^sqrt(2)`. By the
[finite-block amplification criterion](finite-block-crossing-amplification.md),
no decomposition in this full ordered space can cross the strict
`1/sqrt(2)` threshold.

## 4. Receipts, quick validation and a full replay

```sh
make verify-fused-ordered-projection
```

This quick command reconstructs the complete space and the exact binary
audit input, verifies its fingerprint and both frozen native-source
fingerprints, checks the eight shard receipts, recomputes all 84 full
ranks, and verifies the rank-sum argument. It does not repeat the large
threshold enumeration.

```sh
make verify-fused-ordered-projection-full
```

The full command additionally compiles the frozen C++17 source and
repeats both complete audits with four low-priority native processes.
The original combined runs took roughly an hour of wall time. The
affine-zero source is derived from the frozen snapshot by exactly two
literal replacements: the offset becomes zero and the threshold becomes
55. Its separately retained hash matches the source used for that run.
The original eight native computations passed; the packaged full-replay
wrapper was not used to repeat those computations after packaging.

The archived source is intentionally retained byte for byte as an
evidence artifact. The ordinary symmetric checker remains the readable
implementation of the same block-rank operations.

Sources and records:

- [`check_fused_ordered_projection.py`](../tools/check_fused_ordered_projection.py)
- [`fused-ordered-scalar-audit-source.cpp`](../evidence/fused-ordered-scalar-audit-source.cpp)
- [`fused-ordered-scalar-audit.json`](../evidence/fused-ordered-scalar-audit.json)
- [`fused-ordered-scalar-projection.json`](../evidence/fused-ordered-scalar-projection.json)

The seed is the repository's existing finite fusion construction. This
rank-sum argument, native audits and replays were developed with
substantial OpenAI assistant support. No new Lean formalization or
improved factoring exponent is claimed.

## Later piece-count refinement

The [completed F4 subspace analysis](fused-f4-decomposition-obstruction.md) also sharpens the argument for this full binary catalogue. Every decomposition contains a nonzero-functional piece of product at least 162. Three or more nonzero pieces therefore cost at least `sqrt(162)+2sqrt(56)>27`; the catalogue and complement bounds put all two-piece cases above 28, and the single crossing piece costs 81. Thus the balanced finite-block rank-sum parameter here is strictly above `3/4`. This uses the existing exhaustive binary data and fresh quick rank checks; it is not another binary native replay, a full-space F4 result or a factoring improvement.
