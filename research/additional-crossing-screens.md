# Complete control permutations and bipartite Pfaffian codes

This checkpoint preserves three finite searches completed while the
structural obstructions and symplectic codes were being developed. None
gives a seed below the strict threshold
\(S<d^{\sqrt2}\), where \(d\) is the common label dimension.

The integer lower certificates use
\(\lfloor\sqrt{rc}\rfloor\) as each piece's cost. A shortest path in the
binary coefficient space gives a lower bound for every binary
decomposition. The exact comparison \(S_{\rm int}^{70}>d^{99}\) certifies
failure because \(99/70>\sqrt2\). Displayed noninteger costs are numerical
values of recorded decompositions; the exclusions use integer arithmetic.

## 1. Every permutation of the eight D8 control words

Use the ordered binary group word in the generators of \(D_8\), with
parameter factors
\[
 g_j+z^{w_j}(1+g_j).
\]
For each group quotient \(S=xy^{-1}\), choose its word by an arbitrary
permutation of the eight words. This permits controls mixing horizontal
and central coordinates, outside the central-shear hypothesis.

For each of three weight vectors, the native program checks all
\(8!=40320\) permutations and every binary vector in the full coefficient
space. Its row-cut ranks are independent of the control permutation:
the permutation only rearranges the eight convolution blocks. The
column cuts are reconstructed for every control.

| Weights | Coefficient dimension | Controls checked | Minimum integer lower certificate | Best recorded cost |
|---|---:|---:|---:|---:|
| \((1,1,1)\) | 4 | 40320 | 50 | 50.3320209770 |
| \((1,1,2)\) | 5 | 40320 | 50 | 50.3320209770 |
| \((1,2,4)\) | 8 | 40320 | 50 | 50.3320209770 |

Every control in each search has integer lower certificate 50. The
identity control attains the best displayed cost, using the final
coefficient and the sum of the others. This is a search over all
**permutations**, not over all maps from the group to its word set.

The stronger [finite group gate](noncommutative-crossing-finite-gates.md)
now excludes arbitrary orbit kernels on this eight-element group over
every field. The finite computation is retained as reproducible evidence
of the enlarged control search, rather than as the basis of that theorem.

The independent audit constructs \(D_8\) as permutations of the four
vertices of a square. It enumerates all four tensor legs directly,
reproduces every binary rank in each winning profile, checks crossing,
and verifies the integer shortest-path bound. This does not independently
rerun all 120960 native control cases; their complete ranges and rank
digests are recorded in the three native logs.

## 2. Different operator and state occupation codes

Let \(J\) complement every physical mode pair. Partition one full
exterior grade into two disjoint, equally sized, \(J\)-invariant monomial
sets \(A,B\). Take operator labels in \(A\) and state labels in \(B\).
The crossing vertex remains the identity on \(A\times B\).

At the \(p\) vertex, a nonzero entry would require a label in both
\(A\) and \(B\). At the \(q\) vertex it would require a label in both
\(A\) and \(JB=B\). Thus both unwanted endpoint tensors vanish.
The complete three-parameter coefficient space remains available.

The implementation enumerates balanced choices of \(J\)-orbits. It fixes
the first orbit in \(A\), removing the exchange \(A\leftrightarrow B\)
duplication. Dynamic counting and unranking allow sampling without
materializing a large partition space.

| Physical pairs \(P\) | Grade \(L\) | Labels \(d\) | Partitions checked / available | Best recorded cost | Minimum integer lower certificate among checked cases |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 3 | 2 / 2 | 9 | 9 |
| 3 | 3 | 10 | 126 / 126 | 87.9047078536 | 86 |
| 4 | 2 | 14 | 32 / 3564 | 172.6289999049 | 171 |
| 4 | 4 | 35 | 32 / 9405140220 | 1002.6384951654 | 1002 |

The first two partition searches are complete. The four-pair searches
are deterministic samples only. For every listed partition, the binary
coefficient-space calculation is complete.

The independent audit evaluates the local tensor through the separate
[particle-hole exterior formula](pfaffian-particle-hole.md), constructs
the dense cut matrices without using syndrome blocks, and compares
every binary combination with the blocked calculation in both two-pair
cases and the best three-pair case. The archived report retains each
partition, integer certificate, decomposition witness, and rank-table
hash. Its source can regenerate every full rank table.

## 3. Mixed-radix group words

The same endpoint construction works with digits beyond zero and one:
\[
 \Phi_z(w)=\prod_j
 \bigl(g_j+z^{w_j}(1+g_j)\bigr)^{a_j},
 \qquad 0\le a_j<m_j.
\]
Evaluation at \(z=1\) is the identity word, for every digit tuple.
No group-algebra unit assumption is used.

| Group | Labels | Coefficient dimension | Binary integer lower certificate | F4 projective vectors | Affordable span dimension |
|---|---:|---:|---:|---:|---:|
| \(S_3\) | 6 | 4 | 34 | 85 | 1 |
| \(C_3^2\) | 9 | 5 | 81 | 341 | 2 |
| \(C_3^3\) | 27 | 7 | 713 | 5461 | 2 |

Every F4 projective vector is checked. In each case, the vectors below
the integer cost gate fail to span the crossing endpoint. The independent
scalar F4 auditor verifies every affordable vector and additional sampled
vectors. The archive contains the scalar-audit output and the native
logs; the source regenerates the exact bases from the group and weights.

The [involution-free minor theorem](noncommutative-crossing-involution-support.md)
now excludes all three full orbit spaces over every field. These finite
checks therefore close a concrete exploratory branch without requiring
larger searches in it.

## Reproduce

~~~bash
c++ -O3 -std=c++17 -Wall -Wextra -Werror \
  tools/explore_arbitrary_group_controls.cpp -o /tmp/arbitrary-controls
/tmp/arbitrary-controls 1 1 1
/tmp/arbitrary-controls 1 1 2
/tmp/arbitrary-controls 1 2 4

python3 tools/explore_bipartite_pfaffian_codes.py \
  --cases 2:2,3:3 --output /tmp/bipartite-small.json
python3 tools/explore_bipartite_pfaffian_codes.py \
  --cases 4:2,4:4 --limit 32 --seed 20261030 \
  --output /tmp/bipartite-four-pair.json

python3 tools/explore_mixed_radix_group_crossings.py \
  --group S3 --binary-search --output /tmp/mixed-S3.json
python3 tools/check_additional_crossing_screens.py
~~~

For each mixed-radix basis, use the existing native exporter and
`screen_pfaffian_span_f4.cpp`, followed by `check_pfaffian_span_f4.py`,
as documented in the [three-parameter family](pfaffian-ternary-edge-family.md).
The other two group names are `C3xC3` and `C3xC3xC3`.

Evidence:
[control 111](../evidence/arbitrary-D8-controls-111.jsonl),
[control 112](../evidence/arbitrary-D8-controls-112.jsonl),
[control 124](../evidence/arbitrary-D8-controls-124.jsonl),
[bipartite partitions](../evidence/bipartite-pfaffian-partitions.json),
[mixed-radix screens](../evidence/mixed-radix-crossing-screens.json), and
[independent audit](../evidence/additional-crossing-screens-audit.json).
