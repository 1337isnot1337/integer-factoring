# Complete coefficient ranks with the Ward operator retained

For one global form, the Ward factorization uses
\[
 \eta(R_P+\eta O_P)\mathscr M_P,
 \qquad\eta=q+r,
\]
where \(R_P\) is the incidence operator and \(O_P\) inserts an
external full pair. The earlier Euler estimate bounds the derivative
and creation contributions separately. The calculation here combines
them before taking ranks.

The result is a verified improvement of several finite upper bounds.
None of the recorded cases meets the factoring seed criterion.

## 1. Exact type decomposition with one insertion

The local particle-hole decomposition is
\[
 S=2\Delta\oplus2G\oplus H.
\]
Both the incidence and creation operators preserve each local
exterior-degree block. The two copies of \(G\) have identical
operator matrices, and both operators vanish on the scalar blocks.

Center the parameters as \(p=1+x,q=1+y,r=1\), so \(\Delta=x\).
For one site put
\[
 B=(D+\eta O)S,
 \qquad D=\partial_p+\partial_q+\partial_r.
\]
The exact Ward matrix is
\[
 \eta\sum_i S^{\otimes(i-1)}\otimes B\otimes
                       S^{\otimes(P-i)}.
 \tag{1}
\]
For a fixed product of local matrix entries \((s_i,b_i)\), its
ordinary and one-insertion values are calculated by
\[
 (u,v)\longmapsto(us_i,\,vs_i+ub_i),
 \qquad (u,v)_{\rm initial}=(1,0).
 \tag{2}
\]
This is exact arithmetic in the binary polynomial ring. No separate
rank subadditivity is used for the two terms of \(B\).

For a type with \(s\) scalar sites, \(g\) odd sites, and \(h\)
middle sites, the multiplicity is
\[
 \binom Ps\binom{P-s}{g}2^{s+g},\qquad s+g+h=P.
\]
The scalar sites contribute \(x^s\). The implementation enumerates
the full row and column grade slices in each remaining \(G/H\)
type, forms every coefficient matrix from (2), and takes its exact
binary rank. Summing with the type multiplicities gives the exact
ranks of the unrestricted Ward factor matrices.

For the row cut the grades are \((L,L-2;L,L)\), with \(\eta=y\).
For the column cut the exact parameter-swap reshape uses grades
\((L,L;L,L-2)\), with \(\eta=x\), and exchanges the two coefficient
indices afterward. These full matrices bound a chosen homogeneous
code after its one-sided input and output maps are applied.

The implementation retains the minimum with earlier verified bounds
and the logical dimension cap. Its reported cost is
\(\sum_{a,b}\lceil\sqrt{r_{a,b}c_{a,b}}\rceil\), so it is an
integer upper bound on the required coefficient cost.

## 2. Recorded improvements

The following cases use the full homogeneous global-form radical
dimension \(d\).

| \(P\) | \(L\) | \(d\) | Earlier Euler cost bound | Combined-operator cost bound | Largest type dimension |
|---:|---:|---:|---:|---:|---:|
| 3 | 3 | 6 | 110 | 108 | 32 |
| 4 | 4 | 26 | 1,696 | 1,566 | 148 |
| 6 | 4 | 64 | 9,906 | 9,246 | 1,425 |
| 8 | 4 | 118 | 35,770 | 33,478 | 5,936 |

These are bounds on the coefficient decomposition being evaluated.
They need not improve the trivial one-piece decomposition of the
crossing tensor; the recorded values all exceed \(d^2\).
The purpose of these finite controls is to verify that retaining the
operator can remove rank that survives the separate-rectangle bound.

For example, at \(P=L=3\), the full Ward coefficient \((a,b)=(1,2)\)
has row rank 36 and column rank 85. The separate Euler estimates
were 51 and 100. The actual six-label restricted coefficient has
still smaller ranks, 32 and 5, showing that the code maps retain
additional information not captured by these unrestricted factors.

## 3. Verification and scaling limit

The evaluator is
[ward_operator_types.py](../tools/ward_operator_types.py).
The records, including both complete rank profiles, are in
[ward-operator-types.jsonl](../evidence/ward-operator-types.jsonl).

The independent checker
[check_ward_operator_types.py](../tools/check_ward_operator_types.py)
verifies all 512 local operator-block entries. It reconstructs the
complete unrestricted polynomial tensors at \((P,L)=(2,2),(3,3)\)
without type multiplicities and matches both cut-rank profiles.
Those controls contain 216 and 48,000 polynomial entries,
respectively. It also checks four previously recorded restricted
code profiles, including full and half-orientation codes at
\((P,L)=(4,4)\). The evidence is
[ward-operator-types-check.json](../evidence/ward-operator-types-check.json).

This type calculation can still contain exponentially large matrices
when the grades grow proportionally to \(P\). It supplies finite
exact controls for the
[eight-dimensional stable operator interface](ward-stable-operator-interface.md).
Using that interface to carry the operators through further fusion,
while retaining useful grade and code restrictions, remains the
larger-scale rank problem.

