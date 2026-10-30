# A uniform binary barrier for cancelling P6 boundary sections

## Result

Consider the six-pair Pfaffian family with the \(64\) single-occupation
labels as outputs. Modify each input by an arbitrary
\(\Omega\)-boundary, where
\[
\Omega=\sum_{i=1}^{6} e_{2i-1}\wedge e_{2i}.
\]
Suppose the resulting section cancels the pure \(p^6\) coefficient.

For every such section, every decomposition of the crossing using
binary combinations of its universal coefficient tensors has
\[
\boxed{S\ge 834>64^{\sqrt2}=358.3638678\ldots.}
\]
This covers all orthogonal assignments of the active boundary quotient
and all additions from its radical. It replaces the earlier
rank-difference lower bound \(133\), which did not exclude the family.

The result concerns this boundary-section family and binary coefficient
recombinations. It does not address general exterior-space retracts,
higher pair counts, or extension-field coefficient combinations.

The construction, boundary Gram formulas, and exact active operators are
defined in
[noncommutative-crossing-pfaffian-isotropic.md](noncommutative-crossing-pfaffian-isotropic.md)
and
[noncommutative-crossing-pfaffian-boundaries.md](noncommutative-crossing-pfaffian-boundaries.md).

## 1. The common active quotient

With both output labels single-occupied, the boundary tensor splits by
the single-occupation syndrome of its two input boundaries.

For six pairs, only two types of sectors contribute.

- A syndrome of weight two leaves four even pairs. Its even boundary
  Gram is a fixed rank-one form multiplied by
  \(A^2\), where \(A=p(p+q+r)\). The two occupied sites carry the binary
  two-site tensor. There are \(15\) such syndromes and four odd-bit
  orientations, giving a nondegenerate contribution of dimension
  \(60\).
- The zero syndrome has six even pairs. Its Gram is
  \(A^2(A+B)\) times the fixed triangular-graph matrix \(T_6\), where
  \(B=qr\). This matrix has rank four and is alternating.

Every remaining sector vanishes identically. Therefore the radical of
the boundary Gram at \(p=1,q=r=0\) is a common radical of every universal
coefficient tensor.

The full correction consequently factors through one fixed
\(64\)-dimensional quotient with bilinear form
\[
I_{60}\oplus H_4.
\]
The explicit orthonormalization in the isotropic-section construction
identifies this form with \(I_{64}\).

Write \(I_0\) for the original single-occupation injection. Cancellation
of \(p^6\) by a section \(I_0+B_0\) says that the boundary map \(B_0\)
has Gram matrix \(I_{64}\). Its map to the active quotient is therefore
an isometry between two nondegenerate spaces of equal dimension.
It is invertible.

Fix one canonical isometry \(B_{\rm can}\). Every cancelling section has
the form
\[
B_0=B_{\rm can}O+R,\qquad O^{\mathsf T}O=I_{64},
\]
where \(R\) lands in the common radical. The radical term changes no
coefficient tensor. It is enough to control every orthogonal \(O\).

## 2. A nine-dimensional ambient coefficient space

Set
\[
a=r+q,\qquad b=p+q,\qquad A=p(p+q+r).
\]
Let \(E_i\) and \(C_i\) be the crossing and noncrossing tensors at site
\(i\), and define
\[
Q_k=[a^{6-k}b^k]\prod_{i=1}^{6}(aE_i+bC_i),
\qquad 0\le k\le6.
\]
Then \(Q_0=E\), the desired crossing, and
\[
Q_6=C,\qquad
C(x,y,u,v)=[x=u][y=v].
\]

For a cancelling boundary section, the complete tensor has the form
\[
G=
\prod_{i=1}^{6}(aE_i+bC_i)
+A^2\bigl(a^2 Z_0+abZ_1+b^2 C\bigr).
\]
Its coefficients lie in
\[
\mathcal V_O=
\operatorname{span}_{\mathbb F_2}
\{Q_0,\ldots,Q_6,Z_{0,O},Z_{1,O}\}.
\]
Using this entire space permits at least as many decompositions as using
only individual monomial coefficients, so a lower bound for this space
also bounds every original coefficient decomposition.

The two cuts are
\[
(x,y)\mid(v,u),\qquad (x,u)\mid(v,y).
\]
For the canonical section, exact elimination gives:

| Correction \(Z\) | Row rank | Column rank | Row rank of \(Z+C\) | Column rank of \(Z+C\) |
|---|---:|---:|---:|---:|
| \(Z_0\) | 3840 | 154 | 1696 | 154 |
| \(Z_1\) | 1696 | 18 | 3840 | 18 |
| \(Z_0+Z_1\) | 4096 | 154 | 240 | 154 |

All three correction rank pairs differ from the rank pair of every
binary combination of \(Q_0,\ldots,Q_6\). Thus the nine generators are
independent for every orthogonal assignment.

### Why the rank \(240\) appears

Permuting the row-cut columns to group each physical output with its
opposite physical output does not change rank. On a weight-two syndrome
block, the operators become local swaps \(S_i,S_j\). With
\(N_i=I+S_i\), one has
\[
C+Z_0+Z_1=N_iN_j.
\]
Each \(N_i\) has rank one on its four-dimensional local space. Their
product has rank one on the sixteen-dimensional two-site space.

There are sixteen physical spectator states for each syndrome and
fifteen syndromes. The even block cancels. Hence
\[
\operatorname{rank}_R(C+Z_0+Z_1)=15\cdot16=240.
\]
The other row ranks follow from the corresponding local combinations
of \(I,S_i,S_j\). The checker verifies all row and column ranks directly
on the complete \(4096\)-dimensional cuts.

## 3. The shifted rank bound is uniform

Changing the active boundary isometry by \(O\) acts only on the two
input legs \(y,v\). In the row cut,
\[
M_R(Z_O)=
(I_{64}\otimes O^{\mathsf T})\,M_R(Z)\,
(O\otimes I_{64}).
\]
Both factors are invertible. The column cut likewise undergoes an
invertible transformation on its input columns.

The tensor \(C\) is fixed by this transformation because
\(O^{\mathsf T}O=I_{64}\). Therefore the ranks of both \(Z_O\) and
\(Z_O+C\) are exactly those in the table, for every orthogonal \(O\).

For a fixed binary baseline combination \(T\), use
\[
T+Z_O=(T+C)+(Z_O+C).
\]
The rank inequality
\(\operatorname{rank}(U+V)\ge
|\operatorname{rank}U-\operatorname{rank}V|\)
can now be applied in two ways. For either cut \(D\in\{R,C\}\),
\[
\operatorname{rank}_D(T+Z_O)\ge
\max\left\{
\left|\operatorname{rank}_D T-\operatorname{rank}_D Z\right|,
\left|\operatorname{rank}_D(T+C)-
                   \operatorname{rank}_D(Z+C)\right|
\right\}.
\]
When \(Z\ne0\), independence also gives a lower bound of one.

This extra shift is essential. For example, a baseline tensor and a
correction may both have row rank \(4096\), giving the unshifted
difference zero. After adding \(C\) to both, the correction can have
row rank only \(240\), and the new difference is much larger.

## 4. The exact finite certificate

There are \(128\) binary baseline combinations and four correction
choices \(0,Z_0,Z_1,Z_0+Z_1\). For each of the \(512\) resulting
coefficient masks, the checker computes the uniform lower rank pair
\((r_{\rm low},c_{\rm low})\) above and assigns the integer lower cost
\[
w=\left\lfloor\sqrt{r_{\rm low}c_{\rm low}}\right\rfloor.
\]

Use coordinates \(Q_0,\ldots,Q_6,Z_0,Z_1\); the crossing has mask \(1\).
A binary decomposition must have XOR of its masks equal to \(1\).
An integer-only shortest-path calculation on \(\mathbb F_2^9\) gives
minimum lower cost
\[
834.
\]
The same calculation using only the unshifted rank differences gives
\(133\), providing an independent check that the improvement comes from
the invariant \(C\)-shift.

Finally,
\[
834^{70}>64^{99},
\qquad \sqrt2<\frac{99}{70}.
\]
Every actual real-valued decomposition cost is at least the integer
shortest-path lower bound, so every such boundary section fails the
seed criterion.

This proof permits arbitrarily many pieces. It does not depend on a
sample of orthogonal maps or a bound on their matrix support.

## 5. Reproduction

The checker is
[check_pfaffian_p6_boundary_orbit.py](../tools/check_pfaffian_p6_boundary_orbit.py).
It reconstructs the canonical section using the independently audited
isotropic constructor, forms both full cuts, computes all \(128\)
baseline rank pairs, and verifies every rank in the correction table.
It then recomputes all \(512\) uniform lower costs and the exact integer
shortest path.

~~~sh
python3 tools/check_pfaffian_p6_boundary_orbit.py \
  --output evidence/pfaffian-p6-boundary-orbit-barrier.json
~~~

The evidence file
[pfaffian-p6-boundary-orbit-barrier.json](../evidence/pfaffian-p6-boundary-orbit-barrier.json)
contains all baseline ranks, correction ranks, uniform lower rank pairs,
integer weights, shortest-path witness masks, and source hashes.
