# Orthogonal retractions for paired differential codes

## Result and scope

The [paired differential construction](f4-differential-pair-codes.md)
usually uses a coordinate retraction. The same kernel also admits an
orthogonal retraction exactly when a finite list of subset Gram matrices
has the expected rank. This gives a different pair of physical maps,
with the same code dimension and height distribution.

For the four-column F4 MDS example, 18 of the 81 nonzero column scalings
satisfy every subset condition. They give explicit 92-letter orthogonal
codes on eight modes. No improved crossing pressure or factoring exponent
is asserted. The purpose is to construct and test interior changes to
both maps, rather than change only their last height boundary.

## 1. An exact criterion

Work over a characteristic-two field K. Let disjoint mode pairs
`(a_i,b_i)`, for `1 <= i <= m`, have one common total height H. Define

\[
T_j=\sum_i C_{ji}D_{a_i}D_{b_i},\qquad
V=\bigcap_j\ker T_j.
\]

Use the ordinary bilinear pairing in the occupation basis, with no field
conjugation. For an active subset S of pairs, let C_S be the corresponding
column submatrix.

**Theorem.** On the full occupation algebra, the restricted pairing on V
is nondegenerate if and only if

\[
\mathrm{rank}(C_SC_S^{\mathsf T})=\mathrm{rank}(C_S)
\quad\hbox{for every }S\subseteq[m].                 \tag{1}
\]

When (1) holds, every height truncation of V also has nondegenerate pairing.
Unused modes are allowed and do not change the criterion for the full space.
The reverse implication is asserted for the full space, not every truncation.

**Proof.** Split each pair's occupation space into the span of its vacuum
and double occupation, and its two singleton lines. These three summands
are mutually orthogonal. Fix the set S of pairs using the first summand
and fix the singleton choice at each remaining pair. This gives an
orthogonal direct-sum decomposition of the full algebra.

Write z_i for double occupation in an active pair. The sector algebra is
K[z_i : i in S]/(z_i squared), and the operators are contractions by the
rows of C_S. Their common kernel is the exterior algebra of
`U_S = ker C_S`. Here the exterior algebra is commutative, since the
characteristic is two. An invertible linear change of the z variables
reduces the contractions to coordinate derivatives and proves this claim.

The occupation pairing on each exterior degree is the induced determinant
pairing. Thus all degrees are nondegenerate exactly when U_S is nondegenerate.
Necessity already follows from degree one; sufficiency follows by taking
exterior powers of an invertible Gram matrix.

In the standard nondegenerate pairing on K^S, U_S is the orthogonal
complement of the row space of C_S. That row space is nondegenerate exactly
when its Gram matrix has its full row-space rank, which is (1), even if
C_S has redundant rows. A subspace and its orthogonal complement are
nondegenerate simultaneously. This proves the full-space assertion.

Every active generator z_i has the same height H. Within a fixed singleton
sector, height truncation therefore retains whole exterior degrees.
Distinct sectors and exterior degrees are orthogonal, so a direct sum of
these nondegenerate pieces is nondegenerate. Unused modes add orthogonal
one-dimensional occupation factors. This proves the truncation assertion.

## 2. Explicit physical maps

Choose a basis P of the truncated kernel, separately in each height and
occupation cell, using the existing differential constructor. Set

\[
G=P^{\mathsf T}P,\qquad L=G^{-1}P^{\mathsf T}.         \tag{2}
\]

Then LP=I. Both the columns of P and the transposes of the rows of L span
V. They are therefore closed under every physical deletion operator.
The inverse Gram matrix acts within height and occupation cells, so both
maps preserve those grades. All hypotheses of the coherent compound
construction hold for both cuts.

The ambient projector PL is symmetric and idempotent. Equation (2) does
not require an orthonormal basis, and it does not assert L=P transposed.
The two full crossing profiles still have to be computed with these actual
maps. A favorable kernel dimension alone is not a pressure estimate.

### The two natural cut profiles coincide

Let B denote the common ambient compound operator in the two-cut
construction, with the fixed column permutation understood. Its two
compressed matrices have the forms

\[
M_0=(L\otimes L)B(P\otimes P),\qquad
M_\infty=(L\otimes P^{\mathsf T})B(L^{\mathsf T}\otimes P).
\]

Because `P transposed = G L` and `L transposed = P G inverse`,

\[
M_\infty=(I\otimes G)M_0(G^{-1}\otimes I).             \tag{3}
\]

The extra factors are constant invertible matrices. They preserve height
and occupation, so they act inside the same charge blocks and commute
with natural height diagonals. Thus both complete natural Smith profiles
coincide. This is an algebraic consequence for these orthogonal retractions,
not an inference from a numerical equality at small sizes. Arbitrary
additional letter gauges need not commute with G and are not covered.

## 3. An explicit F4 example

Let `omega squared = omega + 1`. Use pairs `(0,7),(1,6),(2,5),(3,4)`,
all of height 16, and

\[
C=\begin{pmatrix}
1&1&\omega&\omega\\
0&1&1+\omega&1
\end{pmatrix}.                                      \tag{4}
\]

This is the matrix with projective columns `(1,t)`, for t in F4,
scaled by `(1,1,omega,omega)`. Every two columns are independent. The
paired Hilbert formula gives dimension 92 on the full eight-mode algebra.
All 16 subset Gram tests pass, including the empty subset.

The [constructor](../tools/orthogonal_differential_pairs_f4.py) implements
(1) and (2). The [checker](../tools/check_orthogonal_differential_pairs_f4.py)
examines all 81 nonzero column scalings of these four projective points.
It compares the criterion with direct Gram inversion on every full kernel.
Exactly 18 pass. For each accepted map it checks LP=I, both deletion
closures and symmetry of the ambient projector. The
[receipt](../evidence/orthogonal-differential-pairs-f4.json) records the
subset ranks and source hashes.

## 4. Column scalings exist over a controlled finite extension

The criterion is not specific to four columns. Start with any fixed
r-by-m coefficient matrix C over a finite characteristic-two field F.
There is a finite extension E/F and a nonzero scaling of every column
such that every subset satisfies (1).

For each S choose a row basis A_S of C_S, of rank r_S. Introduce scaling
variables lambda_i. Cauchy--Binet gives

\[
\det(A_S\,\mathrm{diag}(\lambda_i^2)\,A_S^{\mathsf T})
=\sum_{I\subseteq S,\ |I|=r_S}
 \det(A_{S,I})^2\prod_{i\in I}\lambda_i^2.           \tag{5}
\]

At least one maximal minor is nonzero, and different I give different
monomials. Thus (5) is a nonzero polynomial. For rank zero use the constant
polynomial one. Multiply these polynomials over all S, and multiply by
`product lambda_i`. The result is nonzero and has total degree at most

\[
\Delta=m+2r\,2^m.
\]

Over any extension with more than Delta elements, this polynomial cannot
vanish at every point. One proof applies the usual root bound one variable
at a time; each individual degree is less than the field size. A nonzero
evaluation gives the required nonzero column scalings. The extension can
be chosen with `log2 |E| = O(log2 |F| + m + log(r+1))`.

Exhaustive search is a finite deterministic construction, but its cost
can be exponential in the seed parameters. No polynomial-time search or
uniform pressure improvement follows. The field must be fixed before
applying a word-power transfer, or its growing-field costs must be charged.
