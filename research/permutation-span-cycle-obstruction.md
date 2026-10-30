# Single-cycle blocks obstruct the full permutation coefficient span

## Result and scope

The full coefficient span of the identity-kernel
[split-permutation crossing](split-permutation-seed.md) has a
stronger obstruction than its originally listed pieces.
Let the labels be all \(d=b!\) permutations, and let \(V_b\) be the
linear span of all its coefficient tensors over any field of
characteristic two.

**Theorem.** If

\[
                    E_d=\sum_j B_j,\qquad B_j\in V_b,
\]

then at least one \(B_j\) has

\[
 r(B_j),c(B_j)\ \ge\
 M_b:=\frac{b!(b-1)!}{2}=\frac{d^2}{2b}.                          \tag{1}
\]

There is no restriction on the number of summands or on their
coefficients in the field. In particular,

\[
 \max_j\frac{\log(r(B_j)c(B_j))}{4\log d}
 \ge 1-\frac{\log(2b)}{2\log(b!)}\longrightarrow1.                \tag{2}
\]

Thus arbitrary linear recombination does not give a uniform contraction
below one for this increasing full-label family.

The finite rank-sum criterion also fails for every \(b\ge4\). At
\(b=4,d=24\), every decomposition in the entire 141-dimensional span
satisfies

\[
                   \sum_j\sqrt{r(B_j)c(B_j)}\ge96
                      >24^{\sqrt2}.                            \tag{3}
\]

For \(b\ge5\), (1) alone exceeds the required rank-sum threshold.
Together with the universal exclusions for at most six labels, this
rules out that criterion on every full permutation domain \(b\ge2\).

This is specific to this coefficient span with identity mode kernel.
It does not apply to arbitrary linear restrictions of the external
label space, to a different kernel that breaks mode parity, or to
general crossing tensors. In particular the asymmetric kernels in
[balanced-partition-crossings.md](balanced-partition-crossings.md)
are outside its hypotheses.

## 1. The row and column matrices have fixed blocks

Recall that the polynomial tensor is zero unless

\[
                 A\triangle Y=B\triangle U.
\]

Every coefficient, and hence every linear combination of coefficients,
is therefore block diagonal in the row cut, with blocks indexed by
the symmetric difference \(D=A\triangle Y\).
The same parity identity can be written

\[
                 A\triangle U=B\triangle Y,
\]

so the column cut is block diagonal by its corresponding difference.

Consider a difference graph \(D\) that is one alternating cycle
through all \(b\) row vertices and all \(b\) column vertices. It has
exactly two ordered decompositions into permutation matchings.
There are no remaining vertices on which a common matching could be
chosen. Thus it gives a two-by-two block in each cut.

Put

\[
                  P_*=\prod_{i=1}^b P_i,\qquad
                  Q_+=\prod_{j=1}^b(1+Q_j).
\]

The compact formula in the split-permutation note gives the entire
polynomial row block as

\[
                   \begin{pmatrix}
                     P_*&P_*Q_+\\
                     P_*Q_+&P_*
                    \end{pmatrix},                             \tag{4}
\]

and the entire polynomial column block as

\[
                   \begin{pmatrix}
                     P_*&Q_+\\
                     Q_+&P_*
                    \end{pmatrix}.                             \tag{5}
\]

For (4), an orientation-preserving column pair has \(B=A,U=Y\)
and contributes \(P_*\); the opposite orientation contributes
\(P_*Q_+\). For (5), the opposite orientation instead has
\(A=Y,B=U\), with the two distinct matchings disjoint everywhere,
so its inactive-complement product is \(Q_+\).

The same formulas hold for every full single-cycle graph, irrespective
of its actual edges.

## 2. One linear functional forces large rank in both cuts

An element \(B\in V_b\) is obtained by applying a scalar coefficient
functional \(\mathcal L\) to the polynomial tensor. Write

\[
 a(B)=\mathcal L(P_*),\quad
 b(B)=\mathcal L(P_*Q_+),\quad
 c(B)=\mathcal L(Q_+).
\]

These quantities are well-defined functions of the tensor \(B\),
even if its coefficient-functional representation is not unique:
they are entries of any fixed single-cycle blocks (4)–(5).
Every such row block is
\(\begin{pmatrix}a&b\\b&a\end{pmatrix}\), and every such column
block is
\(\begin{pmatrix}a&c\\c&a\end{pmatrix}\).
If \(a(B)\ne0\), each block has rank at least one in both cuts.

There are \(b!(b-1)!\) ordered pairs of permutation matchings whose
union is one full alternating cycle: choose the first permutation,
then a relative \(b\)-cycle. Each graph is counted twice, once for
each orientation. Therefore the number of distinct blocks is \(M_b\)
from (1). Since the matrices are block diagonal,

\[
                       a(B)\ne0
          \quad\Longrightarrow\quad r(B),c(B)\ge M_b.             \tag{6}
\]

The crossing has \(a(E_d)=1\). In any decomposition of it, linearity
therefore gives \(\sum_j a(B_j)=1\), so at least one summand has
nonzero \(a\). Equation (6) proves the theorem.

The proof uses actual matrix blocks, not a declared factor dimension
or an estimate for a literal source intermediate space.

## 3. Excluding the finite rank-sum criterion

For \(b\ge5\), the elementary inequality
\[
                           b!>4b^2
\]
holds at \(b=5\) and propagates by induction. Indeed, after multiplying
by \(b+1\), the right side is at least \(4(b+1)^2\) because
\(b^2\ge b+1\). Consequently

\[
                 M_b=\frac{d^2}{2b}>d^{3/2}>d^{\sqrt2}.
\]

The distinguished summand alone makes the rank sum too large.

For \(b=4\), put \(N=d^2=576\), \(M=M_4=72\), and suppose for
contradiction that a decomposition has rank sum \(S<96\). Choose
the distinguished summand from (6), with ranks \(r,c\ge72\).
Since \(\sqrt{rc}<96\) and the ranks are integers,

\[
        r,c\le
        \left\lfloor\frac{96^2-1}{72}\right\rfloor=127.
\]

Rank subadditivity in the two cuts forces the remaining summands to
have total row rank and total column rank at least \(576-127=449\).
The elementary rank-budget inequality from Section 3 of
[finite-seed-amplification.md](finite-seed-amplification.md) uses
only those two sums and positivity of both ranks; it does not require
the remaining tensors to sum to a crossing. It gives a remaining
rank sum of at least \(2\sqrt{448}\). Hence

\[
                    S\ge72+2\sqrt{448}>96,
\]

a contradiction. This proves \(S\ge96\).

Finally, \(\sqrt2<10/7\), and
\[
                96^7>24^{10}
       \quad\Longleftrightarrow\quad 4^7>24^3
       \quad\Longleftrightarrow\quad 16\,384>13\,824.
\]
Thus \(96>24^{\sqrt2}\), proving (3).

The domains \(b=2,3\) have two and six labels. Their rank-sum
exclusions follow from
[the small-label bound](finite-seed-amplification.md) and
[the six-label obstruction](six-label-rank-sum-obstruction.md),
respectively.

## 4. What remains open after this obstruction

The 24-label full-span searches in
[seed-subspace-recombination.md](seed-subspace-recombination.md)
are now backed by an all-span exclusion for the half threshold.
The finite bound (3) does not prove that their observed optimum 325
is the optimum in the full span.

Three concrete operations can escape the argument:

- A dual linear label restriction need not preserve the separate
  symmetric-difference blocks. Its new dimension and rank sum require
  their own calculation.
- A nonidentity annihilation kernel can couple modes that are not
  equal. Then the parity condition making (4)–(5) separate matrix
  blocks need not hold.
- Tensor powers of one fixed seed may admit new collective
  recombinations. The finite rank-sum exclusion concerns decompositions
  in the single seed's coefficient span. It does not exclude every
  decomposition in its tensor-power spans. For example, the direct
  tensor-product block bound for the 24-label seed is only `72^t`,
  which is below `24^(sqrt(2)t)`. The single-seed rank-budget argument
  therefore does not establish the desired tensor-power exclusion.

None of these operations is asserted here to achieve the target.

## Verification

[tools/check_permutation_span_cycles.py](../tools/check_permutation_span_cycles.py)
checks every single-cycle row and column block through \(b=4\),
compares arbitrary coefficient recombinations with actual full matrix
ranks, and verifies the exact integer inequalities in Section 3.
The evidence is in
[evidence/permutation-span-cycles.json](../evidence/permutation-span-cycles.json).
