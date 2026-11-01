# Graded rank descent after the Ward factor

The factor \((q+r)^k\) from
[the commuting-form Ward theorem](noncommutative-crossing-radical-ward.md)
also gives a rank factorization through smaller exterior grades. This
note records the factorization explicitly. It improves the available
input to a rank calculation; it does not by itself assert a better
factoring exponent.

The row cut admits a sum of only \(O(k^2)\) rectangular grade blocks,
without a coefficient-channel multiplicity. A separate column bound
has an explicit channel factor. Either bound can be combined with
existing full-grade bounds and the cap \(d^2\).

## 1. Rectangular full-grade matrices

Let
\[
 M_P(R_A,R_Y;C_V,C_U)[a,b]
\]
denote the coefficient of \(x^a y^b\) in the full Pfaffian tensor at
\[
 p=1+x,\qquad q=1+y,\qquad r=1,
 \tag{1}
\]
with row labels \((A,Y)\) of exterior grades \((R_A,R_Y)\), and column
labels \((V,U)\) of grades \((C_V,C_U)\).

The four grades are independent. This is a matrix compression of the
unrestricted \(P\)-pair tensor, before any code is applied. If a grade
is outside \([0,2P]\), or a coefficient index is negative, the
corresponding matrix is zero.

Write \(r(M)\) for its \((A,Y)\mid(V,U)\) rank and \(c(M)\) for its
\((A,U)\mid(V,Y)\) rank. The second rank is a different rectangular
reshaping of the same coefficient tensor.

Let \(I\) be a homogeneous grade-\(L\) input code contained in
\[
 \operatorname{im}(Q_1\cdots Q_k),
\]
where the \(Q_a\) are the commuting weighted forms from the Ward theorem.
Choose any grade-\(L\) output map \(\Pi\). When \(\Pi I=I_d\), this is
the usual exact crossing retract. The rank factorization below holds
without requiring that identity.

Let \(G[a,b]\) be the centered coefficient tensor after applying
\(\Pi\) to \(A,U\) and \(I\) to \(Y,V\).

## 2. Homogeneous preimages of lower grade

**Lemma.** There is a constant matrix \(W\), with all columns of
exterior grade \(L-2k\), such that
\[
 I=Q_1\cdots Q_k W.
 \tag{2}
\]

**Proof.** Work syndrome by syndrome. A nonzero homogeneous product
image has independent restricted forms and zero Gram matrix, by the
dimension theorem in the Ward note. Use its coordinates
\[
 \Omega_a=e_a\wedge,\qquad \Lambda_a=\partial_{f_a},
\]
with \(t-2k\) residual variables \(h_i\). The scalar part of every
\(Q_a\) is zero on that syndrome. A homogeneous image basis consists of
\[
 e_1\cdots e_k h_S.
\]
The vector \(h_S\) is a preimage under the product of the \(Q_a\):
every contraction term is zero, and the creation terms insert the
\(e_a\). Removing \(k\) even-pair variables lowers physical exterior
grade by \(2k\). Linear combinations and different syndromes preserve
the common grade \(L-2k\). \(\square\)

The preimages in (2) are independent of \(p,q,r\). They need not be
orthogonal or coordinate vectors.

## 3. Put the Ward operator on the row-cut side

Let \(A_i\) now denote the external Grassmann variables. Set
\[
 \mathcal L_a=\sum_i w_{ai}
   (A_{\bar i}c_i+A_i c_{\bar i}),\qquad
 \mathcal O_a=\sum_iw_{ai}A_iA_{\bar i},
\]
where the sums have one representative per complementary pair.
Define constant operators
\[
 \mathcal T_j=[z^j]\prod_{a=1}^k(\mathcal L_a+z\mathcal O_a).
 \tag{3}
\]
Every term in \(\mathcal T_j\) has \(k-j\) state Clifford factors and
external \(A\)-degree \(k+j\).

The left-handed form of the Ward commutator identity is
\[
 \operatorname{ad}_{Q_1}\cdots\operatorname{ad}_{Q_k}(\mathscr M)
 =\eta^k\prod_{a=1}^k
     (\mathcal L_a+\eta\mathcal O_a)\mathscr M,
 \qquad \eta=q+r.
 \tag{4}
\]
Indeed, moving the state operator to the left in the first commutator
adds the same \(\eta\mathcal O_a\) term as in the right-handed formula.

Insert (2) into the two input legs, and use the square-zero sandwich
identity. The exact restricted tensor is
\[
 G=\eta^k\sum_{j=0}^k\eta^j
   (\Pi_A\otimes W_Y^{\mathsf T})\,
    \mathcal T_j\,\mathscr M\,
   (I_V\otimes\Pi_U^{\mathsf T}).
 \tag{5}
\]

Both parts of \(\mathcal T_j\) act on the **same row-cut side**:
multiplication by \(A_i\) acts on the \(A\) label, and \(c_i\) acts on
the \(Y\) label. Thus all summands defining \(\mathcal T_j\) can be
combined into one row-side linear map. Expanding its individual
Clifford summands is unnecessary for this cut.

To reach output \(A\)-grade \(L\), the underlying full tensor must have
\[
 R_A=L-k-j.
 \tag{6}
\]
Put \(\ell=k-j\). The \(\ell\) Clifford factors change state degree by
\(\ell-2h\), for some \(0\le h\le\ell\). Since \(W\) has grade
\(L-2k\), the underlying state row grade is
\[
 R_Y=L-3k+j+2h.
 \tag{7}
\]
The column grades remain \(C_V=C_U=L\).

In the centered specialization (1), \(\eta=y\). Therefore (5) gives
a sum of constant two-sided matrix factorizations through the
rectangles (6)--(7), with source coefficient index
\[
 (a,\ b-k-j).
 \tag{8}
\]
This accounts for every coefficient shift.

## 4. The row rank bound

**Theorem 1.** The tensor \(G[a,b]\) vanishes for \(b<k\). Otherwise,
\[
 \boxed{
 \begin{aligned}
 r(G[a,b])\le
 \sum_{j=0}^{\min(k,b-k)}
 \sum_{h=0}^{k-j}
 r\!\left(
 M_P(L-k-j,L-3k+j+2h;L,L)
       [a,b-k-j]
 \right).
 \end{aligned}}
 \tag{9}
\]

There is no factor such as \(\binom{k}{j}\), \(2^{k-j}\), or
\(P^{k-j}\) in this row bound. Each \(\mathcal T_j\) is a single linear
map. The inner sum merely separates its possible source grades.

For a fixed \(j\), one may obtain a better bound by stacking all its
source-\(Y\)-grade matrices vertically and taking their joint rank.
Equation (9) uses the sum of their ranks as a simple sufficient bound.

The largest source \(Y\)-grade in (7) is \(L-k-j\), the same as the
source \(A\)-grade. Thus every nonzero rectangle in (9) has both row
grades at most \(L-k\).

## 5. A column bound with explicit channels

For the column cut, \(A\) and \(Y\) lie on opposite sides. Here a channel
expansion is needed.

For a fixed \(j\), again put \(\ell=k-j\). A product of \(\ell\)
operators \(\mathcal L_a\) has the form
\[
 \sum_T F_T(A)c_{\bar T},
 \tag{10}
\]
where \(|T|=\ell\) and \(T\) contains at most one mode from each
complementary pair. To see the latter restriction, expand
\(\mathcal L_a=\sum_iw_{a,\operatorname{pair}(i)}A_i c_{\bar i}\)
over all modes. A repeated mode gives zero. The coefficient of a
distinct mode subset is the determinant of its weight columns, since
permanent and determinant agree in characteristic two. Choosing both
complementary modes repeats one column and gives zero.

The additional \(\mathcal O_a\) factors act only on the \(A\) side.
After all such terms are collected, (10) still needs at most
\[
 H_\ell=2^\ell\binom P\ell
 \tag{11}
\]
state channels. Each channel is a product of an \(A\)-side map and a
\(Y\)-side map, so it preserves the upper rank bound of the underlying
rectangle in the column cut.

**Theorem 2.** For \(b\ge k\),
\[
 \boxed{
 \begin{aligned}
 c(G[a,b])\le
 \sum_{j=0}^{\min(k,b-k)}
 2^{k-j}\binom P{k-j}
 \sum_{h=0}^{k-j}
 c\!\left(
 M_P(L-k-j,L-3k+j+2h;L,L)
       [a,b-k-j]
 \right).
 \end{aligned}}
 \tag{12}
\]
The actual channel rank may be smaller than (11), depending on the
weight matrix. Formula (12) is a uniform upper bound.

The existing \(p\leftrightarrow q\) symmetry supplies every column
rectangle from a row calculation:
\[
 \boxed{\quad
 c\!\left(M_P(R_A,R_Y;C_V,C_U)[a,b]\right)
 =
 r\!\left(M_P(R_A,C_U;C_V,R_Y)[b,a]\right).
 \quad}
 \tag{13}
\]
The accompanying mode-complement permutations act on full exterior
grades and do not change these ranks.

For either cut, one may take the minimum of (9) or (12), the original
full-grade coefficient bound, and the dimension cap \(d^2\).

## 6. Constructive use and verification

The full local tensor has the same constant grade-preserving stable
decomposition used by the existing Pfaffian rank engine. Its
compression may use independent row and column grades. Thus (9)--(13)
reduce the new code bound to explicit rectangular coefficient matrices;
they do not require constructing the entire restricted tensor first.

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_ward\_rank\_descent.py}
\]
with optional output
\(\texttt{--output evidence/noncommutative-crossing-ward-rank-descent.json}\).

The checker constructs literal full-grade rectangular coefficient
matrices and literal restricted coefficient tensors. It verifies both
bounds for the two-pair radical, the three-pair radical, its orientation
half, and a partial form on two of three physical pairs. All complete
matrix ranks are included; the checks do not use sampled minors.

These finite controls test the grade and coefficient bookkeeping. The
all-size justification is the exact factorization (5), the homogeneous
preimage lemma, and ordinary matrix-rank subadditivity.

