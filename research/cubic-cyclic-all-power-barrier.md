# Cyclic complementary maps: an all-field, all-power barrier

## Theorem

Let \(S\) be the cyclic shift on seven labels, and let \(\mathcal V\) be the full formal coefficient space of the cubic output maps
\(C=\lambda S\). Allow all scalar parameters \(\lambda\), all linear combinations of their polynomial coefficient directions,
and every field of characteristic two.

For every positive integer \(n\), every decomposition of the crossing \(E^{\otimes n}\) into tensors from
\(\mathcal V^{\otimes n}\) has rank cost
\[
 \boxed{\quad
 \sum_i\sqrt{\operatorname{rank}R(T_i)\operatorname{rank}C(T_i)}
 \;\geq\;336^{n/2}\;>\;7^{\sqrt2 n}.
 \quad}
\]
This closes the full cyclic scalar-output family, including recombinations between different output parameters and all tensor powers.
The stronger one-copy estimate is \(\sqrt{546}\).

The theorem follows from a constant row compression and a simultaneous triangular form in the column cut.
The triangular form produces a multiaffine polynomial on a seven-point grid at every tensor power.
No finite-field enumeration is used to establish the field or power quantifiers.

## Space and coefficient notation

Use the cubic construction and raw coefficient order of
[the complementary-output identity](cubic-complementary-output-maps.md):
\[
 (r^3,\;pr^2,\;p^2r,\;qr^2,\;q^3).
\]
The input map is fixed, and the output satisfies \(\Pi\Omega e_i=e_i\) and \(\Pi\Omega f_i=\lambda S e_i\),
with matrix indices interpreted as in the implementation.

Write \(z^3+z+1=0\), so \(z\) has multiplicative order seven.
For the two output maps \(S,zS\), concatenate their five raw directions and select indices
\[
 (0,1,2,3,6,7,8).
\]
These seven directions are an independent basis \(G_0,\ldots,G_6\) with \(G_0=E\).
For \(T=\sum_{i=0}^6w_iG_i\), define the linear forms
\[
 e(T)=w_0,\qquad
 t(T)=(z^2+1)(w_4+w_5),\qquad
 v(T)=w_3+z^2w_6.
\]
In particular, \(e(E)=1\).

Each raw coefficient for \(C=\lambda S\) is quadratic in \(\lambda\), since both output legs depend linearly on \(\lambda\).
The checker reconstructs its three formal parameter coefficients from \(\lambda=0,1,z\).
All linear parameter coefficients vanish, and all remaining formal matrices have entries in F2.
The zero-parameter directions lie in the selected seven-dimensional space.
Thus this space contains every scalar parameter over every field extension, and it is the scalar extension of a space defined over F2.

For a characteristic-two field that does not contain F8, extend it to a field containing a root of \(x^3+x+1\).
Matrix rank is unchanged by field extension. All subsequent calculations can therefore be performed in a field containing the seven distinct roots of unity.

## Fourier blocks

Both flattenings commute with the simultaneous cyclic shift of their two row labels and two column labels.
For a seventh root \(\rho=z^q\), the Fourier basis on ordered label pairs is
\[
 b_{\rho,k}=\sum_{x=0}^6\rho^{-x}e_{(x,x+k)},\qquad k=0,\ldots,6.
\]
Each flattening splits into seven blocks of size seven.
The block entries are computed by
\[
 B_{\rho}[k,l]
 =\sum_{a=0}^6\rho^{-a}
 M\bigl[(0,k),(a,a+l)\bigr].
\]
All label additions are modulo seven.

Every block further splits according to the relative-coordinate sets
\[
 \{0\},\quad\{1,6\},\quad\{2,5\},\quad\{3,4\}.
\]
The checker compares every coefficient of every block entry, rather than checking selected numerical values of the \(w_i\).
The formulas are therefore identities of linear matrix pencils over F8 and all its extensions.

## Constant row compression

For \(k\in\{2,3\}\), the row block on \(\{k,7-k\}\) is
\[
 \begin{pmatrix}
 e&t\rho^{-k}\\
 t\rho^k&e
 \end{pmatrix}.
\]
The relative coordinates 2 and 3 belong to different blocks. Restricting to those two coordinates gives \(eI_2\).
Taking this restriction at all seven charges produces constant linear maps with
\[
 L_R R(T) P_R=e(T)I_{14}.
\]
Consequently, for every \(T\in\mathcal V^{\otimes n}\),
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}R(T)\geq14^n.
\]
The maps used here depend only on the fixed Fourier bases and selected coordinates.

## Column triangular form

The corresponding column block is
\[
 A_{\rho,k}(T)=
 \begin{pmatrix}
 e&v\rho^{-k-1}\\
 v\rho^{k-1}&e
 \end{pmatrix},\qquad k\in\{2,3\}.
\]
The constant change of basis
\[
 P_{\rho,k}=
 \begin{pmatrix}\rho^{-k}&0\\1&1\end{pmatrix}
\]
gives the simultaneous identity
\[
 P_{\rho,k}^{-1}A_{\rho,k}(T)P_{\rho,k}
 =
 \begin{pmatrix}
 e+\rho^{-1}v&\rho^{-1}v\\
 0&e+\rho^{-1}v
 \end{pmatrix}.
\]
Both generic pairs therefore give four diagonal copies of
\(\ell_\rho=e+\rho^{-1}v\) at each charge.
Their direct sum is a constant 28-dimensional restriction of the full column flattening, simultaneously upper triangular for the entire coefficient space.

At tensor power \(n\), the restricted matrices remain upper triangular.
For \(T\in\mathcal V^{\otimes n}\), their diagonal entries are evaluations of
\[
 F_T(x_1,\ldots,x_n)
 =\bigl((e+x_1v)\otimes\cdots\otimes(e+x_nv)\bigr)(T)
\]
on the grid of seventh roots, each repeated \(4^n\) times.
This polynomial has degree at most one in each variable. Its constant coefficient is \(e^{\otimes n}(T)\).

A nonzero polynomial of degree at most one in each variable is nonzero at at least \(6^n\) points of any product of seven-element sets.
For completeness, this follows by induction: write the polynomial as \(A(x_1,\ldots,x_{n-1})x_n+B(x_1,\ldots,x_{n-1})\).
If \(A\ne0\), at least \(6^{n-1}\) choices make \(A\ne0\), and each leaves at most one root in the last coordinate.
If \(A=0\), apply the induction hypothesis to \(B\) and keep all seven choices of the last coordinate.

An upper triangular matrix has rank at least the number of its nonzero diagonal entries: the principal submatrix on those indices is invertible.
It follows that
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}C(T)\geq4^n6^n=24^n.
\]

## Rank cost and the stronger one-copy estimate

If \(\sum_iT_i=E^{\otimes n}\), then
\(\sum_i e^{\otimes n}(T_i)=1\). At least one piece has a nonzero value of this functional.
That piece alone has rank cost at least
\[
 \sqrt{14^n24^n}=336^{n/2}.
\]
The exact integer inequalities
\[
 336^{35}>7^{99},\qquad99^2>8\cdot35^2
\]
give \(336>7^{2\sqrt2}\), proving the strict barrier at every power.

At one copy, the row estimate improves to 21. When \(t\ne e\), each of the two generic row pairs is invertible at every charge.
When \(t=e\ne0\), each generic pair has rank one, and the row block on \(\{1,6\}\) contains the nonzero entry \(t\rho\).
Thus every row Fourier block has rank at least three.

For the column cut, the determinant of each generic pair is \((e+\rho^{-1}v)^2\).
When \(e\ne0\), at most one charge can make this determinant zero.
The two generic pairs contribute rank four at the other six charges and at least two at the exceptional charge.
Hence the column rank is at least 26, giving the one-copy cost bound \(\sqrt{21\cdot26}=\sqrt{546}\).

## Certificate and scope

The [symbolic checker](../tools/check_cubic_cyclic_fourier_bound.py) regenerates the literal seven-label coefficient matrices,
checks all 33,614 cyclic-equivariance entries, all 4,802 Fourier coefficient entries, and all 392 simultaneous-triangularization entries.
It also checks the 72,030 entries of the formal scalar-parameter matrices and their descent to F2.
The [certificate](../evidence/cubic-cyclic-fourier-bound.json) records every Fourier block, the coordinate identities, and source hashes.

```bash
python tools/check_cubic_cyclic_fourier_bound.py \
  --output /tmp/cubic-cyclic-fourier-bound.json
```

This theorem strengthens the cyclic portion of the [mixed F8 search](cubic-mixed-F8-spaces.md).
It concerns the scalar cyclic complementary maps and the full tensor span they generate.
Other complementary matrices, different output selectors, and mixtures with coefficient directions outside this space remain outside the theorem.
