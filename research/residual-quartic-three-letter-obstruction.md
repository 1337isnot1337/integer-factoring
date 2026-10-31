# An all-field obstruction to three-letter residual codes

The residual endpoint from
[`residual-quartic-codes.md`](residual-quartic-codes.md) admits two-letter
physical codes. It admits no three-letter crossing restriction over any field
of characteristic two. This strengthens the earlier exhaustive binary and
four-element-field checks to a field-independent statement.

The obstruction concerns four independent **constant physical linear maps**
applied to one copy of the sixteen-dimensional residual operator. It does not
exclude polynomial codes, larger-block codes, or codes that use bulk modes.

## 1. The residual tensor and a matrix space

Each physical leg has the four states
\(\varnothing,\{0\},\{T\},\{0,T\}\). In exterior variables the residual tensor is

\[
 R=\Gamma(\text{crossing})
 (1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T). \tag{1}
\]

Flatten (1) with rows \((x,v)\) and columns \((y,u)\). This is the cut at which
the desired crossing has rank one. The \(16\times16\) matrix of (1) has rank
eight. Its column space, viewed as a subspace of \(4\times4\) matrices, is

\[
 W=\left\{
 \begin{pmatrix}
 t&b&0&0\\
 0&t&e&f\\
 c&d&t&b\\
 0&c&g&s
 \end{pmatrix}:t,b,c,d,e,f,g,s\in F
 \right\}. \tag{2}
\]

Its row space is \(W^T\). These statements follow by expanding the seven
factors in (1), and are checked directly by the accompanying program. All
entries are binary, so the ranks and spanning identities persist under every
characteristic-two field extension.

## 2. Every three-dimensional compression retains at least five dimensions

### Lemma

Let \(L:F^4\to F^3\) be surjective and \(J:F^3\to F^4\) be injective.
Then

\[
 \dim\{LAJ:A\in W\}\ge5. \tag{3}
\]

**Proof.** Write \(\ker L=Fa\) and \(\mathrm{im}J=\ker b^T\), with
\(a,b\in F^4\setminus\{0\}\). The matrices killed by this compression are
exactly those of the form

\[
 A=ar^T+sb^T,\qquad r,s\in F^4. \tag{4}
\]

Indeed, a killed matrix maps \(\ker b^T\) into \(Fa\), and extending that
one-dimensional image map gives (4). The map \((r,s)\mapsto ar^T+sb^T\)
has a one-dimensional kernel, consisting of \((r,s)=(\lambda b,\lambda a)\).

Membership of (4) in \(W\) imposes the following eight linear equations:

\[
\begin{array}{rlrl}
 a_0r_2+s_0b_2&=0,& a_0r_3+s_0b_3&=0,\\
 a_1r_0+s_1b_0&=0,& a_3r_0+s_3b_0&=0,\\
 a_0r_0+s_0b_0+a_1r_1+s_1b_1&=0,
 &a_0r_0+s_0b_0+a_2r_2+s_2b_2&=0,\\
 a_0r_1+s_0b_1+a_2r_3+s_2b_3&=0,
 &a_2r_0+s_2b_0+a_3r_1+s_3b_1&=0.
\end{array} \tag{5}
\]

We claim that the coefficient matrix in (5) always has rank at least four.

If \(a_0b_0\ne0\), the first four equations already have rank four: the first
two give independent \(r_2,r_3\) pivots and the next two give independent
\(s_1,s_3\) pivots. If \(a_0\ne0,b_0=0\), use the first two equations, the
first diagonal equation, and the first repeated-entry equation. Modulo the
first two pivots, the latter two have independent \(r_0,r_1\) coefficients.

Let \(P\) exchange coordinates one and two. Conjugation by \(P\) sends
\(W^T\) back to \(W\), and fixes coordinate zero. Hence the case
\(a_0=0,b_0\ne0\) follows by transposition and this permutation.

It remains to consider \(a_0=b_0=0\). The first four equations supply an
\(r_0\) pivot when \((a_1,a_3)\ne(0,0)\), and an \(s_0\) pivot when
\((b_2,b_3)\ne(0,0)\).

If both pivots are present, eliminate them. The remaining equations contain
two blocks on disjoint variable sets:

\[
 \begin{pmatrix}a_1&b_1&0\\a_3&0&b_1\end{pmatrix}
 \quad\text{on }(r_1,s_1,s_3),\qquad
 \begin{pmatrix}a_2&0&b_2\\0&a_2&b_3\end{pmatrix}
 \quad\text{on }(r_2,r_3,s_2).
\]

Each block is nonzero, adding at least two ranks to the two pivots.

If only the \(r_0\) pivot is present, then \(b_1\ne0\). The two equations
involving \(s_1\) and \(s_3\), together with
\(b_1s_0+a_2r_3=0\), add three independent ranks. If only the \(s_0\) pivot
is present, then \(a_2\ne0\). The two equations involving \(r_2,r_3\),
together with \(a_2r_0+b_1s_3=0\), again add three independent ranks.

If neither pivot is present, nonzeroness of \(a,b\) gives
\(a_2b_1\ne0\). Four independent equations remain:

\[
 b_1s_1=0,\qquad a_2r_2=0,\qquad
 b_1s_0+a_2r_3=0,\qquad a_2r_0+b_1s_3=0.
\]

This proves the rank claim in every case. The solution space of (5) has
dimension at most four. Quotienting by the one-dimensional redundancy in
(4) shows that the kernel of the compression on \(W\) has dimension at most
three. Since \(\dim W=8\), (3) follows. \(\square\)

The same lemma holds for \(W^T\) by transposition.

## 3. No three-letter crossing restriction

### Theorem

Over every field of characteristic two, no four constant physical linear maps
restrict the tensor (1) to a crossing on an alphabet of size three.

**Proof.** A successful map on each physical leg must have rank three, because
the corresponding one-leg flattening of a three-letter crossing has rank
three. Factor the crossing-cut matrix of \(R\) as

\[
 R_{xv\mid yu}=AB,
 \qquad A\in F^{16\times8},\quad B\in F^{8\times16},
\]

where the columns of \(A\) span \(W\) and the rows of \(B\) span \(W^T\).
After the four physical maps, the two factors each retain rank at least five
by the lemma. The rank inequality for a product through an eight-dimensional
space gives

\[
 \mathrm{rank}(A'B')
 \ge\mathrm{rank}(A')+\mathrm{rank}(B')-8
 \ge2.
\]

A three-letter crossing has rank one across this cut, a contradiction.
\(\square\)

## 4. A stronger obstruction when both crossing pairs use the same code

A diagonal transfer with nonzero entries can sometimes be corrected after
grouping copies. For example, a nonzero \(\mathbb F_4\) phase multiplied by its
Frobenius conjugate equals one. The preceding theorem does not exclude such
phase gates, because their crossing-cut rank can exceed one.

When the two crossing pairs use identical input and output maps, however,
even an invertible three-letter phase gate is impossible.

### Corollary

For arbitrary maps \(B:F^3\to F^4\) and \(C:F^4\to F^3\), the nine-dimensional
operator

\[
 (C\otimes C)R(B\otimes B)
\]

is singular over every characteristic-two field.

**Proof.** For a vector space \(V\), let \(\tau\) exchange the factors of
\(V\otimes V\), and define

\[
 \mathcal S(V)=\ker(\tau-I),\qquad
 \mathcal A(V)=\mathrm{im}(\tau+I).
\]

In characteristic two, \(\mathcal A(V)\subseteq\mathcal S(V)\). For
\(\dim V=d\), their dimensions are \(d(d+1)/2\) and \(d(d-1)/2\), respectively.
The quotient has dimension \(d\), with basis given by the diagonal tensor
classes \([e_i\otimes e_i]\).

The operator \(R\) commutes with \(\tau\). On
\(\mathcal S(F^4)/\mathcal A(F^4)\), its induced matrix is

\[
 N=\begin{pmatrix}
 1&0&0&0\\
 0&1&1&0\\
 0&1&1&0\\
 0&0&0&0
 \end{pmatrix},\qquad \mathrm{rank}N=2. \tag{6}
\]

Tensor-square maps induce entrywise Frobenius maps on these quotients. Thus
the induced three-dimensional quotient map of the proposed coded operator is
\(C^{(2)}NB^{(2)}\), which has rank at most two. The coded operator itself
commutes with factor exchange. If it were invertible, its inverse would also
commute with factor exchange, so its maps on \(\mathcal S\), \(\mathcal A\), and
their quotient would all be invertible. The quotient rank bound is a
contradiction. \(\square\)

This corollary places no left-inverse requirement on \(B,C\). A three-letter
phase-code search must therefore use different codes on the two crossing
pairs, use additional bulk modes, or group copies before applying its code.

## Exact controls

```sh
python tools/check_residual_three_letter_obstruction.py \
  --output evidence/residual-three-letter-obstruction.json
```

The checker reconstructs the tensor, its rank-eight factor spaces, and the
transpose-permutation symmetry. It verifies factor-exchange commutation and
the rank-two quotient matrix (6). It also checks the rank-four equation bound
for every ordered pair of projective nonzero vectors over \(\mathbb F_2\)
and \(\mathbb F_4\). These finite controls accompany the all-field proof;
they are not a substitute for its case argument.
