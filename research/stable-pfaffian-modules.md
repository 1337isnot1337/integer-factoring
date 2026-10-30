# A stable two-factor decomposition of the complete Pfaffian

## Result

The six-dimensional extension in the
[particle-hole description](pfaffian-particle-hole.md) splits after
tensoring with one four-dimensional module. An explicit constant binary
change of basis proves

\[
 H\otimes G\ \cong\
 2\Delta G\ \oplus\
 G\otimes g^{[2]}\ \oplus\
 G\otimes h^{[2]}.
 \tag{1}
\]

Here

\[
 \Delta=p+r,\quad
 g=\begin{pmatrix}r&p\\1&1\end{pmatrix},\quad
 h=\begin{pmatrix}1&1\\q&p+q+r\end{pmatrix},\quad
 G=h\otimes g,\quad H=\Delta^{-1}\Lambda^2G,
\]

and \([2]\) denotes the Frobenius twist, obtained by squaring every matrix
entry. All matrices and changes of basis are over characteristic two.

Let

\[
 S=2\Delta\oplus2G\oplus H
\]

be the exact sixteen-dimensional local Pfaffian matrix. Equation (1) gives
the all-size stable identity

\[
 \boxed{
 S^{\otimes P}\otimes G
 \ \cong\
 G\otimes(g\oplus h)^{\otimes2P}.
 }
 \tag{2}
\]

The changes of basis preserve the two aggregate exterior grades. The
right side admits a polynomial-size decomposition using the elementary
two-dimensional Frobenius algorithm already proved in the repository.
This gives computable bounds for **every complete coefficient matrix** of
the full-grade tensor, including bivariate coefficient bases and arbitrary
characteristic-two fields.

The extra \(G\) is not a heuristic regularization: it has entries
identically equal to one. Fixing such an entry embeds the original
coefficient matrix as a literal submatrix of the corresponding augmented
coefficient matrix. There is no convolution of coefficients in this
embedding.

## 1. Grades

The basis grades for \(G\) are

\[
 (1,0),\ (0,1),\ (2,1),\ (1,2).
 \tag{3}
\]

Those for \(H\), in exterior order \(12,13,14,23,24,34\), are

\[
 (0,0),\ (2,0),\ (1,1),\ (1,1),\ (0,2),\ (2,2).
 \tag{4}
\]

Each scalar \(\Delta\) has grade \((1,1)\).
The two entries of \(g^{[2]}\) have grades \((2,0),(0,2)\), and those of
\(h^{[2]}\) have grades \((0,0),(2,2)\).

These are the original parameter and state grades from the particle-hole
map, not newly imposed gradings.

The identity preserves their sums across sites. It does **not** assert
preservation of each site's separate physical energy. Consequently, the
block decomposition below gives full-grade bounds; any fixed-energy
shell or selected code inherits them by taking submatrices. A sharper
energy-sensitive decomposition would require an additional compatibility
argument.

## 2. An explicit splitting of \(H\otimes G\)

Order the source basis as \(H\)-major and \(G\)-minor. Order the target
blocks as

\[
 \Delta G,\quad \Delta G,\quad G\otimes g^{[2]},\quad G\otimes h^{[2]}.
\]

For the last two blocks use \(G\)-major, Frobenius-factor-minor order.
An intertwining matrix has the following columns, encoded as binary
integers with the first source row in the least significant bit:

\[
\begin{split}
(&4352,8704,17408,34816,4132,66056,1049728,2392064,\\
 &16,74240,4384,131072,64,296960,17536,524288,\\
 &1,1065984,2,2131968,4356,4194304,8712,8388608).
\end{split}
\tag{5}
\]

Its binary rank is \(24\). Every nonzero entry connects equal grades.
For every matrix position the polynomial identity

\[
 (H\otimes G)X
 =X\,\operatorname{diag}
 \bigl(\Delta G,\Delta G,G\otimes g^{[2]},G\otimes h^{[2]}\bigr)
 \tag{6}
\]

holds in \(\mathbb F_2[p,q,r]\).

The checker independently constructs \(H\) from the thirty-six exterior
minors of \(G\), performs exact division by \(p+r\) using sparse
multivariate polynomials, and verifies all 576 equations (6). It does not
take the displayed six-dimensional matrix as an input to this check.
Thus (5), (6) form a finite certificate of the module isomorphism.

## 3. Deriving the stable identity

The previously verified elementary identity

\[
 V^{\otimes3}\cong V\otimes V^{[2]}
 \oplus2(\det V\otimes V)
 \tag{7}
\]

applies to \(g\) and \(h\) separately. Because both determinants equal
\(\Delta\),

\[
 G\otimes g^{\otimes2}
 \cong G\otimes g^{[2]}\oplus2\Delta G,
\]

and the same statement holds with \(h\).
Expanding both sides of

\[
 S\otimes G\cong G\otimes(g\oplus h)^{\otimes2}
 \tag{8}
\]

now yields the same direct summands: \(2G^{\otimes2}\), \(4\Delta G\),
\(G\otimes g^{[2]}\), and \(G\otimes h^{[2]}\).
The left side uses (1); the right side uses (7).
All component intertwiners preserve grades, so (8) does too.

Apply (8) repeatedly, moving the same \(G\) factor through the \(P\)
identical sites. This proves (2). The parameters must be common to those
sites, as they are in the universal complementary-shell tensor.

## 4. Explicit polynomial-size blocks

Expand the right side of (2):

\[
 \bigoplus_{a=0}^{2P}
 \binom{2P}{a}\,
 g^{\otimes(a+1)}\otimes h^{\otimes(2P-a+1)}.
 \tag{9}
\]

Use the explicit modules \(B_m\) from
[the elementary Frobenius decomposition](folded-gram-channels.md).
Their weights are recursively

\[
 \mathcal W_0=(),\quad
 \mathcal W_{2r+1}=(1,2\mathcal W_r),\quad
 \mathcal W_{2r+2}=(1,1,2\mathcal W_r).
\]

If the decomposition of \(V^{\otimes a}\) contains \(\det^eB_m\) with
multiplicity \(\mu_a(m,e)\), then the stable multiplicity is

\[
 M_P(m,n,d)=
 \sum_{a=0}^{2P}\binom{2P}{a}
 \sum_{e+f=d}\mu_{a+1}(m,e)\mu_{2P-a+1}(n,f).
 \tag{10}
\]

Every retained term satisfies

\[
 m+n+2d=2P+2.
 \tag{11}
\]

Thus (2) is a direct sum of
\(\Delta^d B_m(g)\otimes B_n(h)\), with multiplicities (10).
The matrices \(B_m,B_n\) have polynomial dimension in \(P\), and all
multiplicities are computed as integers. No modular character identity is
being used as a substitute for a change of basis: (1) and (7) explicitly
construct the direct sums.

## 5. The ancillary coefficient embedding

The second row of \(G\) is \((1,1,1,1)\).
First fix its entry \(G_{1,0}=1\), using zero-based indices.
For original row and column grades \((L,K)\), the augmented grades are

\[
 \text{row }(L,K+1),\qquad
 \text{column }(L+1,K).
 \tag{12}
\]

Each original polynomial matrix entry is unchanged when tensored with
this constant entry. It follows for every linear coefficient functional
\(\ell\) that

\[
 \operatorname{rank}\ell(T)
 \leq
 \operatorname{rank}
 \ell\bigl((T\otimes G)_{(L,K+1),(L+1,K)}\bigr).
 \tag{13}
\]

This applies equally to univariate coefficients, multivariate jets, and
linear combinations of coefficient tensors.

In a block \(\Delta^dB_m(g)\otimes B_n(h)\), put

\[
 u=(K-L+1+m)/2,\qquad
 v=(L+K-2P-1+n)/2.
 \tag{14}
\]

The \(g\) row slice has weight \(u\), its column slice has weight \(u-1\),
and the \(h\) row and column slices both have weight \(v\). If these are
not valid integral weights, the block is empty.

Let \(w\) be the weighted intersection of the two \(g\)-basis subsets,
and \(z\) that of the \(h\)-basis subsets. The exact polynomial entry of
the rectangular block is

\[
 \boxed{
 \Delta^d
 r^{m-2u+1+w}
 p^{u-1-w}
 q^{v-z}
 (p+q+r)^z.
 }
 \tag{15}
\]

All exponents are nonnegative. Summing the exact ranks of the coefficient
matrices (15), with multiplicities (10), gives the bound (13).

### Other constant entries

One may use any of the four entries \(G_{1,c}=1\).
For general fixed ancillary row and column grades
\((a_r,b_r)\), \((a_c,b_c)\), use

\[
 \begin{aligned}
 u_r&=(K+b_r-L-a_r+m)/2,&
 u_c&=(K+b_c-L-a_c+m)/2,\\
 v_r&=(L+K+a_r+b_r-2P-2+n)/2,&
 v_c&=(L+K+a_c+b_c-2P-2+n)/2.
 \end{aligned}
 \tag{16}
\]

The entry becomes

\[
 \Delta^d
 r^{m-u_r-u_c+w}
 p^{u_c-w}
 q^{v_r-z}
 (p+q+r)^z.
 \tag{17}
\]

Each choice is independently a valid upper bound. Their minimum may be
taken separately for each coefficient.

## 6. Original source coefficients

For \(p=t+t^2,\ q=1+t,\ r=t\), formula (15) is

\[
 t^{\,2d+m-u}
 (1+t)^{u+v-1-w-z}
 (1+t+t^2)^z.
 \tag{18}
\]

The coefficient ranks of these small matrices give the complete row-cut
bounds. On a full fixed-grade state domain, complementation preserves the
domain. Exchanging \(p,q\) reverses the homogeneous source polynomial in
degree \(2L\), so the column bound at degree \(k\) is the row bound at
degree \(2L-k\).

These are complete coefficients: no separate source-angle channels are
added afterward.

## 7. Complete bivariate jets at the common singular point

Set

\[
 p=1+x,\qquad q=1+y,\qquad r=1.
 \tag{19}
\]

Evaluation at \(x=y=1\) is the crossing endpoint \(p=q=0,r=1\).
The sum of **all** bivariate coefficient tensors therefore equals the
crossing tensor.

Now \(\Delta=x\), and (17) becomes

\[
 \boxed{
 x^d(1+x)^{u_c-w}(1+y)^{v_r-z}(1+x+y)^z.
 }
 \tag{20}
\]

For the four ancillary entries in row one, its total degree is at most
\(L\). The two row-zero entries used below can give an augmented degree
up to \(L+1\); only coefficients of total degree at most \(L\) are needed
for the original tensor. The implementation extracts every required
coefficient \(x^a y^b\), not only coefficients on an axis or at one
singular specialization.

Since \(r=1\), the additional entries \(G_{0,0}\) and \(G_{0,2}\) are
also constant one. The centered-jet implementation minimizes each
coefficient's bound over all six available ancillary entries.

The other cut exchanges \(x,y\). If \(U_{a,b}\) is the resulting row
upper bound, valid complete costs are

\[
 R_{a,b}\leq\min(d_0^2,U_{a,b}),\qquad
 C_{a,b}\leq\min(d_0^2,U_{b,a}),
 \tag{21}
\]

where \(d_0\) is the chosen common label dimension.
The exact integer cost upper bound is

\[
 S_{\rm int}
 =\sum_{a+b\leq L}
 \left\lceil
 \sqrt{\min(d_0^2,U_{a,b})\min(d_0^2,U_{b,a})}
 \right\rceil.
 \tag{22}
\]

The reported exponent parameter is
\(\rho=\log_2S_{\rm int}/(2\log_2d_0)\).
A sufficient finite seed needs \(\rho<1/\sqrt2\).
The computed screens do not currently meet this threshold. Failing this
upper-bound test is not a lower bound on the true coefficient costs.

## 8. Exact full-grade types before stabilization

The local decomposition \(S=2\Delta\oplus2G\oplus H\) also gives a useful
exact finite algorithm. A type with \(a\) scalar sites, \(b\) odd sites,
and \(c\) middle sites, \(a+b+c=P\), occurs

\[
 \binom Pa\binom{P-a}b\,2^{a+b}
 \tag{23}
\]

times. All copies have the same grade-restricted coefficient matrices.
Consequently, full-grade actual ranks can be calculated with only
\(O(P^2)\) types, although a pure \(H\) type may still be exponentially
large.

This type algorithm independently reproduces the signed-charge full
profiles. Stabilization then removes the large \(H\) type from the
upper-bound calculation by adding just one ancillary \(G\).

## Verification and scope

The supplied checker verifies:

- the explicit graded \(24\)-dimensional splitting with sparse
  multivariate arithmetic;
- seven exact augmented coefficient profiles, constructed before the
  stable decomposition, including unequal operator and state grades;
- complete actual centered-jet rank profiles against the new bounds;
- exact full-grade type profiles against the independent signed-charge
  calculation.

These checks support the constructive all-size proof. They do not replace
the proof or establish a new factoring exponent.

~~~bash
python3 tools/check_stable_pfaffian_bounds.py \
  --output evidence/stable-pfaffian-bounds.json
~~~

Implementation:
[stable_pfaffian_bounds.py](../tools/stable_pfaffian_bounds.py) and
[pfaffian_grade_types.py](../tools/pfaffian_grade_types.py).
Evidence:
[stable-pfaffian-bounds.json](../evidence/stable-pfaffian-bounds.json).
