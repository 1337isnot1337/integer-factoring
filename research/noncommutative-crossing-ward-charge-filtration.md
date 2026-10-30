# Charge filtration for the Ward rank factors

This note combines the determinant interface from
[the Euler Ward refinement](noncommutative-crossing-ward-euler.md)
with a constant Pfaffian charge gauge. The initial rank bound has
only \(\binom{k+3}{3}\) rectangular terms per cut. Section 6 sharpens
this to \(\binom{k+2}{2}\) terms by an exact descent through either
state leg. Both bounds have no factor
\(P^k\), no independent Clifford-channel count, and no dependence on
the dimension of the quadratic-charge space.

The charge arguments use binary forms over characteristic-two fields.
The direct descent in Section 6 permits arbitrary form weights in
any characteristic-two field.
The code is homogeneous of exterior grade \(L\), with
\(I=Q_1\cdots Q_kW\) and \(\deg W=L-2k\). The output map may be any
grade-\(L\) linear map. Syndrome-preserving outputs are needed only for
the additional joint-order corollary.

The [computed refined bounds](ward-charge-bound-computations.md) record
the current numerical estimates, exact small controls, and reproduction
commands. No improved factoring exponent follows from those estimates.

## 1. Symmetric determinant translation

For symmetric \(n\times n\) matrices \(X,Y\) over a commutative
characteristic-two ring,
\[
 \boxed{
 \det(X+Y)=\sum_{S\subseteq[n]}
       \det(X_{S,S})\det(Y_{S^c,S^c}).
 }
 \tag{1}
\]
The empty determinant is one.

To prove this, pair each determinant permutation containing a cycle
of length at least three with the permutation reversing its first
such cycle. Symmetry makes the terms equal, and characteristic two
cancels the pair. Only involutions remain. A fixed point contributes
a diagonal entry; a transposition contributes the square of an
off-diagonal entry. In \(X+Y\), every such square splits as
\(x_{ij}^2+y_{ij}^2\). Assign each fixed point and transposition to
one of the two matrices. The vertices assigned to \(X\) form \(S\),
giving exactly (1).

Consequently the translation kernel \(\det(C+X+Y)\) has separation
rank at most \(2^n\), regardless of the dimension of the parameter
space describing \(X,Y\).

Let \(V\) be a binary vector space of symmetric charge matrices, and
let \(\mathfrak a\) be the augmentation ideal of \(F[V]\).
A Boolean polynomial of degree at most \(d\) annihilates
\(\mathfrak a^{d+1}\) under the pairing
\[
 \left\langle f,\sum_Q a_Qe_Q\right\rangle
       =\sum_Qa_Qf(Q).
 \tag{2}
\]
Indeed, a product \(\prod_{i=1}^s(1+e_{Q_i})\) evaluates the
\(s\)-fold finite difference of \(f\). Every difference lowers
Boolean degree.

A principal determinant of size \(h\) is a Boolean polynomial of
degree at most \(h\) in the charge coordinates. Restricting the
translation kernel to
\(\mathfrak a^a\times\mathfrak a^b\) therefore gives
\[
 \operatorname{rank}\le
 \sum_{h=a}^{n-b}\binom nh.
 \tag{3}
\]

There is a useful combined version for commuting external operators
\(O_1,\ldots,O_k\):
\[
 \begin{aligned}
 &[z^j]\det(C+X+Y+z\operatorname{diag}(O))\\
 &\quad=\sum_{S\subseteq[k]}
 [z^j]\det((C+X+z\operatorname{diag}(O))_{S,S})
             \det(Y_{S^c,S^c}).
 \end{aligned}
 \tag{4}
\]
The first factor has charge degree at most \(|S|-j\).
Thus the entire \(j\)-th operator-valued interface, without separating
its \(O\)-subsets, has separation rank at most
\[
 \boxed{\quad
 B(k,j,a,b)=\sum_{h=a+j}^{k-b}\binom kh.
 \quad}
 \tag{5}
\]
An empty sum is zero. The bound permits a constant affine offset
\(C\).

## 2. A charge gauge that acts on one side

For one physical pair, let \(e^2=1\) in a commutative
characteristic-two ring. The charged centered parameters are
\[
 p_e=1+xe,\qquad q_e=1+ye,\qquad r=1.
\]
Write
\[
 g=\begin{pmatrix}1&1+x\\1&1\end{pmatrix},\qquad
 h=\begin{pmatrix}1&1\\1+y&1+x+y\end{pmatrix}.
\]
The constant matrices
\[
 T_g(e)=\begin{pmatrix}e&1+e\\0&1\end{pmatrix},\qquad
 T_h(e)=\begin{pmatrix}1&0\\1+e&e\end{pmatrix}
\]
satisfy \(g_e=T_g(e)g\), \(h_e=T_h(e)h\).

The particle-hole representation identifies the local Pfaffian matrix
with its exterior blocks \(2\Delta\oplus2G\oplus H\), where
\(G=h\otimes g\), \(H=\Lambda^2G/\Delta\), and \(\Delta=x\).
Let \(T(e)\) act on exterior degree \(d\) by
\[
 T(e)|_{\Lambda^d}
   = e^{(d+1)\bmod2}\Lambda^d(T_h(e)\otimes T_g(e)).
 \tag{6}
\]
Then
\[
 \boxed{\quad S_e=T(e)S.\quad}
 \tag{7}
\]
The apparent quotient in \(H\) cancels polynomially. Equation (7)
is a polynomial identity, including at \(x=0\).

In the original external/state grades, every nonzero entry of \(T(e)\)
has grade shift
\[
 (0,0),\qquad (1,1),\qquad (1,-1).
 \tag{8}
\]
Every entry with positive first shift is divisible by \(1+e\).
All entries are \(0,1,e,\) or \(1+e\).
More precisely, the entries with shift \((0,0)\) are diagonal and
equal to \(1\) or \(e\); every positive-shift entry equals \(1+e\).
These assertions also follow directly by expanding the exterior
minors in (6). The independent checker verifies every local
polynomial entry and grade shift.

At pair \(i\), take \(e=e_{a_ia_i^{\mathsf T}}\) in the group algebra
of the quadratic-charge span. Tensoring (7) over the \(P\) pairs gives
\[
 \mathscr M_{\rm colored}=T_{\rm all}\mathscr M.
 \tag{9}
\]
A matrix entry of \(T_{\rm all}\) raising external grade by \(s\)
belongs to \(\mathfrak a^s\). Its state-grade change \(t\) satisfies
\[
 |t|\le s,\qquad t\equiv s\pmod2.
 \tag{10}
\]
There is no restriction on the charge-space dimension in this
argument.

## 3. Apply the Ward filter before taking rank

Fix the \(j\)-th Ward term and put \(\ell=k-j\). Its charge filters
are principal determinants of order \(\ell\):
\[
 F_J(Q)=\det(C_A+Q)_{J^c,J^c},
 \qquad |J|=j.
 \tag{11}
\]
Here \(C_A\) is the quadratic charge of the target external monomial.
It may depend on that monomial. This causes no problem: each
row-dependent filter still has Boolean degree at most \(\ell\).

Apply (2) entrywise to the gauge in (9). Every component with
\(s>\ell\) vanishes. For fixed \(s,t\), sum all \(J\), their external
creation maps \(\prod_{a\in J}O_a\), and their charge evaluations.
They form **one constant map on the row side** of the relevant
flattening. Thus rank subadditivity is needed only for the possible
grade shifts, not for the charge values or the \(O\)-subsets.

The external creation maps add full pairs and preserve \(C_A\).
Therefore they may be composed with the row-dependent filters in
(11). A final arbitrary output projection is also a one-sided map.

This proves the following all-size bounds. Put \(q=b-k-j\).
For negative coefficient indices or invalid grades, the corresponding
matrix is zero.

**Theorem.**
\[
 \boxed{
 \begin{aligned}
 r(G[a,b])\le
 \sum_{j=0}^{\min(k,b-k)}
 \sum_{s=0}^{k-j}
 \sum_{\substack{-s\le t\le s\\t\equiv s\ (2)}}
 r\!\left(
 M_P(L-2j-s,L-2k-t;L,L)[a,b-k-j]
 \right).
 \end{aligned}}
 \tag{12}
\]
For the column cut, use the exact \(p\leftrightarrow q\) reshaping
symmetry. The external \(A\) leg and its charge are unchanged, and the
within-pair mode complements preserve all physical grades and
charges. Applying the same one-sided gauge gives
\[
 \boxed{
 \begin{aligned}
 c(G[a,b])\le
 \sum_{j=0}^{\min(k,b-k)}
 \sum_{s=0}^{k-j}
 \sum_{\substack{-s\le t\le s\\t\equiv s\ (2)}}
 r\!\left(
 M_P(L-2j-s,L-t;L,L-2k)[b-k-j,a]
 \right).
 \end{aligned}}
 \tag{13}
\]
No factor \(B(k,j,s,t)\) is needed in these one-sided bounds.
Equation (5) remains available for a two-sided gauge or another
interface decomposition.

The total number of possible \((j,s,t)\) is
\[
 \sum_{j=0}^k\sum_{s=0}^{k-j}(s+1)
       =\binom{k+3}{3}.
\]
One may combine (12)--(13) with every earlier valid bound and the
dimension cap \(d^2\). For one global form, the exact parity-filtered
two-rectangle theorem can be sharper.

## 4. Joint vanishing and constructive entry evaluation

For syndrome-preserving outputs on the admissible homogeneous code,
\(C_A=C_{\rm total}=AA^{\mathsf T}\). The Euler determinant proof gives
\[
 b\ge k,\qquad a+b\ge2k-\operatorname{rank}(AA^{\mathsf T})
 \tag{14}
\]
for every nonzero centered coefficient. The second assertion need not
hold for arbitrary output perturbations.

The constant maps in (12)--(13) need not be evaluated by enumerating
the charge group. A local gauge entry is a unit charge times a product
of factors \(1+e_{uu^{\mathsf T}}\). Pairing a determinant with such a
product is a finite difference along rank-one symmetric updates.

For a symmetric matrix \(B\) and a matrix \(U\) whose columns are
\(u_1,\ldots,u_m\),
\[
 \boxed{
 \sum_{S\subseteq[m]}
 \det\!\left(B+\sum_{i\in S}u_iu_i^{\mathsf T}\right)
 =
 \det\begin{pmatrix}B&U\\U^{\mathsf T}&0\end{pmatrix}.
 }
 \tag{15}
\]
Each rank-one update enters the determinant with degree at most one.
The left side extracts the coefficient of their full product. The
bordered determinant gives that coefficient by the Schur-complement
identity when \(B\) is formally invertible, hence as a polynomial
identity for every \(B\). All signs disappear in characteristic two.

External full-pair insertions fit the same interface. Write
\(Z_i=A_iA_{\bar i}\), so \(Z_i^2=0\). For binary form columns,
\[
 \det(C+z\operatorname{diag}(O_a))
 =
 \det\!\left(C+z\sum_i Z_i a_ia_i^{\mathsf T}\right).
 \tag{16}
\]
The diagonals agree. The off-diagonal difference has square zero,
and the involution expansion used in (1) shows that this does not
change the determinant. Extracting a prescribed set of inserted
pairs is therefore another coefficient of rank-one updates.
Together with the gauge differences, it is evaluated by one bordered
determinant of size at most \(k+j+s\).
There are exactly \(s\) difference factors: the zero-shift gauge
entries contribute only a unit charge. For fixed source and target
monomials, an external-grade increase of two at a physical pair must
come from a full-pair insertion, while an increase of one must come
from the gauge. These choices cannot overlap at a pair, since its
external grade is at most two. Thus a fixed interface entry does not
require a further sum over assignments of insertions to sites.

Thus the interface admits direct polynomial-size determinant
evaluation, rather than an exponential sum over charge values.
The scalar fixed-leg maps and source rectangular factors must still
be supplied by the chosen crossing construction.

## 5. Verification and scope

The local gauge identities, their invertibility, and all grade shifts
are checked in
[check_pfaffian_charge_gauges.py](../tools/check_pfaffian_charge_gauges.py),
with evidence in
[pfaffian-charge-gauges.json](../evidence/pfaffian-charge-gauges.json).

The Euler and typed determinant identities are checked in
[check_noncommutative_crossing_ward_euler.py](../tools/check_noncommutative_crossing_ward_euler.py).
The independent
[charge-filtration checker](../tools/check_ward_charge_filtration.py)
verifies (1) as a formal polynomial identity for sizes one through six,
with repeated variable powers retained. It also checks 800 exact
rank-one bordered-determinant instances, 830 operator-coefficient
augmentation-rank instances, and both-cut bounds against literal
global-form and overlapping-form code tensors. The recorded output is
[ward-charge-filtration.json](../evidence/ward-charge-filtration.json).

The implemented bound evaluators are
[ward_differential_bounds.py](../tools/ward_differential_bounds.py)
for the global-form refinement and
[ward_channel_free_bounds.py](../tools/ward_channel_free_bounds.py)
for the one-sided charge formula. Both retain the minimum with earlier
verified rectangular bounds; the latter also uses (14) when the
syndrome-preserving hypothesis is supplied.

Equations (12)--(13) apply to the homogeneous product-image input
codes and arbitrary grade-preserving output maps. They do not
require that the code is preserved by a stable-module intertwiner.
Their use with an explicit rectangular rank engine is a separate,
checkable step. No score improvement follows until those resulting
bounds meet the seed criterion.

## 6. Exact descent through either state leg

The smaller grade shifts in (12)--(13) vanish identically. There are
two independent ways to prove this refinement. The first does not
need binary weights or charge coordinates.

### Direct dual-incidence identity

Let \(c_i^U=a_i^U+b_{\bar i}^U\) be the same creation/annihilation
Clifford operator, acting on the \(U\) leg. Its first-order identity is
\[
 c_i^U\mathscr M
   =(U_i+p_iA_{\bar i}+Y_{\bar i}+V_i)\mathscr M.
 \tag{17}
\]
Therefore the dual incidence operator
\[
 \widetilde{\mathcal L}_a
   =\sum_i w_{ai}(A_ic_i^U+A_{\bar i}c_{\bar i}^U)
\]
satisfies
\[
 \boxed{\quad
 \widetilde{\mathcal L}_a\mathscr M
   =D_a\mathscr M=\mathcal L_a\mathscr M.
 \quad}
 \tag{18}
\]
Indeed, multiplying the two instances of (17) by the corresponding
external variables cancels the two \(p_iA_iA_{\bar i}\) terms. The
remaining insertion is precisely
\(\partial_{p_i}+\partial_{q_i}+\partial_{r_i}\).

Both incidence operators are parameter-independent and commute with
all parameter derivatives. Induction thus gives, for every ordered
list of form indices,
\[
 \mathcal L_{a_1}\cdots\mathcal L_{a_m}\mathscr M
 =D_{a_1}\cdots D_{a_m}\mathscr M
 =\widetilde{\mathcal L}_{a_1}\cdots
    \widetilde{\mathcal L}_{a_m}\mathscr M.
 \tag{19}
\]
This argument does not assume that the incidence operators commute
as abstract operators. Their commutation with the derivatives is
sufficient. The external creation operators \(\mathcal O_a\) commute
with these transfers and stay on the external leg.

Use the original incidence operators for the row cut and the dual
operators for the column cut. Both act entirely on one side of their
respective cut. In the \(j\)-th Ward term, exactly \(k-j\) incidence
operators raise the external grade by one, while \(j\) full-pair
insertions raise it by two. Put \(\ell=k-j\).

**Exact descent theorem.** For arbitrary form weights in a
characteristic-two field, homogeneous product-image inputs, and any
grade-\(L\) output map,
\[
 \boxed{
 \begin{aligned}
 r(G[a,b])\le{}&
 \sum_{j=0}^{\min(k,b-k)}
 \sum_{\substack{-\ell\le t\le\ell\\t\equiv\ell\ (2)}}
 r\!\left(M_P(L-k-j,L-2k-t;L,L)[a,b-k-j]\right),\\
 c(G[a,b])\le{}&
 \sum_{j=0}^{\min(k,b-k)}
 \sum_{\substack{-\ell\le t\le\ell\\t\equiv\ell\ (2)}}
 r\!\left(M_P(L-k-j,L-t;L,L-2k)[b-k-j,a]\right).
 \end{aligned}}
 \tag{20}
\]
The second line uses the same grade-preserving \(p\leftrightarrow q\)
reshape as (13). Its preimage leg remains of grade \(L-2k\) on the
opposite side. There is no channel multiplier. The total number of
rectangles per cut is at most
\(\sum_{j=0}^k(k-j+1)=\binom{k+2}{2}\).

### Independent charge proof and explicit minors

For binary weights, the gauge gives an independent explanation. Its
diagonal entry is \(e\) exactly when the external occupation is
single, and is \(1\) when it is empty or full. This statement is
independent of the state occupation.

Fix a gauge entry with positive-shift site set \(S\), of size \(s\).
The product of its diagonal unit charges agrees with \(C_A\) outside
\(S\). Consequently every matrix tested by its Ward filter is a sum
of at most \(s\) rank-one charges from these sites. If \(s<\ell\),
every principal minor of order \(\ell\) vanishes by rank. If
\(s>\ell\), the augmentation argument already gives zero.
Thus only \(s=\ell\) survives.

For a surviving row subset \(J^c\), the finite difference is exactly
\[
 \sum_{\delta\in\mathbf F_2^{S}}
 \det\!\left(
 A_{J^c,S}\operatorname{diag}(\epsilon+\delta)
 A_{J^c,S}^{\mathsf T}\right)
 =\det(A_{J^c,S})^2
 =\det(A_{J^c,S}).
 \tag{21}
\]
The equality is independent of the affine parity vector
\(\epsilon\). The last equality uses binary entries.

Let \(T\) be the \(j\) full-pair insertion sites. The sites of \(T\)
and \(S\) are disjoint, since a positive gauge shift followed by a
full-pair insertion would exceed the local external grade two.
Summing the creation coefficient \(\det A_{J,T}\) against (21) is
the Laplace expansion
\[
 \sum_{|J|=j}\det(A_{J,T})\det(A_{J^c,S})
     =\det(A_{[k],T\cup S}).
 \tag{22}
\]
Thus each fixed entry of the combined constant interface is given by
one \(k\times k\) minor. The prescribed source and target occupations
determine the insertion and shift sites; no enumeration of charge
values or subsets of form rows is needed to evaluate that entry.

The independent
[dual-incidence checker](../tools/check_ward_dual_incidence.py)
verifies all 512 generic first-order identities, all 256 local
dual-incidence derivative identities, 2,560 multisite product
identities, the complete gauge support, 1,320 maximal-shift finite
differences, and 2,240 combined-minor identities. The exact evidence is
[ward-dual-incidence.json](../evidence/ward-dual-incidence.json).
The literal code checks in the charge-filtration checker also use
the strengthened bound (20).

