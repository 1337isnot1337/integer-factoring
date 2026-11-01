# Class-two group-word crossings with arbitrary central controls

## Result and scope

All full-domain class-two group-word families studied here satisfy the
binary conservation rank-sum obstruction, including arbitrary central
lookup corrections.

Let
\[
G=\mathbf F_2^n\times\mathbf F_2^t,\qquad
(h,c)(H,C)=(h+H,c+C+\beta(h,H)),
\tag{1}
\]
where \(\beta\) is any bilinear map. Use the full group as the label
domain, any fixed order of its binary normal-form generators, arbitrary
positive integer deformation weights, and any central control
\[
(h,c)\mapsto(h,c+f(h)),
\tag{2}
\]
where \(f:\mathbf F_2^n\to\mathbf F_2^t\) is arbitrary. If tensors
from the complete coefficient span sum exactly to crossing, then
over every field of characteristic two,
\[
\boxed{\sum_j\sqrt{r(T_j)c(T_j)}\ge |G|^{3/2}.}
\tag{3}
\]

The conclusion still holds when coefficient directions from several
central controls, factor orders, or weight assignments are mixed.
The proof has one common degeneration for all of them.

It closes the bilinear, quadratic, and arbitrary central-lookup
Heisenberg and coded-central families for the rank-sum seed criterion.
Arbitrary control permutations mixing horizontal and central coordinates,
arbitrary orbit kernels, arbitrary label restrictions or retracts, and
groups not of the form (1) remain outside the theorem.

The last step is the pinned
[full-conservation theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/full-conservation-rank-sum-barrier.md),
source blob 1b3e426458157f51a3833476d24e9d4bb3d53803.
A related primitive-algebra argument is in
[noncommutative-crossing-pbw-degeneration.md](noncommutative-crossing-pbw-degeneration.md).
No Lean formalization is used.

## 1. Family, endpoint, and Hopf convention

List the horizontal generators \(g_i=(e_i,0)\) first and the central
generators \(c_j=(0,e_j)\) last, for a total of \(\lambda=n+t\).
Every element has a unique binary normal word
\[
g(s)=g_1^{s_1}\cdots g_\lambda^{s_\lambda}.
\tag{4}
\]
Its central group coordinate differs from its central word bits by a
quadratic function of the horizontal bits; adding \(f(h)\) has the
same meaning in either coordinate system.

Let \(U_f\) be the linear extension of (2). Fix a permutation
\(\sigma\) and positive integer weights \(w_i\), and extend linearly
\[
\Phi_z(g(s))=
\prod_{i\ \mathrm{in\ order}\ \sigma}
\bigl(g_i+z^{w_i}(1+g_i)\bigr)^{s_i}.
\tag{5}
\]
Define
\[
K_z(S,a)=[a]\,(\Phi_zU_f)(S),\qquad
T_z(x,y,u,v)=[xy^{-1}=vu^{-1}]K_z(xy^{-1},y^{-1}u).
\tag{6}
\]
At \(z=1\), every factor in (5) is the group identity. Thus
\(K_1(S,a)=[a=e]\) and \(T_1=[x=v][y=u]\).

Write \(\mathcal S\) for the group-algebra antipode. For a linear
map \(\psi\), the cocommutative Hopf expression is
\[
R_\psi(u\otimes v)=\sum(v_1a_1)\otimes(u_1a_2),\qquad
\Delta(a)=\Delta\bigl(\psi(u_2\mathcal S(v_2))\bigr).
\tag{7}
\]
For convention (6), the correct choice is
\[
\psi_z=\mathcal S\,\Phi_z\,U_f\,\mathcal S.
\tag{8}
\]
Indeed, a group-basis term \(b=a^{-1}\) produces
\((va^{-1})\otimes(ua^{-1})\), with coefficient
\([a]\,(\Phi_zU_f)(vu^{-1})\). This is exactly (6).
No Hopf-endomorphism assumption on \(\psi_z\) is needed.

## 2. Make every central correction strictly heavier

Put \(X_i=1+g_i\) for horizontal generators and \(Z_j=1+c_j\)
for central generators. Their ordered squarefree monomials form a
basis by the invertible binary subset transform.

Choose
\[
M>\max(2,n),
\tag{9}
\]
give every \(X_i\) weight one, and every \(Z_j\) weight \(M\).
These are auxiliary filtration weights, separate from \(w_i\).
Let \(F_{\ge a}\) be the span of monomials of weight at least \(a\).

For \(v\in\mathbf F_2^t\), set
\(\zeta_v=\prod_j(1+Z_j)^{v_j}+1\). Every term of \(\zeta_v\)
has weight at least \(M\). The exact relations are
\[
\begin{aligned}
X_i^2&=\zeta_{\beta(e_i,e_i)},\\
X_jX_i+X_iX_j
 &=(1+X_i)(1+X_j)
   \zeta_{\beta(e_i,e_j)+\beta(e_j,e_i)}\quad(i<j),\\
Z_j^2&=0,\qquad [Z_j,X_i]=[Z_j,Z_k]=0.
\end{aligned}
\tag{10}
\]
All commutator and square corrections strictly raise weight because
\(M>2\). Normal-order multiplication never decreases weight.

The associated graded multiplication is therefore
\[
K[x_1,\ldots,x_n,z_1,\ldots,z_t]/(x_i^2,z_j^2).
\tag{11}
\]
The squarefree basis symbols remain independent, so no dimension is
lost and no additional relation is hidden.

Also
\[
\Delta(X_i)=X_i\otimes1+1\otimes X_i+X_i\otimes X_i,
\tag{12}
\]
with the analogous formula for \(Z_j\). Its last term has larger
total weight. Thus all graded generators are primitive.

The antipode is filtered with identity leading map. Its geometric
expansion on \(X_i\) starts with \(X_i\), even when \(g_i\) has
order four; all later terms have higher weight by (10). It reverses
products, whose leading multiplication is commutative.

## 3. Every central lookup disappears in the leading map

For an augmentation monomial with support \(S\), expansion in (5)
gives
\[
\Phi_z(Y^S)=
\left(\prod_{i\in S}(1+z^{w_i})\right)
\prod_{\substack{i\in S\\\mathrm{in\ order}\ \sigma}}Y_i.
\tag{13}
\]
Each selected factor contributes
\(1+(g_i+z^{w_i}(1+g_i))=(1+z^{w_i})(1+g_i)\).
By (10), its leading map is diagonal and independent of order:
\[
\overline\Phi_z(Y^S)
=\left(\prod_{i\in S}(1+z^{w_i})\right)Y^S.
\tag{14}
\]

The dual basis to augmentation monomials consists of the Boolean
coordinate monomials in the word bits. This follows from
\[
\sum_{T\subseteq S}\prod_{i\in A}T_i=[A=S]
\tag{15}
\]
in characteristic two. The dual control is
\[
U_f^*:h_i\mapsto h_i,\qquad c_j\mapsto c_j+f_j(h).
\tag{16}
\]

Every function on \(\mathbf F_2^n\) has a squarefree Boolean
polynomial representation of degree at most \(n\). Every correction
in (16) therefore has weight at most \(n<M\), while \(c_j\) has
weight \(M\). Replacing a central coordinate strictly lowers dual
weight, and Boolean reductions only lower it further. Hence \(U_f\)
never lowers primal weight and
\[
\boxed{\overline U_f=I.}
\tag{17}
\]
The choice (9) works for all central lookup tables simultaneously.

## 4. One four-leg degeneration

Let \(D_s\) scale a basis monomial of weight \(a\) by \(s^a\),
for a fresh indeterminate independent of \(z\) and all coefficient
scalars. Every operation in (7)–(8) preserves or raises total weight.
Thus
\[
(D_s\otimes D_s)T_z(D_s^{-1}\otimes D_s^{-1})
\tag{18}
\]
is polynomial in \(s\): its entry exponent is output weight minus
input weight.

At \(s=0\), (18) is (7) in the primitive commutative algebra
(11), with diagonal map (14). Both antipodes and \(U_f\) have
identity leading map.

For nonzero \(s\), (18) is an invertible change on each of the
four tensor legs, preserving both cut ranks. Specialization at
zero cannot increase either rank. This holds coefficientwise and
for arbitrary scalar recombinations.

The degeneration fixes crossing, since swap commutes with equal
changes on its two tensor factors. Therefore an exact decomposition
gives
\[
\sum_j\widehat T_j=E_d,\qquad
r(T_j)\ge r(\widehat T_j),\quad
c(T_j)\ge c(\widehat T_j).
\tag{19}
\]

## 5. Final binary conservation and its consequence

Put \(q_i=1+z^{w_i}\). The last leading operator fixes \(00,11\)
on coordinate \(i\), and sends
\[
10\mapsto q_i10+(1+q_i)01,\qquad
01\mapsto(1+q_i)10+q_i01.
\tag{20}
\]
Its local kernel, with \(S=x+y=u+v\) and \(b=x+u\), is
\[
\begin{array}{c|cc}
 &b=0&b=1\\ \hline
S=0&1&0\\
S=1&q_i&1+q_i
\end{array}.
\tag{21}
\]

Every leading coefficient, including mixtures from several
admissible families, belongs to the full conservation space
\[
\widehat T(x,y,u,v)=[x+y=u+v]K(x+y,x+u)
\tag{22}
\]
on \(\mathbf F_2^\lambda\). The pinned theorem applied to (19)
gives
\(\sum_j\sqrt{r(\widehat T_j)c(\widehat T_j)}\ge d^{3/2}\),
proving (3).

This permits arbitrarily many pieces, arbitrary coefficient
recombinations, and arbitrary characteristic-two extension fields.
There is no repeated-weight or few-weight-class assumption.

Since \(3/2>\sqrt2\), these families cannot satisfy the finite
seed criterion \(S<d^{\sqrt2}\). With subexponentially many pieces
in a growing family, the normalized width also obeys
\(\rho\ge3/4\).

## 6. Relation to the earlier two-stage proof

The earlier proof used horizontal weight one and central weight
two, then a second ordinary PBW-length filtration. Its first stage
required a degree-two bound on the central control.

The cubic control \(c\mapsto c+h_1h_2h_3\) sends
\(X_1X_2X_3\) to a vector containing \(Z\). This really does
lower weight from three to two in that old filtration. It is not
an obstruction to the present proof: central weight four makes
the same control filtered with identity leading map.

The checker retains the six legacy two-stage cases and this
old-weight-two witness, then adds direct heavier-weight cases.
The all-lookup theorem follows from Sections 2–5, rather than
from extrapolation of the degree-two proof.

## 7. Remaining scope boundaries

Arbitrary control permutations can mix horizontal and central
coordinates. Such a map may send a dual variable of weight one
to one of weight \(M\), violating (16)–(17). Searches over these
noncentral permutations remain outside this theorem.

Arbitrary monomial subsets or linear retracts may not commute
with (18). The full-domain conservation theorem does not apply
to them automatically. Nor does this proof cover arbitrary orbit
kernels \(K(S,a)\) or groups outside (1).

The ordered form (5), a fixed factor order within each family,
the central translation form (2), and a full normal-form domain
are the stated hypotheses.

## 8. Reproduction

The checker
[check_noncommutative_crossing_group_degeneration.py](../tools/check_noncommutative_crossing_group_degeneration.py)
writes
[noncommutative-crossing-group-degeneration.json](../evidence/noncommutative-crossing-group-degeneration.json).

Its direct cases include a cubic control on sixteen labels with
central weight four; a full degree-four lookup on
\(H(\mathbf F_4)\), sixty-four labels, central weight five, and
reordered factors; and the quaternion group with central weight
three, testing nonzero generator squares.

For each direct case it verifies every leading multiplication and
antipode entry, identity leading control, the complete tensor
degeneration, equality with abelian mixing in every coefficient,
and both-cut rank domination. It also checks selected or exhaustive
recombinations.

The degree-four case has 1,583,440 positive-filtration coefficient
entries. Its leading tensor comparison is a full matrix equality.
The six legacy two-stage cases are retained as independent checks.

The ordinary proof establishes all-size scope; the computations
audit conventions and implementation.

The [saved eight-element-field computation](../evidence/noncommutative-crossing-F8-bilinear-checkpoints.json)
contains 32 million completed normalized bilinear cases in the exact interval
`[13800000,45800000)`, with commands, source hash, and per-chunk digests.
All screened ranks are at least 197; larger ranks were stopped at that cap.
This bounded computation is not a complete 67,108,864-case enumeration.
The [resumable runner](../tools/run_group_correction_screen.py) retains the
work, while the stronger theorem above now excludes the entire central
lookup route for the seed criterion.
