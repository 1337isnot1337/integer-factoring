# Two-stage obstruction for class-two group-word crossings

## Result and scope

The full-domain ordered group-word crossings on binary class-two groups
satisfy the binary conservation rank-sum obstruction whenever their central
control correction has Boolean degree at most two.

Let
\[
G=\mathbf F_2^n\times\mathbf F_2^t,\qquad
(h,c)(H,C)=(h+H,c+C+\beta(h,H)),
\tag{1}
\]
where \(\beta\) is any bilinear map. Use the full group as the label
domain, any fixed order of its binary normal-form generators, arbitrary
positive integer deformation weights, and a central control shear of
Boolean degree at most two. If tensors from the complete coefficient span
sum exactly to crossing, then over every field of characteristic two,
\[
\boxed{\quad
\sum_j\sqrt{r(T_j)c(T_j)}\ge |G|^{3/2}.
\quad}
\tag{2}
\]

This includes all bilinear and quadratic central corrections in the
Heisenberg-field and coded-central families, fixed changes of factor order,
and arbitrary scalar recombinations. The bilinear cocycle may have
nonzero diagonal, so horizontal generators may have order four.

The proof uses two common four-leg degenerations. The first replaces the
group algebra by an explicitly presented graded restricted Hopf algebra.
The second replaces that algebra by a commutative square-zero primitive
algebra. Both preserve crossing and can only decrease either cut rank.

Arbitrary monomial subsets, arbitrary linear retracts, higher-degree
control shears, arbitrary orbit kernels on a nonabelian group, and groups
not of the form (1) are outside this theorem.

The last step uses the pinned
[full-conservation theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/full-conservation-rank-sum-barrier.md),
whose source blob is
1b3e426458157f51a3833476d24e9d4bb3d53803.
The separate primitive-algebra degeneration is proved in
[noncommutative-crossing-pbw-degeneration.md](noncommutative-crossing-pbw-degeneration.md).
This note proves the additional group step and precise central-shear
hypotheses directly, without invoking a general Jennings theorem.

## 1. Ordered controls and the exact endpoint

Let \(g_i=(e_i,0)\), \(1\le i\le n\), and let \(c_j=(0,e_j)\),
\(1\le j\le t\). List these generators as
\(g_1,\ldots,g_\lambda\), with the central ones last and
\(\lambda=n+t\). Every group element has a unique word
\[
g(s)=g_1^{s_1}\cdots g_\lambda^{s_\lambda},
\qquad s\in\mathbf F_2^\lambda.
\tag{3}
\]

Its central coordinate in (1) differs from its last \(t\) word bits
by a quadratic function of the horizontal bits. Adding a fixed function
of the horizontal bits to the central word bits therefore adds that same
function to the central coordinate in (1).

Let \(f:\mathbf F_2^n\to\mathbf F_2^t\) have Boolean algebraic
normal form of degree at most two. Define the permutation
\[
\theta_f(g(h,c))=g(h,c+f(h)),
\tag{4}
\]
where \(h,c\) now denote word bits. Extend it linearly to \(U_f\).

Fix a permutation \(\sigma\) of the generator positions and positive
integer weights \(w_i\). Define, and extend linearly,
\[
\Phi_z(g(s))=
\prod_{i\ \mathrm{in\ order}\ \sigma}
  \bigl(g_i+z^{w_i}(1+g_i)\bigr)^{s_i}.
\tag{5}
\]
Set
\[
K_z(S,a)=[a]\,(\Phi_zU_f)(S),
\qquad
T_z(x,y,u,v)
=[xy^{-1}=vu^{-1}]K_z(xy^{-1},y^{-1}u).
\tag{6}
\]

At \(z=1\), each selected factor in (5) is the identity. Hence
\[
K_1(S,a)=[a=e],\qquad
T_1(x,y,u,v)=[x=v][y=u].
\tag{7}
\]
This remains true for affine terms in \(f\), any factor order, and
generators of order four. There is no division by the group order.

## 2. The correct Hopf formula

Let \(\mathcal S\) be the group-algebra antipode. For a linear map
\(\psi\) on a cocommutative Hopf algebra define
\[
R_\psi(u\otimes v)=
\sum(v_1a_1)\otimes(u_1a_2),
\qquad
\Delta(a)=\Delta\bigl(\psi(u_2\mathcal S(v_2))\bigr).
\tag{8}
\]
The map \(\psi\) need not be a Hopf endomorphism.

For the orbit convention (6), the correct map is
\[
\boxed{\psi_z=\mathcal S\,\Phi_z\,U_f\,\mathcal S.}
\tag{9}
\]
Indeed, on a group-basis input, (8) sends a term \(b\) of
\(\psi_z(uv^{-1})\) to \((vb)\otimes(ub)\). Writing \(b=a^{-1}\),
its coefficient is \([a]\,(\Phi_zU_f)(vu^{-1})\), and its output is
\(x=va^{-1},\,y=ua^{-1}\), exactly as in (6).

Omitting either antipode in (9) gives a different convention on a
nonabelian group. The checker independently constructs both formulas
and compares their complete polynomial matrices.

## 3. Augmentation weight and the first leading algebra

Put
\[
X_i=1+g_i\quad(1\le i\le n),\qquad
Z_j=1+c_j\quad(1\le j\le t).
\tag{10}
\]
Their ordered squarefree monomials form a basis: the change from (3)
is the invertible binary subset transform.

Give \(X_i\) weight one and \(Z_j\) weight two. Let \(F_{\ge r}\)
be the span of monomials of weight at least \(r\). For
\(v\in\mathbf F_2^t\), write
\[
\zeta_v=\prod_j(1+Z_j)^{v_j}+1.
\]
Its lowest-weight part is \(\sum_jv_jZ_j\). The exact relations are
\[
\begin{aligned}
X_i^2&=\zeta_{\beta(e_i,e_i)},\\
X_jX_i+X_iX_j
&=(1+X_i)(1+X_j)
  \zeta_{\beta(e_i,e_j)+\beta(e_j,e_i)}\quad(i<j),\\
Z_j^2&=0,\qquad [Z_j,X_i]=[Z_j,Z_k]=0.
\end{aligned}
\tag{11}
\]
Every right-hand term has at least the weight of the left-hand side.
Thus multiplication never decreases total weight.

The coproduct is
\[
\Delta(X_i)=X_i\otimes1+1\otimes X_i+X_i\otimes X_i,
\tag{12}
\]
with the analogous formula for \(Z_j\). Its last term has higher
total weight. The antipode is filtered as well: it reverses products,
and \(\mathcal S(X_i)\) has leading term \(X_i\), even if \(g_i\)
has order four. For example, the finite geometric expansion of
\((1+X_i)^{-1}\) gives this assertion immediately.

Keeping only weight-preserving terms gives the Hopf algebra \(H\)
with primitive generators \(x_i,z_j\) and relations
\[
\begin{aligned}
x_i^2&=\sum_k\beta(e_i,e_i)_kz_k,\\
x_jx_i+x_ix_j
&=\sum_k\bigl(\beta(e_i,e_j)+\beta(e_j,e_i)\bigr)_kz_k,\\
z_j^2&=0,\qquad z_j\ \mathrm{central}.
\end{aligned}
\tag{13}
\]
The ordered squarefree monomials remain a basis. There are
\(2^\lambda\) independent such symbols in the associated graded
vector space; relations (13) span them by normal-order rewriting.
This proves the presentation without assuming additional relations
or appealing to a general structure theorem.

The checker verifies all leading multiplication and antipode tables
against independent rewriting using (13), including the quaternion
case with nonzero generator squares.

## 4. Ordered deformations and quadratic controls are filtered

For an augmentation monomial with support \(S\), expansion of its
selected subsets in (5) gives
\[
\Phi_z(Y^S)=
\left(\prod_{i\in S}(1+z^{w_i})\right)
\prod_{\substack{i\in S\\\mathrm{in\ order}\ \sigma}}Y_i,
\tag{14}
\]
where \(Y_i\) denotes \(X_i\) or a central \(Z_j\).
Each selected generator contributes
\(1+(g_i+z^{w_i}(1+g_i))=(1+z^{w_i})(1+g_i)\).
The formula is diagonal in canonical order. In any other fixed order
its normal-order rewrite still never decreases augmentation weight.

For the control map, the dual basis to augmentation monomials consists
of the Boolean coordinate monomials in the word bits \(h_i,c_j\).
Indeed, in characteristic two,
\[
\sum_{T\subseteq S}\prod_{i\in A}T_i=[A=S].
\tag{15}
\]
Use dual weights one for \(h_i\) and two for \(c_j\). The dual map is
\[
U_f^*:h_i\mapsto h_i,\qquad c_j\mapsto c_j+f_j(h).
\tag{16}
\]
Degree at most two means every term replacing \(c_j\) has weight
at most two. Boolean reductions \(h_i^2=h_i\) only lower weight.
Thus \(U_f^*\) never increases dual weight, and \(U_f\) never
decreases primal augmentation weight.

In the first associated graded dual algebra, coordinate squares
become zero, and only the quadratic part remains:
\[
\overline U_f^*:h_i\mapsto h_i,\qquad
c_j\mapsto c_j+f_{j,2}(h).
\tag{17}
\]
The affine terms disappear at this first stage.

The checker obtains (17) independently by a group-basis permutation
followed by subset transforms, and by direct substitution in the
square-zero dual coordinate algebra.

## 5. First four-leg degeneration

Let \(D_s\) multiply an augmentation basis monomial of weight \(a\)
by \(s^a\), where \(s\) is independent of \(z\) and all coefficient
scalars. Every operation in (8)–(9) preserves or increases total
weight. Hence a nonzero tensor entry from total input weight \(a\)
to total output weight \(b\) has \(b\ge a\).

Therefore
\[
(D_s\otimes D_s)T_z(D_s^{-1}\otimes D_s^{-1})
\tag{18}
\]
is polynomial in \(s\). Its value at \(s=0\) is exactly (8)
computed with (13), the leading deformation (14), and control (17).

For nonzero \(s\), this is an invertible change on each tensor leg.
It preserves both cut ranks. Specialization at zero cannot increase
either rank, because a nonzero minor at zero was already a nonzero
polynomial minor. This degeneration is common to every coefficient
and every scalar recombination.

Beyond the operation-table comparisons, the checker compares the
entire first leading tensor with a separately assembled primitive
Hopf calculation in every recorded case.

## 6. Ordinary PBW length and the second degeneration

Now give every primitive generator of \(H\), including each \(z_j\),
ordinary length one. This is different from augmentation weight.

Each commutator or square relation in (13) replaces two generators
by one. Multiplication and antipode preserve or lower ordinary
length; the primitive coproduct preserves total length. The next
associated graded algebra is
\[
K[u_1,\ldots,u_\lambda]/(u_1^2,\ldots,u_\lambda^2),
\tag{19}
\]
with primitive generators.

Reordering in (14) only lowers length. Its length-preserving part
is diagonal, independent of \(\sigma\), and multiplies support \(S\)
by
\[
q_S(z)=\prod_{i\in S}(1+z^{w_i}).
\tag{20}
\]
The dual control (17) preserves or raises ordinary monomial length:
using a quadratic replacement adds one factor. Its primal map
therefore preserves or lowers length, with identity leading map.
This removes the quadratic control.

Both antipodes in (9) also have identity leading map in (19).
Thus the last leading \(\psi_z\) is exactly (20).

Let \(E_t\) scale a PBW monomial of length \(a\) by \(t^a\), for
a second fresh variable \(t\). On the first leading tensor, use
\[
(E_t^{-1}\otimes E_t^{-1})\overline T_z(E_t\otimes E_t).
\tag{21}
\]
The exponents are input length minus output length and are
nonnegative. The same four-leg rank argument applies.

The scalings in (18) and (21) point in opposite directions because
group augmentation products can raise weight, whereas primitive
PBW rewrites can lower length. They are two separate filtrations.

## 7. Binary conservation and the rank-sum transfer

In (19), the final Hopf tensor factors independently over its
primitive coordinates. Put \(q_i=1+z^{w_i}\). Its local two-bit
operator fixes \(00,11\) and sends
\[
10\mapsto q_i10+(1+q_i)01,\qquad
01\mapsto(1+q_i)10+q_i01.
\tag{22}
\]
With \(S=x+y=u+v\), \(b=x+u\), its kernel is
\[
\begin{array}{c|cc}
 &b=0&b=1\\ \hline
S=0&1&0\\
S=1&q_i&1+q_i
\end{array}.
\tag{23}
\]
Consequently every last leading coefficient and every recombination
has the full-domain form
\[
\widehat T(x,y,u,v)
=[x+y=u+v]K(x+y,x+u)
\tag{24}
\]
on \(\mathbf F_2^\lambda\).

Both degenerations fix crossing, since swap commutes with equal
changes on its two tensor factors. Thus an exact decomposition gives
\[
\sum_j\widehat T_j=E_d,\qquad
r(T_j)\ge r(\widehat T_j),\qquad
c(T_j)\ge c(\widehat T_j),\qquad d=2^\lambda.
\tag{25}
\]
The pinned full-conservation theorem applies to the entire
decomposition in (25), giving
\[
\sum_j\sqrt{r(\widehat T_j)c(\widehat T_j)}\ge d^{3/2}.
\]
This proves (2).

No repeated-weight or small-weight-class assumption is needed.
Arbitrarily many pieces, arbitrary coefficient recombinations, and
arbitrary characteristic-two extension fields are allowed.
The two new scaling variables can always be chosen independently
of the existing coefficient field.

## 8. Consequences and explicit limits

Since \(3/2>\sqrt2\), no full-domain seed covered here satisfies
\[
\sum_j\sqrt{r(T_j)c(T_j)}<d^{\sqrt2}.
\tag{26}
\]
For a growing family with subexponentially many pieces, the usual
normalized width parameter also satisfies \(\rho\ge3/4\).

Thus the bilinear correction screen on \(H(\mathbf F_8)\) is
not an unresolved route to a higher factoring score. Its finite
rank optima remain valid data, but its entire coefficient family
is covered. The same applies to the quadratic central corrections
already identified on \(H(\mathbf F_4)\).

The [saved eight-element-field screen](../evidence/noncommutative-crossing-F8-bilinear-checkpoints.json)
records 32 million completed normalized bilinear cases, the exact interval
`[13800000,45800000)`, per-chunk digests, commands, and the native source
hash. Every screened rank is at least the baseline 197; larger ranks were
stopped at that cap. This is a bounded screen, not an exhaustive result for
the 67,108,864-case family. It was paused when the all-size obstruction
above closed this score route. The [resumable runner](../tools/run_group_correction_screen.py)
saves each completed interval atomically and checks command and source
identity before reusing a checkpoint.

The degree restriction is essential to this proof. For three
horizontal coordinates and one central coordinate, take
\[
c\mapsto c+h_1h_2h_3.
\tag{27}
\]
Then \(U_f(X_1X_2X_3)\) has a nonzero \(Z\) term. The input
weight is three and that output weight is two. The first filtered
map hypothesis actually fails. The checker records this witness.

Arbitrary lookup corrections may have such higher-degree terms.
They remain outside the theorem. Likewise, arbitrary restrictions
or linear retracts need not commute with either scaling and cannot
be passed automatically to the full-domain conservation theorem.

Finally, the proof does not concern arbitrary kernels \(K(S,a)\)
on a nonabelian group. The ordered polynomial form, the low-degree
central control, and the common degenerations are essential.

## 9. Exact reproduction

Run the checker
[check_noncommutative_crossing_group_degeneration.py](../tools/check_noncommutative_crossing_group_degeneration.py)
with the option
\(\texttt{--output evidence/noncommutative-crossing-group-degeneration.json}\).

Its cases include an unweighted dihedral group of order eight;
affine and quadratic controls; reordered distinct weights;
the quaternion group testing nonzero restricted squares;
a direct product with an extra horizontal coordinate; and an
extraspecial group of order thirty-two.

For each case it checks the orbit/Hopf convention, all leading
multiplication and antipode tables, the central-control substitution,
the complete first and second tensor degenerations, and equality
with direct abelian mixing. It also checks both-cut rank domination
for every coefficient and selected or exhaustive recombinations.

The evidence includes the cubic counterexample (27), auditing
the scope as well as the positive identities. The all-size claim
is the ordinary proof above; the finite computations check its
conventions and implementation.
