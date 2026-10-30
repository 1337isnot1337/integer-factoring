# Euler refinement and a small Ward interface

The Ward factorization can be evaluated without expanding the Clifford
channels separately. Its incidence operator is a derivative in a common
parameter direction. After fixing suitable grades and coefficient
degrees, that derivative becomes a scalar.

For one global symplectic form, this gives **two rectangular rank
bounds for either cut**, with no channel multiplier. For several binary
forms, the remaining interface is a \(k\times k\) determinant of
external creation operators. Syndrome-preserving outputs also give a
joint vanishing order
\[
 \boxed{\qquad
 G(x,y)\in (x,y)^{\,2k-\operatorname{rank}C_{\rm total}},
 \qquad C_{\rm total}=AA^{\mathsf T}.
 \qquad}
 \tag{1}
\]
In particular, even multiplicities of every form-column type give
order at least \(2k\).

These are positive rank factorizations and coefficient-support
restrictions. They do not assert an improved factoring exponent.
The [recorded rank computations](ward-charge-bound-computations.md)
give the current numerical results and reproduction commands.

## 1. Setup

Use the universal Pfaffian kernel and operators from
[the Ward theorem](noncommutative-crossing-radical-ward.md). Let
\[
 c_i=a_i+b_{\bar i},\qquad
 Q_a=\sum_s w_{as}c_sc_{\bar s}.
\]
The \(k\times P\) form matrix is \(A=(w_{as})\).
Let the input code \(I\) have physical exterior grade \(L\) and lie in
\(\operatorname{im}(Q_1\cdots Q_k)\). Choose the homogeneous preimage
\[
 I=Q_1\cdots Q_k W,\qquad \deg W=L-2k.
 \tag{2}
\]
The existence proof is in
[the graded Ward descent note](noncommutative-crossing-ward-rank-descent.md).
Let \(\Pi\) be any grade-\(L\) output map. The statements through
Section 4 do not require a coordinate output map or
\(\Pi I=I_d\).

Write
\[
 \mathcal L_a=\sum_s w_{as}
       (A_{\bar s}c_s+A_sc_{\bar s}),\qquad
 \mathcal O_a=\sum_s w_{as}A_sA_{\bar s}.
\]
The symbol \(A_s\) here denotes an external Grassmann variable, not an
entry of the form matrix. Define
\[
 \mathcal T_j=[z^j]\prod_{a=1}^k(\mathcal L_a+z\mathcal O_a).
\]
With \(\eta=q+r\), the exact left-sided Ward factorization is
\[
 G=\eta^k\sum_{j=0}^k\eta^j
  (\Pi_A\otimes W_Y^{\mathsf T})\,
  \mathcal T_j\mathscr M\,
  (I_V\otimes\Pi_U^{\mathsf T}).
 \tag{3}
\]

## 2. The incidence operator is a parameter derivative

Temporarily give every physical pair its own parameters
\(p_s,q_s,r_s\). The three parameter-bearing edges are
\(p_sU_sA_s\), \(q_sA_sY_{\bar s}\), and \(r_sA_sV_s\),
including both modes of that pair. Put
\[
 D_a=\sum_s w_{as}
   (\partial_{p_s}+\partial_{q_s}+\partial_{r_s}).
 \tag{4}
\]

**Lemma 1.**
\[
 \boxed{\qquad \mathcal L_a\mathscr M=D_a\mathscr M.\qquad}
 \tag{5}
\]

**Proof.** The first-order Pfaffian identity gives
\[
 c_s\mathscr M
   =(Y_s+q_sA_s+U_{\bar s}+V_{\bar s})\mathscr M.
\]
Multiply by \(A_{\bar s}\), add the complementary-mode equation
multiplied by \(A_s\), and use characteristic two. The two
\(q_sA_sA_{\bar s}\) terms cancel. The result is insertion of
\[
 A_sY_{\bar s}+A_{\bar s}Y_s+
 A_sU_s+A_{\bar s}U_{\bar s}+
 A_sV_s+A_{\bar s}V_{\bar s}.
\]
This is exactly the sum of the three parameter derivatives at the
pair. Sum with the weights \(w_{as}\). \(\square\)

The operators \(\mathcal L_a,\mathcal O_a\) are parameter-independent.
They commute with the \(D_b\). Therefore
\[
 \boxed{
 \mathcal T_j\mathscr M
 =
 \sum_{\substack{J\subseteq[k]\\|J|=j}}
  \left(\prod_{a\in J}\mathcal O_a\right)
  \left(\prod_{a\notin J}D_a\right)\mathscr M .
 }
 \tag{6}
\]
Every \(\mathcal O_a\) acts only on the external \(A\) leg. It is thus a
one-sided linear map for **both** tensor cuts. This is the useful
difference from expanding \(\mathcal L_a\) into separate state
channels.

The independent local parameters are an auxiliary device for (5)--(6).
The Ward factor \(\eta^{k+j}\) in (3) is used after identifying all local
parameters with the original \(p,q,r\).

## 3. One global form: two rectangles for both cuts

Suppose \(k=1\) and \(w_{1s}=1\) for every pair. Write
\[
 D=\partial_p+\partial_q+\partial_r,\qquad
 \mathcal O=\sum_sA_sA_{\bar s}.
\]
Every fixed-\(A\)-grade-\(R\) entry of the unrestricted Pfaffian is
homogeneous of parameter degree \(R\): every parameter-bearing edge
uses exactly one external \(A\) mode.

In the centered coordinates \(p=1+x,q=1+y,r=1\), Euler's homogeneous
identity gives
\[
 (D\mathscr M)(1+x,1+y,1)
 =
 (R+x\partial_x+y\partial_y)
      \mathscr M(1+x,1+y,1).
 \tag{7}
\]
Thus the coefficient \([x^ay^b]\) is multiplied by
\(R+a+b\), interpreted in characteristic two.

Use \(M_P(R_A,R_Y;C_V,C_U)[a,b]\) for the full rectangular coefficient
tensor defined in the graded Ward note. Put
\[
 \epsilon_{L,a,b}=(L+a+b-1)\bmod2.
\]
Equation (3) becomes the following exact coefficient factorization:
\[
 \begin{aligned}
 G[a,b]={}&
 \epsilon_{L,a,b}\,
 (\Pi_A\otimes W_Y^{\mathsf T})
 M_P(L,L-2;L,L)[a,b-1]
 (I_V\otimes\Pi_U^{\mathsf T})\\
 &+
 (\Pi_A\otimes W_Y^{\mathsf T})
 \mathcal O_A
 M_P(L-2,L-2;L,L)[a,b-2]
 (I_V\otimes\Pi_U^{\mathsf T}).
 \end{aligned}
 \tag{8}
\]
Here the displayed grouping is the row cut. After reshaping to the
column cut, all four leg maps and \(\mathcal O_A\) remain one-sided
maps. No cross-cut channel is introduced.

**Theorem 2.** For either cut rank \(\rho\in\{r,c\}\),
\[
 \boxed{
 \begin{aligned}
 \rho(G[a,b])\le{}&
 \epsilon_{L,a,b}\,
 \rho\!\left(M_P(L,L-2;L,L)[a,b-1]\right)\\
 &+\rho\!\left(M_P(L-2,L-2;L,L)[a,b-2]\right).
 \end{aligned}}
 \tag{9}
\]
Negative coefficient indices mean zero. As usual, one may take the
minimum with other verified bounds and \(d^2\).

For example, for the full \(P=3,L=3\) radical, the literal coefficient
\((a,b)=(0,1)\) has cut ranks \((14,14)\). The new bound is
\((29,29)\); the earlier separate-channel column bound was \(312\).
At \((0,2)\), the new bounds are \((14,28)\), against actual
\((14,14)\). These are finite controls, not asymptotic claims.

## 4. Repeated column types and the determinant interface

Group the physical pairs by equal columns of the form matrix. Let
\(\mathcal V\) be the set of column types, including a zero type when
present. Give type \(v\) common parameters \(p_v,q_v,r_v\), and put
\[
 D_v=\partial_{p_v}+\partial_{q_v}+\partial_{r_v},\qquad
 D_a=\sum_{v\in\mathcal V}v_aD_v.
\]
Ordinary derivatives satisfy \(D_v^2=0\) in characteristic two.

Fix a full Pfaffian component with external \(A\)-grade \(R_v\) in
each type, and centered coefficient degrees \((\alpha_v,\beta_v)\)
in that type. Define
\[
 \delta_v=(R_v+\alpha_v+\beta_v)\bmod2.
 \tag{10}
\]
For distinct types, the product of the corresponding \(D_v\) acts on
this component by the product of the \(\delta_v\). For a row subset
\(R\subseteq[k]\), the exact scalar is therefore
\[
 \chi_R(\delta)=
 \sum_{\substack{S\subseteq\mathcal V\\|S|=|R|}}
  \det(A^{\rm type}_{R,S})\prod_{v\in S}\delta_v.
 \tag{11}
\]
Repeated derivatives must be discarded. In particular, (11) is
generally **not** the product of the individual row scalars
\(\sum_vv_a\delta_v\).

For binary forms, Cauchy--Binet gives
\[
 \chi_R(\delta)=\det C(\delta)_{R,R},\qquad
 C(\delta)=
 A^{\rm type}\operatorname{diag}(\delta)
       (A^{\rm type})^{\mathsf T}.
 \tag{12}
\]
Binary minors equal their squares, which is the reason the Gram
determinant has exactly the unsquared minors in (11).

For the fixed component, all terms in (6) combine into the single
external creation map
\[
 \begin{aligned}
 \mathcal B_j(\delta)
 &=\sum_{\substack{J\subseteq[k]\\|J|=j}}
      \det C(\delta)_{J^c,J^c}\prod_{a\in J}\mathcal O_a\\
 &= [z^j]\det\!\left(
       C(\delta)+z\operatorname{diag}
          (\mathcal O_1,\ldots,\mathcal O_k)\right).
 \end{aligned}
 \tag{13}
\]
The entries belong to a commutative exterior algebra, so this
determinant is an ordinary finite polynomial expression. It supplies a
\(k\times k\) interface rather than a separate state-channel sum.

For each \(j\), the source full tensor has total grades
\[
 (R_A,R_Y;C_V,C_U)=(L-2j,L-2k;L,L),
 \tag{14}
\]
and its centered coefficient totals are
\[
 \sum_v\alpha_v=a,\qquad
 \sum_v\beta_v=b-k-j.
 \tag{15}
\]
The map (13) acts only on \(A\). Consequently each typed component
factors for both cuts through the corresponding full rectangular
coefficient tensor with **no channel multiplier**. Typed components
with the same map (13) can be combined before taking rank.

This statement does not identify the rank of a charge-filtered
rectangle with the rank of the unfiltered rectangle. Such a rank
estimate requires a separate construction or computation.

## 5. Quadratic jet charges and the joint vanishing order

Assume now that the forms are binary and the output map preserves the
admissible single-mode syndromes of the input code. Coordinate
left inverses from the overlapping-form construction have this
property. Arbitrary left inverses need not have it.

Let
\[
 C_{\rm total}=\sum_{s=1}^P a_sa_s^{\mathsf T}=AA^{\mathsf T},
\]
where \(a_s\) is the form column at physical pair \(s\). For an
admissible syndrome \(\sigma\), the even positions \(E\) satisfy
\(A_EA_E^{\mathsf T}=0\). Hence
\[
 \sum_{s\in\sigma}a_sa_s^{\mathsf T}=C_{\rm total}.
 \tag{16}
\]
The parity of the physical \(A\)-grade in a type is the parity of its
number of single sites. The operators \(\mathcal O_a\) add full pairs
and preserve that parity. Therefore every component that survives the
output map has
\[
 C(\delta)=C_{\rm total}+Q_{\rm jet},\qquad
 Q_{\rm jet}=\sum_v(\alpha_v+\beta_v)vv^{\mathsf T}.
 \tag{17}
\]
Thus all typed components of the same quadratic jet charge share
the same external map (13).

A nonzero coefficient in (13) requires
\[
 k-j\le \operatorname{rank}C(\delta).
 \tag{18}
\]
Indeed, its scalar coefficients are principal minors of order \(k-j\).
On the other hand, (15) and rank subadditivity give
\[
 \begin{aligned}
 \operatorname{rank}C(\delta)
 &\le \operatorname{rank}C_{\rm total}
       +\operatorname{rank}Q_{\rm jet}\\
 &\le \operatorname{rank}C_{\rm total}
       +\sum_v(\alpha_v+\beta_v)\\
 &=\operatorname{rank}C_{\rm total}+a+b-k-j.
 \end{aligned}
\]
Together with (18), this proves \(a+b\ge
2k-\operatorname{rank}C_{\rm total}\), which is (1).

The original Ward factor still gives \(b\ge k\). Both restrictions
hold simultaneously. If every nonzero column type occurs an even
number of times, then \(C_{\rm total}=0\), and all centered
coefficients of total degree less than \(2k\) vanish.

One can retain the exact charge rather than only its rank. Introduce
the group algebra of the binary span of the symmetric matrices
\(a_sa_s^{\mathsf T}\), with basis elements \(e_Qe_{Q'}=e_{Q+Q'}\).
The substitutions
\[
 p_s=1+x e_{a_sa_s^{\mathsf T}},\qquad
 q_s=1+y e_{a_sa_s^{\mathsf T}},\qquad r_s=1
 \tag{19}
\]
produce precisely the charge-filtered rectangular coefficients.
Equation (13) then applies once to each charge component. This gives
a constructive interface for a further rank calculation.

## 6. Exact controls and limitations

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_ward\_euler.py}.
\]
The checker solves for genuine homogeneous preimages, verifies the
complete polynomial identity (8), and checks both cut bounds against
literal full rectangular ranks. It includes coordinate and
noncoordinate output maps, so the global-form theorem is not relying
on an output-syndrome assumption.

It also checks all 256 local polynomial derivative identities and all
256 local Euler identities, 1,096 exact typed Gram/minor cases, and
complete polynomial tensors for the six-label balanced two-form code
and a three-label subcode of the 34-label unbalanced example. Their
joint multiplicities satisfy (1); the latter attains the predicted
order three.

The evidence is in
[noncommutative-crossing-ward-euler.json](../evidence/noncommutative-crossing-ward-euler.json).
The general determinant identity follows from (5), ordinary
characteristic-two differentiation, Cauchy--Binet, and the
principal-minor expansion of a determinant.

The joint multiplicity theorem in Section 5 does require
syndrome-preserving outputs on the stated homogeneous product-image
code. The two-rectangle theorem in Section 3 does not. Neither theorem
covers a raw mixture of different global exterior grades without an
independent crossing-endpoint argument.

