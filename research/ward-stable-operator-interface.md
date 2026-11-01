# An eight-dimensional interface for the actual Ward operators

The complete-rectangle Ward bound discards the incidence and
full-pair creation maps before calculating ranks. This note retains
those maps through the stable Pfaffian splitting. It gives an exact
operator recurrence with an eight-dimensional ancillary algebra and
a matrix-product bond bound of nine, independent of the number of
physical pairs.

This is a constructive representation theorem. It does not yet give
a polynomial-time coefficient-rank algorithm or an improved factoring
exponent. In particular, a small matrix-product bond does not by
itself bound the rank of a large coefficient matrix.

## 1. Local operators and the stable splitting

Use the exact local Pfaffian matrix and stable identity
\[
 S=2\Delta\oplus2G\oplus H,
 \qquad S\otimes G\cong G\otimes F^{\otimes2},
 \qquad F=g\oplus h,
\]
from [the stable module construction](stable-pfaffian-modules.md).
Here \(G=h\otimes g\), and \(H\) has dimension six.

Let \(R\) be the local incidence operator on the external and state
legs, and let \(O\) multiply the external monomial by its full pair.
They have grade shifts \((1,\pm1)\) and \((2,0)\), respectively.
Both preserve each local particle-hole exterior degree. Their two
four-dimensional copies agree.

On \(G\), write
\[
 e=\begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
 f=\begin{pmatrix}0&0\\1&0\end{pmatrix}.
\]
Then
\[
 R_G=f\otimes I+I\otimes e,
 \qquad O_G=f\otimes e.
 \tag{1}
\]
On \(H\), in the basis used by the stable module note,
\(O_H=E_{2,1}+E_{6,5}\). Matrix units here use one-based indices.

The checker constructs a constant binary matrix \(X\) of size 64
such that
\[
 (S\otimes G)X=X(G\otimes F\otimes F).
 \tag{2}
\]
It derives \(X\) from the existing 24-dimensional \(H\otimes G\)
intertwiner and the eight-dimensional natural-cube splitting. It
verifies all 4,096 identities in the polynomial ring
\(\mathbf F_2[p,q,r]\), proves invertibility by exact elimination,
and checks both exterior grades of every nonzero basis entry.
The columns are recorded in the evidence file.

## 2. Why the ancillary correction matters

An intertwiner for a representation preserves the total infinitesimal
action, including the ancillary factor. It need not preserve the
physical operator alone.

For the 24-dimensional splitting, use target order
\[
 \Delta G,\quad\Delta G,\quad G\otimes g^{[2]},\quad G\otimes h^{[2]}.
\]
The total incidence
\(R_H\otimes I+I\otimes R_G\) becomes block diagonal: it is \(R_G\)
on the first two blocks and \(R_G\otimes I_2\) on the last two.

The total creation operator
\(O_H\otimes I+I\otimes O_G\) has the following exact form. It is
\(O_G\) on each \(\Delta G\). On the merged last two blocks it is
\[
 O_G\otimes I_4+I_G\otimes K,
 \qquad K=E_{1,3}+E_{4,2},\qquad K^2=0,
 \tag{3}
\]
where the four coordinates are
\(g^{[2]}_0,g^{[2]}_1,h^{[2]}_0,h^{[2]}_1\).
Thus the Frobenius blocks are coupled, but their coupling is independent
of \(G\). The physical-only operator has additional couplings; the
certificate records them separately.

For the full splitting (2), let \(R_F=e\oplus f\). The total
incidence transports to
\[
 X^{-1}(R_S\otimes I+I\otimes R_G)X
 =R_G\otimes I\otimes I
  +I\otimes R_F\otimes I+I\otimes I\otimes R_F.
 \tag{4}
\]
The creation operator requires a larger interface. After subtracting
\(O_G\otimes I_{16}\), its ancillary operator Schmidt rank is four.
Here that rank means the dimension of the span of its ancillary
coefficient matrices, with the two \(F\) factors regarded as the
other tensor factor. Equation (3) and this rank statement are checked
over exact binary matrices.

## 3. A closed ancillary algebra

Define the unital algebra
\[
 \mathcal A=\operatorname{Mat}_2(\mathbf F_2)
                 \otimes\mathbf F_2[e]/(e^2)
       \subseteq\operatorname{End}(G).
 \tag{5}
\]
The first factor acts on \(h\); the second acts on \(g\), with the
matrix \(e\) from (1). Its dimension is eight.

For \(A\in\operatorname{End}(G)\), define the transport map
\[
 \Phi(A)=X^{-1}(I_S\otimes A)X.
\]
Also set
\[
 \Psi_R=X^{-1}(R_S\otimes I_G)X,
 \qquad
 \Psi_O=X^{-1}(O_S\otimes I_G)X.
\]
The exact finite certificate proves
\[
 \boxed{\quad
 \Phi(\mathcal A)\subseteq
       \mathcal A\otimes\operatorname{End}(F^{\otimes2}),
 \qquad
 \Psi_R,\Psi_O\in
       \mathcal A\otimes\operatorname{End}(F^{\otimes2}).
 \quad}
 \tag{6}
\]

The checker starts with the ancillary coefficients of the two physical
operators and repeatedly closes their span under \(\Phi\). It obtains
exactly the eight-dimensional space (5). It reconstructs every
transported matrix from the recorded coefficients and checks all
64 products of basis matrices. The chosen basis is homogeneous in
both exterior grades, with shifts
\[
 (0,0),\ (-1,-1),\ (0,0),\ (1,1),\ (0,-2),\
 (1,-1),\ (1,-1),\ (2,0).
\]

Thus (6) is an algebraic closure statement, not a dimension estimate
inferred from a few numerical specializations. Its matrices are
constant and apply over every characteristic-two scalar field.

## 4. An explicit recurrence for all numbers of pairs

Choose the recorded basis \(A_1,\ldots,A_8\) of \(\mathcal A\), and
write the exact finite expansions
\[
 \Phi(A_\alpha)=\sum_{\beta=1}^8 A_\beta\otimes B_{\alpha\beta},
 \qquad
 \Psi_R=\sum_{\beta=1}^8A_\beta\otimes C^R_\beta,
 \qquad
 \Psi_O=\sum_{\beta=1}^8A_\beta\otimes C^O_\beta.
 \tag{7}
\]
All \(B,C\) are explicit \(16\times16\) binary matrices in the
evidence file.

Let \(U_P:G\otimes F^{\otimes2P}\to S^{\otimes P}\otimes G\)
be the sequential stable isomorphism obtained by applying (2) at
each site. Define
\[
 T_P(A)=U_P^{-1}(I_{S^{\otimes P}}\otimes A)U_P,
\]
and let \(\widehat R_P,\widehat O_P\) be the conjugates of the sums
of the physical local operators. Then
\[
 \begin{aligned}
 T_0(A)&=A,\qquad \widehat R_0=\widehat O_0=0,\\
 T_P(A_\alpha)&=
   \sum_\beta T_{P-1}(A_\beta)\otimes B_{\alpha\beta},\\
 \widehat R_P&=\widehat R_{P-1}\otimes I_{16}
       +\sum_\beta T_{P-1}(A_\beta)\otimes C^R_\beta,\\
 \widehat O_P&=\widehat O_{P-1}\otimes I_{16}
       +\sum_\beta T_{P-1}(A_\beta)\otimes C^O_\beta.
 \end{aligned}
 \tag{8}
\]
These identities follow directly by factoring the last local
isomorphism out of \(U_P\). They prove the recurrence for every \(P\)
from the finite identities (6)--(7).

The combination \(\widehat R_P+\eta\widehat O_P\) uses the same
eight transported states and one accumulated insertion state, with
emission \(C^R_\beta+\eta C^O_\beta\). Therefore
\[
 \boxed{\quad\text{matrix-product bond dimension}\ \le9.\quad}
 \tag{9}
\]
For one global Ward form, the desired operator is
\(\eta(R_P+\eta O_P)\), so (9) retains the full combination before
taking coefficient ranks. The outside factor \(\eta\) adds no
operator state.

The checker independently constructs the complete two-site stable
isomorphism of dimension 1,024 and compares 2,048 columns of the
directly conjugated physical operators with (8). This checks the
ordering of tensor factors and the ancillary correction in the
recurrence.

## 5. Verification and remaining rank problem

Run
\[
 \texttt{python tools/check\_ward\_stable\_operators.py}.
\]
The checker is
[check_ward_stable_operators.py](../tools/check_ward_stable_operators.py),
and its complete finite matrices are saved in
[ward-stable-operators.json](../evidence/ward-stable-operators.json).

The recurrence keeps the constant operators that the unrestricted
rectangle bound loses. The remaining task is to calculate, or bound
sharply, the ranks of its coefficient matrices after the grade
projections and the chosen code maps. Further natural-module fusion
can mix multiplicity spaces, so the eight-dimensional ancillary
closure is not a proof that all such ranks reduce to independent
polynomial-size blocks.

