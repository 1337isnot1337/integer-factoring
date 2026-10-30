# Several Ward forms in the same eight-dimensional operator algebra

The [stable operator interface](ward-stable-operator-interface.md)
keeps eight ancillary operator channels. Multiplying several Ward
forms does not require taking a tensor product of those ancillary
spaces: the eight-dimensional space is an algebra. The form-selection
information remains a separate coefficient problem.

This note gives the exact distinction. There is an eight-state
matrix-product recurrence over a ring of form variables. Expanding
that ring into scalar coordinates still uses \(2^k\) coordinates.
No polynomial-time coefficient-rank algorithm is asserted.

## 1. Form products as one generating coefficient

At physical pair \(i\), let
\[
 T_i(\eta)=R_i+\eta O_i,
 \qquad\eta=q+r.
\]
The local operators satisfy
\[
 R_i^2=O_i^2=R_iO_i=O_iR_i=0,
 \qquad T_i(\eta)^2=0.
 \tag{1}
\]
Operators on different sites commute. Let \(A=(w_{ai})\) be any
\(k\times P\) weight matrix over a characteristic-two field. Put
\[
 L_a=\sum_iw_{ai}R_i,\qquad
 O_a=\sum_iw_{ai}O_i,
 \qquad\ell_i(z)=\sum_{a=1}^kw_{ai}z_a.
\]
Then
\[
 \boxed{\quad
 \prod_{a=1}^k(L_a+\eta O_a)
 =[z_1\cdots z_k]\prod_{i=1}^P
                 (I+\ell_i(z)T_i(\eta)).
 \quad}
 \tag{2}
\]
Indeed, expanding the left side selects a site for every form.
Repeated sites vanish by (1). On the right, extracting the indicated
coefficient makes exactly the same injective assignments. This proof
works with arbitrary characteristic-two weights, not only binary
ones.

One may calculate in the commutative ring
\[
 \mathcal R_k=F[z_1,\ldots,z_k]/(z_1^2,\ldots,z_k^2).
 \tag{3}
\]
Discarding the squares does not change the requested coefficient.

The coefficient for a fixed active site set \(S\), of size \(k\),
is
\[
 \det(A_{[k],S})\prod_{i\in S}T_i(\eta).
 \tag{4}
\]
In characteristic two the determinant signs disappear. In the
original occupation basis, the source and target determine the
active sites and which of them use full-pair insertion. Thus an
operator entry is either zero or
\(\eta^j\det(A_{[k],S})\), where \(j\) is the number of such
insertions. This gives direct determinant evaluation before the
stable change of basis.

## 2. Insertion transitions stay in the same algebra

Use the local stable isomorphism \(X\) and basis
\(A_1,\ldots,A_8\) of
\[
 \mathcal A=\operatorname{Mat}_2(F)\otimes F[e]/(e^2)
       \subseteq\operatorname{End}(G).
\]
The prior certificate gives
\[
 \Phi(A_\alpha)=X^{-1}(I_S\otimes A_\alpha)X
        =\sum_\beta A_\beta\otimes B_{\alpha\beta}.
\]
Since \(\mathcal A\) is an algebra, the inserted maps also have
eight-channel expansions:
\[
 \begin{aligned}
 X^{-1}(R_S\otimes A_\alpha)X
      &=\sum_\beta A_\beta\otimes D^R_{\alpha\beta},\\
 X^{-1}(O_S\otimes A_\alpha)X
      &=\sum_\beta A_\beta\otimes D^O_{\alpha\beta}.
 \end{aligned}
 \tag{5}
\]
Equivalently, these maps are the products \(\Psi_R\Phi(A_\alpha)\)
and \(\Psi_O\Phi(A_\alpha)\). The coefficient matrices in (5) are
explicit \(16\times16\) binary matrices. The new checker constructs
them, reconstructs all 16 inserted operator matrices, and compares
the resulting two-site recurrences with 4,096 directly conjugated
operator columns.

## 3. The all-size recurrence

Let \(U_P\) be the sequential stable isomorphism, and define
\[
 V_P(A)=U_P^{-1}
 \left(\prod_{i=1}^P(I+\ell_i(z)T_i(\eta))\otimes A\right)U_P.
\]
The last local isomorphism and (5) give
\[
 \boxed{
 \begin{aligned}
 V_0(A)&=A,\\
 V_P(A_\alpha)
 &=\sum_{\beta=1}^8 V_{P-1}(A_\beta)\otimes
 \left[B_{\alpha\beta}
  +\ell_P(z)\bigl(D^R_{\alpha\beta}
                         +\eta D^O_{\alpha\beta}\bigr)\right].
 \end{aligned}}
 \tag{6}
\]
This is a matrix-product recurrence with eight states over
\(\mathcal R_k[\eta]\), regardless of \(P\) and \(k\).
Extracting \(z_1\cdots z_k\) from \(V_P(I_G)\) gives the stable
conjugate of the product in (2). The outside Ward factor \(\eta^k\)
is then applied as usual.

The state count in (6) refers to ancillary operator channels.
It does not remove the \(2^k\) square-free monomials of
\(\mathcal R_k\). A scalar expansion of the recurrence can still
require that form-subset memory. The determinant formula (4)
evaluates entries in the original occupation basis; it is not a
proof of polynomial-time coefficient extraction after the stable
basis change. Nor does it give a sharp rank bound for the resulting
large matrices.

This separates two issues that a direct product of one-form
matrix-product representations would conflate. The ancillary algebra
remains eight-dimensional. The correlations among form selections
must still be retained or compressed by an additional argument.

## 4. Exact certificate

The checker is
[check_ward_form_operator_transport.py](../tools/check_ward_form_operator_transport.py).
Its output is
[ward-form-operator-transport.json](../evidence/ward-form-operator-transport.json).
It verifies the 64 local square-zero and mixed-zero operator columns,
all insertion transitions in (5), the two-site ordered-insertion
recurrence, and 1,280 form-product versus site-minor polynomial
controls, including rank-deficient form matrices.

The proof of (2) and the local closure (5) make (6) valid for every
number of sites and forms. The finite controls verify the concrete
matrices, tensor ordering, and scalar-ring conventions used to
implement that proof.

