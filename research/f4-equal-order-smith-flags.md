# Reusing coefficient fibers while changing F4 Smith row flags

## Result and scope

Equal-order Smith blocks admit a useful exact change of search coordinates.
A constant invertible row change within one such block preserves its common
row floor. The complete coefficient fiber can therefore be eliminated once
and reused while the requested higher row conditions change.

This enlarges the row conditions explored by the
[nonbinary reference-frame constructor](f4-recentered-crossing-lattices.md).
It does not change the crossing endpoint, the polynomial degree, the physical
heights, or the set of coefficients allowed by the common floor. It does not
establish an exponent improvement or a uniform obstruction.

The implementation is
[`f4_smith_flag_fibers.py`](../tools/f4_smith_flag_fibers.py), with the finite
search driver
[`search_f4_smith_flags.py`](../tools/search_f4_smith_flags.py). These files use
the predecessor's complete coefficient space and exact F4 arithmetic. The
construction and checks used substantial OpenAI assistant support.

## 1. The common-floor identity

Fix either crossing cut and a polynomial reference matrix \(A(B)\).
Let \(L,R\) be regular polynomial frames whose complete Smith orders are
\(\sigma_1,\ldots,\sigma_n\). The reference need not have binary
coefficients. For another crossing tensor \(T\), write

\[
F(T)=L A(T)R.
\]

Choose common row floors \(f_i\) and finite row caps \(c_i\). Assume
that both are constant on each block of indices having equal \(\sigma_i\).
The searches use

\[
f_i=\max(0,\sigma_i-a),\qquad c_i=\sigma_i+b
\]

for fixed positive integers \(a,b\).

Let \(U\) be a constant invertible matrix over \(\mathbf F_4\), block
diagonal on these equal-order blocks. Define

\[
L'=UL,\qquad F'(T)=UF(T).
\]

Then \(L'\) is a regular frame. Moreover,

\[
F(T)\in\mathrm{diag}(z^{f_i})\mathrm{Mat}_n(\mathbf F_4[[z]])
\quad\Longleftrightarrow\quad
F'(T)\in\mathrm{diag}(z^{f_i})\mathrm{Mat}_n(\mathbf F_4[[z]]).
\tag{1}
\]

To prove (1), note that \(U\) and \(U^{-1}\) commute with
\(\mathrm{diag}(z^{f_i})\). Multiplication by either matrix preserves the
displayed module. The argument applies independently at both places.

Consequently the two-place common-floor affine fiber is exactly the same
before and after either frame change. Its free-coordinate parameterization
can be reused. No earlier floor equation is discarded or weakened.

For higher targets \(\tau_i\ge f_i\), the new conditions are

\[
[z^k]\sum_j U_{ij}F(T)_{j\ell}=0
\qquad(k<\tau_i,\ 1\le\ell\le n).
\tag{2}
\]

The targets need not be equal inside an equal-order block. Equation (2)
therefore selects different affine subspaces as \(U\) changes. These
conditions remain linear in the complete crossing coefficients.

## 2. Why a fixed row frame can miss a row lift

A two-by-two local example makes the distinction explicit. Take reference
matrix \(B=I\), whose two Smith orders and common row floors are both zero.
For another candidate matrix and an allowed flag, put

\[
A(z)=\begin{pmatrix}1&0\\ \omega&z\end{pmatrix},\qquad
U=\begin{pmatrix}1&0\\ \omega&1\end{pmatrix},
\qquad \omega^2+\omega+1=0.
\]

Both rows of \(A\) have order zero. The candidate has complete Smith orders
\((0,1)\), and

\[
UA=\begin{pmatrix}1&0\\0&z\end{pmatrix}.
\]

Thus a row lift can require a nonbinary row combination even when the
old coordinate rows show no lift. This is a local matrix illustration,
not a new binary crossing construction or an exception to an existing
small-alphabet bound.

The search samples invertible matrices by elementary row transvections
inside each equal-order block. Every retained matrix is stored literally.
Its invertibility and block support are checked before its equations are
used. Singular changes and mixtures across unequal reference orders are
rejected.

## 3. Exact reuse of the affine coefficient coordinates

The underlying tensor still has the form

\[
T=B+(1+w)Q,\qquad\deg Q_{xyuv}<D.
\]

The common floors give an affine system in all \(Dq^4\) field variables.
The preceding implementation eliminates that system and expresses the
original variables in free affine coordinates.

The flag constructor retains the matrix-column label of every coefficient
equation. For each row, coefficient power, and matrix column, it stores the
corresponding affine expression in those same free coordinates. The column
labels matter: row changes must combine coefficients from matching columns.
Combining an unordered list of nonzero equations would be incorrect.

After independently reconstructing these column-resolved expressions, the
implementation removes their column labels and zero equations. It requires
exact equality, in order, with every equation in the existing compressed
fiber. This check covers the reversed polynomial at infinity and all
coefficient powers in the requested windows.

For a new flag, equation (2) is then a finite F4 linear combination of the
stored expressions. The implementation does not repeat elimination of the
original \(Dq^4\)-variable floor system. It samples literal tensors only
after the new target system is solved. All reported objective values use
both complete realized profiles and both physical charge lists.

## 4. Independent finite controls

Run:

```sh
python3 tools/check_f4_smith_flag_fibers.py \
  --output /tmp/f4-smith-flag-controls.json
```

The [recorded controls](../evidence/f4-smith-flag-controls.json) cover
three complete coefficient-channel comparisons, twelve nonbinary row flags,
2,832 literal mixed polynomial coefficients, ninety feasible tensor
substitutions, and six rejected invalid frames.

The coefficient controls choose a free affine assignment, reconstruct its
full literal tensor, and multiply its polynomial matrices by the changed
frames without truncation. Each resulting coefficient is compared with the
stored affine expression evaluated on that assignment. The cases include
two- and three-symbol tensors and polynomial degrees two and three.

These are finite implementation controls. They verify the stated coordinate
identity and the cases recorded in the evidence; they do not establish that
the chosen flags exhaust all useful row subspaces.
