# Terminal identity operators just above a power of two

## Result and scope

Let \(n=N+1\), where \(N\ge2\) is a power of two. An arbitrary
occupation-preserving operator filter on the last \(r\) polynomial
identity directions has a bounded zero-cut core with \(2(r+1)\)
Gaussian modes. Its complete profiles follow from a Weyl matrix of
order \(2^{r+1}\), independently of the interval length.

The extra identity direction is the existing physical direction
\(j=0\), whose filter remains the identity. It adds no physical mode
or polynomial degree. Its ordinary Gaussian contribution is included
in the zero core and remains a spectator at infinity. This gives a
complementary family to the
[reduction just below a power of two](near-power-identity-operator-reduction.md).
No improvement of the factoring exponent is asserted here.

## 1. The old reflected basis isolates all required forms

Work in characteristic two. Set

\[
q=z^2,\qquad \tau=q^P,\qquad c=1+\tau,
\qquad n=N+1,
\tag{1}
\]

where \(P\) is a positive odd integer. At aggregate gauge \(J=1\),
the diagonal Gaussian has identity-mode coefficients \(\tau^j\)
and auxiliary-mode coefficients \(\tau^{i+1}\). The polynomial
Pascal forms satisfy

\[
f_j(i;c)=\sum_h\binom ih\binom jh c^h
        =\tau^{i\mathbin\&j}.
\tag{2}
\]

Here the bitwise expression in (2) follows directly from Lucas'
theorem and the binary factorization of \(1+c\).

Restrict first to the old rows and columns \(0,\ldots,N-1\).
Let \(H_{ij}=f_j(i;c)\), let \(R\) reverse their order, and put
\(B=HR\). The
[prime-power reflected basis](prime-power-joint-identity-filters.md)
gives

\[
B^{-1}H=R,\qquad
H^{\mathsf T}D^{-1}BD=c^{N-1}R,
\qquad D=\operatorname{diag}(\tau^{i+1}).
\tag{3}
\]

Both this change and its Gaussian companion are regular, with regular
inverses. Explicitly,
\(B_{ij}=\tau^{i\mathbin\&(N-1-j)}\),
\(B^2=c^{N-1}I\), and

\[
\operatorname{ord}_\tau(B_{ij})+j-i
 =j\mathbin\&(N-1-i)\ge0.
\tag{4}
\]

Extend this change by the identity on the new auxiliary index \(N\).
The terminal selected indices, augmented identity indices, and
retained auxiliary indices are

\[
S=\{N-r+1,\ldots,N\},\qquad
\widehat S=(S,0),\qquad
T=\{0,\ldots,r-2,N-1,N\}.
\tag{5}
\]

For \(r=1\), the first set in the union defining \(T\) is empty.
Require \(1\le r\le N\), so the augmented indices are distinct.

On an old row, the new column \(j=N\) equals column zero. On the
new row, every old column equals one, and the new column equals
\(\tau^N\). Thus the exact concentrated form matrices have entries

\[
\begin{aligned}
(F_X)_{ij}&=\mathbf1_{i=N-1-j},
& (F_V)_{ij}&=c^{N-1}\mathbf1_{i=N-1-j}
&& (i<N,\ j<N),\\
(F_X)_{iN}&=\mathbf1_{i=N-1},
& (F_V)_{iN}&=c^{N-1}\mathbf1_{i=N-1}
&& (i<N),\\
(F_X)_{Nj}&=(F_V)_{Nj}=1 &&&& (j<N),\\
(F_X)_{NN}&=(F_V)_{NN}=\tau^N.
\end{aligned}
\tag{6}
\]

Rows are restricted to \(T\), and columns to \(\widehat S\).
All other auxiliary coordinates vanish from every selected form.
After eliminating the individual old-column pivots, the only special
two-column blocks, with columns ordered \((0,N)\), are

\[
\begin{pmatrix}1&1\\1&\tau^N\end{pmatrix},
\qquad
\begin{pmatrix}c^{N-1}&c^{N-1}\\1&\tau^N\end{pmatrix}.
\tag{7}
\]

Their determinants are respectively \(c^N\) and \(c^{2N-1}\),
both units at zero. Consequently \(F_X\) and \(F_V\) are square
unit matrices. This concentration uses only the already proved old
reflected basis and an unchanged final coordinate; it introduces no
new interpolation assumption.

## 2. The square active tensor and its Weyl parameter

Retain the existing identity direction \(j=0\) with the identity
filter. The active Gaussian and forms are

\[
\begin{aligned}
G_{\rm act}
 &=\prod_{j\in\widehat S}(1+\tau^jY_jU_j)
   \prod_{i\in T}(1+\tau^{i+1}V_iX_i),\\
L_j&=\tau^jY_j+c\sum_{i\in T}(F_V)_{ij}V_i,\\
R_j&=U_j+\tau^{-j}\sum_{i\in T}(F_X)_{ij}X_i.
\end{aligned}
\tag{8}
\]

If \(F\) is the original \(r\)-direction occupation kernel, use
\(F\otimes I_2\) on \(\widehat S\). The removed factor for
\(Q_0=L_0R_0\) cancels its restored identity factor exactly:
\((1-Q_0)(1+Q_0)=1\). Therefore (8) is precisely a regular local
basis representation of the original physical tensor. No local
Gaussian factor is divided out afterward.

Put \(A=\operatorname{diag}(\tau^j)_{j\in\widehat S}\) and
\(D_T=\operatorname{diag}(\tau^{i+1})_{i\in T}\). The inverse
contraction parameter in the
[weighted Weyl theorem](identity-operator-weyl-reduction.md) is

\[
C^{-1}=I+c^{-1}A F_X^{-1}D_TF_V^{-\mathsf T}.
\tag{9}
\]

All inverse form matrices are regular. Every correction in (9) has
positive order, as does every correction in
\(A^{-1}C^{-1}A\). Indeed the first correction includes a row factor
\(\tau^j\), a factor \(\tau^{i+1}\), and only regular remaining
factors. In the conjugated matrix, the same \(A\) factor occurs on
the column instead. Hence both parameters and their inverses are
regular. This verifies the complete weighted hypothesis of that
theorem for every interval and selected core in (5).

The implementation clears only a local zero unit:

\[
\Delta=c\det(F_X)\det(F_V),\qquad
\mathcal N=\Delta I+
A\operatorname{adj}(F_X)D_T\operatorname{adj}(F_V)^{\mathsf T}.
\tag{10}
\]

It constructs \(\Delta^{r+1}\Psi_C(F\otimes I_2)\) from complete
polynomials. Multiplication by this unit changes no zero Smith order.
The denominator clearing in (10) is local profile arithmetic and is
never inserted into the global polynomial crossing.

## 3. Both complete profiles and physical costs

Let \(\sigma\) range over the \(q\)-Smith orders of

\[
\operatorname{diag}_{\mathbf a}
  (q^{P\sum_{j\in\widehat S}j a_j})
\Delta^{r+1}\Psi_C(F\otimes I_2)
\tag{11}
\]

in a fixed occupation block \(|\mathbf a|=h\). For every binary
auxiliary state \(\boldsymbol\epsilon\) on \(\widehat S\), the
corresponding zero-cut order is

\[
2\sigma-2P\sum_{j\in\widehat S}j\epsilon_j
+P(1-J)(h-|\boldsymbol\epsilon|).
\tag{12}
\]

This gives exactly \(2^{2(r+1)}\) active orders. The zero spectators
have orders

\[
P(2j+1-J)\quad(j\notin\widehat S),\qquad
P(2i+1+J)\quad(i\notin T).
\tag{13}
\]

At infinity, only the original selected set \(S\) is replaced. The
active matrix is

\[
\left(F_{\mathbf a\mathbf b}(z^{-2})
 z^{P\sum_{j\in S}(2j+1-2n)b_j}\right)_{\mathbf a,\mathbf b}.
\tag{14}
\]

Every ordinary Gaussian direction outside \(S\), including \(j=0\),
remains an infinity spectator. The paired physical energies are still
\(P(4i-2(n-1)+J)\) and \(P(4i-2(n-1)-J)\), for all
\(0\le i<n\). These full lists are included by the evaluator.

The actual polynomial seed uses only its original \(r\) identity
forms. Occupation-preserving kernels with invertible endpoints that
preserve the physical energy filtration have the same constant
physical-Y endpoint repair as the
[operator construction](coupled-identity-operator-filters.md).
For common Laurent clearing \(q^e\), its ordinary degree is bounded
by

\[
P\bigl(\deg G+2nr\bigr)
+\max_{\mathbf a,\mathbf b}\deg(q^eF_{\mathbf a\mathbf b}).
\tag{15}
\]

The common clearing contributes opposite scalar shifts at the two
cuts. The extra identity direction in the local calculation adds
nothing to (15).

## 4. Implementation and verification

The constructor and complete affine compiler are
[`above_power_terminal_operators.py`](../tools/above_power_terminal_operators.py).
Its `geometry_data`, `zero_profiles`, `infinity_profiles`,
`compile_family`, and `polynomial_seed` interfaces keep the actual
selected directions distinct from the augmented zero directions.

The independent
[`checker`](../tools/check_above_power_terminal_operators.py)
checks the old reflected changes coefficientwise, compares full and
reduced Gram matrices, expands the complete square active tensor
before applying Weyl reduction, and replays the original physical
tensor at both cuts. It also compares complete compiled sources with
literal F4 profiles and includes a degree-600 coefficient in the
second F4 bit plane. The pressure controls independently expand all
spectator and physical partitions.

The saved
[`evidence`](../evidence/above-power-terminal-operators.json)
records six complete active-tensor versus Weyl comparisons, 44
compiled two-cut comparisons, eight original physical cut profiles,
3,062 endpoint coefficient equalities, and 54 independently expanded
pressure evaluations. It also contains 10,960 reflected Pascal
identities, 99 full-versus-reduced Gram identities, and 11,059
weighted-regularity checks. The structural cases reach \(n=65\)
with four filtered directions; the complete active-tensor comparisons
reach \(n=33\) with three. Executed source hashes are retained.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/check_above_power_terminal_operators.py --extended --output /tmp/above-power-terminal-replay.json
```

The construction gives a new finite-dimensional search family. Exact
profiles and numerical pressure searches remain distinct from a
proved improvement of the asymptotic factoring bound.
