# Mixed-end identity cores of arbitrary fixed size

## Result and scope

Let the interval length \(n\) be a power of the characteristic
\(p\). Select \(a\) initial and \(b\) terminal polynomial
identity directions, with \(r=a+b\le n\). Their arbitrary
joint occupation filters and occupation-preserving operator kernels
have a regular zero reduction to \(2r\) active modes. Retain
the \(b\) lowest and \(a\) highest auxiliary \(X,V\) modes.
The [Weyl reduction](identity-operator-weyl-reduction.md) then gives
a matrix of order \(p^r\), together with explicit digit shifts.

For \((a,b)=(1,2)\), the selected directions are
\(Q_0,Q_{n-2},Q_{n-1}\), and the retained auxiliary indices
are \(0,1,n-1\). In characteristic two its zero profile therefore
needs an \(8\)-by-\(8\) reduced matrix. This supplies two low
auxiliary energies and the initial direction's larger infinity
degree budget. No pressure improvement is asserted.

## 1. The form matrix and retained minor

Use \(\tau=z^{2P}\), \(c=1-\tau\), and the polynomial
forms

\[
f_j(i;c)=\sum_h\binom jh\binom ih c^h.
\tag{1}
\]

The selected columns are

\[
S=\{0,\ldots,a-1\}\cup\{n-b,\ldots,n-1\},
\tag{2}
\]

and the retained rows are

\[
T=L\cup H,\qquad
L=\{0,\ldots,b-1\},\quad H=\{n-a,\ldots,n-1\}.
\tag{3}
\]

Write the full form matrix as \(E=(I\mid J)\), separating
the initial and terminal column sets. Let \(F=E_T\).
The terminal reversal identity, valid because \(n\) is a
power of \(p\), is

\[
f_{n-1-h}(i;c)=\tau^i f_h(i;-c/\tau).
\tag{4}
\]

For \(0\le i,h<b\), the constant coefficient of (4) is
\((-1)^i\binom hi\). Thus \(J_L\), with its terminal
columns in reversed order if needed, is triangular modulo
\(\tau\) with nonzero diagonal. It is a regular unit.

At \(\tau=0\), the initial columns satisfy

\[
f_j(i;1)=\binom{i+j}{j}.
\tag{5}
\]

The consecutive \(a\)-by-\(a\) matrix \(I_H(0)\) has
determinant one over the integers. For example, it factors into
a consecutive Pascal evaluation matrix and a triangular Pascal
matrix, both of determinant one. It is consequently a unit in
every characteristic, including \(a\ge p\).

Every row in \(H\) has index at least \(b\), so all entries
of \(J_H\) vanish modulo \(\tau\). These two unit diagonal
blocks prove that \(F\) is a regular unit. Empty initial or
terminal blocks recover the already-proved one-end cases.

## 2. The interpolation has the required weighted orders

Terminal-only interpolation on the first \(b\) rows gives

\[
J_i=r_iJ_L,\qquad
(r_i)_t=\tau^{i-t}(B_0)_{it},\quad t\in L,
\tag{6}
\]

where \(B_0\) is the constant consecutive-Pascal interpolation
matrix. The same formula gives \(J_H=R J_L\), with
\(R_{ut}=\tau^{u-t}(B_0)_{ut}\) for \(u\in H\).
This is the terminal interpolation identity from Section 6 of the
[coupled reduction](coupled-identity-filter-reduction.md).

For a discarded middle index \(b\le i<n-a\), write its
interpolation row as \(B_i=(B_{i,L},B_{i,H})\). The terminal
columns give

\[
B_{i,L}=r_i-B_{i,H}R.
\tag{7}
\]

The initial columns then give

\[
B_{i,H}(I_H-RI_L)=I_i-r_iI_L.
\tag{8}
\]

Every entry of \(R\) has positive order, so
\(I_H-RI_L\) reduces to the invertible matrix \(I_H(0)\).
Equation (8) proves that \(B_{i,H}\) is regular. Equation (7)
now gives, for every low retained index \(t\),

\[
\operatorname{ord}_\tau B_{i,t}\ge i-t.
\tag{9}
\]

Indeed, the first term of (7) has this order bound, while every
term from \(B_{i,H}R\) has order at least \(u-t>i-t\)
because each \(u\in H\) exceeds the discarded index \(i\).
For high retained indices, ordinary regularity is enough.
Therefore, for all retained \(t\),

\[
\operatorname{ord}_\tau B_{i,t}\ge\max(0,i-t).
\tag{10}
\]

The auxiliary diagonal Gaussian orders differ by \(2P(i-t)\).
Thus (10) proves that both \(B\) and
\(D_I^{-1}BD_T\) are regular. The weighted concentration
criterion applies and leaves exactly \(2r\) active modes.
Its basis changes and inverse Gaussian companions are regular.

## 3. Polynomial active matrices

The implementation computes the retained minor and interpolation
without using an unproved symbolic simplification. Put

\[
\delta=\det F,\qquad
N_i=E_i\operatorname{adj}(F),\qquad B_i=N_i/\delta.
\tag{11}
\]

The determinant \(\delta\) is a zero unit. After the regular
concentration, scale every retained \(V\) variable by
\(\delta\). The complete active coefficient matrices become

\[
F_X=F,\qquad
(F_V)_{tj}=\delta F_{tj}
 +\sum_{i\notin T}\tau^{t-i}(N_i)_t E_{ij}.
\tag{12}
\]

Every power in (12) is nonnegative by (10). The auxiliary Gaussian
diagonal is simultaneously multiplied by the unit \(\delta\).
These polynomial data retain all local unit factors; no global
endpoint denominator is removed.

The original crossing still uses the polynomial identity summands
\(Q_j\), for \(j\in S\), and the same Laurent kernel
\(F(q)\) and endpoint repair as in the
[operator construction](coupled-identity-operator-filters.md).
Thus its endpoint, degree bound, and infinity matrix are unchanged
by this local reduction. The full infinity active orders come from
the selected energies \(P(2j+1-2n)\), for \(j\in S\).

## 4. Exact constructor and scope of its checks

[split_identity_mode_filters.py](../tools/split_identity_mode_filters.py)
constructs (11)--(12), checks that the retained determinant is a
unit, and checks (10) coefficient by coefficient for every discarded
row. The existing operator and Weyl interfaces accept geometry names
such as `split:1:2` and `split:2:2`.

For prime-power lengths, Sections 1--2 prove these conditions
uniformly. The constructor also accepts another length when its
finite exact checks pass; in that case the weighted concentration
criterion still supplies a valid certificate for that particular
length. It does not infer the reversal identity (4) outside its
stated prime-power scope.

The independent checker compares the complete original physical
tensors, the earlier full active matrices, and the smaller Weyl
profiles. It also checks Gaussian recovery and the weighted
interpolation certificates at larger lengths. The resulting finite
screening family is broader than the two-direction extreme pair,
but it is not itself an improved complexity bound.
