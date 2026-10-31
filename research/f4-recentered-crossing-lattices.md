# Complete crossing synthesis from nonbinary reference frames

## Result and scope

The complete-entry affine constructor now accepts a reference tensor over
the four-element field, rather than requiring a binary reference. A retained
nonbinary tensor can therefore supply the regular frames for the next exact
coefficient search. Every entry coefficient remains free, subject to the
strict crossing endpoint, the stated polynomial degree, and the selected
local row conditions.

The new implementation is
[`f4_recentered_lattices.py`](../tools/f4_recentered_lattices.py). Its independent
finite controls include binary restriction of scalars, complete Toeplitz
ranks, direct witness substitution, and agreement with the previous binary
constructor. This is a construction tool and an ordinary algebraic argument.
It does not establish a factoring exponent above one half.

The predecessor is the
[complete-entry affine construction](full-alphabet-two-cut-affine-synthesis.md),
with the [finite-precision frame method](finite-precision-smith-frames.md).
The previous coefficient/height searches could reuse selected binary planes
when those planes happened to retain the same profiles. That condition is
not required here. The construction and checks used substantial OpenAI
assistant support and retain the predecessor's source attribution.

## 1. Full coefficient space over the four-element field

Write

\[
F=\mathbf F_2[\omega]/(\omega^2+\omega+1).
\]

Let \(B(w;x,y,u,v)\) be any \(F\)-polynomial tensor of degree at most
\(D\), on an alphabet of size \(q\), with

\[
B(1;x,y,u,v)=[x=v][y=u].
\tag{1}
\]

Every other tensor of degree at most \(D\) with the same endpoint has
the unique form

\[
T=B+(1+w)Q,\qquad \deg Q_{xyuv}<D.
\tag{2}
\]

Indeed, each entry of \(T-B\) vanishes at one. Division by the monic
polynomial \(w+1\) is exact and lowers degree by one. Thus (2) has exactly
\(Dq^4\) independent unknowns over \(F\), irrespective of the
reference tensor's coefficient colors. No occupation, parity, or Gaussian
condition is imposed.

At infinity use the cut \((x,u)\mid(y,v)\); at zero use
\((x,y)\mid(u,v)\). Eight independent integer physical height lists are
retained. If \(h_e\) denotes the sum of the four heights on one entry,
the matrices used for the finite equations are

\[
A_0(T)_{xyuv}=z^{c_0+h_0(x,y,u,v)}T(z;x,y,u,v),
\]

\[
A_\infty(T)_{xyuv}
 =z^{c_\infty+h_\infty(x,y,u,v)}T(z^{-1};x,y,u,v),
\]

where

\[
c_0=-\sum_l\min_a h_{0,l}(a),\qquad
c_\infty=D-\sum_l\min_a h_{\infty,l}(a).
\tag{3}
\]

These are polynomial matrices for every choice of the unknown coefficients.
Both scalar shifts are removed from the reported profiles. The original
degree \(D\) remains explicit throughout the search.

## 2. Exact finite polynomial frames certify complete orders

Let \(A\) be an \(n\)-square polynomial matrix, where \(n=q^2\).
The frame routine computes actual polynomial matrices \(L,R\), with
invertible constant terms, such that

\[
LAR\equiv\mathrm{diag}(d_1,\ldots,d_n)\pmod {z^P},
\qquad \sigma_i=\mathrm{ord}(d_i)<P,
\tag{4}
\]

and the \(\sigma_i\) are sorted.

Every elimination uses a least-order pivot. Writing its value as
\(z^a u\), with \(u(0)\ne0\), an off-diagonal entry is eliminated by
multiplying its row or column by \(u\) and subtracting the pivot row or
column times the corresponding polynomial quotient. These operations are
regular units. Truncating their polynomial representatives modulo \(z^P\)
retains their constant terms and all congruences below \(P\). The routine
checks (4) and the ranks of \(L(0),R(0)\) explicitly.

Equation (4) certifies the complete profile of the untruncated matrix. In
fact, divide row \(i\) of the actual product \(LAR\) by
\(z^{\sigma_i}\). Every error term has order at least
\(P-\sigma_i>0\), and the remaining diagonal entries have nonzero
constant terms. The resulting matrix is a regular unit. Hence

\[
LAR=\mathrm{diag}(z^{\sigma_1},\ldots,z^{\sigma_n})U,
\qquad U(0)\text{ invertible}.
\tag{5}
\]

The \(\sigma_i\) are therefore the full Smith orders. A proposed precision
that misses any of the \(n\) pivots is rejected. The native profile engine
supplies a useful precision hint; the exact polynomial identity and (5)
are the independent complete-profile check.

For a requested row vector \(\tau\), the constructor chooses
\(P\ge\max_i\tau_i\). It then requires

\[
[z^k]\bigl(L A_e(T)R\bigr)_{ij}=0
