# Two unequal Boolean intervals split into scalar Gaussian sectors

## Result and scope

The first two-block extension of the offset-interval construction has
non-scalar holonomy, but that holonomy splits exactly over
\(\mathbb F_2(z)\). After local rational changes of basis, the Gaussian
kernel is a direct sum of two scalar offset kernels. One of their
parameters need not be a power of two.

This gives an explicit extension of the
[offset-interval construction](schur-offset-interval-reduction.md), and
explains why unequal dyadic intervals do not automatically give a new
coupled Gaussian construction. The statement concerns the displayed
family. It does not exclude other down-sets, arbitrary coefficient
recombinations, or an improvement after Schur energy projection. It
establishes no new factoring exponent.

The calculation and exact checks were developed with substantial OpenAI
assistant support. The complete Schur contraction formula and its source
references are in [the fermion query note](schur-fermion-queries.md).

## 1. The two-block family

Let \(N,L,M\) be powers of two with

~~~math
2\le N,\qquad 2N\le L,\qquad 2L\le M.
~~~

Put \(q=1+z\), \(d=L+2N\), and

~~~math
S=[0,2N)\ \cup\ \bigl(L+[0,N)\bigr).
~~~

This is a Boolean down-set of size \(3N\). In the order
\[
S=(0+[0,N))\ \sqcup\ (N+[0,N))\ \sqcup\ (L+[0,N)),
\]
write

~~~math
G_{a,b}=z^{a\mathbin{\&}b},\qquad
H_{a,b}=z^{M-1-(a\mathbin{|}b)}.                           \tag{1}
~~~

Let \(G_N,H_N\) denote the full-interval matrices on \([0,N)\).
Separating the low bits gives

~~~math
G=G_{\mathrm c}\otimes G_N,\qquad
H=H_{\mathrm c}\otimes H_N,\qquad
G_NH_N=q^{N-1}I_N,                                      \tag{2}
~~~

where the three core indices are \(0,N,L\), and

~~~math
(G_{\mathrm c})_{a,b}=z^{a\mathbin{\&}b},\qquad
(H_{\mathrm c})_{a,b}=z^{M-N-(a\mathbin{|}b)}.             \tag{3}
~~~

## 2. Exact holonomy decomposition

Set

~~~math
u=(1,z^N,z^L)^T,\quad
\lambda=1+z^N+z^L,\quad
f=z^{M-L-2N}.
~~~

The full four-element subcube \(\{0,N,L,N+L\}\) has scalar
holonomy. Removing its last row and column gives the rank-one identity

~~~math
G_{\mathrm c}H_{\mathrm c}
 = f\bigl(q^{N+L}I_3+u(1,1,1)\bigr).                    \tag{4}
~~~

The two-dimensional kernel of \((1,1,1)\) has eigenvalue
\(f q^{N+L}\). The vector \(u\) has eigenvalue

~~~math
f(q^{N+L}+\lambda)=f z^{N+L}=z^{M-N}.                    \tag{5}
~~~

Here \(q^{N+L}=(1+z^N)(1+z^L)\). In particular, these two
eigenvalues are distinct rational functions.

An explicit diagonalizing matrix is

~~~math
C=
\begin{pmatrix}
1&1&1\\
1&0&z^N\\
0&1&z^L
\end{pmatrix},
\qquad \det C=\lambda\ne0.                              \tag{6}
~~~

Combining (2)--(6), \(C\otimes I_N\) diagonalizes \(GH\) with

~~~math
\boxed{\quad
c_d=z^{M-d}q^{d-1}\quad\text{with multiplicity }2N,\qquad
c_N=z^{M-N}q^{N-1}\quad\text{with multiplicity }N.
\quad}                                                  \tag{7}
~~~

This is an exact polynomial matrix identity after multiplying by
\(C\otimes I_N\); it does not rely on numerical eigenvalue fitting.
The associated core projector is \(u(1,1,1)/\lambda\).
Its entries are regular at both zero and infinity. This last observation
concerns the projector only; it is not a claim that every boundary
change used below is regular at those two places.

## 3. The whole Gaussian graph reduces with the holonomy

The same argument applies to a general Boolean down-set whenever its
holonomy admits a useful decomposition. Define

~~~math
P_{a,b}=\mathbf1_{b\subseteq a},\qquad
Q=\mathrm{diag}(q^a).
~~~

The relabelled raw Schur contraction matrix is

~~~math
\mathcal R=
\begin{pmatrix}
0&H&P^T&I\\
G&0&q^{-M}PQ&z^{-M}G\\
P^T&q^{-1}Q^{-1}P&0&P^T\\
I&z^{-M}H&P^T&0
\end{pmatrix}.                                          \tag{8}
~~~

Since \(G=PQP^T\) and \(P^2=I\), apply the row and column
changes from the offset reduction, first with its \(q^{N/2}\)
normalization replaced by \(1\). All nonzero blocks become scalar
matrices except

~~~math
(\text{row }V,\text{column }X)
  =z^M q^{-1}(HG)^{-1}.                                 \tag{9}
~~~

Every other edge is one of \(1,q^{M/2},z^M,z^M q^{M/2}\).
Consequently a common similarity change at all four boundaries
diagonalizes the entire graph when it diagonalizes \(HG\). Diagonalizing
one product matrix is enough; no simultaneous diagonalization of
\(G,H,P\) is asserted.

For (7), use the transpose inverse of \(C\otimes I_N\), which
diagonalizes \(HG=(GH)^T\). In the \(c_e\) sector, for
\(e=d\) or \(e=N\), multiply the \(V\) rows and \(X\) columns by
\(q^{e/2}\). Then (9) becomes

~~~math
z^M q^{e-1}/c_e=z^e.                                   \tag{10}
~~~

Both \(d\) and \(N\) are even. The resulting one-particle matrix is a
direct sum of \(2N\) copies of \(K_{M,d}\) and \(N\) copies of
\(K_{M,N}\), where

~~~math
K_{M,e}=
\begin{pmatrix}
0&q^{e/2}&q^{M/2}&1\\
q^{M/2}&0&z^M&q^{M/2}\\
q^{e/2}&z^e&0&z^M q^{e/2}\\
1&q^{e/2}&z^M q^{M/2}&0
\end{pmatrix}.                                          \tag{11}
~~~

The corresponding Gaussian tensor factors over these sectors. The
changes are rational, local to the four boundaries, and preserve total
species occupation. They need not preserve individual index sums.
Therefore this identity alone does not identify the ranks after the
Schur energy projection.

## 4. Direct seeds beyond power-of-two interval lengths

Formula (11) defines a polynomial crossing for every pair of positive
even integers \(M,e\), regardless of whether they arise from intervals.
At \(z=1\), its only nonzero edges are \(B,Y\) and \(X,V\),
both directions with coefficient one.

Encode a physical label by \(0=\varnothing\), \(1=+\), \(2=-\),
\(3=\{+,-\}\), as in the offset note. The selected minors in boundary
states \((\mathrm{flip}(u),\mathrm{flip}(x),v,y)\)
therefore satisfy

~~~math
T_{M,e}(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.           \tag{12}
~~~

Their maximum degree is \(3M+e\). For non-power-of-two \(e\), one
must use the actual binomial polynomial \(q^{e/2}\); replacing it
by \(1+z^{e/2}\) would change the construction.

More generally, both \(q\)-exponents can be varied independently:
replace \(q^{M/2}\) by \(q^a\), and \(q^{e/2}\) by \(q^b\),
for arbitrary positive integers \(a,b\). The endpoint in (12) remains
valid. This gives a four-parameter direct family, without a claimed
Schur interpretation for all parameter choices.

For this last family, writing \(Z=z^M\), \(W=z^e\),
\(A=q^{2a}\), \(B=q^{2b}\), the four useful two-by-two
determinants are

~~~math
\begin{array}{c|c}
\text{cross block}&\text{determinant}\\\hline
R^+&Z(1+A)\\
R^-&W+ZB\\
C^+&Z+A\\
C^-&W+B .
\end{array}                                             \tag{13}
~~~

These formulas expose the cancellation conditions directly and provide
compact inputs for further gauge searches. They are not pressure
certificates.

## 5. Exact reproduction

Run

~~~bash
python tools/check_schur_two_block_holonomy.py \
  --output evidence/schur-two-block-holonomy.json
~~~

The checker verifies the polynomial eigenvector identities in (7),
the twelve raw block transformations over a finite field, direct
crossing endpoints for non-power-of-two parameters, and the determinant
formulas in (13). Finite-field evaluations are supplementary checks of
the displayed rational changes; the holonomy identity itself is checked
as a polynomial identity over \(\mathbb F_2[z]\).
