# Prime-power Gram involutions and streamed Schur factors

## Result

At every prime-power weight \(h=p^r\), the homogeneous vacuum Gram
matrix of the full integral Schur lattice is an involution modulo \(p\):

\[
B_h^2=I\quad\text{on }\Lambda_{\mathbb F_p,h}.
\tag{1}
\]

The homogeneous unit dual is therefore explicitly \(a_x=B_hs_x\).
There is no need to invert a partition-sized matrix.

The construction also admits factor queries in
\(\exp(O(\sqrt h))\) time and \(\operatorname{poly}(h)\) bit space.
The proof supplies streamed integer algorithms for symmetric-group
characters and source mode coefficients. They avoid a table of all
partitions and division by a nonunit in the coefficient field.

With the prime-power capacity construction in
[integral-schur-crossing.md](integral-schur-crossing.md), this gives
full-Schur factor queries in \(2^{O(\lambda)}\) time and polynomial
space in \(\lambda\). It preserves the existing half-rate factoring
bound and polynomial-space implementation. Neither inversion alone
nor these coefficient algorithms improve the joint width exponent.

The operators come from
[Anthropic Appendix B at the pinned revision](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The integral lattice, homogeneous Gram identity, and characteristic-zero
primitive decomposition are proved in
[integral-schur-crossing.md](integral-schur-crossing.md),
[schur-fusion-gram.md](schur-fusion-gram.md), and
[homogeneous-duals.md](homogeneous-duals.md).
The additional square identity below is integral before prime reduction.

## 1. Operators and integral divided powers

On the power sums, put

\[
D(p_m)=
\begin{cases}m p_{m-1},&m>1,\\0,&m=1,\end{cases}
\qquad U(p_m)=mp_{m+1},\qquad H(p_\mu)=|\mu|p_\mu.
\]

Both \(D,U\) are derivations, and

\[
[D,U]=2H,\qquad [H,U]=U,\qquad [H,D]=-D.
\tag{2}
\]

Write \(D^{(k)}=D^k/k!\) and \(U^{(k)}=U^k/k!\). These are
integral operators on the symmetric-function lattice. With
\(\mathcal H(u)=\sum_mh_mu^m\), their generating transformations are

\[
e^{zD}\mathcal H(u)=\mathcal H\!\left(\frac{u}{1-zu}\right),
\qquad
e^{zU}\mathcal H(u)=\frac{\mathcal H(u+z)}{\mathcal H(z)}.
\tag{3}
\]

The coefficients are integral polynomials in the complete and elementary
symmetric functions. Thus the divided powers can be reduced modulo any
prime without inverting a factorial.

Let \(\omega(p_m)=(-1)^{m-1}p_m\). For \(h>0\), the full homogeneous
Gram operator is

\[
G_h=\sum_{k=0}^{h-1}U^{(k)}D^{(k)},\qquad B_h=G_h\omega.
\tag{4}
\]

Every \(D^{(k)}\) with \(k\ge h\) vanishes on degree \(h\):
downward translation of a positive-degree function has no constant
component. The truncation is exact, including degree one.
The Gram identity means
\([1]a[-h]b=\langle B_ha,b\rangle\) for the Hall form.

The involution \(\omega\) anticommutes with \(U,D\), so it commutes
with every \(U^{(k)}D^{(k)}\). Consequently \(B_h^2=G_h^2\).

## 2. An integral square identity

The divided-power form of (2) gives

\[
D^{(a)}U^{(b)}
=\sum_{j=0}^{\min(a,b)}
U^{(b-j)}
\binom{2H+a+b-j-1}{j}
D^{(a-j)}.
\tag{5}
\]

The polynomial in \(H\) is placed between the displayed operators.
Equation (5) follows by induction from (2); its first-degree recurrence
is \(DU^{(b)}=U^{(b)}D+U^{(b-1)}(2H+b-1)\).
One can first work over the rationals and then use the equality of
integral operators.

In \(U^{(a)}D^{(a)}U^{(b)}D^{(b)}\) on degree \(h\), put
\(k=a+b-j\). The term from (5) is \(U^{(k)}D^{(k)}\) times

\[
\binom ka\binom kb\binom{2h-k-1}{a+b-k}.
\]

The shift of \(H\) matters: the rightmost \(D^{(k)}\) first lowers
the degree to \(h-k\). Terms with \(k\ge h\) vanish. For \(k<h\),
two Vandermonde sums give

\[
\begin{aligned}
&\sum_{\substack{0\le a,b\le k\\a+b\ge k}}
\binom ka\binom kb\binom{2h-k-1}{a+b-k}\\
&\quad=\sum_{u=0}^k
\binom{2k}{k+u}\binom{2h-k-1}{u}
=\binom{2h+k-1}{k}.
\end{aligned}
\tag{6}
\]

Therefore the exact integral identity is

\[
\boxed{G_h^2=
\sum_{k=0}^{h-1}\binom{2h+k-1}{k}U^{(k)}D^{(k)}.}
\tag{7}
\]

As a separate scalar check, on a primitive chain of lowest degree \(j\),
put \(t=h-j\). Equation (7) becomes

\[
\binom{2h-1}{t}^2
=\sum_k\binom{2h+k-1}{k}\binom tk\binom{2h-t-1}{k}.
\tag{8}
\]

The proof of (7) does not require an eigenbasis to survive reduction.

## 3. Prime-power reduction

If \(h=p^r\), then in characteristic \(p\),

\[
(1-z)^{-2h}=(1-z^h)^{-2}.
\]

The coefficient \(\binom{2h+k-1}{k}\) therefore vanishes modulo \(p\)
for \(0<k<h\). Equation (7) reduces to \(G_h^2=I\), giving (1).
The integral Schur matrix \(B_h\) is symmetric. Hence

\[
a_x=B_hs_x,\quad b_v=s_v
\quad\Longrightarrow\quad
[1]a_x[-h]b_v
=\langle B_h^2s_x,s_v\rangle=[x=v].
\tag{9}
\]

This applies to the **full** degree space. Restricting external labels
afterwards retains duality. An arbitrary principal submatrix of \(B_h\)
is not asserted to be an involution.

## 4. A streamed character algorithm

Let \(\lambda,\mu\) partition \(n\), and write \(m_k=m_k(\mu)\).
The Frobenius expansion is

\[
s_\lambda=\sum_{\nu\vdash n}\frac{\chi^\lambda_\nu}{z_\nu}p_\nu,
\qquad z_\mu=\prod_k k^{m_k}m_k!.
\tag{10}
\]

Set every power sum outside the support of \(\mu\) to zero. Apply the
forward difference of order \(m_k\) in each remaining \(p_k\), at zero.
A monomial survives only when all its exponents are at least \(m_k\).
Weighted homogeneity forces equality, since its weighted degree is
\(\sum k m_k=n\). Thus

\[
\chi^\lambda_\mu=
\left(\prod_k k^{m_k}\right)
\sum_{0\le t_k\le m_k}
\left(\prod_k(-1)^{m_k-t_k}\binom{m_k}{t_k}\right)
s_\lambda(p_k=t_k).
\tag{11}
\]

The factorials from the forward differences cancel those in \(z_\mu\).
The remaining rational sum is taken exactly before obtaining the integer.

There are \(\prod_k(m_k+1)\) evaluations. Each tuple describes a
subpartition of \(\mu\), so their number is at most
\(\sum_{j\le n}p(j)=\exp(O(\sqrt n))\).

At each point, compute the complete symmetric functions by

\[
h_0=1,\qquad jh_j=\sum_{k=1}^j p_kh_{j-k},
\]

then the Jacobi--Trudi determinant. It has at most \(n\) rows and uses
only \(h_j\) with \(j\le n\). Each denominator divides \(j!\).
Multiplying the matrix by \(n!\) gives integer entries. Their values,
determinants, and the running sum have polynomial-in-\(n\) bit lengths.
Fraction-free elimination takes polynomial bit time and space per point.

The grid and accumulator are streamed. No partition-indexed table or
rim-hook recursion cache is required. One character therefore takes
\(\exp(O(\sqrt n))\) time and \(\operatorname{poly}(n)\) space.

## 5. A streamed source-mode algorithm

For partitions \(\mu,\nu,\rho\) of weights at most \(W\), consider
the coefficient of \(p_\rho\) in the source operator
\((p_\mu)[r]p_\nu\), where \(r=|\rho|-|\nu|\).

Enumerate common subpartitions \(\tau\) of \(\nu,\rho\). These are
the particles left untouched by the normally ordered operator. Put
\(x_k=m_k(\rho)-m_k(\tau)\) and
\(y_k=m_k(\nu)-m_k(\tau)\).
Unless \(\sum_k(x_k+y_k)=\ell(\mu)\), the tail contributes zero.
Otherwise form

\[
P(X,Y)=\prod_{m\in\mu}
\left(\sum_k C(m,k)X_k+\sum_k A(m,k)Y_k\right),
\tag{12}
\]

where the exact source coefficients are

\[
C(m,k)=\binom{k-1}{m-1},\qquad
A(m,k)=(-1)^{m-1}k\binom{k+m-1}{m-1}.
\tag{13}
\]

A nonnegative upper argument smaller than the lower argument gives zero.
The contribution of \(\tau\) is

\[
[X^xY^y]P(X,Y)
\prod_k\frac{m_k(\nu)!}{m_k(\tau)!}.
\tag{14}
\]

This is normal ordering: the \(X\) choices create particles, the
\(Y\) choices differentiate the input, and the falling factorial
includes every annihilation multiplicity.

Extract the coefficient by finite differences over
\(0\le t_k\le x_k\), \(0\le u_k\le y_k\), dividing the final
integer numerator by \(\prod_k x_k!y_k!\).
Ordinary homogeneity of (12) makes this exact. Each grid has at most

\[
\prod_k(x_k+1)(y_k+1)=\exp(O(\sqrt W))
\]

points, since the two factors count subpartitions of partitions of
weight at most \(W\). The number of tails has the same bound.

Evaluate (12) as a product of linear forms. The binomials in (13)
have \(O(W+\log(W+2))\) bits, and all products, factorials, and
accumulators have polynomially many bits. Stream the tails and grids.
This proves \(\exp(O(\sqrt W))\) time and polynomial space.
It does not store every state of the source normal-order recursion.

## 6. Full-Schur factors without stored inverses

Combining (10) with the two coefficient algorithms, a Schur-basis source
entry of degrees \(O(h)\) is a sum over a fixed number of partitions
of degrees \(O(h)\). Enumerate them sequentially. Their number is
\(\exp(O(\sqrt h))\), and a common denominator divides a product
of factorials of \(O(h)\), with polynomial bit length.

In particular, each integral entry
\(B_h(x,y)=[1]s_x[-h]s_y\) is computed in the same time and
polynomial space, then reduced modulo \(p\).
Equation (9) supplies each dual coefficient on demand. A source entry
of the dual is one further sum over \(y\), keeping a fixed number
of nested partition enumerations. Its time is still
\(\exp(O(\sqrt h))\) and space polynomial in \(h\).

Every homogeneous half-piece factor is such a mode coefficient times
an integral angle binomial. Channel partitions have degree at most
\(4h\). Their counts and rank/unrank metadata use a polynomial-in-\(h\)
partition dynamic-programming table, not a table indexed by every
partition.

For arbitrary binary label lengths, use the capacity-one block and
greedy prime-power tensor blocks from the integral-Schur note.
Their capacities sum to \(\lambda\), there are \(O(\log(\lambda+2))\)
blocks, and the sum of their degree square roots is \(O(\lambda+1)\).
Multiply block factor values sequentially. The result takes
\(2^{O(\lambda)}\) time and polynomial space in \(\lambda\).
The full tuple of pieces and the channel-width slack both have
\(O(\log^2(\lambda+2))\) overhead in the relevant logarithmic units.

Use the streamed evaluator's volume cutoff with a fixed-power top cap.
The relaxed piece transfer preserves the exact half endpoint; streamed
execution and bit-tagged witness recovery preserve polynomial space.
These consequences use the explicit factor queries above, not generic
Gaussian factors or a unit-cost unbounded coefficient representation.

## 7. Verification

Run:

    python3 tools/check_prime_power_gram_involution.py --output evidence/prime-power-gram-involution.json

The checker verifies (7) as an equality of integer operators, compares
the streamed character algorithm with independent rim-hook characters,
and compares the streamed source algorithm with the literal recursion.
It also constructs Schur Gram matrices over several prime fields and
checks their squares, including weight eight in characteristic two,
where the matrix is not the identity.

The finite matrices in the checker are for independent comparison.
The algorithms whose space is bounded above stream their summands and
do not retain those matrices.
