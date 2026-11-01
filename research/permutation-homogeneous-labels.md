# Large homogeneous permutation label spaces

## Result and scope

For every even integer \(a\ge4\), consider the \(a!\) bipartite
perfect matchings between two sets of \(a\) vertices. Encode each edge
\(e\subseteq\{0,\ldots,2a-1\}\) by the odd mode

\[
n_e=1+2\sum_{j\in e}2^j.
\]

The matching monomials have a common weight
\(h=a+2(2^{2a}-1)\). Their homogeneous vacuum Gram matrix \(G_a\)
over \(\mathbb F_2\) satisfies

\[
\operatorname{rank}G_a\ \ge\
\frac{a!}{147^a}=2^{a\log_2a-O(a)}.
\tag{1}
\]

Consequently one can construct homogeneous unit duals for
\(2^\lambda\) labels using \(a=O(\lambda/\log(\lambda+2))\).
The construction takes \(2^{O(\lambda)}\) time and produces a
scale-one crossing family with

\[
\log_2|\mathcal I|=O(a)
 =O\!\left(\frac{\lambda}{\log(\lambda+2)}\right).
\tag{2}
\]

Every coefficient and exact rank factorization can be generated in
\(2^{O(\lambda)}\) time. This meets the piece-count and effectivity
conditions of [subexponential-piece-transfer.md](subexponential-piece-transfer.md).
It does not prove a joint-width bound below the existing half-rate
threshold. Gaussian generation here is a time statement; no polynomial
space bound is claimed for storing these matrices.

The full permutation Gram matrix need not be invertible. It is already
singular at \(a=10\). Its square also need not be the identity, even when
it is invertible. Both failures are recorded as negative controls.

The source operators are from
[Anthropic Appendix B, revision 16a5ce41304fcdaea4e0be7c132e1c361dc26da5](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The finite matching formula is proved in
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md).
The representation-theoretic input is the integral Jucys--Murphy
center theorem, stated explicitly in the introduction of
[Francis--Jones, arXiv:0709.0326](https://arxiv.org/pdf/0709.0326),
and its content-evaluation theorem, stated as Theorem 3.1 of
[Novak, arXiv:0811.3595](https://arxiv.org/pdf/0811.3595).
These ordinary mathematical inputs are separate from the pinned Lean
source. No new Lean formalization is asserted.

## 1. The Gram matrix is a central convolution

Write the matching of \(\sigma\in S_a\) as
\(\{\{i,a+\sigma(i)\}:0\le i<a\}\), and its monomial as \(b_\sigma\).
Only annihilation contributes to the vacuum coefficient, so

\[
G_a(\sigma,\tau)
 =[1]\,b_\sigma[-h]b_\tau
 =\det\bigl([e\cap f=\varnothing]\bigr)_{e\in\sigma,f\in\tau}
 \quad\text{over }\mathbb F_2.
\tag{3}
\]

The determinant equals the permanent in characteristic two.
Relabel the two vertex sets so that the relative permutation is
\(\pi=\sigma^{-1}\tau\). Let \(P_\pi\) be its permutation matrix and
let \(F_\pi\) be diagonal with a one at each fixed point. Then the
disjointness matrix is

\[
J+C,\qquad C=I+P_\pi+F_\pi,
\tag{4}
\]

where \(J\) is the all-ones matrix. The extra \(F_\pi\) is essential:
a common matching edge intersects in two vertices, but it is still
not disjoint.

Each fixed point contributes a \(1\times1\) identity block to \(C\).
A nontrivial cycle of length \(r\) contributes \(I+P_r\), of rank
\(r-1\), with both left and right kernels spanned by the all-ones
vector. Therefore:

* two or more nontrivial cycles give nullity at least two, which the
  rank-one update \(J\) cannot remove;
* one nontrivial \(r\)-cycle gives determinant \(r\bmod2\), by the
  adjugate form of the determinant lemma;
* no nontrivial cycle gives \(\det(I+J)=1+a\bmod2\).

Thus for even \(a\), \(G_a\) is convolution by the reduction of

\[
\mathcal G_a
 =1+\sum_{\substack{3\le r\le a\\r\text{ odd}}}
       C_{(r,1^{a-r})}\ \in Z(\mathbb Z[S_a]),
\tag{5}
\]

where \(C_\mu\) denotes the conjugacy class sum.

## 2. Balanced content parities force an odd central eigenvalue

The integral center theorem says that every element of
\(Z(\mathbb Z[S_a])\) is an integral symmetric polynomial in the
Jucys--Murphy elements. Write \(\mathcal G_a=F(J_1,\ldots,J_a)\),
with \(F\in\mathbb Z[X_1,\ldots,X_a]\) symmetric.

On the ordinary irreducible representation of shape \(\nu\vdash a\),
the scalar is

\[
\kappa_\nu=F\bigl(\{j-i:(i,j)\in\nu\}\bigr)\in\mathbb Z.
\tag{6}
\]

If the Young diagram has \(a/2\) even contents and \(a/2\) odd contents,
its content multiset modulo two equals that of the one-row diagram
\((a)\). Symmetry and integrality of \(F\) therefore give
\(\kappa_\nu\equiv\kappa_{(a)}\pmod2\).
Every row of even length has equally many even and odd contents, so
all diagrams whose row lengths are even satisfy this condition.

The one-row representation is trivial, and hence

\[
\kappa_{(a)}
 =1+\sum_{\substack{3\le r\le a\\r\text{ odd}}}
       \binom ar(r-1)!
 \equiv1\pmod2.
\tag{7}
\]

This proof uses content residues directly and requires no assertion
about modular semisimplicity or lifting an eigenbasis modulo two.

## 3. A large rank survives reduction

Let \(f^\nu\) denote the ordinary dimension of shape \(\nu\).
The integral regular-convolution matrix has characteristic polynomial

\[
\det(tI-\mathcal G_a)
 =\prod_{\nu\vdash a}(t-\kappa_\nu)^{(f^\nu)^2}.
\tag{8}
\]

Reduce this identity modulo two. A matrix of rank \(r\) on a
space of dimension \(D\) has characteristic polynomial divisible by
\(t^{D-r}\). Therefore

\[
\operatorname{rank}_{\mathbb F_2}G_a
 \ge\sum_{\kappa_\nu\text{ odd}}(f^\nu)^2.
\tag{9}
\]

In particular, take \(b\) to be the largest even integer at most
\(\sqrt a\), and put

\[
\nu=(b+a-b^2,b,\ldots,b),
\tag{10}
\]

with \(b\) rows. All rows have even length. Also
\(a-b^2\le4b+2\), so every hook is at most
\(a-b^2+2b-1\le6b+1\le7\sqrt a\).
The hook-length formula gives

\[
(f^\nu)^2
 \ge\frac{(a!)^2}{(49a)^a}
 \ge\frac{a!}{147^a},
\tag{11}
\]

where \(a!\ge(a/e)^a\) and \(e<3\).
Equations (9)--(11) prove (1). The explicit hook dimension
\((f^\nu)^2\), rather than the weaker last expression, is available
for selecting practical parameters.

## 4. Select labels and compute homogeneous duals

Form \(G_a\) and find a rank-sized nonsingular minor \(G_a[R,C]\).
All its rows and columns are permutation monomials of weight \(h\).
If \(H=G_a[R,C]^{-1}\), define

\[
a_x=\sum_{r\in R}H_{x,r}b_r,\qquad b_v=b_v\ (v\in C).
\tag{12}
\]

Then \([1]a_x[-h]b_v=[x=v]\).
Retain any \(2^\lambda\) columns of \(C\) and the corresponding
already-computed dual rows. The retained duals can still use every
monomial in \(R\); taking an arbitrary smaller principal Gram submatrix
and assuming it remains invertible would be incorrect.

For each \(\lambda\), choose an even \(a\) with the explicit
\((f^\nu)^2\ge2^\lambda\), treating finitely many small lengths
directly. Equation (11) permits \(a\log a=O(\lambda)\).
Enumeration of \(a!\) permutations, the Gram matrix, and ordinary
finite-field Gaussian elimination all take \(2^{O(\lambda)}\) time.
Every mode index has \(O(a)\) bits, and every dual has at most \(a!\)
monomials. All duals remain homogeneous of weight \(h\).

## 5. An explicit multigraded crossing lift

For an edge \(e\), write
\(Z^e=\prod_{j\in e}Z_j\) and
\((1+Z)^e=\prod_{j\in e}(1+Z_j)\).
Use the four vertex groups \(U,B,A,Y\), one vertex for each matching
edge. Define a symmetric matrix with zero diagonal and zero entries
within each group, and the following cross-group entries:

| Pair | Entry |
| --- | --- |
| \(U_u,B_b\) | \([u=b] (1+L)(1+Z)^b\) |
| \(U_u,A_e\) | \([u=e]L(1+L)Z^u(1+Z)^e\) |
| \(U_u,Y_y\) | \([u=y]LZ^y\) |
| \(B_b,A_e\) | \([b\cap e=\varnothing]LZ^e\) |
| \(B_b,Y_y\) | \([b\cap y=\varnothing]L(1+L)Z^y(1+Z)^b\) |
| \(A_e,Y_y\) | \([e\cap y=\varnothing] (1+L)(1+Z)^e\) |

Its Pfaffian is a polynomial \(\mathcal P(L,Z)\) over \(\mathbb F_2\).
Under

\[
L=z,\qquad Z_j=z^{2^{j+1}},
\tag{13}
\]

it equals the source polynomial
\(z^{2h}\sum_s[t^u]Q_s(t^e,t^b)t^y z^s\).
This follows term by term from the six-block table in
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md):
equal-size edge masks turn the creation test into equality, while
annihilation is disjointness.

At \(L=Z_j=1\), only the \(UY\) and \(BA\) blocks survive.
Consequently

\[
\mathcal P(1,1,\ldots,1)
 =[u=y]\,G_a(e,b).
\tag{14}
\]

Extend the \(A\) argument linearly using (12). The sum of all
multigraded coefficient tensors is then
\([u=y][x=v]\), the scale-one crossing.

Every physical vertex belongs once to each of the four matching labels.
The table assigns each occurrence at most one factor \(Z_j\) or
\((1+Z_j)\). Hence each variable \(Z_j\) has degree at most four.
There are \(2a\) matched edges in a Pfaffian term, each of \(L\)-degree
at most two. Thus

\[
|\mathcal I|\le(4a+1)5^{2a}.
\tag{15}
\]

Subset matching uses \(2^{O(a)}\) states, and the bounded polynomial
box in (15) has \(2^{O(a)}\) monomials. Straight polynomial arithmetic
therefore computes a monomial-label entry in \(2^{O(a)}\) time.
Summing the at most \(a!\) dual monomials, generating all external
entries, and exact Gaussian factorization remain \(2^{O(\lambda)}\).
This proves (2) and the claimed effectivity.

The univariate specialization is also usable: \(h=2^{O(a)}\), so its
whole shift interval has only \(2^{O(a)}\) positions. Splitting into
independent multidegrees can change ranks, so no rank comparison
between these two coefficient families is presumed.

## 6. An exact particle-count constraint on the coefficient profile

In any nonzero perfect-matching term, the numbers of the six edge
types satisfy

\[
\#UB=\#AY=b,\qquad
\#UA=\#BY=c,\qquad
\#UY=\#BA=d,\qquad b+c+d=a.
\]

This is just the equality of the four vertex-group sizes. Its complete
\(L\) factor is

\[
L^{2(c+d)}(1+L)^{2(b+c)}
 =T^{a-b}(1+T)^{a-d},\qquad T=L^2.
\tag{16}
\]

In particular, no odd \(L\) degree occurs. A coefficient of \(T^k\)
can receive a term only if

\[
a-b\le k\le a+c.
\tag{17}
\]

For \(k<a\), this forces at least \(a-k\) internal \(AY\) and
\(UB\) contractions in one flattening. For \(k>a\), it forces at
least \(k-a\) internal \(UA\) and \(BY\) contractions in the other.
At the central coefficient \(k=a\), the binomial multiplier is
\(\binom{b+c}{b}\bmod2\), equal to one precisely when \(b\mathbin{\&}c=0\).
This describes a structural target for rank compression. It does not
on its own bound the rank of the sum of all surviving patterns.

## 7. Finite verification and failed stronger claims

The checker compares the Gram formula against direct disjointness
determinants, computes small exact Gram ranks, verifies central
eigenvalue parity by the independent Murnaghan--Nakayama character
rule, and checks the multigraded lift against the direct source
polynomial. It also checks balanced hook shapes and the particle-count
constraint.

Two negative controls prevent extrapolating the small invertible cases:

* at \(a=8\), the coefficient of cycle type \((7,1)\) in
  \(\mathcal G_a^2\) is one modulo two, so \(G_a^2\ne I\);
* at \(a=10\), shape \((4,3,2,1)\) has central eigenvalue \(-176\).
  Equation (8) then makes the determinant even, so \(G_{10}\) is
  singular modulo two.

The output is saved in
[evidence/permutation-homogeneous-labels.json](../evidence/permutation-homogeneous-labels.json).
