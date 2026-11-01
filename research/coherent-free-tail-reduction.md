# A fixed active matrix for an arbitrarily long free tail

This note gives an exact reduction for a coherent core followed by a full
exterior-algebra tail. It preserves the complete Smith profiles of both
crossing cuts. It also identifies an order-of-limits issue: a favorable
determinant moment of a finite core does not by itself give a favorable
infinite-tail pressure.

The reduction leaves a fixed active matrix and an explicit Gaussian
spectator product. It therefore makes the remaining infinite-tail search a
finite algebraic problem. No factoring exponent improvement is claimed.

## 1. Setup

Work over the characteristic-two series ring and use the deflated variable
\(t=w^2\). For a dyadic mode count \(s\), write

\[
B_s=\begin{pmatrix}D_s&J_s\\J_s&D_s\end{pmatrix},\qquad
(D_s)_{ab}=(1+t)t^{2(a\mathbin{|}b)+1},\qquad
(J_s)_{ab}=\delta_{ab}t^{2a+1}.
\]

Here the two blocks are the two physical legs of a crossing cut. The
complete one-particle Smith orders of \(B_s\) are
\(1,3,\ldots,4s-1\). Exterior powers turn a direct sum of one-particle
matrices into a tensor product of complete compound matrices.

Fix a dyadic core size \(k\) and a matched pair of maps \(P,L\) with
\(LP=I_q\). Assume their paired labels have the same height and occupation,
both represented spaces are closed under each core deletion, and the
vacuum is represented by its ordinary coordinate vector. These are the
hypotheses checked by `validate_maps` in the coherent compound engine.

For a larger dyadic \(m=kN\), attach every occupation of modes
\(k,\ldots,m-1\) to each core label. The physical alphabet has size
\(q2^{m-k}\). The two complete cut matrices use the row/column spaces

\[
\bigl(L\otimes L,P\otimes P\bigr),\qquad
\bigl(L\otimes P^T,L^T\otimes P\bigr),
\]

with the free tail on every factor. Column permutations used to restore
the endpoint identity do not change any Smith order.

## 2. A regular change of the tail coordinates

Assume \(N\ge2\), and set

\[
z=t^{2k},\qquad
G=1+z+\cdots+z^{N-2},\qquad
\alpha_N=zG.
\]

Index a tail mode by \(c+kh\), where \(0\le c<k\) and
\(1\le h<N\). On the \(N-1\) tail indices with fixed residue \(c\),
put

\[
e=(1,\ldots,1)^T,\qquad v=(1,z,\ldots,z^{N-2})^T.
\]

Let the first row of \(R\) be \(e^T\). For \(i=1,\ldots,N-2\),
let row \(i\) be \(e_i^T+z^i e_0^T\). Then

\[
Rv=(G,0,\ldots,0)^T,\qquad \det R=G.
\]

Since \(G(0)=1\), this is a regular invertible transformation at
\(t=0\). Apply it independently to every residue and both physical
legs. It acts only on the free tail. Its complete exterior transformation
therefore preserves every core-code tensor free-tail space, for both cut
orientations.

Here are the identities responsible for the reduction. The core-to-tail
coupling at tail index \(h\) is \(z^h\,\mathrm{diag}(D_k,D_k)\).
For a positive dyadic index \(h<N\),

\[
\sum_{\ell=0}^{N-1}z^{h\mathbin{|}\ell}=0.
\]

Indeed, a set bit of \(h\) pairs equal summands. Consequently the tail
block \(E\) satisfies

\[
E(e\otimes x)=(zv)\otimes B_kx.
\]

The regular congruence by the tail transformation thus makes the first
tail coordinate interact with the core and makes all remaining tail
coordinates orthogonal to both. The active one-particle block is exactly

\[
\boxed{\quad
M_N=\begin{pmatrix}
B_k&\alpha_N\mathrm{diag}(D_k,D_k)\\
\alpha_N\mathrm{diag}(D_k,D_k)&\alpha_NB_k
\end{pmatrix}.
\quad}                                                    \tag{1}
\]

This matrix has size \(4k\), independently of \(m\). Its core carries
the original maps; its aggregate tail is a full \(2k\)-mode exterior
algebra across the two legs. All other modes are spectators.

## 3. The complete spectator profile

Let \(T=\mathrm{diag}(D_k,D_k)\). The dyadic kernel identities give

\[
TB_k^{-1}T=(z^{-1}+1)B_k.                                  \tag{2}
\]

Both \(\alpha_N B_k^{-1}T\) and
\(\alpha_N TB_k^{-1}\) are regular. For example, after taking out
\(z^{-1}\), their nonconstant entries have exponents
\(2(i\mathbin{\&}\neg j)\) or
\(2(j\mathbin{\&}\neg i)\), which are nonnegative. Regular block
elimination of (1) leaves

\[
B_k\quad\text{and}\quad
\alpha_N\bigl[1+(1+z)G\bigr]B_k
=\alpha_Nz^{N-1}B_k.
\]

The latter scalar has valuation \(2m\). Hence the active full-Fock
one-particle orders are

\[
\{1,3,\ldots,4k-1\}\ \cup\
\{2m+1,2m+3,\ldots,2m+4k-1\}.
\]

Subtract these orders from those of \(B_m\). The spectator multiset is

\[
\boxed{\quad
\Omega_{k,m}=
\{4k+1,4k+3,\ldots,2m-1\}
\cup
\{2m+4k+1,2m+4k+3,\ldots,4m-1\}.
\quad}                                                    \tag{3}
\]

Both intervals are empty when \(m=2k\). For each of the two compressed
cuts, let \(Z^{\mathrm{act}}_{c,N}(r)\) count the **complete** Smith
orders of the compressed compound matrix of (1). The active matrix has
\(q^2 4^k\) channels. Regularity of the tail-only transformations proves
the exact identity

\[
\boxed{\quad
Z^{\mathrm{full}}_{c,m}(r)=
Z^{\mathrm{act}}_{c,N}(r)
\prod_{\omega\in\Omega_{k,m}}(1+r^\omega),
\qquad c\in\{0,\infty\}.
\quad}                                                    \tag{4}
\]

This is a complete profile identity, not merely a determinant identity.
The physical height partition remains that of the original core and
tail; (4) is not a replacement of its physical cost by an active-alphabet
cost.

## 4. An effective infinite-tail limit

In the series ring,

\[
\alpha_\infty=\frac{z}{1+z},\qquad
\alpha_N-\alpha_\infty=\frac{z^N}{1+z}.
\]

Thus every entry of the finite active compound matrix agrees with its
limiting entry modulo \(t^{2m}\). Put
\(S=TB_k^{-1}\). Since
\(\alpha_\infty TB_k^{-1}T=B_k\), the limiting matrix factors as

\[
M_\infty=
\begin{pmatrix}I&S\\0&I\end{pmatrix}
\begin{pmatrix}0&0\\0&\alpha_\infty B_k\end{pmatrix}
\begin{pmatrix}I&0\\S^T&I\end{pmatrix}.                    \tag{5}
\]

The induced shears replace core particles by free-tail particles. Core
deletion closure makes both compressed transformations invertible over
the Laurent series field. Therefore each limiting compressed cut has
exact rank

\[
R_k=4^k,
\]

although its size is \(D_{\mathrm{act}}=q^2 4^k\). Equivalently, the sole surviving
core sector in the middle matrix of (5) is the core vacuum.

All entries of \(S\) have valuation at least \(-2k\). Each exterior
shear and its inverse therefore has entries of valuation at least
\(-4k^2\). The largest full exterior Smith order of
\(\alpha_\infty B_k\) is \(8k^2\). Applying the two bounded-pole
shears gives the conservative uniform bound

\[
0\le a\le16k^2                                             \tag{6}
\]

for every nonzero limiting Smith order. One way to justify this bound is
to multiply each shear and its inverse by \(t^{4k^2}\); the resulting
regular matrices give the usual componentwise Smith-order inequalities
under row and column multiplication. The original limiting compressed
matrix is regular, which gives the lower bound.

Agreement modulo \(t^{2m}\) fixes every Smith order strictly below
\(2m\). Consequently, once \(m>8k^2\), exactly \(4^k\) finite active
orders equal their limiting values and all other active orders are at
least \(2m\). Write the finite limiting counting polynomials as
\(Z^{\mathrm{lim}}_0,Z^{\mathrm{lim}}_\infty\). They have degree at
most \(16k^2\), coefficient sum \(4^k\), and constant term at least
one. For every real \(u>0\),

\[
0\le \log Z^{\mathrm{act}}_{c,N}(e^{-u})
       -\log Z^{\mathrm{lim}}_c(e^{-u})
\le(D_{\mathrm{act}}-R_k)e^{-2mu}.                          \tag{7}
\]

This gives a finite algebraic target: compute two fixed polynomials,
then combine them with the explicit infinite Gaussian product. A strict
negative limiting pressure at a fixed positive scale transfers to a
sufficiently large finite tail with an explicit error bound.

For completeness, let

\[
F(u)=\sum_{j\ge0}\log(1+e^{-(2j+1)u}),\qquad
A(v)=\sum_{j=1}^q e^{-v h_j}.
\]

The limiting complete two-cut pressure is

\[
\begin{split}
\mathcal G_\infty(u,\rho)={}&
\log Z^{\mathrm{lim}}_0(e^{-u})+
\log Z^{\mathrm{lim}}_\infty(e^{-u})\\
&+2F(u)-2\sum_{j=0}^{2k-1}\log(1+e^{-(2j+1)u})\\
&-4\rho\left[
\log A(u/\rho)+F(u/\rho)
-\sum_{a=0}^{k-1}\log(1+e^{-(2a+1)u/\rho})
\right].                                                  \tag{8}
\end{split}
\]

For \(m>8k^2\), the absolute difference between the finite and limiting
pressures is at most

\[
\begin{split}
&2(D_{\mathrm{act}}-R_k)e^{-2mu}\\
&\quad+4k e^{-(2m+1)u}
 +\frac{2e^{-(4m+1)u}}{1-e^{-2u}}\\
&\quad+\frac{4\rho e^{-(2m+1)u/\rho}}{1-e^{-2u/\rho}}.
                                                               \tag{9}
\end{split}
\]

The three lines bound, respectively, the escaping active channels, the
missing spectator modes, and the missing physical tail modes. Thus a
strict inequality in (8) is an effective finite-tail witness, subject to
the existing whole-machine crossing transfer conditions.

## 5. What survives at the determinant level

There is a simpler identity that holds before any limiting operation.
Let \(r=m-k\), and let \(\sigma_c\) be the sum of all \(q^2\) core
Smith orders for cut \(c\). If label \(j\) has occupation \(n_j\),
then the sum of the extended cut orders is

\[
\boxed{\quad
\sigma^{\mathrm{full}}_c
=4^r\left[\sigma_c+4rq\sum_j n_j+2mrq^2\right].
\quad}                                                    \tag{10}
\]

To prove this, eliminate the tail over the Laurent field. The dyadic
inverse identity gives the core Schur complement
\(t^{2r}B_k\). The two induced shears strictly lower core occupation
and have determinant one on the represented spaces. Their entries may
have poles, so this step proves a determinant identity without claiming
Smith equivalence. The tail one-particle determinant has valuation
\(4mr\); its full exterior determinant has valuation
\(2mr4^r\). Taking the tensor-product determinant gives (10).

Two cores with the same occupation multiset therefore preserve their
normalized difference in summed Smith orders under every dyadic free-tail
extension. For the published coherent and coordinate \(q=32,k=8\)
maps, the sum of the differences over both cuts is \(540\). Its
normalized value remains

\[
\frac{540}{32^2}=\frac{135}{256}.
\]

This statement concerns dyadic total mode counts. Appending only mode
eight to that eight-mode core gives nine total modes and is outside the
hypothesis.

## 6. Why the first moment does not settle the threshold

The unrestricted infinite Gaussian profile has the midpoint expansion

\[
F(u)=\frac{\pi^2}{24u}-\frac{u}{24}+O(u^3).
\]

At \(\rho_*=1/\sqrt2\), its critical pressure is
\(2F(u)-4\rho_*F(u/\rho_*)=u/12+O(u^3)\).
It is tempting to compare \(135/256\) directly with \(1/12\).
That comparison interchanges limits which are not uniform.

Indeed, the limiting active polynomials each have coefficient sum
\(4^k\), independent of \(q\), whereas the physical core has
\(q\) labels. Substituting these values into (8) gives

\[
\boxed{\quad
\lim_{u\downarrow0}\mathcal G_\infty(u,\rho_*)
=4\rho_*\log\frac{2^k}{q}.
\quad}                                                    \tag{11}
\]

This is strictly positive for every proper fixed core \(q<2^k\).
For the \(q=32,k=8\) example, it is
\(4\rho_*\log8\). The missing active channels escape to orders of
size at least \(2m\), as (7) quantifies; the finite determinant moment
records those channels, while a fixed positive Laplace scale suppresses
them.

Equation (11) rules out the proposed small-scale argument in this order
of limits. It does not rule out a negative value of (8) at a finite
positive scale, nor a family whose core itself grows. Those are concrete
remaining construction problems.

## 7. Retained checks and scope

`tools/check_coherent_free_tail.py` checks the regular transformation and
active/spectator separation by exact binary polynomial arithmetic. It
also replays small complete compressed profiles for the determinant
identity and checks the published \(q32\) moment. Its finite checks
support the implementation; the identities above supply the arbitrary
tail-length proof.

The tail-only reduction and the Laurent determinant argument were also
independently reviewed by the parallel cube-construction researcher.
All complete channels are retained. No score gain or all-scale pressure
improvement is inferred from a determinant moment.
