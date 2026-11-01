# Coherent quotient codes and an exact zero-cut invariant

## Result and scope

A homogeneous quotient of a Boolean algebra gives a physical retraction
whose labels can mix occupation numbers. This extends the homogeneous
inverse-system construction with an explicit algebraic parameterization.
Both map spaces are closed under every bit deletion.

When the number of virtual modes is a power of two, the complete zero-cut
Smith profile is unchanged by every nonlinear quotient in this family.
This is an exact theorem for arbitrary finite support and characteristic
two. It is not a conclusion drawn from a finite search. The infinity cut
can change, so the theorem alone does not exclude a pressure improvement.

The complete four-virtual-mode binary family contains sixteen maps. Every
nontrivial map worsens the numerical balanced pressure in the retained
scan. A separate five-virtual-mode sample checks 24 of 1,024 maps, with the
coordinate map again best. These are finite search results. They do not
give a factoring exponent above one half or prove a general infinity-cut
inequality.

## 1. Dualizing a graded quotient

Let

\[
A_m=F[x_0,\ldots,x_{m-1}]/(x_0^2,\ldots,x_{m-1}^2)
\]

over a field of characteristic two. Give `x_a` height `2a+1`. Define the
analogous algebra `A_k` on virtual generators `z_0,...,z_(k-1)`, with
`k <= m`. Choose a height-preserving homomorphism

\[
\phi:A_m\longrightarrow A_k,\qquad
\phi(x_a)=g_a(z),\qquad g_a=z_a\quad(a<k).
\tag{1}
\]

For `a >= k`, the image `g_a` can be any homogeneous polynomial of height
`2a+1`. Every positive-height polynomial squares to zero in this algebra,
so these choices define a homomorphism. They have odd occupation parity.
Moreover, their monomials have degree at least three: a linear virtual
monomial has height at most `2k-1`.

Use the occupation bases to set `P=phi^T`. Let `L` select the physical
monomials using only the first `k` generators. Then

\[
LP=I_{2^k}.
\tag{2}
\]

The columns of `P` and the corresponding rows of `L` share one height.
For physical bit deletion `D_a`, duality with multiplication gives

\[
D_aP=P\,g_a(D_0,\ldots,D_{k-1}).
\tag{3}
\]

Thus the column space of `P` is deletion invariant. The row space of `L`
is the full coordinate downset on the first `k` modes, so it is also
deletion invariant. Equations (2) and (3) verify the exact physical-map
hypotheses of the [mixed-occupation construction](coherent-hot-mixed-boundary.md).

The [constructor](../tools/coherent_graded_quotients.py) multiplies the
generator images in the square-zero algebra, transposes the complete
monomial table, and checks the physical map identities. It retains every
generator image, physical vector, and height. Its finite inputs therefore
determine the complete physical contraction without an existence oracle.

## 2. Integral interpolation of the hot cross block

Assume now that `k` is a power of two. Let `B_0` be the hot one-particle
cross block on the first `k` modes of each of its two sides. Thus

\[
B_0=\begin{pmatrix}D&J\\J&D\end{pmatrix},\quad
D_{ab}=(1+t)t^{2(a\mathbin{|}b)+1},\quad
J_{ab}=[a=b]t^{2a+1}.
\tag{4}
\]

Write `H_aa=t^(a+1/2)` formally, and

\[
B_0=\mathrm{diag}(H,H)K\,\mathrm{diag}(H,H),\qquad
K=\begin{pmatrix}C&I\\I&C\end{pmatrix},\quad
C_{ab}=(1+t)t^{a\mathbin{\mathrm{xor}}b}.
\]

The half-powers cancel in all expressions below; equivalently this step
can be performed over `F(t^(1/2))`. If `k=2^s`, the matrix with entries
`t^(a xor b)` is the tensor product of the matrices

\[
\begin{pmatrix}1&t^{2^j}\\t^{2^j}&1\end{pmatrix},
\qquad 0\le j<s.
\]

Squaring in characteristic two gives

\[
C^2=(1+t^{2k})I,\qquad K^2=t^{2k}I,
\qquad K^{-1}=t^{-2k}K.
\tag{5}
\]

Take any additional physical mode `b=c+h`, where `0 <= c < k` and
`h` is a positive multiple of `k`. Its cross-block column, restricted to
the low rows, is `t^(2h)(D e_c,0)^T` for a column on the first side.
Using (5), its coordinates in the `B_0` columns are

\[
R_{(0,i),(0,b)}
 =[i=c]t^{2(h-k)}(1+t^{2k}),
\]

\[
R_{(1,i),(0,b)}
 =(1+t)t^{2(h-k)+2(c\mathbin{\&}\neg i)}.
\tag{6}
\]

Swap the two sides for a column on the second side. Low columns have
identity coordinates. Every exponent in (6) is nonnegative because
`h >= k`. Consequently the entire restricted one-particle cross block
factors as

\[
B_{\mathrm{low},\mathrm{all}}=B_0R,
\qquad R\in F[t]^{2k\times2m},
\qquad R_{\mathrm{low}}=I.
\tag{7}
\]

This integrality is the essential point. A factorization over `F(t)`
alone would not preserve Smith valuations.

## 3. The complete zero-cut profile is unchanged

Use the complete compound operator `E_B`, whose equal-occupation entries
are determinants of selected one-particle minors. Cauchy–Binet applied
to (7) gives

\[
(L\otimes L)E_B(P\otimes P)
 =E_{B_0}\,T,\qquad
T=E_R(P\otimes P).
\tag{8}
\]

Order the virtual occupation basis by total particle number. A term in
`P` that uses a high physical generator has physical degree at least two
less than its virtual degree, by section 1. The compound operator `E_R`
preserves the physical degree. The unique degree-preserving term uses only
low physical generators, where both the quotient and `R` are identities.
Therefore

\[
T=I+N,
\]

where `N` strictly lowers total occupation by at least two. Its entries
are polynomials by (7). The matrix `N` is nilpotent, so the finite series
`I+N+N^2+...` is a polynomial inverse for `T` in characteristic two.
Thus `T` is unimodular over `F[t]`.

Equation (8) proves that all `2^(2k)` zero-cut Smith orders equal those of
the bare full Fock alphabet on the first `k` modes. This is a complete
profile identity, including every occupation block and channel. The
within-cut Gaussian factors already descend by deletion invariance, so
the conclusion applies to the actual compressed crossing as well as the
compound representative.

The infinity cut uses the distinct compression
`(L tensor P^T) E_B (L^T tensor P)`. Both its row and column supports can
contain high modes. Equation (7) does not factor this entire matrix, and
no unchanged-profile claim is made for that cut.

## 4. Exact controls and the retained scans

For four virtual modes, the only nonlinear generator images are scalar
multiples of the four triples of heights 9, 11, 13, and 15. Over `F_2`,
this gives exactly sixteen choices. The
[complete scan](../evidence/coherent-graded-quotients-4-8.json) retains both
full profiles for each map. The coordinate map has numerical balanced
rate `rho = 0.7166536597281816`; every nonzero parameter choice has the
larger value `0.7168366165249694`.

The [independent checker](../tools/check_coherent_graded_quotients.py)
verifies the polynomial column identity (7), the exact differential action
(3), and a fresh replay of both complete profiles for all sixteen maps.
Its [control record](../evidence/coherent-graded-quotient-controls.json)
contains the precise counts and source hashes.

A second [scan](../evidence/coherent-graded-quotients-5-13-screen.json)
evaluates 24 explicitly recorded members of the 1,024-map binary family
on five virtual modes. It includes the coordinate map, the all-parameter
map, every single-parameter map, and a deterministic sample of further
maps. The coordinate rate `0.7140695939909403` is best among those tested.
This sample is not exhaustive, and the dyadic proof above does not cover
five virtual modes.

Finally, the [shared-gauge checker](../tools/check_coherent_shared_gauges.py)
includes the nonlinear eight-letter map with
`phi(x_4)=z_0 z_1 z_2`. It reconstructs the literal ambient Pfaffian tensor
and compares both complete profiles after a nontrivial rational gauge.
This supplies a direct tensor control for the mixed-degree quotient
construction, independently of the one-particle interpolation argument.
