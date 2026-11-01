# Sparse symplectic radical codes with two vanishing endpoints

This note gives an explicit family of linear retracts for the
[ternary Pfaffian crossing](pfaffian-ternary-edge-family.md).
The input code is a radical of the symplectic boundary pairing. A
consistent choice of odd-mode orientation makes the \(p^L\) and \(q^L\)
endpoint tensors vanish, while the \(r^L\) endpoint remains the crossing.

The construction retains a constant fraction of the full-grade label
dimension at every fixed density below the middle grade. Its vectors
have at most \(P^2/4+1\) nonzero monomial terms on \(P\) complementary
mode pairs. These are exact construction and dimension statements; a
better complete rank-sum exponent is not claimed.

## 1. Boundary and contraction operators

Work over a field \(F\) of characteristic two. Let \(V\) have basis
\[
 e_1,f_1,\ldots,e_P,f_P,
\]
with complementary modes \(e_i,f_i\). Give exterior monomials the
orthonormal bilinear form. Let
\[
 \Omega:\Lambda^jV\to\Lambda^{j+2}V,\qquad
 v\longmapsto\left(\sum_{i=1}^P e_i f_i\right)v,
\]
and let \(\Lambda=\Omega^{\mathsf T}\) be its adjoint contraction.
Thus \(\Lambda\) removes every fully occupied pair and sums the results.
Both operators square to zero.

A syndrome specifies the pairs that are singly occupied, together with
the choice of their occupied mode. Fix such a syndrome with \(s\) singles.
The remaining \(t=P-s\) pairs are empty or fully occupied. In total
exterior grade \(L\), exactly
\[
 k=(L-s)/2
\]
of these pairs are full.

Within this syndrome, use commuting square-zero variables \(z_i=e_if_i\)
for the even pairs. The relevant space is the degree-\(k\) part of
\[
 A_t=F[z_1,\ldots,z_t]/(z_1^2,\ldots,z_t^2).
\]
On it,
\[
 \Omega=\sum_i z_i,\qquad
 \Lambda=\sum_i\partial_i,\qquad
 \Lambda\Omega+\Omega\Lambda=t\,I.
 \tag{1}
\]
The last identity follows one variable at a time from
\(\partial_i z_i+z_i\partial_i=I\); mixed terms cancel in characteristic
two.

Let
\[
 B_{t,k}=\operatorname{im}\bigl(\Omega:A_{t,k-1}\to A_{t,k}\bigr).
\]
For \(t\ge1\),
\[
 \dim B_{t,k}=\binom{t-1}{k-1}.
 \tag{2}
\]
Indeed, a linear change of square-zero variables can make \(\Omega\) the
first variable, after which its image consists of the monomials
containing that variable. Binomial coefficients outside their natural
range are interpreted as zero. For \(t=0\), there are no boundaries.

## 2. Exact rank and radical of the boundary Gram form

**Theorem 1.** For \(t\ge1\):

- if \(t\) is odd, the restricted Gram form on \(B_{t,k}\) is
  nondegenerate;
- if \(t\) is even, its rank and radical dimension are
  \[
   \operatorname{rank}\operatorname{Gram}(B_{t,k})
      =\binom{t-2}{k-2},\qquad
   \dim\operatorname{Rad}B_{t,k}
      =\binom{t-2}{k-1},
   \tag{3}
  \]
  and
  \[
   \operatorname{Rad}B_{t,k}
       =\operatorname{im}\bigl(\Omega\Lambda:A_{t,k}\to A_{t,k}\bigr).
   \tag{4}
  \]

**Proof.** If \(t\) is odd, (1) makes \(\Omega\Lambda\) a self-adjoint
projection onto \(\operatorname{im}\Omega\). The image of a self-adjoint
projection in a nondegenerate ambient bilinear space is nondegenerate.

For the even case, distinguish one variable \(z_t\) and write
\(\Omega=\Omega'+z_t\). Every boundary has a unique expression
\[
 \Omega v=(\Omega'v,v),\qquad v\in A_{t-1,k-1},
 \tag{5}
\]
in the splitting according to whether \(z_t\) occurs. The Gram matrix in
this basis is
\[
 I+\Lambda'\Omega'=\Omega'\Lambda',
 \tag{6}
\]
since \(t-1\) is odd. This is the self-adjoint projection onto
\(\operatorname{im}\Omega'\) in degree \(k-1\). Its rank is
\(\binom{t-2}{k-2}\), and the complementary dimension gives (3).

For (4), choose linear coordinates \(w_0,w_1,\ldots,w_{t-1}\) such that
\(\Omega\) is multiplication by \(w_0\) and \(\Lambda=\partial_{w_1}\).
They exist because the defining vector and covector pair to
\(t=0\) in the field. Then
\[
 \Omega\Lambda=w_0\partial_{w_1}
\]
has rank \(\binom{t-2}{k-1}\). Its image is contained in
\(\operatorname{im}\Omega\cap\ker\Lambda=\operatorname{Rad}B_{t,k}\).
The dimensions agree. \(\square\)

On the full grade, \(t=P-s\equiv P-L\pmod2\). Therefore every boundary
syndrome is nondegenerate when \(P-L\) is odd. The radical construction
uses \(P\equiv L\pmod2\).

## 3. A closed sparse radical basis and its coordinate projection

Assume \(t\ge2\) is even. Order the even pairs, choose
\[
 p=\min E,\qquad q=\max E,
\]
and let \(S\subseteq E\setminus\{p,q\}\) have size \(k-1\). Suppress the
fixed single-mode factor in the notation. Define
\[
 R_S=\Omega\prod_{i\in S}(z_i+z_p)
     =\Omega\Lambda(z_pz_S).
 \tag{7}
\]
These are a basis of \(\operatorname{Rad}B_{t,k}\).

To see this directly, take
\[
 w_0=\Omega,\quad w_1=z_p,\quad
 w_i=z_i+z_p\quad(i\in E\setminus\{p,q\}).
\]
These \(t\) forms are independent, \(\Lambda(w_1)=1\), and \(\Lambda\)
annihilates all other displayed forms. Equation (7) is precisely the
monomial basis \(w_0w_S\) for the image of \(w_0\partial_{w_1}\).

There is also a particularly simple left inverse. The coefficient of
\(z_qz_S\) in \(R_{S'}\) is \([S=S']\). Thus the output map can be the
coordinate projection
\[
 \Pi_S(v)=[z_qz_S]v.
 \tag{8}
\]
No matrix inversion is needed.

For an explicit support formula, put \(U=\{p\}\cup S\), so \(|U|=k\).
Then
\[
 R_S=(k\bmod2)z_U+
      \sum_{\substack{i\in U\\j\in E\setminus U}}
       z_{(U\setminus\{i\})\cup\{j\}}.
 \tag{9}
\]
Every term in the double sum is different. Hence
\[
 |\operatorname{supp}R_S|
   =k(t-k)+(k\bmod2)\le t^2/4+1.
 \tag{10}
\]
This is polynomial in the number of pairs, even though the number of
basis vectors is exponential.

With the ambient monomials ordered by increasing occupied-mode subsets,
the pivot \(z_qz_S\) is the last term of \(R_S\). No other radical basis
vector contains that pivot. Thus (7) is already a reduced coordinate
basis. The checker compares it exactly, vector for vector, with the
independent syndrome Gram-kernel and row-reduction constructor in
[explore_symplectic_boundary_codes.py](../tools/explore_symplectic_boundary_codes.py).

The input coefficients are explicit Boolean predicates: a target
\(k\)-subset is present if it equals \(U\) with \(k\) odd, or differs from
\(U\) by one removal and one insertion. A basis label and one coefficient
can therefore be described and evaluated with \(O(P)\) bits of working
storage. Enumerating an input vector uses \(O(P^2)\) output terms.

## 4. Select half the orientations and preserve the crossing

Let \(J\) swap the two modes in every complementary pair. It fixes empty
and fully occupied pairs, and exchanges the two choices at every
singly occupied pair.

Take all radical basis vectors (7) with a nonempty single-occupation
syndrome, but retain only those whose least singly occupied pair uses
its low mode \(e_i\). Omit the zero-syndrome sector. Let \(I\) be the
matrix of these input vectors, and use the coordinate projections (8)
as the rows of \(\Pi\).

Then
\[
 \Pi I=I_d,\qquad I^{\mathsf T}I=0,\qquad \Pi JI=0.
 \tag{11}
\]
The first identity is (8). The second holds because every retained input
belongs to a Gram radical, and different syndrome/orientation sectors
have disjoint support. For the third, \(J\) changes the retained low mode
at the least single site to the omitted high mode. None of the retained
coordinate pivots can occur.

For the full-grade ternary Pfaffian tensor, the three vertex coefficients
are
\[
 p^L:\ \delta_{x,u}\delta_{y,v},\qquad
 q^L:\ J_{x,y}J_{u,v},\qquad
 r^L:\ \delta_{x,v}\delta_{u,y}.
\]
Apply \(\Pi\) on the \(x,u\) legs and \(I\) on the \(y,v\) legs. The three
coefficients become, respectively,
\[
 (\Pi\Pi^{\mathsf T})_{x,u}(I^{\mathsf T}I)_{y,v}=0,\qquad
 (\Pi JI)_{x,y}(\Pi JI)_{u,v}=0,\qquad
 (\Pi I)_{x,v}(\Pi I)_{u,y}=E.
 \tag{12}
\]
This proves the exact two-vertex cancellation and crossing endpoint.

This construction uses a common input/output retract on both logical
label pairs. All input vectors remain homogeneous of exterior grade
\(L\). Their physical mode energies need not coincide.

## 5. Exact dimensions

Let \(R(P,L)\) be the radical dimension in the full grade. For
\(P\equiv L\pmod2\),
\[
 R(P,L)=
 \sum_{\substack{0\le s\le P-2\\s\equiv L\ (2)}}
 \binom Ps\,2^s
 \binom{P-s-2}{(L-s)/2-1}.
 \tag{13}
\]
For opposite parity, \(R(P,L)=0\). The zero-syndrome contribution is
\[
 R_0(P,L)=
 \begin{cases}
 \binom{P-2}{L/2-1},&P,L\text{ even},\\
 0,&\text{otherwise}.
 \end{cases}
\]
The code in Section 4 has dimension
\[
 \boxed{\quad d(P,L)=\frac{R(P,L)-R_0(P,L)}2.\quad}
 \tag{14}
\]

Some exact values are:

| \(P\) | \(L\) | Ambient \(\binom{2P}{L}\) | Radical \(R(P,L)\) | Selected code \(d(P,L)\) |
|---:|---:|---:|---:|---:|
| 3 | 3 | 20 | 6 | 3 |
| 4 | 4 | 70 | 26 | 12 |
| 5 | 3 | 120 | 10 | 5 |
| 5 | 5 | 252 | 100 | 50 |
| 6 | 4 | 495 | 64 | 30 |
| 6 | 6 | 924 | 366 | 180 |
| 7 | 3 | 364 | 14 | 7 |
| 7 | 5 | 2002 | 336 | 168 |
| 7 | 7 | 3432 | 1316 | 658 |

For odd \(P\) and \(L=3\), (14) simplifies to \(d=P\). Its input vectors
are especially simple:
\[
 I_i=e_i\sum_{j\ne i}e_jf_j,
\]
and \(\Pi_i\) selects the term with the largest \(j\ne i\).

## 6. The full leading entropy survives

Fix \(0<\alpha<1/2\), choose \(L/(2P)\to\alpha\), and enforce
\(P\equiv L\pmod2\). Then
\[
 \frac{R(P,L)}{\binom{2P}{L}}
 \longrightarrow
 c(\alpha):=
 \frac{\alpha^2(1-\alpha)^2}
      {(\alpha^2+(1-\alpha)^2)^2}>0,
 \tag{15}
\]
and
\[
 \boxed{\quad
 d(P,L)\sim\frac{c(\alpha)}2\binom{2P}{L}.
 \quad}
 \tag{16}
\]

Here is an elementary proof of (15). Choose a uniformly random
\(L\)-element subset of the \(2P\) modes. Let \(k\) count its fully
occupied pairs, and let \(t-k\) count its empty pairs. Inside a syndrome,
the radical-to-ambient dimension ratio is
\[
 \frac{\binom{t-2}{k-1}}{\binom tk}
 =\frac{k(t-k)}{t(t-1)}
\]
for \(t\ge2\), with zero in the exceptional \(t=0\) sector. Therefore
the left side of (15) is the expectation of this ratio.

Pair-occupancy indicators give
\[
 \frac{k}{P}\longrightarrow\alpha^2,\qquad
 \frac{t-k}{P}\longrightarrow(1-\alpha)^2
\]
in probability. For example, their means follow from
\((L)_2/(2P)_2\) and \((2P-L)_2/(2P)_2\); disjoint-pair covariances are
\(O(P^{-1})\), so the variances of the totals are \(O(P)\).
The displayed dimension ratio is bounded by \(1/2\). Bounded convergence
in probability therefore proves (15).

The omitted zero-syndrome dimension has exponential order
\(2^{P H_2(\alpha)+o(P)}\), whereas
\(\binom{2P}{L}=2^{2P H_2(\alpha)+o(P)}\), so its relative contribution
vanishes. This proves (16). In particular,
\[
 \log_2 d(P,L)=2P H_2\bigl(L/(2P)\bigr)+O(\log P).
 \tag{17}
\]
The same leading entropy is available as in the original full grade.
If \(L=2\alpha P+O(1)\), obtained by rounding with the required parity,
the right side is \(2P H_2(\alpha)+O(\log P)\).

For example, near \(\alpha=0.2255\), the selected code keeps about
\(c(\alpha)/2\approx0.036\) of the full-grade dimension. This costs a
constant number of label bits asymptotically, not a positive fraction
of the entropy.

## 7. Verification and the remaining rank question

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_symplectic\_radicals.py}
\]
with optional output
\(\texttt{--output evidence/noncommutative-crossing-symplectic-radicals.json}\).

The checker verifies the explicit supports, both annihilation identities
\(\Omega R_S=\Lambda R_S=0\), Gram isotropy, coordinate left inverses,
and the \(J\)-orientation selector through \(P=7\). Through \(P=6\), it
compares every basis vector and pivot with the independent Gram-kernel
constructor. It also records exact dimension and entropy controls up
to \(P=512\).

The particle-hole formulation in
[pfaffian-particle-hole.md](pfaffian-particle-hole.md) sends the physical
operators to the summed roots
\[
 \Omega_X=\sum_i E_{32}^{(i)},\quad
 \Lambda_X=\sum_i E_{23}^{(i)},\qquad
 \Omega_Y=\sum_i E_{41}^{(i)},\quad
 \Lambda_Y=\sum_i E_{14}^{(i)}.
\]
These constraints give a concrete target for refining the stable
\(GL_2\times GL_2\) coefficient bounds. Their compatibility with the
stable block decomposition requires an additional rank argument.
Equations (11)--(17) establish the code and its size independently of
that remaining question.

