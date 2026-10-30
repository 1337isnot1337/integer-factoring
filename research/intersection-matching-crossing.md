# Matching crossings with monomial duals

## Result and scope

There is an explicit characteristic-two crossing family whose labels
are arbitrary perfect matchings on \(2b\) physical vertices. Its vacuum
pairing is exactly the identity, including at odd particle count.
Each label therefore has a single monomial dual.

The family has at most

\[
(4b+1)5^{2b}
\tag{1}
\]

coefficient pieces and \((2b-1)!!\) available labels. A prefix of
\(2^\lambda\) labels uses \(b=O(\lambda/\log(\lambda+2))\).
All coefficients and exact rank factorizations can be generated in
\(2^{O(\lambda)}\) time. The piece entropy is
\(O(\lambda/\log(\lambda+2))\), within the relaxed evaluator interface.

This is an effective construction with a simpler dual pairing. Its
unmodified finite coefficient profiles do not establish an improvement
over the half-rate bound. A scalar deformation and a low-rank correction
formula below give precise further rank-analysis targets.

The doubled construction uses the pinned
[Anthropic Appendix B operators](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth)
and the matching identity proved in
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md).
The smaller construction can also be defined directly by its Pfaffian
and checked at the crossing evaluation point. No new Lean work is used.

## 1. Intersection gives an identity pairing

Let \(\mathcal M_b\) be the perfect matchings of
\(\{0,\ldots,2b-1\}\). For edges \(e,f\), put

\[
K(e,f)=[e\cap f\ne\varnothing].
\]

For \(A,B\in\mathcal M_b\), the matrix
\((K(e,f))_{e\in A,f\in B}\) splits according to the alternating
components of \(A\cup B\). A common edge gives a \(1\times1\) identity
block. Every nontrivial alternating cycle gives \(I+P_r\), where
\(P_r\) is a cyclic permutation matrix of size \(r\ge2\). Such a
block has the all-ones vector in its kernel. Hence

\[
\det K[A,B]=[A=B]\quad\text{over }\mathbb F_2.
\tag{2}
\]

This applies to every subset or prefix of the matching labels.
No Gram inversion or principal-minor selection is necessary.

## 2. The direct polynomial crossing

Write \(Z^e=\prod_{j\in e}Z_j\), and define

\[
z_e=LZ^e,\qquad
t_e=(1+L)\prod_{j\in e}(1+Z_j).
\]

For labels \(A,B,U,Y\), form an alternating matrix with vertex groups
\(U,B,A,Y\). Within-group entries vanish. Its cross-group entries are

| Pair | Entry |
| --- | --- |
| \(U_u,B_e\) | \([u=e]t_e\) |
| \(U_u,A_e\) | \([u=e]z_u t_e\) |
| \(U_u,Y_y\) | \([u=y]z_y\) |
| \(B_e,A_f\) | \(K(e,f)z_f\) |
| \(B_e,Y_y\) | \(K(e,y)z_y t_e\) |
| \(A_e,Y_y\) | \(K(e,y)t_e\) |

Let \(\mathcal P_K(A,B,U,Y;L,Z)\) be its Pfaffian.
At \(L=Z_j=1\), every \(t_e\) vanishes and every \(z_e\) is one.
Only \(UY\) and \(BA\) pairings remain. By (2),

\[
\mathcal P_K(A,B,U,Y;1,\mathbf1)
 =[U=Y][A=B].
\tag{3}
\]

Thus the coefficient tensors of \(\mathcal P_K\) sum to the crossing
with multiplier one. This is an identity over every commutative
\(\mathbb F_2\)-algebra, by base change.

Each physical vertex occurs once in each external label. The table
uses each occurrence at most once in a \(Z_j\) or \(1+Z_j\) factor.
Therefore \(\deg_{Z_j}\mathcal P_K\le4\).
There are \(2b\) paired edges in a Pfaffian term, each of \(L\)-degree
at most two, so \(\deg_L\mathcal P_K\le4b\). This proves (1).

The opposite edge-type counts satisfy the same relations as in
[permutation-homogeneous-labels.md](permutation-homogeneous-labels.md),
so every nonzero \(L\) degree is even. No rank estimate follows merely
from this parity.

## 3. Derivation by duplicating and folding source labels

Replace each matching label by two disjoint physical copies of itself.
Each new monomial has \(2b\) distinct odd parts and common weight

\[
H=2b+2(2^{4b}-1).
\]

Any nontrivial alternating component between two original labels now
occurs twice. The source disjointness Gram formula therefore vanishes
between distinct duplicated labels; on equal labels its even particle
count gives one. These are homogeneous monomial unit duals in the
original source construction.

Use the source multigraded six-block Pfaffian, then identify the
variables of corresponding physical vertices in the two copies.
After reordering particles by copy, the matrix has the form

\[
\begin{pmatrix}M&N\\N&M\end{pmatrix}.
\tag{4}
\]

Both \(M\) and \(N\) are alternating. In particular, a particle and
its copy belong to the same operator or state group, so the diagonal
of the cross-copy block \(N\) is zero.

For alternating \(M,N\) over a characteristic-two polynomial ring,

\[
\operatorname{Pf}
\begin{pmatrix}M&N\\N&M\end{pmatrix}
=\operatorname{Pf}(M+N)^2.
\tag{5}
\]

One proof pairs perfect matchings under exchange of the two copies.
Nonfixed matchings cancel. A fixed matching consists of two-edge
orbits between different base vertices; an edge cannot be fixed by
the exchange because the corresponding entry of \(N\) is zero.
For each paired base-vertex pair, the parallel and crossed lifts
contribute \(M_{ij}^2+N_{ij}^2=(M_{ij}+N_{ij})^2\).
Summing the quotient matchings gives (5).

In \(M\), annihilation uses disjointness. Between distinct physical
copies, annihilation is always allowed; creation and state equality
are never allowed. Therefore \(M+N\) has creation and equality
kernels unchanged, and annihilation kernel
\(1+[e\cap f=\varnothing]=K(e,f)\).
It is exactly the matrix of Section 2.

Squaring a polynomial over \(\mathbb F_2\) doubles every exponent and
leaves its coefficient tensors unchanged. The direct crossing family
is therefore a reindexing of this folded source family.

The alternating hypothesis on \(N\) cannot be dropped. For example,
\(M=0\) and \(N=I_2\) give a nonzero four-dimensional Pfaffian,
whereas discarding the diagonal of \(M+N\) gives zero.

## 4. A scalar family of exact identity pairings

Let \(B\) be the vertex-edge incidence matrix of the complete graph
on the \(2b\) physical vertices. In characteristic two,

\[
K=I+B^{\mathsf T}B.
\tag{6}
\]

The incidence product is zero on a common edge, one on two distinct
edges meeting in one vertex, and zero on disjoint edges.
More generally, for any scalar \(s\) in a commutative
\(\mathbb F_2\)-algebra, set

\[
K_s=I+sB^{\mathsf T}B.
\tag{7}
\]

Between matching labels, a common edge still gives a one. A
nontrivial alternating-cycle block is \(s(I+P_r)\), whose determinant
vanishes identically as a polynomial in \(s\). Consequently

\[
\det K_s[A,B]=[A=B]
\tag{8}
\]

for every \(s\). Replacing \(K\) by \(K_s\) in Section 2 produces a
scale-one crossing for every such parameter. Fixed finite extension
fields are allowed by the algebraic evaluator interface.

The correction has rank \(2b-2\) for \(b\ge1\). Indeed, incidence
has image the even-weight vertex vectors, of dimension \(2b-1\);
its transpose has kernel spanned by the all-ones vector, which lies
in that image. Also \(BB^{\mathsf T}=J\) for the complete graph on
an even number of vertices, and \(JB=0\). Hence

\[
(B^{\mathsf T}B)^2=0,\qquad K_s^2=I.
\tag{9}
\]

Equation (9) concerns the complete edge universe. A restriction to
only bipartite edges need not retain the same square identity.
The matching pairing (8) remains valid on bipartite labels.

## 5. Effectivity and the precise low-rank expansion

The number of perfect matching labels is
\[
d_b=(2b-1)!!=2^{b\log_2b+O(b)},\qquad d_{b+1}/d_b=2b+1.
\]
Taking the least \(b\) with \(d_b\ge2^\lambda\) gives
\(b=O(\lambda/\log(\lambda+2))\) and \(d_b=2^{O(\lambda)}\).
Canonical recursive generation, subset Pfaffians with \(2^{O(b)}\)
states, and the coefficient box in (1) give \(2^{O(\lambda)}\)
generation time. Complete finite-field rank factorizations have the
same time form.

There is also a concrete reduction to the identity annihilation
kernel on partial-matching domains. The three corrected pair blocks
are \(BA,BY,AY\). For each pair block and each physical vertex,
introduce one auxiliary channel. There are \(6b\) channels.
Each channel is a rank-two alternating update \(uv^{\mathsf T}
+vu^{\mathsf T}\), with \(u,v\) supported on their respective
external groups. Its entries are incidence bits times the table's
one-particle polynomial weights.

The Pfaffian minor expansion selects a subset of these channels.
For a selected subset, the coefficient determinant factors by
external group, since each auxiliary column is supported on one
group. Removing the selected particles leaves the identity-kernel
Pfaffian on four partial matchings. Thus the whole polynomial is
a sum of at most \(2^{6b}\) local linear images of that
partial-matching tensor, with polynomial-valued local maps.

After coefficient extraction, keep the channel subset and the local
polynomial charge splits as separate piece indices. The additional
number of charge splits is \(\exp(O(b))\): every physical variable
has total target degree at most four, and the \(L\) degree is \(O(b)\).
For each resulting piece, both flattenings are local linear images
of one identity-kernel coefficient tensor. Their ranks therefore
cannot exceed the corresponding two ranks of that tensor.

This formulation preserves a joint rank bound piece by piece.
Bounding the row rank and column rank of an aggregated sum by two
unrelated maxima would not preserve such a bound. The reduction
also does not assert that the identity-kernel partial-matching
coefficients themselves have sufficiently small ranks.

## 6. Checks and current finite limits

The checker verifies the identity pairing on small full matching
domains, the complete-edge incidence rank and square identities,
the generic alternating copy-fold identity, and independent
Pfaffian evaluations. A nonalternating cross block is retained as
a negative control.

The exploratory bipartite permutation label domains of sizes
\(2,6,24\) have uniform central coefficient ranks \(2,16,364\),
respectively, in each flattening. These are the coefficients
\([L^{2b}\prod_jZ_j^2]\) at \(s=1\). The one-variable specialization
\(L=z,\ Z_j=z^{2^{j+1}}\) has the following complete finite profiles:

| \(b\) | Labels | Nonzero pieces | Maximum rank product |
| ---: | ---: | ---: | ---: |
| 2 | 2 | 39 | \(2\cdot4=8\) |
| 3 | 6 | 221 | \(24\cdot33=792\) |
| 4 | 24 | 957 | \(360\cdot576=207360\) |

These finite profiles motivate testing the parameter \(s\), different
coefficient recombinations, and the identity-kernel reduction. They
do not support an asymptotic improvement for the unmodified family.

Run the constructor checks with:

    python3 tools/check_intersection_matching_lift.py \
      --output evidence/intersection-matching-crossing.json

The complete profile replay is a separate, larger calculation:

    python3 tools/check_intersection_rank_profiles.py \
      --maximum-particles 4 --output evidence/intersection-rank-profiles.json

The profile checker uses NumPy only to accelerate exact binary
substitutions. Its rank calculation uses integer bit operations.
It verifies the complete crossing sums and independently checks
selected polynomial substitutions. The standard verification target
runs the constructor checks; it does not repeat the larger profile sweep.
