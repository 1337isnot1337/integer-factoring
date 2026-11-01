# Optimal common retracts that cancel two Pfaffian vertices

## Result and scope

The equal-length Pfaffian family has three distinguished coefficient
tensors, at its pure \(p^L,q^L,r^L\) vertices. One can cancel the
first two while retaining the exact crossing at the third.

For a single common retract used on all four label roles, the largest
possible dimension has an exact answer. Let \(J\) be the mode-complement
permutation on the ambient monomial labels. Suppose it has \(m\) two-cycles
and \(f\) fixed labels. Then the optimal dimension is
\[
\boxed{
d_{\max}=
\begin{cases}
m,&f>0,\\
\max(0,m-1),&f=0.
\end{cases}}
\tag{1}
\]
This holds over every field of characteristic two, and explicit
constructions attain the bound.

The constructions admit a dimension-versus-sparsity tradeoff. For example,
the middle grade on eight modes has \(D=70,m=32,f=6\). It supports
\(32\) labels with at most twelve terms per input lift and eleven per
output row, while both pure vertex coefficients vanish.

This is an exact endpoint construction, not a successful factoring seed.
Its intermediate coefficient ranks still need to be bounded.

**Scope.** Formula (1) concerns one common pair of maps \(I,P\), with
\(PI=I_d\), used for both types of labels. Independent operator and state
spaces can do better and are explicitly outside this optimality theorem.
In particular, disjoint \(J\)-stable monomial domains can have dimension
\(D/2\) in each role.

The underlying universal tensor is defined in
[complement-shell-crossings.md](complement-shell-crossings.md).
The boundary-only retracts in
[noncommutative-crossing-pfaffian-isotropic.md](noncommutative-crossing-pfaffian-isotropic.md)
are a separate subclass; the present maps need not be cycles.

## 1. The three exact vertex tensors

Let \(V\) be the space of degree-\(L\) exterior monomials on \(P\)
complementary mode pairs. Its dimension is
\[
D=\binom{2P}{L}.
\]
The involution \(J\) swaps the two physical modes within each pair:
occupations \(0,3\) remain fixed and \(1,2\) are exchanged.
It is not the operation that takes the complement of an entire subset.

For the homogeneous universal tensor \(G(p,q,r)\), the three vertex
coefficients are
\[
\begin{aligned}
[p^L]G(x,y,u,v)&=[x=u][y=v],\\
[q^L]G(x,y,u,v)&=J_{x,y}J_{u,v},\\
[r^L]G(x,y,u,v)&=[x=v][y=u].
\end{aligned}
\tag{2}
\]
Equal length on all four legs forces opposite matching counts to agree.
Thus a pure vertex uses only its two opposite edge types, and (2) follows
from their coordinate or mode-complement matchings.

Let \(I:K^d\to V\) and \(P:V\to K^d\) satisfy \(PI=I_d\).
Apply \(P\) to the \(x,u\) legs and \(I\) to the \(y,v\) legs.
The resulting vertices are
\[
\begin{aligned}
[p^L]G_{P,I}
  &=(PP^{\mathsf T})_{x,u}(I^{\mathsf T}I)_{y,v},\\
[q^L]G_{P,I}
  &=(PJI)_{x,y}(PJI)_{u,v},\\
[r^L]G_{P,I}
  &=[x=v][y=u].
\end{aligned}
\tag{3}
\]

It is therefore enough to arrange
\[
PI=I_d,\qquad I^{\mathsf T}I=0,\qquad PJI=0.
\tag{4}
\]
The complete homogeneous coefficient span then has at most
\(\binom{L+2}{2}-2\) generators.

## 2. The upper bound for a common retract

Write \(U=\operatorname{im}I\). Conditions \(PI=I_d\) and \(PJI=0\)
imply
\[
U\cap JU=0.
\tag{5}
\]
Indeed, \(Ia=JIb\) gives \(a=0\) after applying \(P\), and then \(b=0\).

Let \(N=J+I_D\). Since \(J\) has \(m\) two-cycles,
\(\operatorname{rank}N=m\). Its restriction to \(U\) is injective:
if \(Nu=0\), then \(Ju=u\), so \(u\in U\cap JU\). Thus
\[
d\le m.
\tag{6}
\]

Suppose \(f=0\) and \(d=m\). Then \(D=2m\), and the isotropic
subspace \(U\) has half the ambient dimension. The ordinary dot product
is nondegenerate, so \(U=U^\perp\).

Let \(\mathbf1\) be the all-ones vector. For every \(u\in U\),
\[
0=\langle u,u\rangle
=\left(\sum_i u_i\right)^2.
\]
A field has no nonzero element with square zero, so
\(\langle u,\mathbf1\rangle=0\). Hence
\(\mathbf1\in U^\perp=U\). But \(J\mathbf1=\mathbf1\), contradicting
(5). Therefore \(d\le m-1\).

This also bounds retracts whose \(p^L\) vertex vanishes through
\(PP^{\mathsf T}=0\) instead of \(I^{\mathsf T}I=0\).
By (3), at least one of those two factors must vanish. In the former
case replace \(I,P\) by \(P^{\mathsf T},I^{\mathsf T}\); the same
argument applies because \(J=J^{\mathsf T}\).

Likewise, the pure \(q^L\) tensor vanishes only if \(PJI=0\).
Thus the bound concerns all common retracts that annihilate both
vertices, not just the constructions below.

## 3. Attaining \(m-1\) when there is no fixed label

Write the ambient space as \(K^m\oplus K^m\), with
\(J(x,y)=(y,x)\). For \(0\le i\le m-2\), define
\[
z_i=e_i+e_{m-1},\qquad
a_i=\sum_{j=i+1}^{m-2}e_j,
\qquad
I_i=(a_i,a_i+z_i).
\tag{7}
\]
These are the input columns. For \(0\le k\le m-2\), define the
output rows
\[
b_k=\sum_{j<k}e_j,\qquad
P_k=(b_k,b_k+e_k).
\tag{8}
\]

A direct calculation gives
\[
\begin{aligned}
P_k I_i
&=b_k\cdot z_i+e_k\cdot a_i+e_k\cdot z_i
=[k=i],\\
P_k JI_i
&=b_k\cdot z_i+e_k\cdot a_i=0.
\end{aligned}
\]
Also
\[
I_i^{\mathsf T}I_j
=a_i\cdot z_j+a_j\cdot z_i+z_i\cdot z_j=0.
\tag{9}
\]
For \(i\ne j\), exactly one of the first two terms is one and the
last term is one; for \(i=j\), all contributions cancel.
Thus (4) holds with \(d=m-1\).

For \(m=0,1\), the zero-dimensional construction matches (1).
Every formula uses only binary coefficients and therefore works
over any extension field of characteristic two.

## 4. Attaining \(m\) with fixed labels

Now write
\[
J(x,y,z)=(y,x,z),\qquad
x,y\in K^m,\ z\in K^f.
\]
Choose an \(f\times m\) matrix \(C\) whose columns have norm one,
and choose an \(m\times m\) matrix \(A\) satisfying
\[
A+A^{\mathsf T}=I_m+C^{\mathsf T}C.
\tag{10}
\]
The right side is alternating, so a strictly triangular choice always
exists. Put
\[
\boxed{
I=
\begin{pmatrix}I_m+A\\ A\\ C\end{pmatrix},
\qquad
P=\begin{pmatrix}I_m+A&A&0\end{pmatrix}.
}
\tag{11}
\]
Then
\[
PI=(I_m+A)^2+A^2=I_m,\qquad
PJI=(I_m+A)A+A(I_m+A)=0,
\]
and
\[
I^{\mathsf T}I
=I_m+A+A^{\mathsf T}+C^{\mathsf T}C=0.
\tag{12}
\]
This attains the upper bound \(d=m\).

One fixed coordinate suffices: take every column of \(C\) equal to
the same fixed unit vector, and let \(A\) be the strictly triangular
all-ones matrix.

### A sparse colored form

Using more fixed coordinates reduces support sizes. Partition the
\(m\) swap orbits into color classes, and assign each class a distinct
fixed monomial. Within a class, order the orbit representatives as
\(a_i,b_i=Ja_i\), and let \(e_c\) be its fixed monomial. Set
\[
\begin{aligned}
I_i&=e_{a_i}
 +\sum_{\substack{j<i\\j\ {\rm in\ class}\ c}}
  (e_{a_j}+e_{b_j})+e_c,\\
P_i&=e_{a_i}^{\mathsf T}
 +\sum_{\substack{j>i\\j\ {\rm in\ class}\ c}}
  (e_{a_j}^{\mathsf T}+e_{b_j}^{\mathsf T}).
\end{aligned}
\tag{13}
\]
This is the block-diagonal triangular case of (11).

If a class has \(s\) orbits, its input columns have at most \(2s\)
terms and its output rows at most \(2s-1\). Balanced use of all
available fixed labels gives \(s\le\lceil m/f\rceil\).

For \(m=32,f=6\), classes of sizes \(6,6,5,5,5,5\) give the
dimension-thirty-two example announced above. The complete map uses
204 input terms and 172 output terms.

The choices of orbit representatives, classes, within-class order,
and the matrix \(C\) can all change intermediate ranks while keeping
the vertex identities exact. They are finite optimization parameters,
not additional hypotheses.

## 5. A uniform sparsity tradeoff

Partition the \(m\) orbits into
\[
g=\lceil m/s\rceil
\]
nonempty blocks of size at most \(s\). Give up to \(f\) blocks
distinct fixed anchors and use (13) on them. On each unanchored
block use (7)–(8).

The direct sum has
\[
\boxed{
d=m-\max(0,\lceil m/s\rceil-f),
}
\tag{14}
\]
with at most \(2s\) terms per input column and \(2s-1\) per output
row. In particular, sparse maps can retain any prescribed constant
fraction below the maximum when the number of orbit pairs grows.

For the seventy-dimensional middle grade, \(s=3,4,6\) give
\(d=27,30,32\), respectively. This supplies intermediate choices
between very sparse pairing maps and the maximum-dimensional map.

## 6. Exterior dimensions and entropy

On degree-\(L\) monomials of \(2P\) modes, a \(J\)-fixed monomial
is a union of complete mode pairs. Thus
\[
f=
\begin{cases}
\binom P{L/2},&L\ {\rm even},\\
0,&L\ {\rm odd},
\end{cases}
\qquad
m=\frac{\binom{2P}{L}-f}{2}.
\tag{15}
\]

Some explicit cases are:

| Mode pairs \(P\) | Degree \(L\) | Ambient \(D\) | Fixed \(f\) | Optimal common labels |
|---:|---:|---:|---:|---:|
| 2 | 2 | 6 | 2 | 2 |
| 3 | 3 | 20 | 0 | 9 |
| 4 | 2 | 28 | 4 | 12 |
| 4 | 3 | 56 | 0 | 27 |
| 4 | 4 | 70 | 6 | 32 |
| 5 | 5 | 252 | 0 | 125 |
| 6 | 6 | 924 | 20 | 452 |

For the middle grade \(L=P\), the elementary bounds
\[
\frac{2^{2P}}{2P+1}\le\binom{2P}{P}\le2^{2P},
\qquad f\le2^P
\]
show that the construction retains \(2^{2P-O(\log P)}\) labels.
Its dimension is asymptotic to half the ambient middle grade.

Any subset of the constructed label coordinates still satisfies (4).
Thus a power-of-two label count can be obtained by restriction without
changing the two canceled vertices.

The construction works on any \(J\)-invariant monomial domain of a
fixed degree, not just the full grade. A fixed-energy shell can be used
when mode complement preserves that shell; its actual orbit and fixed
counts determine the dimension.

## 7. Independent operator and state domains are different

The bound (1) does not apply when different spaces serve the two roles.
For example, take disjoint \(J\)-stable monomial sets \(A,B\) of equal
size \(d\), use \(x,v\in A\), and \(u,y\in B\). Both first vertices
in (2) vanish, and the third remains crossing.

This variant may achieve \(d=D/2\). It has no common isotropic image
\(U\), so the obstruction in Section 2 is inapplicable.

It also has an exact cut symmetry:
\[
G_{p,q,r}(x,y,u,v)
=G_{q,p,r}(Jx,u,y,Jv).
\tag{16}
\]
Swapping the two \(B\) legs and applying mode complement on the two
\(A\) legs exchanges the \(p,q\) edge types and preserves \(r\).
Consequently, for \(J\)-stable \(A,B\),
\[
c([p^a q^b r^c]G)
=r([p^b q^a r^c]G).
\tag{17}
\]

If \(A,B\) are fibers of a linear parity observable on pair
occupations, the original parity-syndrome block decomposition is
preserved. This can make complete finite rank computations smaller.

These are exact structural identities. Neither the common-retract
dimension theorem nor the independent-domain construction establishes
a favorable rank-sum bound by itself.

## 8. Reproduction

The checker
[check_noncommutative_crossing_optimal_retracts.py](../tools/check_noncommutative_crossing_optimal_retracts.py)
audits \(PI=I\), \(PJI=0\), input isotropy, full rank, and
transversality on 287 abstract involution cases. It then constructs
the actual exterior mode-complement matrices and checks the cases
in the table above.

The evidence file is
[noncommutative-crossing-optimal-retracts.json](../evidence/noncommutative-crossing-optimal-retracts.json).
For small exterior cases it contains complete sparse input columns
and output rows, with the ambient monomial ordering.

The maximal dimension is established by the ordinary proof in
Section 2. The computations independently verify the attaining
constructions and their conventions; complete intermediate-rank
searches are separate.
