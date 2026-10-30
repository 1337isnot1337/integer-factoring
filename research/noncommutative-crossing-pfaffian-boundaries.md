# Pfaffian boundary sections

## Result and scope

The equal-length Pfaffian crossing admits an explicit family of linear
retracts obtained by changing the standard single-occupation labels by
symplectic boundaries. These changes preserve the exact crossing endpoint.

The complete tensor of every such section has an ordinary, all-size
formula. In particular:

* Boundary sections on at most three complementary pairs leave the
  crossing family unchanged.
* On four pairs, every boundary section gives exactly one additional
  tensor direction. Its two cut ranks are \(16\) and \(1\).
* This first nontrivial family cannot satisfy the finite seed criterion:
  every exact decomposition from its full coefficient span has
  \(\sum_i\sqrt{r_i c_i}\ge120>16^{\sqrt2}\), over every field of
  characteristic two.
* At larger sizes the correction splits into explicit parity-syndrome
  channels. Every correction is divisible by
  \(p^2(p+q+r)^2\). The formula permits genuinely different families;
  the four-pair obstruction does not exclude all larger sections.

This note concerns the standard coordinate projection onto single
occupations and arbitrary lifts into the cycles of the symplectic
differential. It does not classify arbitrary projections, arbitrary
noncycle lifts, or all linear retracts of the exterior algebra.

The starting tensor is the universal six-edge Pfaffian of
[complement-shell-crossings.md](complement-shell-crossings.md), with
the coefficient representation in
[explore_complement_edge_pencils.py](../tools/explore_complement_edge_pencils.py).
The proofs below use exact polynomial identities, without numerical
rank assumptions or new Lean formalization.

## 1. The exterior complex and exact sections

Let \(K\) be any field of characteristic two. Let \(W\) have the
ordered basis

\[
e_1,\bar e_1,\ldots,e_P,\bar e_P,
\qquad
\Omega=\sum_{j=1}^P e_j\bar e_j.
\tag{1}
\]

Products are exterior products. They commute in characteristic two,
and every basis generator squares to zero. Since
\(\Omega^2=0\), multiplication by \(\Omega\) is a differential of
degree two.

On each complementary pair use the occupation labels \(0,1,2,3\):
empty, \(e_j\), \(\bar e_j\), and \(e_j\bar e_j\). The local
differential sends \(0\) to \(3\) and kills the other three states.
Its two-dimensional homology has representatives \(1,2\).

Let \(H\subset\Lambda^P W\) be the space spanned by the \(2^P\)
words with exactly one occupation on each pair. Let \(I_0\) be its
coordinate inclusion and \(P_0\) its coordinate projection. Then

\[
\ker(\Omega:\Lambda^P W\longrightarrow\Lambda^{P+2}W)
=
H\ \oplus\
\operatorname{im}(\Omega:\Lambda^{P-2}W\longrightarrow\Lambda^P W).
\tag{2}
\]

To see this without a homological theorem, split the local space as
the direct sum of the single occupations and the two-state complex
\(0\longmapsto3\). Every global summand containing an even pair is
contractible: on one chosen even pair, the map \(3\mapsto0\) is a
contracting homotopy. The only summand without an even pair is \(H\).

Write the boundary space on the right of (2) as \(B\). For any linear
map \(L:H\to B\),

\[
I=I_0+L,\qquad P_0I=\operatorname{id}_H.
\tag{3}
\]

Thus \(I\) and \(P_0\) are an exact linear retract.

### The universal crossing tensor

Use four exterior blocks \(U,A,V,Y\). The six edge weights, in the
order \(UA,AY,AV,UV,VY,UY\), are

\[
(p,q,r,1,1,1).
\tag{4}
\]

The \(UA,AV,VY,UY\) edges join equal local modes. The \(AY,UV\)
edges join complementary local modes. The tensor entry is the
Pfaffian matching polynomial on the occupied vertices.

When all four legs have degree \(P\), this is a homogeneous
polynomial \(G(p,q,r)\) of degree \(P\). It has the exact endpoint

\[
G(0,0,1)_{x,y,u,v}=[x=v][y=u].
\tag{5}
\]

Apply \(P_0\) on the \(x,u\) legs and \(I\) on the \(v,y\) legs.
Equation (3) gives the same endpoint on \(H\):

\[
G_L(0,0,1)_{x,y,u,v}=[x=v][y=u].
\tag{6}
\]

Every construction in this note is therefore an exact crossing,
before any rank estimate is made.

The coefficient span has at most \(\binom{P+2}{2}\) dimensions.
For example, substituting \(p=z,\ q=z^{P+1},\ r=1\) separates
the homogeneous coefficients and gives \(O(P^2)\) univariate
positions. Replacing \(z\) by \(1+z\) moves the crossing endpoint
to \(z=1\).

## 2. The complete local calculation

All output labels \(x,u\) are single occupations. There are only two
nonzero parity cases for the input labels \(v,y\).

If both inputs are single occupations, the local tensor is

\[
T=(r+q)E+(p+q)C,
\tag{7}
\]

where

\[
E_{x,y,u,v}=[x=v][y=u],
\qquad
C_{x,y,u,v}=[x=u][y=v].
\tag{8}
\]

If both inputs are even, index their states \(0,3\) by \(0,1\).
The local tensor is

\[
[x=u]\,M_{v,y},\qquad
M=
\begin{pmatrix}
p&q\\
r&p+q+r
\end{pmatrix}.
\tag{9}
\]

If one input is odd and the other is even, the local tensor vanishes.

These identities follow by listing the perfect matchings on at most
six occupied vertices. For example, the \((3,3)\) input entry has
three surviving matchings, with weights \(p,q,r\). The mixed
\((3,0)\) and \((0,3)\) entries have weights \(r\) and \(q\).

The checker compares the literal local Pfaffian contractions, not
only the asserted rank consequences.

Define the two scalar products

\[
A=p(p+q+r),\qquad B=qr.
\tag{10}
\]

Their sum is

\[
A+B=(p+q)(p+r)=\det M.
\tag{11}
\]

The notation \(A\) in (10) denotes a scalar polynomial; it is
distinct from the exterior block named \(A\).

## 3. Parity syndromes separate every boundary section

For a word, its parity syndrome is the set \(\sigma\subseteq[P]\)
of pairs carrying a single occupation. Equation (9) and the
vanishing mixed-parity case imply that the two input syndromes
must coincide.

Suppose \(|\sigma|=k\), and put \(m=P-k\). Since the input degree
is \(P\), \(m=2h\) is even. On the even pairs an input word must
have exactly \(h\) filled pairs and \(h\) empty pairs.

The input sector with syndrome \(\sigma\) is

\[
K^{\{1,2\}^{\sigma}}\otimes
K^{\{S\subseteq[P]\setminus\sigma:\ |S|=h\}}.
\tag{12}
\]

The boundary subspace in this sector is the first factor tensored
with the image of

\[
D e_S=\sum_{i\notin S}e_{S\cup\{i\}},
\qquad |S|=h-1.
\tag{13}
\]

Here \(D\) is multiplication by the sum of the even-pair variables.
Its image has dimension \(\binom{2h-1}{h-1}\). This follows from
the contracting homotopy for any chosen even pair.

Consequently

\[
\dim B=
\sum_{\substack{2\le m\le P\\m\ {\rm even}}}
\binom Pm\,2^{P-m}\binom{m-1}{m/2-1}.
\tag{14}
\]

In particular, \(\dim B=6\) for \(P=3\), and \(\dim B=27\)
for \(P=4\).

There are no terms linear in the section map \(L\). A single input
has the full syndrome \([P]\), whereas every boundary lies in
sectors with at least two even pairs. Thus

\[
G_L=T^{\otimes P}
+\text{the contraction with \(L\) on both input legs}.
\tag{15}
\]

The quadratic term is a direct sum over input syndromes. On
\(\sigma\), its odd-pair factor is \(T^{\otimes k}\); its output
factor on the remaining pairs is
\(\prod_{j\notin\sigma}[x_j=u_j]\).

It remains to compute one explicit bilinear form on the even
boundaries.

## 4. The all-size even-boundary Gram formula

For two \(h\)-subsets \(R,T\) of \(2h\) even pairs, (9) gives

\[
H_{R,T}
=A^{|R\cap T|}B^{h-|R\cap T|}.
\tag{16}
\]

Indeed, the number of simultaneously empty pairs equals the
number of simultaneously filled pairs, namely \(|R\cap T|\).
The two kinds of mismatched pairs both occur
\(h-|R\cap T|\) times.

Let \(S,T\) now have size \(h-1\). Write

\[
k=|S\cap T|,\qquad a=h-1-k.
\]

**Theorem 1.** The form induced on boundary generators is

\[
\boxed{
\begin{aligned}
H(D e_S,D e_T)
={}&(k\bmod2) A^{k+1}B^{h-k-1}\\
&+(a\bmod2) A^{k+2}B^{h-k-2}.
\end{aligned}}
\tag{17}
\]

A summand with zero scalar coefficient is omitted; in every
remaining summand its exponents are nonnegative.

**Proof.** Expand the two sums in (13). Their common outside set
has size \(k+2\).

If the two newly filled positions coincide in the common outside
set, the intersection size becomes \(k+1\), with multiplicity
\(k+2\equiv k\pmod2\).

If the first new position is in \(T\setminus S\) and the second
is in \(S\setminus T\), the intersection size becomes \(k+2\),
with multiplicity \(a^2\equiv a\pmod2\).

The two mixed cases, with one new position in a symmetric
difference and the other outside both sets, have equal
multiplicities and cancel.

The remaining case chooses two distinct positions outside both
sets. Its multiplicity is \((k+2)(k+1)\), which is even.
These are all possibilities. Substituting (16) proves (17).
\(\square\)

Equations (7), (13), (15), and (17) give every entry of every
boundary-section crossing, at every size. Redundant boundary
generators cause no ambiguity: (17) was obtained by applying the
literal form to their actual images under \(D\).

### A universal factor

Every nonzero term in (17) contains \(A^2\). The first term can
survive only when \(k\ge1\), and the second already has exponent
\(k+2\). Therefore

\[
\boxed{
G_L=T^{\otimes P}
+p^2(p+q+r)^2 R_L(p,q,r),
}
\tag{18}
\]

where \(R_L\) is homogeneous of degree \(P-4\). For \(P<4\),
the correction is zero.

This proves the exact boundary-section invariance on two and
three pairs. It also explains why a test confined to those
sizes would incorrectly suggest invariance at every size.

Equation (18) concerns a polynomial tensor, including every
coefficient. It is stronger than equality at a few field points.

## 5. Four pairs: every section collapses to one functional

Take \(P=4\).

A boundary with two odd pairs has only two even pairs, so \(h=1\).
Both terms in (17) vanish. These \(24\) boundary directions
contribute zero to the complete crossing.

The fully even boundary space has dimension three. It is generated
by the four stars

\[
b_i=D e_{\{i\}}
=\sum_{j\ne i}e_{\{i,j\}},
\qquad \sum_{i=1}^4 b_i=0.
\tag{19}
\]

For any \(i,j\), (17) gives

\[
H(b_i,b_j)=A^2.
\tag{20}
\]

Thus this form has a single effective functional
\(\ell(b_i)=1\). The functional is well defined, since its value
on the only relation in (19) is \(4=0\).

**Theorem 2.** For every boundary section on four pairs there is
a functional \(f:H\to K\) such that

\[
\boxed{
G_L=T^{\otimes4}
+p^2(p+q+r)^2 Z_f,
\qquad
Z_f(x,y,u,v)=[x=u]\,f(y)f(v).
}
\tag{21}
\]

Conversely, every \(f\) occurs: add \(f(v)b_1\) to the lift of
each single label \(v\).

For \(f\ne0\), the row cut \((x,y)\mid(u,v)\) of \(Z_f\) is

\[
I_{16}\otimes(f f^{\mathsf T}),
\]

up to harmless index order, and the column cut
\((x,u)\mid(v,y)\) is an outer product. Hence

\[
r(Z_f)=16,\qquad c(Z_f)=1.
\tag{22}
\]

The correction is nonzero. For instance, both inputs \(b_1\)
and equal single outputs give \(p^2(p+q+r)^2\).

Put \(a=r+q,\ b=p+q\). The binary part has the five coefficient
tensors

\[
H_j=[a^{4-j}b^j]T^{\otimes4},\qquad 0\le j\le4.
\tag{23}
\]

The additional polynomial is

\[
p^2(p+q+r)^2
=a^2b^2+b^4+a^2q^2+q^4.
\tag{24}
\]

Therefore the full coefficient span for \(f\ne0\) is exactly
\(\operatorname{span}(H_0,\ldots,H_4,Z_f)\). The independent
\(q^4\) coefficient isolates \(Z_f\), and then recovers every
\(H_j\). The crossing is \(H_0=E^{\otimes4}\).

## 6. A rank law over every field of characteristic two

A better basis for the five-dimensional binary span is

\[
B_i=[z^i]\bigl(C+z(E+C)\bigr)^{\otimes4},
\qquad 0\le i\le4.
\tag{25}
\]

The crossing is

\[
E^{\otimes4}=\sum_{i=0}^4 B_i.
\tag{26}
\]

Write

\[
v=(256,120,46,10,1).
\tag{27}
\]

**Lemma 3.** Let
\(F=\sum_{i=0}^4 a_iB_i\ne0\), and let \(\ell\) and \(h\)
be the least and greatest indices with \(a_i\ne0\). Over every
field \(K\) of characteristic two,

\[
r(F)=v_\ell,\qquad c(F)=v_{4-h}.
\tag{28}
\]

**Proof.** In the row cut, \(C\) is the identity on a
four-dimensional space and \(N=E+C\) is a rank-one square-zero
operator. Its module is a two-dimensional Jordan block plus
two trivial one-dimensional blocks.

In the fourfold tensor product, choose \(k\) factors carrying the
Jordan block. There are \(\binom4k2^{4-k}\) copies of the regular
module

\[
K[x_1,\ldots,x_k]/(x_1^2,\ldots,x_k^2).
\]

On this copy, \(F\) is multiplication by

\[
\sum_i a_i e_i(x_1,\ldots,x_k),
\tag{29}
\]

where \(e_i\) is the elementary symmetric polynomial.

If \(a_0\ne0\), the multiplier is a unit. If \(a_0=0\) and
\(a_1\ne0\), its nonzero linear term can be completed to an
invertible linear change of variables. Sending one variable to
the whole multiplier, and the others to a complement of its
linear term, is an automorphism of the truncated polynomial
algebra: every image squares to zero, and its linear part is
invertible. Its multiplication rank is consequently \(2^{k-1}\).

For least degree two, use \(e_3=e_1e_2\). On \(k=4\),
also \(e_4=x_1x_2e_2\). Thus every higher term factors from
\(e_2\) by multiplication with a unit. The ranks of \(e_2\)
for \(k=2,3,4\) are respectively \(1,2,6\). For \(k=4\),
these are the ranks \(1,4,1\) on source degrees zero, one,
and two. The middle matrix is \(I+J\) of order four, which
is invertible in characteristic two.

For least degree three, \(e_4=x_1e_3\) on \(k=4\).
The corresponding ranks for \(k=3,4\) are \(1,2\).
The top monomial \(e_4\) has rank one.

Summing the module multiplicities gives exactly (27).

In the column cut, \(C\) is a rank-one square-zero operator
\(N\), while \(E+C=I+N\). Thus the local pencil in (25)
is \(zI+(1+z)N\). A combination with largest index \(h\)
has least nonzero elementary-symmetric degree \(4-h\),
with coefficient \(a_h\). All other terms have larger
degree. The same calculation proves the second equality
in (28). \(\square\)

In particular, \(Z_f\) is never in the binary span when \(f\ne0\).
Its column rank one would force \(h=0\) in (28), but then its
row rank would be \(256\), contrary to (22).

## 7. An exact all-field obstruction for the first boundary family

**Theorem 4.** Let \(f\ne0\), and decompose the crossing exactly as

\[
E^{\otimes4}=\sum_j F_j,\qquad
F_j\in\operatorname{span}_K(B_0,\ldots,B_4,Z_f).
\tag{30}
\]

Then, over every field \(K\) of characteristic two,

\[
\boxed{\ \sum_j\sqrt{r(F_j)c(F_j)}\ \ge120.\ }
\tag{31}
\]

**Proof.** Since the six directions are independent, write each
piece uniquely as

\[
F_j=A_j+\alpha_j Z_f,\qquad A_j=\sum_i a_{j,i}B_i.
\]

A piece with \(A_j=0\) receives no charge. Otherwise let
\([\ell,h]\) be the smallest interval containing its nonzero
binary coordinates.

If \(\alpha_j=0\), Lemma 3 gives the exact rank product
\(v_\ell v_{4-h}\). If \(\alpha_j\ne0\), the rank inequality
\(\operatorname{rank}(X+Y)\ge
|\operatorname{rank}X-\operatorname{rank}Y|\), together with
independence and (22), gives

\[
r(F_j)c(F_j)\ge
\max(1,|v_\ell-16|)
\max(1,|v_{4-h}-1|).
\tag{32}
\]

Assign the five binary coordinates the integer weights

\[
w=(15,30,36,23,16).
\tag{33}
\]

For all fifteen intervals \(0\le\ell\le h\le4\), direct integer
arithmetic gives

\[
\left(\sum_{i=\ell}^h w_i\right)^2
\le
\min\!\left\{
v_\ell v_{4-h},
\max(1,|v_\ell-16|)\max(1,|v_{4-h}-1|)
\right\}.
\tag{34}
\]

Thus a piece's rank cost pays for all its nonzero binary
coordinates, and even for the entire containing interval.

By (26) and (30), for each \(i\) the sum of the scalars
\(a_{j,i}\) is one. At least one piece must therefore have
a nonzero coefficient on that coordinate. Charge weight \(w_i\)
to one such piece. Equation (34) pays every charge, so the
total cost is at least \(\sum_iw_i=120\). \(\square\)

The proof permits arbitrarily many pieces, arbitrary coefficient
recombinations, and arbitrary extensions of the coefficient field.
If \(f=0\), the same charging argument with the unperturbed
rank products still gives the lower bound \(120\).

Since \(16^{\sqrt2}<51\), this family cannot produce a finite
seed below the required rank-sum threshold. The proof is stronger
than checking binary or four-element-field coefficient combinations.

## 8. What changes at the next sizes

The all-size formula gives compact descriptions before any dense
matrix is built.

For \(P=5\), only sectors with four even pairs contribute. Each
choice of the remaining odd pair has a two-dimensional effective
boundary quotient. Thus the entire section is described, for
rank purposes, by five linear maps

\[
\gamma_j:H\longrightarrow K^2.
\]

The correction has only two polynomial directions:

\[
A^2\bigl((r+q)Z_E+(p+q)Z_C\bigr).
\tag{35}
\]

Explicitly,

\[
\begin{aligned}
Z_E(x,y,u,v)
&=\sum_{j=1}^5
 \prod_{i\ne j}[x_i=u_i]\,
 \gamma_j(v)_{x_j}\gamma_j(y)_{u_j},\\
Z_C(x,y,u,v)
&=[x=u]\sum_{j=1}^5
 \sum_{a=1}^2\gamma_j(v)_a\gamma_j(y)_a.
\end{aligned}
\tag{36}
\]

For example, \(c(Z_C)\le1\), because its column cut is one
outer product. The obvious bounds are

\[
r(Z_E)\le320,\quad c(Z_E)\le20,\qquad
r(Z_C)\le320,\quad c(Z_C)\le1.
\tag{37}
\]

They are upper bounds, not a favorable complete profile.

For a fully even sector with six pairs, (17) becomes

\[
A^2(A+B)\,J,
\tag{38}
\]

where \(J\) is indexed by two-subsets of six positions and has
entry one exactly when their intersection has size one.
Its rank is four in characteristic two. Indeed,
\(J=W^{\mathsf T}W\), where \(W\) is the vertex-edge incidence
matrix of the complete graph on six vertices. The image of
\(W\) is the five-dimensional even-sum subspace, whose radical
for the ordinary dot product is the all-ones line.

Thus higher sections have additional explicit active quotients.
The four-pair rank certificate does not automatically extend to
these coupled channels. Their complete coefficient ranks, rather
than only one favorable coefficient, determine whether they can
improve a seed.

## 9. Reproduction and evidence

Run

\[
\texttt{python tools/check\_noncommutative\_crossing\_pfaffian\_boundaries.py}
\]

from the repository root. The checker:

* contracts the literal local Pfaffians on every generated
  boundary pair for \(P=2,3,4\), including \(28{,}672\)
  complete contractions for \(P=4\);
* checks the all-size Gram identity (17) by direct expansion
  through eight even pairs;
* checks the fifteen exact integer inequalities (34);
* can export the five binary coefficient tensors and one chosen
  \(Z_f\) as complete dense matrices in both cuts.

For example:

\[
\begin{aligned}
\texttt{python tools/check\_noncommutative\_crossing\_pfaffian\_boundaries.py}\\
\texttt{ --functional 1 --export-prefix /tmp/pfaffian\_boundary\_P4\_f1}\\
\texttt{ --output evidence/noncommutative-crossing-pfaffian-boundaries.json}.
\end{aligned}
\]

In exported label \(i\), pair \(k\) has occupation
\(1+\operatorname{bit}_k(i)\). The row cut is
\((x,y)\mid(u,v)\), and the column cut is
\((x,u)\mid(v,y)\). The first basis vector is the exact crossing.

The evidence is a finite audit of the stated formulas. The
all-size claims are proved in Sections 1–4, and the all-field
obstruction is proved in Sections 6–7.
