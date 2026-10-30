# Sparse isotropic Pfaffian sections

## Result and scope

For every number \(P\ge6\) of complementary mode pairs, there is an
explicit sparse linear retract of the middle exterior grade that:

* has \(d=2^P\) labels and exact crossing endpoint;
* changes only the input lifts, retaining the standard single-label
  coordinate projection;
* makes the entire coefficient of \(p^P\) vanish.

For \(P\ge7\), every lifted label is the sum of just four exterior
monomials. The construction is uniform and does not require matrix
inversion. The first possible size, \(P=6\), has a separate sparse
construction with at most ten monomials per lifted label.

The six-pair family has an explicit coefficient space of dimension at
most nine. Its two added directions have exact ranks independent of the
chosen boundary isometry. Constant-coefficient cancellation is therefore
a proved construction, not an inferred numerical effect.

No improved factoring exponent is asserted. The remaining coefficient
ranks must be controlled together. In particular, the next full
derivative cannot be canceled completely, although partial cancellation
is possible and admits a concrete combinatorial design.

The starting family and the all-size boundary Gram formula are proved in
[noncommutative-crossing-pfaffian-boundaries.md](noncommutative-crossing-pfaffian-boundaries.md).
This note records the new positive sections, their exact compression, and
the limit on first-order matching.

## 1. The pairing to cancel

Use the notation of the boundary-section note: \(W\) has \(P\)
complementary pairs, \(\Omega=\sum_j e_j\bar e_j\), and \(H\)
is the \(2^P\)-dimensional space of single occupations. Write \(I_0\)
for its coordinate inclusion and \(P_0\) for its projection.

On single output labels, the local Pfaffian tensor at \((p,q,r)=(1,0,0)\)
is
\[
[x=u][y=v]
\tag{1}
\]
for arbitrary input occupations. Hence, for any lift map \(I\), the
complete \(p^P\) coefficient after projection is
\[
[x=u]\,(I^{\mathsf T}I)_{y,v}.
\tag{2}
\]

Suppose a boundary map \(B:H\to\operatorname{im}\Omega\) has
\[
B^{\mathsf T}B=I_d.
\tag{3}
\]
Boundary words have at least two even pairs, so their coordinate supports
are disjoint from the single words. Therefore
\[
I=I_0+B,\qquad P_0I=I_d,\qquad
I^{\mathsf T}I=I_d+B^{\mathsf T}B=0.
\tag{4}
\]
The endpoint remains crossing, and (2) vanishes.

All pairings here are ordinary bilinear dot products over a field of
characteristic two. No conjugation, positivity, or characteristic-zero
orthogonality is being assumed.

## 2. Exact rank of the boundary pairing at every size

Consider an input parity sector with \(m=2h\) even pairs and \(P-m\)
odd pairs. On the even pairs, let
\[
D e_S=\sum_{i\notin S}e_{S\cup\{i\}},
\qquad |S|=h-1.
\tag{5}
\]
Its image is the boundary space in that sector. The Gram matrix on these
generators at \((1,0,0)\) is \(D^{\mathsf T}D\).

Let \(V\) be the \(m\)-dimensional space whose basis indexes the even
pairs. Then \(D\) is exterior multiplication by the all-ones vector,
and \(D^{\mathsf T}\) is contraction by the all-ones covector.
Their pairing is \(m=0\) in characteristic two.

Choose a basis in which the vector is the first basis vector and the
covector is the second dual basis vector. This is possible because
both are nonzero and their pairing is zero. The operator
\(D^{\mathsf T}D\) becomes
\[
e_1\wedge\partial_{e_2}
\]
on exterior degree \(h-1\). Its nonzero independent images correspond
exactly to monomials containing \(e_2\) and not \(e_1\). Thus
\[
\boxed{\operatorname{rank}(D^{\mathsf T}D)
=\binom{m-2}{h-2}.}
\tag{6}
\]
The value is zero when \(h=1\).

Summing the parity sectors gives the exact rank of the full middle-grade
boundary pairing:
\[
\boxed{
R_P=
\sum_{\substack{4\le m\le P\\m\ {\rm even}}}
\binom Pm\,2^{P-m}\binom{m-2}{m/2-2}.
}
\tag{7}
\]
The factor \(2^{P-m}\) counts the odd occupations.

At \(P=4,5,6\), this gives respectively \(1,10,64\). Hence six pairs
are the first size at which a \(2^P\)-dimensional boundary isometry
can exist.

## 3. A four-term lift for every \(P\ge7\)

Choose sixteen distinct four-element subsets
\[
J_0,\ldots,J_{15}\subseteq[P].
\]
This is possible since \(\binom74=35\). For each \(J_a\), let \(j_a\)
be its least element.

Label \(H\) by a four-bit block index \(a\) and a word
\(\eta\in\{1,2\}^{[P]\setminus J_a}\). The latter word uses
\(P-4\) bits, so there are exactly \(16\cdot2^{P-4}=2^P\)
labels. Let \(s_\eta\) be its single-occupation exterior monomial, and set
\[
B(a,\eta)=
\Omega\wedge(e_{j_a}\bar e_{j_a})\wedge s_\eta.
\tag{8}
\]
Exactly three terms survive: the additional filled pair can be any
element of \(J_a\setminus\{j_a\}\).

Different \(J_a\) have different parity syndromes. Within one sector,
different \(\eta\) have disjoint coordinate supports. Each vector in
(8) has three nonzero coordinates, so its self-pairing is one.
Thus these \(2^P\) boundary vectors are pairwise orthonormal and
satisfy (3).

Adding each three-term boundary to its standard single label gives
the promised four-term lift. Taking the first sixteen subsets in
lexicographic order makes the construction fully explicit. A label
can be constructed with space polynomial in \(P\).

This theorem gives an all-size exact cancellation. A particular fixed
choice of sectors may be a tensor extension of a smaller construction;
the theorem alone supplies no asymptotic rank improvement. Additional
sector choices and input identifications are available for optimization.

## 4. The first size: sixty orthonormal directions and two hyperbolic pairs

Take \(P=6\), with pair indices \(0,\ldots,5\).

For each two-element odd syndrome \(\sigma\), choose the least even
pair and form the three-term boundary star as in (8). There are fifteen
such syndromes and four odd words in each. These give sixty mutually
orthonormal boundary vectors, denoted \(e_0,\ldots,e_{59}\), in
lexicographic syndrome order followed by binary odd-word order.

For the fully even sector write
\[
b_{ij}=\Omega\wedge
(e_i\bar e_i)(e_j\bar e_j).
\tag{9}
\]
Its Gram matrix is the adjacency matrix of the triangular graph on
the fifteen two-subsets: two generators pair to one exactly when their
index sets intersect in one element. Its rank is four by (6).

The following four vectors form a hyperbolic basis of a nondegenerate
four-dimensional subspace of this sector:
\[
\begin{aligned}
u_1&=b_{01}, &v_1&=b_{02},\\
u_2&=b_{03}+b_{01}+b_{02},
&v_2&=b_{04}+b_{01}+b_{02}.
\end{aligned}
\tag{10}
\]
They satisfy
\(\langle u_i,u_j\rangle=\langle v_i,v_j\rangle=0\) and
\(\langle u_i,v_j\rangle=[i=j]\).

For a norm-one vector \(e\) orthogonal to a hyperbolic pair \(u,v\),
the three vectors
\[
e+u,\qquad e+v,\qquad e+u+v
\tag{11}
\]
form an orthonormal basis of their three-dimensional space.

Keep \(e_2,\ldots,e_{59}\), and replace \(e_0,u_1,v_1\) and
\(e_1,u_2,v_2\) by their triples (11). This is an explicit
orthonormal set of sixty-four boundary vectors.

The canonical logical assignment is:

* labels \(0,\ldots,57\) receive \(e_2,\ldots,e_{59}\);
* labels \(58,59,60\) receive the first triple in (11);
* labels \(61,62,63\) receive the second triple.

The checker expands these actual exterior vectors, checks their complete
Gram matrix, checks that they are cycles with no single coordinates,
and checks that the lifted Gram matrix vanishes.

## 5. Exact compression of every six-pair boundary isometry

Write
\[
a=r+q,\qquad b=p+q,\qquad A=p(p+q+r),
\qquad T=aE+bC.
\tag{12}
\]

In each two-odd-pair sector \(\sigma\), the three-dimensional even
boundary space has only one effective functional. The section therefore
induces a map
\[
\Gamma_\sigma:H\longrightarrow K^4.
\]
The fully even sector induces an alternating Gram matrix \(G_0\) of
rank at most four. For an isometry these maps satisfy
\[
\sum_\sigma\Gamma_\sigma^{\mathsf T}\Gamma_\sigma+G_0=I_{64}.
\tag{13}
\]

Define \(Z_0\) by using \(E\otimes E\) on the two odd output sites,
with \(\Gamma_\sigma\) on each input, and the identity condition
\([x_i=u_i]\) on the other four output sites.

Define \(Z_1\) similarly using
\(E\otimes C+C\otimes E\), and add
\([x=u]\,(G_0)_{y,v}\). Then the complete section tensor is
\[
\boxed{
G_B(p,q,r)
=T^{\otimes6}
+A^2\bigl(a^2 Z_0+ab Z_1+b^2 C^{\otimes6}\bigr).
}
\tag{14}
\]
The last direction is exactly \(C^{\otimes6}\) by (13).

This follows directly from the all-size boundary Gram formula. Four-even
sectors contribute \(A^2T^{\otimes2}\); the six-even sector contributes
\[
A^2(A+qr)=A^2(ab+b^2).
\]
No inactive boundary direction affects (14).

Thus seven binary coefficient directions and \(Z_0,Z_1\) contain the
entire coefficient space.

### Isometry-invariant ranks

For every such boundary isometry, over every characteristic-two field,
the exact ranks are
\[
\begin{array}{c|rr}
 &r&c\\ \hline
Z_0&3840&154\\
Z_1&1696&18\\
Z_0+Z_1&4096&154\\
C^{\otimes6}&4096&1
\end{array}.
\tag{15}
\]

For the row cut, change the input basis to the active boundary quotient.
The fifteen two-odd sectors become separate blocks. Their local row
ranks are \(16\) for \(E\otimes E\), \(6\) for
\(E\otimes C+C\otimes E\), and \(16\) for their sum.
Each has an additional identity factor of rank sixteen.
The fully even sector adds \(64\cdot4=256\) to \(Z_1\) and
\(Z_0+Z_1\). This gives the first column of (15).

For the column cut, \(Z_0\) spans the operators supported on at most
two output sites. Their dimension is
\[
1+6\cdot3+\binom62\,3^2=154.
\]
The local input-to-output map for \(Z_0+Z_1\) differs from the identity
by square-zero partial-trace maps, so it has the same image.

For \(Z_1\), there are twelve independent off-diagonal one-site
operators, five independent diagonal projection coefficients, and one
identity coefficient. The six diagonal projections sum to zero because
every two-site sector is counted twice. There are no further relations,
giving \(12+5+1=18\).

The checker also verifies the two-site operator dimension and the
compressed nineteen-vector calculation of rank eighteen directly.

## 6. Exact first-jet cancellation is impossible

At \((p,q,r)=(1,0,0)\), both \(T^{\otimes6}\) and the correction
have the constant term \(C^{\otimes6}\), which cancels.

The first nonzero direction before the correction is
\[
J_1=\sum_{j=1}^6
C^{\otimes(j-1)}\otimes E\otimes C^{\otimes(6-j)}.
\tag{16}
\]
It has row rank \(2016\). This is the rank of the nonzero linear
elementary-symmetric multiplier in six rank-one square-zero local
operators. Explicitly it is \((4^6-2^6)/2\).

The corrected first derivative is \(J_1+Z_1\). By (15),
\[
r(J_1+Z_1)\ge2016-1696=320.
\tag{17}
\]
Hence no boundary isometry cancels this derivative completely.

A more local comparison explains the missing dimensions. Fix output
site \(j\) and an off-diagonal matrix unit there. Its input coefficient
in \(Z_1\) comes from the five sectors containing \(j\), with rank
two in each. A boundary isometry changes this by invertible congruence,
so its rank is ten. The corresponding logical input matrix unit in
\(J_1\) has rank thirty-two. Their difference has rank at least
twenty-two. The fully even sector affects only the diagonal output.

These statements exclude exact jet matching, while permitting substantial
partial cancellation.

### An explicit six-piece crossing

Set \(p=1+z,\ q=0,\ r=z\), and let
\[
B_i=[z^i]\bigl(C+z(E+C)\bigr)^{\otimes6}.
\]
Then (14) has zero constant term. Its six nonzero positions are
\[
\begin{array}{c|l}
z^1&B_1+Z_1\\
z^2&B_2+Z_0+Z_1\\
z^3&B_3+Z_1\\
z^4&B_4+Z_0+Z_1+C^{\otimes6}\\
z^5&B_5\\
z^6&B_6.
\end{array}
\tag{18}
\]
Their sum is the exact crossing. This is a concrete complete decomposition,
but its rank sum still requires optimization; fewer coefficient positions
alone do not establish a better seed.

## 7. A structured partial-matching assignment

Assigning a two-site boundary module to a logical coordinate face with
the same two varying bits aligns its matrix units with actual logical bit
flips. The following thirteen disjoint faces provide an explicit assignment.
A face \((i,j;b)\) consists of
\(b,\ b+2^i,\ b+2^j,\ b+2^i+2^j\), with bits \(i,j\) absent from \(b\).

| Varying bits | Base |
|---|---:|
| \(0,1\) | 0 |
| \(2,3\) | 16 |
| \(0,5\) | 4 |
| \(1,3\) | 17 |
| \(1,5\) | 21 |
| \(0,3\) | 6 |
| \(2,5\) | 26 |
| \(3,4\) | 35 |
| \(1,2\) | 40 |
| \(0,4\) | 32 |
| \(1,4\) | 45 |
| \(0,2\) | 8 |
| \(2,4\) | 34 |

The unused orientations are \((3,5)\) and \((4,5)\). The twelve
remaining logical labels are
\[
10,11,18,22,29,31,39,41,52,56,57,60.
\]
Assign six untouched stars from the two remaining modules to the first
six labels, then assign the two triples (11) to the final six labels.

This gives a second explicit isometry and aligns fifty-two logical
flip edges. Its complete compressed polynomial identity is checked
independently against literal local Pfaffians. It is an optimization
candidate, not a claim that all remaining cuts improve.

## 8. Reproduction

The checker is
[check_noncommutative_crossing_pfaffian_isotropic.py](../tools/check_noncommutative_crossing_pfaffian_isotropic.py).
Its output is
[noncommutative-crossing-pfaffian-isotropic.json](../evidence/noncommutative-crossing-pfaffian-isotropic.json).

It verifies the actual canonical and thirteen-face boundary vectors,
their Gram identities, all active operator dimensions used above,
and 512 polynomial entries for each construction. The sample includes
deliberately nonzero entries, rather than only random zeros.
It also checks (6) through twelve even pairs and the sparse all-size
construction at \(P=7,8,9,10\).

The ordinary proofs establish the all-size statements. Complete finite
rank searches for a chosen isometry are a separate task and are not
replaced by the correction-only rank table (15).
