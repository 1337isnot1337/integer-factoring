# Exact one-particle reduction of the four-state Pfaffian

## Result

The complete local Pfaffian tensor has an exact exterior-power description.
It follows from a **permutation of the monomial basis**, with no quotient,
rank relaxation, or loss of the homogeneous energy constraints.

Put

\[
 \Delta=p+r,\qquad
 g=\begin{pmatrix}r&p\\1&1\end{pmatrix},\qquad
 h=\begin{pmatrix}1&1\\q&p+q+r\end{pmatrix},\qquad G=h\otimes g.
 \tag{1}
\]

Both two-dimensional matrices have determinant \(\Delta\). The local
six-edge Pfaffian uses

\[
 (UA,AY,AV,UV,VY,UY)=(p,q,r,1,1,1),
 \tag{2}
\]

with the usual equality/complementation rules for physical mode indices.
As a row-cut matrix, it is exactly

\[
 \boxed{\quad
 \Delta\,\Lambda^\bullet(G/\Delta).
 \quad}
 \tag{3}
\]

The apparent denominators cancel. Formula (3) is a polynomial identity,
including at \(\Delta=0\).

Across \(P\) complementary physical mode pairs, the entire tensor is thus
the appropriate grade and energy restriction of

\[
 \boxed{\quad
 \Delta^P\,
 \Lambda^\bullet\bigl((G/\Delta)\otimes I_P\bigr).
 \quad}
 \tag{4}
\]

The other flattening has the same formula after exchanging \(p,q\) and
complementing the two state legs. This identifies both **actual** complete
cut matrices, rather than merely ambient channel spaces.

The result also gives a stronger finite block decomposition. Local signed
charges, taking five possible sum values, replace the previous binary
parity syndrome. The supplied implementation computes exact complete ranks
with this smaller collection of entries. It does not claim a better
factoring exponent.

## 1. The particle-hole basis

For one physical pair, let its low and high modes be \(a\) and \(M-a\).
Write the parameter occupation as
\(x=(x_{\rm low},x_{\rm high})\), and the state occupation as
\(y=(y_{\rm low},y_{\rm high})\).

The new four occupation bits are

\[
 \boxed{
 (n_1,n_2,n_3,n_4)
 =(1-y_{\rm low},\ 1-x_{\rm low},\ x_{\rm high},\ y_{\rm high}).
 }
 \tag{5}
\]

This simply exchanges particles and holes in the two low-mode coordinates
and reorders the coordinates. It is a bijection on all sixteen monomials.
The new vacuum is the old state in which both legs occupy the low mode.

The row indices use \((x,y)\); the column indices use \((v,u)\).
The same map (5) is used on both.
In occupation-mask order \(x,y=0,1,2,3\), the new bit masks are

\[
 (3,2,11,10,\ 1,0,9,8,\ 7,6,15,14,\ 5,4,13,12).
 \tag{6}
\]

## 2. The five exterior-degree blocks

The matrix has no entry between different new exterior degrees. Its block
sizes are \(1,4,6,4,1\).

### Degrees zero and four

The old states \((x,y)=(1,1)\) and \((2,2)\) are isolated. Both carry the
scalar \(p+r=\Delta\).

### Degree one

In the old-state order

\[
 (1,0),\ (0,1),\ (3,1),\ (1,3),
\]

the block is

\[
 G=
 \begin{pmatrix}
 r&p&r&p\\
 1&1&1&1\\
 qr&pq&r(p+q+r)&p(p+q+r)\\
 q&q&p+q+r&p+q+r
 \end{pmatrix}
 =h\otimes g.
 \tag{7}
\]

The new degree-three states, in increasing bit-mask order
\(7,11,13,14\), give an identical block.

### Degree two

In the old-state order

\[
 (0,0),\ (3,0),\ (1,2),\ (2,1),\ (0,3),\ (3,3),
\]

the block is

\[
 H=
 \begin{pmatrix}
 1&0&1&1&0&1\\
 0&r^2&pr&pr&p^2&0\\
 q&r&r+q&p+q&p&p+q+r\\
 q&r&p+q&r+q&p&p+q+r\\
 0&1&1&1&1&0\\
 q^2&0&q(p+q+r)&q(p+q+r)&0&p^2+q^2+r^2
 \end{pmatrix}.
 \tag{8}
\]

With wedge basis \(12,13,14,23,24,34\), its identity is exactly

\[
 \Delta H=\Lambda^2G.
 \tag{9}
\]

Thus every two-by-two minor of \(G\) is divisible by \(\Delta\), and (8)
gives the polynomial quotient.

### Completion of the proof

The two-dimensional symplectic matrix
\(J_2=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)\) satisfies
\(g^{\mathsf T}J_2g=h^{\mathsf T}J_2h=\Delta J_2\).
Consequently \(G\) preserves \(J_2\otimes J_2\) up to the multiplier
\(\Delta^2\), and \(\det G=\Delta^4\).
The degree-three exterior matrix is therefore \(\Delta^2G\), after the
displayed complementary-index identification. Its factor
\(\Delta^{1-3}\) in (3) leaves \(G\). Degree four gives
\(\Delta^{-3}\det G=\Delta\). Together with (7) and (9), this proves all
five blocks of (3).

Alternatively, the accompanying checker verifies all 256 generic
polynomial entries directly against the independent Pfaffian-matching
definition. For degree \(k\geq1\), it checks

\[
 \Delta^{k-1}T_{I,J}=\det G[I,J].
 \tag{10}
\]

This direct polynomial verification also covers the specialization
\(\Delta=0\), without dividing in a finite field.

## 3. Original grades and energies remain exact

Let \(N_b=\sum_i n_{i,b}\), summing new bits over the \(P\) physical pairs.
The old parameter and state grades are

\[
 \boxed{
 L=P+N_3-N_2,\qquad K=P+N_4-N_1.
 }
 \tag{11}
\]

Thus a restriction to any original pair of grades is an explicit weight
restriction of the exterior representation.

For a description using the two tensor factors of \(G=h\otimes g\), let
\(A=N_1+N_2+N_3+N_4\) be the new exterior degree, and put

\[
 G_c=N_2+N_4,\qquad H_c=N_3+N_4.
\]

Then

\[
 L=P+H_c-G_c,\qquad
 K=P+H_c+G_c-A.
 \tag{12}
\]

Equivalently, in exterior degree \(A\), the selected torus weights are

\[
 G_c=(K-L+A)/2,\qquad
 H_c=(L+K-2P+A)/2.
 \tag{13}
\]

These formulas include both even and odd new exterior degrees.

For low modes \(a_1,\ldots,a_P\), write \(E_0=\sum_i a_i\). The old
parameter energy \(j\) and state energy \(E\) become

\[
 \boxed{
 \begin{aligned}
 j&=E_0+\sum_i\bigl((M-a_i)n_{i,3}-a_i n_{i,2}\bigr),\\
 E&=E_0+\sum_i\bigl((M-a_i)n_{i,4}-a_i n_{i,1}\bigr).
 \end{aligned}}
 \tag{14}
\]

For the original power-of-two construction,
\(P=M/4\), \(a_i=2i-1\), and \(E_0=P^2\).
No averaging or loss of shell homogeneity is involved in (14).

Taking exterior powers of a direct sum gives (4). More explicitly, in
new exterior degree \(A\), the operator is

\[
 \Delta^{P-A}\Lambda^A(G\otimes I_P),
 \tag{15}
\]

restricted to the affine grade and energy conditions (11), (14).
The local polynomial cancellations prove that the resulting operator is
polynomial, including when the displayed exponent of \(\Delta\) is
negative.

## 4. Both cuts and the source pencil

Let \(\bar y\) complement every physical mode of a state label. The local
tensor obeys the exact identity

\[
 T_{p,q,r}(x,y,u,v)
 =T_{q,p,r}(x,\bar u,\bar y,v).
 \tag{16}
\]

The right-hand side is the row-cut formula for the other flattening.
Physical complementation preserves the state grade and sends its energy
from \(E\) to \(KM-E\). This supplies the precise grade and energy map for
the column cut.

For the original complete polynomial, with \(t=z^M\), use

\[
 p=t+t^2,\qquad q=1+t,\qquad r=t.
 \tag{17}
\]

The row-cut determinant is then \(\Delta_R=t^2\).
For the column cut, (16) instead gives \(\Delta_C=1\). In particular,
the latter is an ordinary exterior-power operator with no scalar
denominator:

\[
 \Lambda^\bullet(G_C(t)\otimes I_P),\qquad
 G_C=h_C\otimes g_C,
\]

where

\[
 g_C=\begin{pmatrix}t&1+t\\1&1\end{pmatrix},\qquad
 h_C=\begin{pmatrix}1&1\\t+t^2&1+t+t^2\end{pmatrix}.
 \tag{18}
\]

Both matrices in (18) have determinant one.
Every coefficient of these restricted exterior matrices is an **actual**
complete coefficient matrix. The only remaining difficulty in rank
calculation is the coefficient extraction and weight restriction.

## 5. Signed-charge blocks

Define the local signed charge

\[
 \operatorname{ch}(x)=x_{\rm high}-x_{\rm low}\in\{-1,0,1\}.
\]

The new local exterior degree is

\[
 2+\operatorname{ch}(x)+\operatorname{ch}(y).
\]

Hence a row-cut entry vanishes unless, at every physical pair,

\[
 \operatorname{ch}(x_i)+\operatorname{ch}(y_i)
 =\operatorname{ch}(v_i)+\operatorname{ch}(u_i).
 \tag{19}
\]

The corresponding column block label is

\[
 \operatorname{ch}(x_i)-\operatorname{ch}(u_i)
 =\operatorname{ch}(v_i)-\operatorname{ch}(y_i).
 \tag{20}
\]

Reducing these labels modulo two recovers the older parity syndrome.
The integer labels (19), (20) give a strict refinement.

The exact implementation groups row and column pairs by these labels.
It evaluates only entries within a row block, then places the same
coefficient into its independently indexed column block. Rank is the sum
of the ranks of these diagonal blocks.

For the full middle grade on ten physical modes, \(P=L=5\), the label
dimension is \(\binom{10}{5}=252\). A dense four-leg enumeration has
\(252^4=4,032,758,016\) entries. Signed charges reduce this to
9,343,296 possibly nonzero entries, in 1,563 blocks. The largest block has
712 rows. The reduction is exact.

## 6. Structure that a further rank algorithm must preserve

As a graded polynomial matrix, the local operator is

\[
 2\Delta\ \oplus\ 2(h\otimes g)\ \oplus\
 \Delta^{-1}\Lambda^2(h\otimes g).
 \tag{21}
\]

The odd blocks admit the elementary two-dimensional Frobenius
decomposition independently in \(g\) and \(h\). The six-dimensional block
requires more care.

Its character resembles the two quadratic modules \(g\otimes g\) and
\(h\otimes h\), but it is not their direct sum. In the six-dimensional
block, the map from the two central-grade vectors to the four endpoint
grades has rank one. In a direct sum of the two quadratic modules, it
would have rank two. A character-only replacement would therefore be
invalid.

Formula (9) identifies the six-dimensional extension exactly. The
[stable two-factor decomposition](stable-pfaffian-modules.md) now resolves
it after tensoring with one four-dimensional factor. Its explicit
intertwiner preserves the aggregate grades and gives polynomial-size
blocks for complete coefficient-rank upper bounds. The particle-hole
description remains exact for the finer energy conditions as well.

## Reproduce

~~~bash
python3 tools/check_pfaffian_particle_hole.py --include-large \
  --output evidence/pfaffian-particle-hole.json
~~~

The checker verifies generic multivariate identities, affine grade and
energy maps, and exact complete profiles. Small full-grade profiles and a
homogeneous shell are compared with the separate parity-block Pfaffian
implementation.

Implementation:
[pfaffian_particle_hole.py](../tools/pfaffian_particle_hole.py).
Evidence:
[pfaffian-particle-hole.json](../evidence/pfaffian-particle-hole.json).
