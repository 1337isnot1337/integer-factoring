# Finite rank-sum gates for crossing constructions

This note gives three exact exclusions before any coefficient search:

1. Every crossing decomposition on at most five labels misses the strict finite-seed target, over every field.
2. Every orbit-kernel decomposition on a group of order six, seven, or eight that is not elementary abelian of exponent two also misses that target, over every field.
3. On a finite two-group in characteristic two, constant-augmentation pieces satisfy a stronger bound. This includes every coefficient recombination of the word deformation, with an arbitrary control map.

These statements concern the rank-sum target
\[
 S=\sum_i\sqrt{r_i c_i}<d^{\sqrt2}.
 \tag{1}
\]
The finite-seed transfer using this target is recorded in
[simultaneous-seed-endpoint.md](simultaneous-seed-endpoint.md).
No Lean formalization is used here.

## 1. Conventions

Let
\[
 E(x,y,u,v)=[x=v][y=u]
\]
be the crossing on \(d\) labels. For a four-tensor \(T\), write
\[
 r(T)=\operatorname{rank}T_{(x,y),(u,v)},\qquad
 c(T)=\operatorname{rank}T_{(x,u),(v,y)}.
\]
Both cuts of \(E\) have rank \(n=d^2\). We consider exact decompositions
\(E=\sum_{i=1}^m T_i\), with zero pieces removed. Thus
\[
 1\le r_i,c_i\le n,\qquad
 \sum_i r_i\ge n,\quad \sum_i c_i\ge n.
 \tag{2}
\]

All ranks in this note are ordinary matrix ranks over the coefficient field.
The first two results impose no characteristic restriction.

## 2. A universal gate

**Theorem 1.** For every \(d\ge2\) and every exact crossing decomposition,
\[
 \boxed{\quad S\ge2\sqrt{d^2-1}.\quad}
 \tag{3}
\]

**Proof.** For \(r,c\ge1\),
\[
 \sqrt{rc}\ge\sqrt r+\sqrt c-1
\]
because \((\sqrt r-1)(\sqrt c-1)\ge0\).

If \(2\le m\le n\), concavity of the square root and the constraints
\(r_i\ge1,\ \sum r_i\ge n\) give
\[
 \sum_i\sqrt{r_i}\ge\sqrt{n-m+1}+m-1.
 \tag{4}
\]
One can prove (4) without an optimization theorem: lower the sum of the
\(r_i\) to \(n\), then transfer mass from a smaller coordinate to a larger
coordinate until all but one equal one. Every transfer weakly decreases
the sum of square roots. The same argument applies to \(c_i\).

Consequently
\[
 S\ge2\sqrt{n-m+1}+m-2.
\]
The expression is increasing for real \(m\in[2,n]\), since its derivative is
\(1-1/\sqrt{n-m+1}\ge0\). Its minimum is \(2\sqrt{n-1}\).
If \(m>n\), then \(S\ge m>n\ge2\sqrt{n-1}\).
If \(m=1\), the only piece is \(E\), so \(S=n\).
This proves (3). \(\square\)

**Corollary 1.** No crossing decomposition on \(d=2,3,4,5\) satisfies (1).

This comparison needs no floating-point arithmetic. Since
\[
 \sqrt2<\frac{99}{70},
\]
the four integer inequalities
\[
 2^{140}(d^2-1)^{70}>d^{198},
 \qquad 2\le d\le5,
 \tag{5}
\]
imply \(2\sqrt{d^2-1}>d^{99/70}>d^{\sqrt2}\).
The accompanying checker verifies (5) exactly.

## 3. A column-rank-one orbit tensor forces an elementary binary group

Let \(G\) be any finite group of order \(d\), with identity \(e\). For an
arbitrary kernel \(K:G\times G\to F\), define
\[
 T_K(x,y,u,v)
 =[xy^{-1}=vu^{-1}]K(xy^{-1},y^{-1}u).
 \tag{6}
\]
A nonzero kernel entry \(K(S,a)\) appears at
\[
 x=Sy,\qquad u=ya,\qquad v=Sya
 \tag{7}
\]
for every \(y\in G\). In particular, \(K\mapsto T_K\) is injective.

**Theorem 2.** If \(T_K\ne0\) and \(c(T_K)=1\), then every element of \(G\)
has square \(e\). Moreover, \(r(T_K)=d^2\).

**Proof.** Factor the rank-one column cut as
\[
 T_K(x,y,u,v)=A(x,u)B(v,y).
 \tag{8}
\]
Choose a nonzero entry \(K(S,a)\). As \(y\) ranges over \(G\), (7) shows that
every row of \(A\) and every column of \(B\) is nonzero.

In particular, choose \(c\in G\) with \(B(c,e)\ne0\). Whenever \(A(x,u)\ne0\),
the support condition in (6), with \(y=e,v=c\), forces
\[
 u=x^{-1}c.
 \tag{9}
\]
Thus \(A\) has exactly one nonzero entry in each row, on the graph of a
permutation; in particular \(A(e,c)\ne0\). Fixing this entry in (8) forces
every nonzero entry of \(B\) to satisfy
\[
 v=y^{-1}c.
 \tag{10}
\]
Because each column is nonzero, \(B\) is also an invertible monomial matrix.

Now \(K(e,c)=T_K(e,e,c,c)\ne0\). Formula (7), with this kernel entry, gives
\(T_K(y,y,yc,yc)\ne0\) for every \(y\). Equation (9) requires
\(yc=y^{-1}c\), hence \(y^2=e\).

Finally, the row cut in (8) is \(A\otimes B^{\mathsf T}\), so it has rank
\(d^2\). \(\square\)

A group in which every element has square \(e\) is abelian:
\(xy=(xy)^{-1}=y^{-1}x^{-1}=yx\). It is therefore an elementary abelian
two-group. Conversely, on such a group, the kernel
\[
 K(S,a)=[a=S^{-1}c]
\]
has column rank one and row rank \(d^2\). Thus the group condition in
Theorem 2 is exact.

**Corollary 2.** If \(G\) is not elementary abelian of exponent two, then
every nonzero orbit tensor (6) has \(c(T_K)\ge2\).

There is no augmentation hypothesis in this corollary. The kernel may be
arbitrary, and the coefficient field may have any characteristic.

## 4. The resulting group gate

**Theorem 3.** Suppose every nonzero piece of a crossing decomposition has
column rank at least two. Then
\[
 \boxed{\quad
 S\ge\sqrt{2(d^2-1)}+\sqrt{d^2-2}.
 \quad}
 \tag{11}
\]

**Proof.** Put \(n=d^2\). Since \(r\ge1,c\ge2\),
\[
 \sqrt{rc}\ge\sqrt2\,\sqrt r+\sqrt c-\sqrt2.
 \tag{12}
\]
For \(2\le m\le n/2\), the mass-transfer argument used in (4) gives
\[
 S\ge F(m):=
 \sqrt2\sqrt{n-m+1}+\sqrt{n-2m+2}+(m-2)\sqrt2.
 \tag{13}
\]
The derivative is nonnegative throughout this interval:
\[
 F'(m)=\sqrt2-\frac{\sqrt2}{2\sqrt{n-m+1}}
                    -\frac1{\sqrt{n-2m+2}}
 \ge\frac1{\sqrt2}
       -\frac{\sqrt2}{2\sqrt{n-m+1}}\ge0.
\]
For \(n/2\le m\le n\), use \(\sum_i\sqrt{c_i}\ge m\sqrt2\) instead. This
gives
\[
 S\ge G(m):=\sqrt2\bigl(\sqrt{n-m+1}+m-1\bigr),
 \tag{14}
\]
which is increasing. At \(m=n/2\), (13) and (14) agree. Thus the minimum
over \(2\le m\le n\) occurs at \(m=2\), giving (11).

For \(m>n\), \(S\ge m\sqrt2>n\), which is larger than the proposed bound.
For \(m=1\), \(S=n\), and
\[
 \sqrt{2(n-1)}+\sqrt{n-2}\le n\qquad(n\ge4).
\]
Indeed, the function \(n-\sqrt{2(n-1)}-\sqrt{n-2}\) is increasing for
\(n\ge4\): each of the two subtracted derivative terms is less than
\(1/2\). Its value at \(n=4\) is \(4-\sqrt6-\sqrt2>0\), since
\((\sqrt6+\sqrt2)^2=8+4\sqrt3<16\). \(\square\)

**Corollary 3.** For every field, no orbit-kernel decomposition on a
non-elementary group of order \(d=6,7,8\) meets (1).

Here are exact rational certificates. Each entry in the second column is
a strict lower bound for (11); each entry in the third follows from an
integer power comparison and \(\sqrt2<99/70\).

| \(d\) | Rank sum exceeds | Target is below |
|---:|---:|---:|
| 6 | \(\sqrt{70}+\sqrt{34}>8+5=13\) | \(6^{\sqrt2}<6^{99/70}<13\) |
| 7 | \(\sqrt{96}+\sqrt{47}>48/5+32/5=16\) | \(7^{\sqrt2}<7^{99/70}<16\) |
| 8 | \(\sqrt{126}+\sqrt{62}>56/5+39/5=19\) | \(8^{\sqrt2}<8^{99/70}<19\) |

The checker verifies \(13^{70}>6^{99}\), \(16^{70}>7^{99}\), and
\(19^{70}>8^{99}\). Together with Corollary 1, this excludes every
non-elementary group of order at most eight.

In particular, the conclusion applies to \(S_3\) in characteristic two
even though its group algebra is not a two-group algebra. No unit
argument is needed for that case.

## 5. A stronger gate for constant augmentation on two-groups

This section assumes \(\operatorname{char}F=2\) and \(G\) is a finite
two-group.

### 5.1. Augmentation one means a convolution unit

Let \(I=\ker(\varepsilon:F[G]\to F)\) be the augmentation ideal. It is
nilpotent. Here is an elementary induction.

A nontrivial finite two-group has a central involution \(z\). For example,
the class equation shows that the center has even order, and repeated
squaring of a nonidentity central element produces an involution.
The kernel of
\[
 F[G]\longrightarrow F[G/\langle z\rangle]
\]
is \(J=(1+z)F[G]\). It satisfies \(J^2=0\), since \(z\) is central and
\((1+z)^2=0\). If the augmentation ideal of the quotient has nilpotence
index at most \(d/2\), then \(I^{d/2}\subseteq J\), so \(I^d=0\).
The trivial group starts the induction.

Consequently, an element of nonzero augmentation \(\alpha\) can be written
\(\alpha(1+b)\), with \(b\in I\), and has inverse
\[
 \alpha^{-1}\sum_{j=0}^{d-1}b^j.
 \tag{15}
\]

### 5.2. Constant augmentation forces a full row cut

For the orbit tensor (6), reindex row pairs by \((S,y)\), where \(x=Sy\),
and column pairs by \((S,u)\), where \(v=Su\). The row cut is block diagonal:
\[
 r(T_K)=\sum_{S\in G}
 \operatorname{rank}\bigl[K(S,y^{-1}u)\bigr]_{y,u\in G}.
 \tag{16}
\]
Each block is a convolution matrix.

Suppose each piece has augmentation
\[
 \sum_a K_i(S,a)=\alpha_i
\]
independent of \(S\). Since the crossing kernel is \(K_E(S,a)=[a=e]\),
an exact decomposition gives \(\sum_i\alpha_i=1\). At least one \(\alpha_i\)
is nonzero. By (15) and (16), that piece has row rank \(d^2\).

**Theorem 4.** Let \(G\) be a finite two-group that is not elementary
abelian, and let all pieces have augmentation constant in \(S\). Then
\[
 \boxed{\quad S\ge d\sqrt2+\sqrt{d^2-2}.\quad}
 \tag{17}
\]

**Proof.** Choose a piece with nonzero augmentation and column rank \(c\).
It has row rank \(d^2\), and Corollary 2 gives \(2\le c\le d^2\).
Column-rank subadditivity requires the remaining pieces to have total
column rank at least \(d^2-c\). Their total cost is at least
\(\sqrt{d^2-c}\), since every nonzero row rank is at least one and
\(\sum_j\sqrt{c_j}\ge\sqrt{\sum_jc_j}\). Therefore
\[
 S\ge d\sqrt c+\sqrt{d^2-c}.
\]
This is a concave function of \(c\), so its minimum on \([2,d^2]\) occurs
at an endpoint. For \(d\ge3\), the endpoint \(c=2\) is smaller than
\(c=d^2\), giving (17). Every non-elementary two-group has \(d\ge4\).
\(\square\)

For \(d=8\), the bound in (17) is strictly larger than \(19\), because
\(\sqrt2>7/5\) and \(\sqrt{62}>39/5\). The target is strictly smaller
than \(19\), as above.

### 5.3. Arbitrary control maps are included

Fix group generators \(g_1,\ldots,g_\ell\), positive integers
\(w_1,\ldots,w_\ell\), and any map
\(\theta:G\to\{0,1\}^{\ell}\). No bijectivity or triangularity is required.
Form the kernel in the group algebra
\[
 K_z(S,\cdot)=
 \prod_{j=1}^{\ell}
 \bigl(g_j+z^{w_j}(1+g_j)\bigr)^{\theta_j(S)}
 \tag{18}
\]
in a fixed order. Every factor has augmentation one, so
\[
 \varepsilon K_z(S,\cdot)=1
 \quad\text{for all }S.
\]
At \(z=1\), every factor is the group-algebra identity, hence
\(K_1(S,a)=[a=e]\) and \(T_{K_1}=E\).

The constant coefficient of (18) has augmentation one; every positive
coefficient has augmentation zero. Thus every scalar recombination of
the coefficients has augmentation constant in \(S\). The same remains
true for a space spanned by coefficients from several control maps,
weight assignments, or factor orders. Theorem 4 applies to every exact
decomposition drawn from any such space.

This result is separate from the filtered group-algebra argument in
[noncommutative-crossing-group-degeneration.md](noncommutative-crossing-group-degeneration.md).
It includes noncentral control maps at order eight, but it does not give
an all-size binary conservation theorem for arbitrary controls.

## 6. Reproducible checks and scope

Run:
\[
 \texttt{python tools/check\_noncommutative\_crossing\_finite\_gates.py}
\]
or pass \(\texttt{--output evidence/noncommutative-crossing-finite-gates.json}\).

The checker verifies:

- the integer target comparisons in Sections 2 and 4;
- all forced rank-one permutation-support candidates for
  \(D_8,Q_8,C_4\times C_2,C_8,C_2^3,S_3,D_{10},C_3,C_9\);
- every binary augmentation-one convolution kernel for the five
  order-eight groups above;
- sixteen arbitrary control maps on each of \(D_8,Q_8\), with three
  weight assignments, checking every coefficient recombination that
  has nonzero augmentation.

These finite checks audit the formulas; the theorems are the all-size
proofs above. The bounds are lower bounds on rank sums. They do not
supply a new factoring upper bound. The universal theorem applies to
all tensors, the orbit theorem applies to arbitrary kernels on the
specified group, and the stronger unit theorem requires the stated
constant-augmentation and characteristic-two hypotheses.
