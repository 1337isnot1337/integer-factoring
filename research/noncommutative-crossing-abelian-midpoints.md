# Abelian orbit kernels: midpoint blocks and an all-cyclic exclusion

For every finite abelian group \(G\), the column cut of an orbit tensor
has an exact block multiplicity controlled by the doubling map
\(x\mapsto2x\). The proof is a change of indices, and works over every
coefficient field.

The consequences are:

- for odd-order \(G\), every exact orbit decomposition has rank sum at
  least \(|G|^{3/2}\);
- more generally, the lower bound is
  \(|G|^{3/2}/\sqrt{|G[2]|}\);
- combining the block multiplicity with the finite rank-sum gate excludes
  **every cyclic group, of every order**, from the strict target
  \(S<d^{\sqrt2}\), even with arbitrary orbit kernels.

The last statement also holds for any abelian group with at most two
elements in \(G[2]=\{x:2x=0\}\). No group-algebra unit, Fourier, control-map,
or characteristic hypothesis is required.

## 1. The exact column blocks

Write the group additively and put \(d=|G|\). For an arbitrary kernel
\(K:G\times G\to F\), let
\[
 T_K(x,y,u,v)=[x-y=v-u]K(x-y,u-y).
 \tag{1}
\]
The row and column cut conventions are
\[
 r(T)=\operatorname{rank}T_{(x,y),(u,v)},\qquad
 c(T)=\operatorname{rank}T_{(x,u),(v,y)}.
\]
The support condition is equivalent to
\[
 t=x+u=v+y.
\]
Thus the column cut is block diagonal with one block for each \(t\in G\).
Using \(x\) and \(y\) as free coordinates inside the block gives
\[
 M_t(x,y)=K(x-y,t-x-y),\qquad
 c(T_K)=\sum_{t\in G}\operatorname{rank}M_t.
 \tag{2}
\]

For every \(h\in G\),
\[
 M_{t+2h}(x+h,y+h)=M_t(x,y).
 \tag{3}
\]
Therefore the blocks indexed by one coset of \(2G\) are equivalent by row
and column permutations. If \(R\) is a set of coset representatives,
\[
 \boxed{\quad
 c(T_K)=|2G|\sum_{t\in R}\operatorname{rank}M_t.
 \quad}
 \tag{4}
\]
This is an equality, not just a divisibility bound.

The map \(K\mapsto T_K\) is injective: every value \(K(S,a)\) occurs, for
example, at \(x=S,y=0,u=a,v=S+a\). Hence a nonzero kernel has a nonzero
column cut. Formula (4) implies
\[
 K\ne0\quad\Longrightarrow\quad c(T_K)\ge|2G|.
 \tag{5}
\]

For an exact crossing decomposition \(E=\sum_i T_{K_i}\), remove zero
pieces. Since both cuts of \(E\) have rank \(d^2\),
\[
 \sum_i\sqrt{r_i c_i}
 \ge\sqrt{|2G|}\sum_i\sqrt{r_i}
 \ge\sqrt{|2G|}\sqrt{\sum_i r_i}
 \ge d\sqrt{|2G|}.
 \tag{6}
\]
The kernel of doubling has size \(|G[2]|\), so
\[
 \boxed{\quad S\ge\frac{d^{3/2}}{\sqrt{|G[2]|}}.\quad}
 \tag{7}
\]

## 2. Odd order: one common midpoint block

If \(G\) has odd order, doubling is invertible. Define
\[
 t=\frac{x+u}{2}=\frac{v+y}{2},\qquad
 b=x-u,\quad c=v-y.
\]
Then
\[
 x=t+b/2,\quad u=t-b/2,\quad
 v=t+c/2,\quad y=t-c/2,
\]
and
\[
 x-y=(b+c)/2,\qquad u-y=(c-b)/2.
\]
Every column block is consequently the same matrix
\[
 B(b,c)=K\bigl((b+c)/2,(c-b)/2\bigr).
\]
Thus
\[
 \boxed{\quad c(T_K)=d\,\operatorname{rank}B.\quad}
 \tag{8}
\]
The change from \((b,c)\) to the two kernel arguments is bijective, so a
nonzero kernel has \(\operatorname{rank}B\ge1\). In particular,
\[
 S\ge d^{3/2}>d^{\sqrt2}\qquad(d>1).
 \tag{9}
\]

This closes the odd abelian orbit lane over every field, including cases
where the coefficient characteristic divides the group order. The
argument is independent of semisimplicity.

## 3. A rank floor gives a sharper finite bound

The following variant of the gate in
[noncommutative-crossing-finite-gates.md](noncommutative-crossing-finite-gates.md)
is useful when \(|2G|\) is moderately large.

**Lemma.** Suppose an exact crossing decomposition on \(d\ge2\) labels has
\(c_i\ge q\) for every nonzero piece, where \(1\le q\le d^2/2\) is real.
Put \(n=d^2\). Then
\[
 S\ge
 \min\left\{n,\sqrt{q(n-1)}+\sqrt{n-q}\right\}.
 \tag{10}
\]

**Proof.** For \(r\ge1,c\ge q\),
\[
 \sqrt{rc}\ge\sqrt q\,\sqrt r+\sqrt c-\sqrt q,
 \tag{11}
\]
because \((\sqrt r-1)(\sqrt c-\sqrt q)\ge0\).

If the number \(m\) of nonzero pieces satisfies \(2\le m\le n/q\),
concavity of the square root and the two rank-sum constraints give
\[
 S\ge F_q(m):=
 \sqrt q\sqrt{n-m+1}+\sqrt{n-q(m-1)}+(m-2)\sqrt q.
 \tag{12}
\]
Throughout this interval,
\[
 F_q'(m)
 =\sqrt q-\frac{\sqrt q}{2\sqrt{n-m+1}}
             -\frac{q}{2\sqrt{n-q(m-1)}}
 \ge\frac{\sqrt q}{2}
       -\frac{\sqrt q}{2\sqrt{n-m+1}}\ge0.
\]
If \(n/q\le m\le n\), use
\(\sum_i\sqrt{c_i}\ge m\sqrt q\), obtaining
\[
 S\ge G_q(m):=\sqrt q\bigl(\sqrt{n-m+1}+m-1\bigr).
 \tag{13}
\]
This is increasing, and agrees with \(F_q\) at \(m=n/q\).
For \(m>n\), \(S\ge m\sqrt q>n\). For \(m=1\), \(S=n\).
Thus the minimum is either the one-piece value \(n\), or
\(F_q(2)=\sqrt{q(n-1)}+\sqrt{n-q}\). \(\square\)

The use of a real \(q\) is legitimate: it is only a lower bound on the
integer ranks. This lets us use \(q=d/2\) uniformly below.

## 4. Every cyclic group is excluded

**Theorem.** Let \(G\) be a finite abelian group of order \(d\ge2\) with
\(|G[2]|\le2\). Then every exact orbit-kernel decomposition satisfies
\[
 \boxed{\quad S>d^{17/12}>d^{\sqrt2}.\quad}
 \tag{14}
\]

**Proof.** Equation (5) gives \(c_i\ge |2G|\ge d/2\). Apply (10) with
\(q=d/2\). The second entry of its minimum is at most \(n\), because
\[
 \sqrt{(d/2)(d^2-1)}+\sqrt{d^2-d/2}
 \le d\bigl(\sqrt{d/2}+1\bigr)\le d^2
\]
for \(d\ge2\). Hence
\[
 S\ge\sqrt{(d/2)(d^2-1)}+\sqrt{d^2-d/2}
 \ge\frac{\sqrt3}{2}d\bigl(\sqrt{d/2}+1\bigr).
 \tag{15}
\]

Put \(t=\sqrt d\). Weighted arithmetic-geometric mean, applied to five
copies of \(t/(5\sqrt2)\) and one copy of \(1\), gives
\[
 \left(\frac{t}{\sqrt2}+1\right)^6
 \ge\frac{6^6}{5^5\,2^{5/2}}t^5.
\]
Consequently
\[
 \left[
 \frac{(\sqrt3/2)(t/\sqrt2+1)}{t^{5/6}}
 \right]^{12}
 \ge\frac{3^{18}}{2^5\,5^{10}}>1.
 \tag{16}
\]
The strict inequality is the integer comparison
\(387420489>312500000\). Thus the right side of (15) is strictly larger
than \(d\,t^{5/6}=d^{17/12}\). Finally,
\(17^2=289>288=2\cdot12^2\), so \(17/12>\sqrt2\). \(\square\)

A cyclic group has one element killed by doubling when its order is odd
and two when its order is even. Therefore (14) excludes every cyclic
group. More generally, it covers every finite abelian group whose
two-primary subgroup is cyclic.

The bound does not exclude every abelian group with larger two-torsion.
For example, \(|2G|=1\) for an elementary abelian two-group, so (7) alone
gives only \(S\ge d\). The separate full-domain binary conservation
theorems require their own hypotheses.

## 5. Exact checks and scope

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_abelian\_midpoints.py}
\]
with optional output
\(\texttt{--output evidence/noncommutative-crossing-abelian-midpoints.json}\).

The checker uses twelve abelian groups through order sixteen, including
cyclic groups, \(C_2^2,C_2\times C_4,C_3^2\), and
\(C_2^2\times C_3\). Over both \(\mathbb F_2\) and \(\mathbb F_3\), it checks
the exact entry identity (3), rank divisibility (4), the odd-order midpoint
formula (8), and direct dense versions of both cuts against the block
formulas. It also verifies the integer constants in (16).

The proof applies to arbitrary scalar combinations and arbitrary
kernels, not only polynomial word coefficients. It requires that each
piece be an orbit tensor on the same abelian group. Arbitrary tensors
outside this orbit space, or restrictions that change the full group
label domain, are not covered.
