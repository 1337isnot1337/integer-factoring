# Involution-free minors of group orbit tensors

This addendum extends the abelian midpoint obstruction to many
nonabelian groups. Every nonzero orbit tensor contains a diagonal minor
indexed by any set whose pairwise quotients are not involutions.

As consequences, **all odd-order groups** and **all groups with an
odd-order subgroup of index two** miss the strict rank-sum seed target,
over every coefficient field. The latter includes every dihedral group
of order \(2m\) with \(m\) odd and \(S_3\times C_3\). Groups with a unique
nonidentity involution, including generalized quaternion groups, are
also excluded.

## 1. Normalize one nonzero kernel entry

Let \(G\) be a finite group of order \(d\), with identity \(e\). As in
[the finite rank-sum gates](noncommutative-crossing-finite-gates.md), put
\[
 T_K(x,y,u,v)
 =[xy^{-1}=vu^{-1}]K(xy^{-1},y^{-1}u).
 \tag{1}
\]
The column cut has rows \((x,u)\) and columns \((v,y)\).

Choose \(K(S_0,a_0)\ne0\). Define
\[
 \widetilde K(S,a)=K(S_0S,aa_0).
\]
A direct substitution gives
\[
 T_{\widetilde K}(x,y,u,v)
 =T_K(S_0x,y,ua_0,S_0va_0).
 \tag{2}
\]
These are independent permutations of the four label sets. They
preserve both cut ranks, and
\(\widetilde K(e,e)=K(S_0,a_0)\ne0\).
Thus one may normalize a nonzero kernel entry to \((e,e)\) without any
assumption on the kernel, group, or field.

## 2. The diagonal minor

Call \(A\subseteq G\) *involution-free in its quotients* if
\[
 x\ne y,\quad x,y\in A
 \quad\Longrightarrow\quad (xy^{-1})^2\ne e.
 \tag{3}
\]
Let \(\alpha_2(G)\) be the maximum size of such a set.

**Theorem.** Every nonzero orbit tensor satisfies
\[
 \boxed{\quad c(T_K)\ge\alpha_2(G).\quad}
 \tag{4}
\]

**Proof.** Normalize as in (2). In the column cut, retain the rows
\((x,x)\), \(x\in A\), and columns \((y,y)\), \(y\in A\). A retained entry
can be nonzero only if
\[
 xy^{-1}=yx^{-1}=(xy^{-1})^{-1}.
\]
By (3), this forces \(x=y\). Every diagonal entry equals
\(\widetilde K(e,e)\), which is nonzero. The retained matrix is therefore
\(\widetilde K(e,e)I_{|A|}\). \(\square\)

Equivalently, \(\alpha_2(G)\) is the independence number of the Cayley
graph whose edges join \(x\ne y\) when \(xy^{-1}\) is an involution.
The generating set of nonidentity involutions is closed under inversion
and conjugation, so this is an undirected graph. This graph description
is useful for computing the parameter; the proof needs only the
diagonal minor.

For every exact crossing decomposition into nonzero orbit tensors,
\[
 S=\sum_i\sqrt{r_i c_i}
 \ge d\sqrt{\alpha_2(G)}.
 \tag{5}
\]
When \(\alpha_2(G)\) is at least \(d/2\), the sharper rank-floor lemma in
[the abelian midpoint note](noncommutative-crossing-abelian-midpoints.md)
applies without any abelian hypothesis. Its proof uses only the lower
bound on each column rank. The same exact arithmetic-geometric-mean
certificate gives
\[
 \boxed{\quad
 \alpha_2(G)\ge d/2
 \ \Longrightarrow\
 S>d^{17/12}>d^{\sqrt2}
 \qquad(d\ge2).
 \quad}
 \tag{6}
\]

## 3. Families excluded at every size

### Odd-order groups

A finite group of odd order has no nonidentity involution. Thus
\(A=G\) is admissible, \(\alpha_2(G)=d\), and
\[
 S\ge d^{3/2}>d^{\sqrt2}\qquad(d>1).
 \tag{7}
\]
The group need not be abelian, and the field characteristic may divide
its order.

### An odd-order subgroup of index two

Let \(H\le G\) have odd order. For distinct \(x,y\in H\), the quotient
\(xy^{-1}\) belongs to \(H\) and is not an involution. Hence
\[
 \alpha_2(G)\ge |H|.
 \tag{8}
\]
If \([G:H]=2\), this is \(d/2\), so (6) excludes the seed.

For a dihedral group of order \(2m\) with \(m\) odd, use the cyclic
rotation subgroup. For \(S_3\times C_3\), use \(C_3\times C_3\).
More generally, the argument applies to any group admitting such a
subgroup; no chosen group presentation or control map is involved.

### A unique nonidentity involution

Suppose \(G\) has exactly one nonidentity involution \(z\). It is central,
because every conjugate is again an involution. The forbidden-quotient
graph consists of the disjoint pairs \(\{x,zx\}\). Choosing one element
from every pair gives
\[
 \alpha_2(G)=d/2.
 \tag{9}
\]
Therefore (6) also applies.

For clarity, the generalized quaternion group
\[
 Q_{4m}=\langle a,b:
 a^{2m}=e,\ b^2=a^m,\ bab^{-1}=a^{-1}\rangle
\]
has a unique involution whenever \(m\) is a power of two. Indeed, among
the rotations only \(a^m\) has order two, while every element \(a^jb\)
has square \(a^m\ne e\). Thus every generalized quaternion group is
covered. The same calculation covers the usual dicyclic presentations
with arbitrary \(m\ge2\).

## 4. Exact audits and limits

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_involution\_support.py}
\]
with optional output
\(\texttt{--output evidence/noncommutative-crossing-involution-support.json}\).

The checker constructs ten groups, verifies their multiplication laws,
and computes \(\alpha_2(G)\) by exact backtracking. It checks the
normalization identity, the full diagonal minors, and dense column ranks
for nonzero kernels over both \(\mathbb F_2\) and \(\mathbb F_3\).

| Group | Order | \(\alpha_2(G)\) |
|---|---:|---:|
| \(C_4\) | 4 | 2 |
| \(C_9\) | 9 | 9 |
| \(C_2^3\) | 8 | 1 |
| \(D_8\) | 8 | 2 |
| \(Q_8\) | 8 | 4 |
| \(Q_{16}\) | 16 | 8 |
| \(D_{10}\) | 10 | 5 |
| \(D_{12}\) | 12 | 3 |
| \(S_3\times C_3\) | 18 | 9 |
| \(C_7\rtimes C_3\), with nontrivial action | 21 | 21 |

The lower bound applies to arbitrary kernels, arbitrary coefficient
recombinations, and arbitrary word controls, provided every piece
remains in the full orbit-tensor space (1). It does not assert a
conservation theorem for all groups. For example, the displayed
parameter is only \(d/4\) for \(D_8,D_{12}\), and is one for an elementary
abelian two-group. Restricted or retracted label domains require a
separate argument: the entire diagonal minor used above may no longer
be present.
