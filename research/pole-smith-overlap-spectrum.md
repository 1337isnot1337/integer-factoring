# Balanced pole Smith spectra are overlap-count partitions

## Exact result and scope

The [residue-flow formula](pole-smith-residue-flow.md) computes every
Smith order of an explicitly listed pole Cauchy matrix. When the full
row and transformed-column counts differ by at most one in every binary
residue class, the convex optimization has an explicit solution. A
single postorder traversal of the residue tree gives the complete Smith
histogram.

This hypothesis holds after the
[balanced interval padding](pole-cauchy-conditioning.md). The resulting
histogram has at most `k` distinct orders, even when the matrix dimension
is exponential in its interval description. Its implicit computation
still requires reciprocal-residue overlap counts. The final overlap is
the desired factor count. No improved factoring exponent follows.

## 1. Relative balance and the overlap sequence

Let `N` and every node of the distinct positive sets `A,B` be odd. Assume
`|A|=|B|=m`, and choose `h>=1` such that `|N-ab|<2^h` throughout their
rectangle. For `k>=h`, put

\[
 M_k(a,b)=\frac{2^k w}{N+2^k w-ab},\qquad w^2+w+1=0.
 \tag{1}
\]

For each odd residue `c modulo 2^ell`, define the full counts

\[
 r_{\ell,c}=\#\{a\in A:a\equiv c\pmod{2^\ell}\},\qquad
 s_{\ell,c}=\#\{b\in B:Nb^{-1}\equiv c\pmod{2^\ell}\}.
 \tag{2}
\]

The hypothesis in this note is **relative balance**:

\[
 |r_{\ell,c}-s_{\ell,c}|\le 1
 \quad\text{for every }1\le\ell\le k\text{ and }c.
 \tag{3}
\]

Define the overlap sequence

\[
 H_\ell=\sum_c\min(r_{\ell,c},s_{\ell,c}).
 \tag{4}
\]

Refining a residue class cannot increase the sum of its minimum counts,
so

\[
 m=H_1\ge H_2\ge\cdots\ge H_k=D,
 \tag{5}
\]

where `D` is the number of pairs with `ab=N`. The last identity follows
because a congruence `ab=N mod 2^k` is an integer equality under the
strict difference bound. Such pairs form a matching on the distinct
nodes. In fact `H_ell=D` for all `ell>=h`.

Let `Gamma_j` be the minimum valuation of a size-`j` minor of (1), with
`Gamma_0=0`. Then

\[
 \boxed{\displaystyle
 \Gamma_j=\sum_{\ell=1}^k\max(0,j-H_\ell).
 }
 \tag{6}
\]

Consequently, the Smith orders in increasing order are

\[
 \boxed{\displaystyle
 \sigma_j=\Gamma_j-\Gamma_{j-1}
       =\#\{\ell:H_\ell<j\},\qquad 1\le j\le m.
 }
 \tag{7}
\]

The multiplicity of order zero is `H_k`. For `1<=t<=k-1`, the
multiplicity of order `t` is

\[
 H_{k-t}-H_{k-t+1}.
 \tag{8}
\]

Thus the Smith profile is the conjugate partition of the overlap
sequence.

### Proof of the lower bound

Choose any `j` rows and `j` columns, and write `r'`, `s'` for their
counts. The minor-energy identity gives

\[
 v_2(\det M_k[R,S])
   =\frac12\sum_{\ell,c}(r'_{\ell,c}-s'_{\ell,c})^2.
 \tag{9}
\]

An integer square is at least its absolute value. At each level,

\[
 \frac12\sum_c|r'_{\ell,c}-s'_{\ell,c}|
 =j-\sum_c\min(r'_{\ell,c},s'_{\ell,c})
 \ge j-\min(j,H_\ell).
 \tag{10}
\]

Summing proves the lower bound in (6). This direction does not require
relative balance.

### A matching that attains every cardinality

Build the binary residue tree, with row leaves at `a mod 2^k` and
column leaves at `N/b mod 2^k`. Match rows and columns at the leaves,
then proceed from children to parents. Inside each subtree match all
possible opposite-color nodes before passing leftovers to its parent.

Condition (3) says that every subtree has at most one leftover. The
paths of the resulting full matching therefore use each tree edge at
most once. The root has no leftover since both sides have size `m`.
Every depth-`ell` subtree contains exactly the minimum of its row and
column counts in internal matches. Hence exactly `H_ell` of the full
matching pairs have their least common ancestor at depth at least
`ell`.

Select the `j` matching pairs with deepest least common ancestors. At
each level exactly `min(j,H_ell)` selected pairs have both ends in the
same residue class. Their paths are still edge-disjoint. Thus every
selected edge flow is zero or has absolute value one. Squared energy
equals absolute-flow energy, and its contribution at level `ell` is
exactly `j-min(j,H_ell)`. This attains (6) simultaneously for every `j`.

## 2. One postorder traversal gives the histogram

Suppress each unbranched path of the residue tree and retain its length.
The direct compressed construction has at most `4m` vertices. At a
vertex `v`, sum its row and column capacities over the subtree. The
number of matching pairs first completed at that vertex is

\[
 \min(r_v,s_v)-\sum_{u\text{ child of }v}\min(r_u,s_u).
 \tag{11}
\]

If `v` has depth `d`, each such pair contributes one Smith order `k-d`.
All counts and histogram entries follow in one postorder traversal.
The number of integer operations after constructing the tree is
`O(m+k)`, or `O(m)` if only nonzero histogram entries are emitted.
Residue generation and compression cost `m poly(b)` bit operations for
`b`-bit input nodes and `k=O(b)`.

For a very large requested `k`, build only to depth `h`. The orders
zero remain zero, and every positive order increases by `k-h`, by the
[uniform stabilization theorem](pole-smith-stabilization.md). The code
returns `H_1,...,H_h`; the omitted overlap entries are all `D`. Its
dependence on a larger `k` is polynomial in `log(k+1)`.

The implementation is
[pole_smith_balanced.py](../tools/pole_smith_balanced.py). It rejects a
residue imbalance larger than one. This check is necessary. For

\[
 N=101,\quad A=(1,3,5,7),\quad B=(1,9,17,25),\quad k=7,
\]

the true Smith orders are `(2,5,7,11)`. The overlap sequence is
`(4,2,1,1,1,0,0)`, which would incorrectly predict `(2,5,6,6)` without
hypothesis (3).

## 3. Power-of-two padding gives direct diagonal equivalence

The interval-padding construction can pad both domains to

\[
 M=2^s\ge m,
 \qquad M<2m\text{ unless }m\text{ is already a power of two}.
 \tag{12}
\]

Use modulus `Q=2M`. For each original odd progression, retain its
original nodes. Continue its missing tail to length `M`, then shift
that entire tail by a sufficiently large multiple of `Q` so that every
added node exceeds `N`. Each padded domain is represented by at most
two odd progressions. Each contains every odd residue modulo `Q`
exactly once. The construction retains all original factor pairs and
adds none.

Pair each padded row `a` with the unique padded column `b(a)` satisfying

\[
 a b(a)\equiv N\pmod{2M}.
 \tag{13}
\]

After this column permutation, every diagonal entry of (1) has valuation
at most `k-s-1`. Every off-diagonal entry has valuation at least `k-s`.
Let `T` be the diagonal matrix of the diagonal entries. Then

\[
 M_k=T(I+E),\qquad E\equiv0\pmod2.
 \tag{14}
\]

Here take `k>=max(h,s+1)`, enlarging the displacement if needed. The
matrix `I+E` is unimodular. Thus the Smith multiset is exactly

\[
 \left\{\,k-\min\bigl(k,v_2(N-a b(a))\bigr):a\in A\,\right\},
 \tag{15}
\]

with `v_2(0)=+infinity`. This also follows from (7), since `H_ell=M`
through `ell=s+1`. The improved conditioning bound is

\[
 2^{k-s-1}M_k^{-1}\text{ is integral}.
 \tag{16}
\]

The pairing in (13) has a polynomial-bit-cost query for a single row:
invert its residue modulo `2M`, then locate that residue in one of the
two stored column progressions. Computing or counting all zero entries
of (15) remains an implicit modular-inverse problem. Formula (14) is
an exact normal form, not a polynomial-time global evaluator.

## 4. What an implicit evaluator must still compute

Suppose each full domain is balanced at every level, meaning its counts
are either `u` or `u+1` when

\[
 m=u\,2^{\ell-1}+r,\qquad 0\le r<2^{\ell-1}.
\]

Let `E_A` be the `r` odd classes with the larger row count, and `E_B`
the `r` classes with the larger transformed-column count. Then

\[
 H_\ell=m-r+|E_A\cap E_B|.
 \tag{17}
\]

For ordinary odd intervals, the exceptional classes before transforming
the columns are cyclic intervals of odd residues. The unresolved part
of (17) is therefore a modular-inverse overlap count. Complete cycles
contribute the known term `m-r`. When `m=2^s`, all levels through
`s+1` are complete; higher levels retain the entire domain as their
remainder. In particular, `H_k=D` is the original factor count.

The preceding results remove the matrix optimization and large output
size from the balanced problem. They leave a precise scalar arithmetic
task: evaluate (17), or the count of true zeros in (15), on binary
interval descriptions at an improved uniform cost.

## 5. Exact computational checks

[check_pole_smith_balanced.py](../tools/check_pole_smith_balanced.py)
compares the histogram with the general convex-flow profile on 1,000
balanced padded matrices. The first 100 also use independent exact
Gaussian elimination over `Q(w)`. Every case checks (7), and the
unbalanced control above must be rejected.

There are also 200 power-of-two padding controls checking the direct
diagonal formula (15); 100 use independent exact Gaussian elimination.
One displacement of `10^80` checks that the implementation stores the
large positive orders without constructing that many residue levels.

A separate 201-bit control has dimension `65,536`. The postorder
algorithm visits `262,142` compressed vertices and returns 14 nonzero
histogram entries. Its complete expanded profile has the same SHA-256
digest as the earlier convex-flow calculation:

```text
5119a2e4c7efcd01771df1906c0608e5a40d88465147b05081747ee83d4a6b5f
```

The numerical node list is explicit in that control. It verifies the
algorithmic simplification and makes no claim about an implicit
exponential-dimension computation. Results are recorded in
[pole-smith-balanced.json](../evidence/pole-smith-balanced.json).
