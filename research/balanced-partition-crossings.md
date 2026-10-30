# Crossings from partitions into blocks of total weight zero

## Results and scope

There is a general monomial label construction behind the
[intersection-matching crossing](intersection-matching-crossing.md).
Assign weights to physical vertices. Use blocks whose total vertex
weight is zero, and label a partition of the vertices into a fixed
number of such blocks. The kernel

\[
 K(e,f)=[e=f]+\sum_{v\in e\cap f}w_v                              \tag{1}
\]

has determinant Gram matrix exactly the identity on these labels.
This holds over every commutative ring of characteristic two.
It gives scale-one coefficient crossings without selecting or
inverting a Gram minor.

Two families retain factorial label entropy while using only a
subexponential number of physical coefficient pieces:

- Partitions into blocks of a fixed even size \(k\), with all vertex
  weights one.
- Perfect \(k\)-partite hypermatchings, with constant weights on the
  parts that sum to zero. Three parts are already possible over
  \(\mathbb F_4\), using weights \(1,\omega,\omega^2\), where
  \(\omega^2+\omega+1=0\).

The constructor is unconditional ordinary mathematics. A rank profile
that yields an exponent above one half has not been proved for either
family. The finite profiles below are stated at their actual scope.

## 1. A determinant identity for balanced blocks

Let \(V\) be a finite set with weights \(w_v\) in a commutative
characteristic-two ring \(R\). Call a nonempty block \(e\subseteq V\)
balanced if \(\sum_{v\in e}w_v=0\). Let \(\mathcal L\) be any family
of partitions of \(V\) into exactly \(q\) balanced blocks.

For \(A,B\in\mathcal L\), form the \(q\)-by-\(q\) matrix
\(K[A,B]\) using (1). Then

\[
                         \det K[A,B]=[A=B].                      \tag{2}
\]

Every block common to \(A\) and \(B\) contributes an isolated diagonal
entry one: its own total weight is zero, and it is disjoint from all
other blocks of either partition. Remove these common blocks.
If nothing remains, the determinant is one. Otherwise, the remaining
blocks still partition the same remaining vertex set, and no row
block equals a column block. The sum of a row indexed by \(e\) is

\[
 \sum_{f\in B}\sum_{v\in e\cap f}w_v=\sum_{v\in e}w_v=0.
\]

Thus the sum of the matrix columns is zero. Multilinearity of the
determinant makes its determinant zero: replace one column by the
sum of all columns; all other summands have a repeated column.
This argument works over rings with zero divisors as well.

The proof also applies to a restricted label family. Its conclusion
is exact identity, not merely nonsingularity or a dimension bound.

### A nonsymmetric rank-one alternative

There is a second dual-kernel construction that requires only a
condition on the total vertex weight. Any two distinct partitions
into \(q\) blocks have at most \(q-2\) common blocks: \(q-1\) common
blocks would force the remaining blocks to be equal as well.

Let functions \(u,v\) on the allowed modes satisfy
\(\sum_{e\in A}u(e)v(e)=0\) for every label \(A\). Then the possibly
nonsymmetric kernel

\[
                   K_1(e,f)=[e=f]+u(e)v(f)                       \tag{2a}
\]

also has \(\det K_1[A,B]=[A=B]\). Off the diagonal, its equality
matrix has at most \(q-2\) nonzero entries. Adding a rank-one matrix
cannot give full rank. Equivalently, in the determinant expansion
at most one column can come from the outer product, while at least
two are needed. This proves vanishing over arbitrary commutative
rings, without interpreting rank over a ring with zero divisors.
On the diagonal the determinant lemma gives
\[
 \det(I+u_Av_A^{\mathsf T})=1+v_A^{\mathsf T}u_A=1.
 \]

For example, set \(u(e)=1\) and
\(v(e)=\sum_{x\in e}w_x\), assuming only \(\sum_{x\in V}w_x=0\).
Each label is a partition, so the required sum is zero. Individual
blocks need not be balanced. For bipartite permutation labels this
allows \(v(i,j)=x_i+y_j\), with
\(\sum_i x_i+\sum_j y_j=0\).

This alternative fits the same direct Pfaffian endpoint argument.
Only three rank-two alternating updates are needed to pass from its
identity kernel to \(K_1\). It supplies asymmetric parameters that
the symmetric incidence family does not have. No rank improvement
is implied just by the update count.

## 2. A direct coefficient crossing

Introduce \(L\) and a variable \(Z_v\) for each physical vertex. For
every mode block \(e\), put

\[
 z_e=L\prod_{v\in e}Z_v,\qquad
 t_e=(1+L)\prod_{v\in e}(1+Z_v).                                  \tag{3}
\]

Use four groups of \(q\) particles, ordered \(U,B,A,Y\). The six
off-diagonal blocks of an alternating matrix are

| Group pair | Entry at modes \(e,f\) |
| --- | --- |
| \(U,B\) | \([e=f]t_f\) |
| \(U,A\) | \([e=f]z_et_f\) |
| \(U,Y\) | \([e=f]z_f\) |
| \(B,A\) | \(K(e,f)z_f\) |
| \(B,Y\) | \(K(e,f)z_ft_e\) |
| \(A,Y\) | \(K(e,f)t_e\) |

Let \(\mathcal F(A,B,U,Y;L,Z)\) be its Pfaffian. At
\(L=Z_v=1\), all \(t\)'s vanish and all \(z\)'s equal one. Its
remaining \(U,Y\) block is an equality matrix, and its \(B,A\) block
is \(K[B,A]\). Equation (2) therefore gives

\[
 \mathcal F(A,B,U,Y;1,\mathbf1)=[A=B][U=Y].                        \tag{4}
\]

Expanding \(\mathcal F\) into coefficients gives a scale-one
crossing. Every physical vertex appears in exactly four external
labels. Its degree is therefore at most four. The \(L\)-degree is at
most \(4q\). There are at most

\[
                         (4q+1)5^{|V|}                          \tag{5}
\]

coefficient pieces.

The endpoint proof also permits other polynomial weights
\(z_e,t_e\) satisfying \(z_e(\mathbf1)=1,t_e(\mathbf1)=0\).
Those alternative weights need their own degree and rank accounting.
Some finite profiles use that flexibility explicitly.

## 3. Even-uniform set partitions

Take \(a=|V|\), let \(k\ge2\) be a fixed even integer dividing \(a\),
and assign every vertex weight one. Every \(k\)-subset is balanced.
Labels are the partitions into \(q=a/k\) blocks of size \(k\). Their
number is

\[
 d_{a,k}=\frac{a!}{(k!)^{a/k}(a/k)!},
 \qquad
 \log d_{a,k}=(1-1/k)a\log a+O_k(a).                              \tag{6}
\]

For a requested label width \(\lambda\), choose the smallest admissible
\(a\) with \(d_{a,k}\ge2^\lambda\), and retain the first \(2^\lambda\)
labels in a fixed canonical enumeration. The ratio

\[
 \frac{d_{a+k,k}}{d_{a,k}}
     =\frac{(a+1)\cdots(a+k)}{k!(a/k+1)}
\]

has logarithm \(O_k(\log(a+2))\). Hence

\[
 a=\Theta_k(\lambda/\log(\lambda+2)),\qquad
 \log d_{a,k}=\lambda+O_k(\log(\lambda+2)).                        \tag{7}
\]

The logarithm of the piece count (5) is
\(O_k(\lambda/\log(\lambda+2))\). This is within the \(\delta=1\)
entropy allowance of
[subexponential-piece-transfer.md](subexponential-piece-transfer.md),
provided the additional joint rank hypotheses of that transfer can
be proved.

The kernel in this family is
\[
 K(e,f)=[e=f]+|e\cap f|\pmod2.
\]
For \(k=2\), this is the intersection kernel. For \(k>2\), it uses
the parity of the intersection; it is not the indicator of a
nonempty intersection.

### An incidence involution

Let \(B\) be the vertex-by-mode incidence matrix of all \(k\)-subsets,
and \(N=B^{\mathsf T}B\), over \(\mathbb F_2\). If \(a\ge2k\),
then

\[
                  N^2=0,\qquad \operatorname{rank}N=a-2.          \tag{8}
\]

The diagonal entries of \(BB^{\mathsf T}\) are
\(\binom{a-1}{k-1}\), and its off-diagonal entries are
\(\binom{a-2}{k-2}\). Their sum is
\(\binom{a-2}{k-1}=0\pmod2\), since the upper index is even and
the lower index is odd. Thus \(BB^{\mathsf T}\) is a scalar
multiple of the all-ones matrix. Also \(B^{\mathsf T}\mathbf1=0\),
because \(k\) is even. This proves \(N^2=0\).

The column space of \(B\) is the even-weight subspace \(H\) of
\(\mathbb F_2^a\). It contains every \(e_i+e_j\): subtract the
incidence vectors of two \(k\)-subsets that differ only at \(i,j\).
These vectors span \(H\). Hence \(\ker B^{\mathsf T}\) is the line
spanned by \(\mathbf1\), and this line belongs to \(H\) because \(a\)
is even. Therefore the rank of \(B^{\mathsf T}\) restricted to \(H\)
is \((a-1)-1=a-2\), proving the second claim.

For every scalar \(s\) in every commutative characteristic-two
algebra,

\[
                         (I+sN)^2=I.                            \tag{9}
\]

Its label Gram matrix remains exactly identity: the row-sum proof
of (2) works with all vertex weights equal to \(s\).

## 4. Weighted multipartite hypermatchings

Let \(V\) consist of \(k\) parts, each of size \(b\). A mode selects
one vertex from each part. A label partitions \(V\) into \(b\) such
modes. Choose constants \(w_1,\ldots,w_k\) in a fixed finite field
of characteristic two with \(\sum_jw_j=0\); every mode is balanced.

After the first part is used to index the blocks, a label is a tuple
of \(k-1\) permutations. Thus

\[
 d=(b!)^{k-1},\qquad |V|=kb,\qquad q=b,\qquad
 \log d=(k-1)b\log b+O_k(b).                                     \tag{10}
\]

The same minimal-size choice and restriction to binary labels give
the bounds (7), with constants depending on fixed \(k\). The
piece count (5) again has logarithm
\(O_k(\lambda/\log(\lambda+2))\).

For modes \(e=(e_1,\ldots,e_k)\) and \(f\),

\[
 K(e,f)=[e=f]+\sum_{j:e_j=f_j}w_j.                               \tag{11}
\]

In particular, \(k=3\) over \(\mathbb F_4\) allows
\((w_1,w_2,w_3)=(1,\omega,\omega^2)\). All three coordinates enter
the kernel, and the two freely chosen weights give a larger parameter
space than the one-scalar bipartite matching kernel. The incidence
involution (8) is not asserted for arbitrary multipartite weights.

## 5. Effectivity and an incidence reduction

For either fixed-\(k\) family, all labels can be enumerated in
\(2^{O(\lambda)}\) time. Every label uses \(O(|V|\log(|V|+2))\)
bits. The Pfaffian has \(4q=O(|V|)\) particles. Its coefficient
polynomial has at most (5) positions. A subset dynamic program for
the Pfaffian, using dense bounded-degree coefficient vectors,
computes it in \(2^{O_k(|V|)}\) field operations. Ordinary bit costs
remain \(2^{O(\lambda)}\) over a fixed field.

Generating every entry of every coefficient tensor therefore costs
\(2^{O(\lambda)}\). Exact Gaussian elimination can produce rank
factorizations within the same time bound. This supplies an
effective family; the ranks themselves still need useful bounds.

There is also a structural reduction. Write (1) as
\(K=I+B^{\mathsf T}DB\), where \(D\) contains the vertex weights.
Only the \(B,A\), \(B,Y\), and \(A,Y\) blocks differ from the
identity-kernel Pfaffian. Each correction is a sum of \(|V|\)
rank-two alternating updates. Expand the Pfaffian by the selected
update columns, as in the intersection-matching note. There are at
most \(2^{3|V|}\) channel subsets.

After fixing such a subset, its determinant splits into four local
maps, one on each external label. The remaining particles are
partial partitions, and their Pfaffian uses the identity kernel.
Refining the polynomial charges multiplies the piece count by at
most

\[
            70^{|V|}\binom{4q+4}{4}.
\]

Here each physical degree at most four is split among five factors,
and the \(L\)-degree is split in the same way. Thus the refinement
still has only \(2^{O_k(|V|)}\) pieces. Each refined tensor is the
image of one identity-kernel partial-partition coefficient tensor
under four separate linear maps, so neither flattening rank
increases.

This does not bound the rank after arbitrary channel pieces are
summed. It identifies a concrete partial-partition family whose
coefficient ranks would suffice for a constructive advance.

## 6. First finite profiles

For the smallest four-uniform example, \(a=8,k=4,q=2\), there are
35 labels. The exact calculation uses all 152 quartet orbits under
physical-vertex permutations and all \(35^4\) entries.

For the physical-central coefficient
\([L^4\prod_vZ_v^2]\), both flattening ranks are

\[
                       596=\binom{35}{2}+1.
\]

The coefficient tensor agrees entrywise with

\[
 [A=B][U=Y]+[A=U][B=Y]+[A=Y][B=U].
\]

For this finite instance the claim follows from complete exact
enumeration; no all-size identity is inferred from it.

The specialization \(L=Z_v=z\) has seven nonzero pieces:

| Degree | Row rank | Column rank |
| --- | ---: | ---: |
| 0 | 1 | 1,225 |
| 4 | 1 | 1,225 |
| 16 | 1 | 1,225 |
| 20 | 596 | 596 |
| 24 | 1,225 | 1 |
| 36 | 1,225 | 1 |
| 40 | 1,225 | 1 |

The worst rank product is 355,216. This curve does not beat the half
threshold.

For three-partite labels over \(\mathbb F_4\), use
\((w_1,w_2,w_3)=(1,\omega,\omega^2)\), and the alternative weights
\(z_e=X_{e_1}\), \(t_e=1+X_{e_1}\). They satisfy the endpoint
conditions but do not use the \(L,Z\) grading (3). The first two
complete profiles are:

| Part size \(b\) | Labels \((b!)^2\) | Nonzero pieces | Maximum rank product | Worst ranks |
| --- | ---: | ---: | ---: | --- |
| 2 | 4 | 9 | 49 | \(7,7\) |
| 3 | 36 | 99 | 675,684 | \(822,822\) |

The four-label value \(49\) is a finite gain below the half threshold.
As with the six-label seed, it is not an asymptotic amplification
theorem. The next size already has a much worse profile.

## Verification files

The generic even-uniform implementation is
[tools/even_partition_crossing.py](../tools/even_partition_crossing.py).
It provides canonical label generation, the parity kernel, direct
coefficient extraction, and an independent subset-Pfaffian polynomial
routine.

- [tools/check_balanced_partition_crossings.py](../tools/check_balanced_partition_crossings.py)
  checks the Gram identities, incidence identities, independent
  polynomial expansions, and rank-one alternative.
- [evidence/balanced-partition-crossings.json](../evidence/balanced-partition-crossings.json)
  records 4,459 binary Gram checks, 288 scalar checks over
  \(\mathbb F_4\), six incidence rank and square identities, 144
  independent polynomial comparisons, 96 comparisons with the
  matching constructor, 2,304 coefficient extractions, 2,312
  multipartite Gram checks, and 192 rank-one deformation checks.
- [tools/check_balanced_partition_profiles.py](../tools/check_balanced_partition_profiles.py)
  reproduces the complete finite profiles specified in Section 6.
- [evidence/balanced-partition-profiles.json](../evidence/balanced-partition-profiles.json)
  contains those profiles, including every nonzero coefficient and
  both flattening ranks. The bit-packed \(\mathbb F_4\) elimination
  is also compared with independent entrywise elimination on 192
  random matrices.

Run each checker with its evidence path as the argument to
\(\texttt{--output}\). The tests are finite checks of the identities
and implementations; the all-size statements are proved above.
