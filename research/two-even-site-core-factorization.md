# All-size assembly from four-site cores

The two-even-site construction has an exact decomposition into cores on
at most four physical sites and elementary factors on all remaining
sites. This gives an all-size polynomial description with at most
\(2P-2\) coefficient directions. It also gives a practical method to
assemble the full matrices before computing rank.

The implementation reaches **six sites and 240 logical labels**. Its
complete ten-dimensional binary coefficient space fails the finite-seed
criterion: all 1,023 nonzero combinations were tested, and the eight
affordable combinations span a five-dimensional space that misses the
crossing. A separate
[constant-compression certificate](two-even-site-P6-character-barrier.md)
extends the exclusion for this seed to every characteristic-two field
and every tensor power. That extension uses explicit constant maps,
not an extrapolation from the finite scan.

The construction itself and the four-site factorization hold at every
size. This advances the available exact evaluator beyond the fixed
24-, 26-, and 80-label spaces covered by earlier all-power theorems.

## 1. Logical labels and the syndrome blocks

Use \(P\) physical mode pairs. A label consists of a two-element set
\(e=\{a,b\}\), with \(a<b\), and one single-mode orientation at every
site outside \(e\). Its input is

\[
I_{e,\omega}=m_\omega(D_a+D_b),
\]

where \(D_i\) is the full pair at site \(i\). Its output selects
\(m_\omega D_b\). The label count is

\[
d_P=\binom P2 2^{P-2}.
\tag{1}
\]

All vectors have physical grade \(P\), each input has two terms, and
the selected outputs form a left inverse. These are homogeneous
product-image codes of the kind used in the
[Ward theorem](noncommutative-crossing-radical-ward.md).

For a fixed tensor entry, let \(e_x,e_y,e_u,e_v\) be the four even-site
sets. A local Pfaffian entry has even total occupation. Therefore a
nonzero entry requires

\[
e_x\mathbin\triangle e_y=e_v\mathbin\triangle e_u.
\tag{2}
\]

The same condition gives the column-cut blocks after regrouping the
four labels. The possible block syndromes have weights zero, two,
and four.

For a syndrome \(S\), let \(m_S\) be the number of ordered pairs of
two-element sets with symmetric difference \(S\). Then

\[
m_S=
\begin{cases}
\binom P2,&S=\varnothing,\\
2(P-2),&|S|=2,\\
6,&|S|=4.
\end{cases}
\]

Each ordered pair has \(4^{P-2}\) row orientations, so the cut block
has dimension \(m_S4^{P-2}\).

## 2. At most four exceptional sites

Let
\[
U=e_x\cup e_y\cup e_u\cup e_v.
\]
Equation (2) says that every site in this union occurs an even number
of times among the four two-element sets. There are eight incidences
in total, and every site occurs at least twice. Hence

\[
|U|\le4.
\tag{3}
\]

Outside \(U\), all four logical occupations are single modes. All
two-term input expansions occur inside \(U\). Thus, after fixing
the four even-site sets, the corresponding orientation submatrix is
an exact tensor product of:

- a two-, three-, or four-site two-even-site core, restricted to those
  even-site sets; and
- the ordinary two-label single-mode local tensor at every site
  outside \(U\).

There are only independent row and column coordinate permutations
between this tensor-product ordering and the global ordering. The
upper-site output selector is preserved when the exceptional sites
are placed in their original increasing order.

This factorization applies to the fixed-edge orientation submatrices.
In general, different submatrices in one syndrome block have different
exceptional sets. Their polynomial entries must be assembled into
the full block before rank is taken. The implementation does this;
it does not replace the full rank by a sum of submatrix ranks.

For weight-four syndromes, the exceptional set is the syndrome itself
for every submatrix. In that case the entire block is the
96-dimensional four-site block tensored with the free-site factors.

## 3. The free-site factors

Put

\[
u=q+r,\qquad v=p+q,\qquad s=qr,
\qquad t=v/u,\qquad z=s/u^2.
\]

These ratios describe a homogeneous polynomial identity; one may
clear powers of \(u\) and specialize even when \(u=0\).

In the row cut, the local single-mode factor is

\[
H_r=uI_4+vF,
\qquad
F=\begin{pmatrix}
1&0&0&0\\0&0&1&0\\0&1&0&0\\0&0&0&1
\end{pmatrix}.
\tag{4}
\]

In the column cut, it is

\[
H_c=uI_4+vN',
\qquad
N'=\begin{pmatrix}
1&0&0&1\\0&0&0&0\\0&0&0&0\\1&0&0&1
\end{pmatrix}.
\tag{5}
\]

Thus the normalized factors are \(I+tF\) and \(I+tN'\).
Both \(F+I\) and \(N'\) have rank one and square zero.

The two-site core is simply \(u^2\). The three-site normalized core
has directions \(1,t,t^2,z\). The four-site normalized core has
directions \(1,t,t^2,t^3,z,zt\), as independently verified in the
[24-label certificate](symplectic-seed-character-barrier.md).

Multiplying a core on \(h\le4\) sites by \(P-h\) free factors gives
the following all-size identity:

\[
\boxed{
G_P(p,q,r)=
\sum_{i=0}^{P-1}A_i u^{P-i}v^i
+\sum_{i=0}^{P-3}B_i u^{P-i-2}v^i s.
}
\tag{6}
\]

The matrices \(A_i,B_i\) are constant tensors over \(\mathbb F_2\).
Their span has dimension at most \(P+(P-2)=2P-2\). The first
coefficient \(A_0\) is the crossing tensor. Equation (6) includes all
scalar polynomial coefficients of the seed; it is not a selected
one-parameter curve. It also displays the common factor \(q+r\).

The implementation uses the ordered normalized directions

\[
1,t,\ldots,t^{P-1},z,zt,\ldots,zt^{P-3}.
\tag{7}
\]

For comparison with the raw literal basis, the \(A_i\) are the
coefficients of \(p^ir^{P-i}\). If \(C_i\) is the raw coefficient
of \(p^iqr^{P-i-1}\), then

\[
B_i=C_i+((P-i)\bmod2)A_i+((i+1)\bmod2)A_{i+1}.
\tag{8}
\]

This relation is checked during the exhaustive five-site comparison.

## 4. Exact six-site result

For \(P=6\), equation (1) gives \(d=240\). There are 31 syndrome
blocks:

| Syndrome weight | Number of blocks | Block dimension |
|---:|---:|---:|
| 0 | 1 | 3,840 |
| 2 | 15 | 2,048 |
| 4 | 15 | 1,536 |

Their dimensions sum to \(240^2=57600\). Direct literal evaluation
of every block entry would require 113,049,600 contracted entries.
The small-core construction avoids that contraction while producing
the exact same polynomial matrices.

Ten independent coordinate witnesses certify that all ten directions
in (7) are independent for this seed. The complete binary scan gives:

| Labels | Coefficient dimension | Nonzero combinations | Integer gate | Affordable combinations | Affordable span dimension | Crossing in span |
|---:|---:|---:|---:|---:|---:|:---:|
| 240 | 10 | 1,023 | 2,325 | 8 | 5 | No |

The gate is \(g=\lceil240^{99/70}\rceil=2325\), with
\(99^2>2\cdot70^2\), hence \(g>240^{\sqrt2}\).
A combination is affordable when its rank product is less than
\(g^2\). A crossing decomposition of total cost less than \(g\)
could use only affordable pieces. Since their span misses the
crossing, no such decomposition exists over \(\mathbb F_2\).

In the normalized coordinate order (7), all eight affordable masks
also annihilate the functional selecting coordinates zero and four.
Its polynomial is
\((q+r)^2(p+r)^4\). This finite observation is not an all-field rank
certificate, and is not extrapolated to tensor powers.

## 5. Verification

The five-site implementation was compared with the independently
generated 80-label literal basis. Every one of the 102,400 normalized
bit rows, covering all eight generators and both cuts, matched.
This checks the core selection, orientation permutations, free-site
factors, and normalization together.

For six sites, 346 additional complete polynomial entries were
checked directly from the local Pfaffian kernel and the four input
support contributions. These controls cover exceptional unions of
sizes two, three, and four on both cuts. They compare formal binary
polynomials in \(p,q\) at \(r=1\), which suffices for the homogeneous
degree-six identity.

All eight affordable rank pairs and twelve additional rejection
controls were independently recomputed by scalar elimination. The
complete native log covers all 1,023 nonzero binary vectors. Source
and input hashes bind the certificate, matrix export, and scan.

The files are:

- [assemble_two_even_site_cores.py](../tools/assemble_two_even_site_cores.py);
- [two-even-site-core-P5.json](../evidence/two-even-site-core-P5.json);
- [two-even-site-core-P6.json](../evidence/two-even-site-core-P6.json);
- [two-even-site-core-P6-F2.jsonl](../evidence/two-even-site-core-P6-F2.jsonl).

Reproduce the six-site input and exact checks:

```sh
python tools/assemble_two_even_site_cores.py --pairs 6 \
  --native /tmp/two-even-P6.in \
  --output /tmp/two-even-P6-assembly.json
g++ -O3 -DNDEBUG -std=c++20 tools/screen_weighted_pfaffian_span.cpp \
  -o /tmp/two-even-screen
/tmp/two-even-screen --field2 < /tmp/two-even-P6.in \
  > /tmp/two-even-P6-screen.jsonl
python tools/assemble_two_even_site_cores.py --pairs 6 \
  --screen /tmp/two-even-P6-screen.jsonl \
  --output /tmp/two-even-P6-audit.json
```

## 6. What the local factors do and do not prove

The exact assembly provides a route to larger seeds, but the full
matrices still grow rapidly. A useful next step is to combine the
rank-one nilpotent decompositions in (4)--(5) with weighted small
modules before assembling large blocks. Weight-two and weight-zero
blocks require retaining the marked even-site positions as well.

The simplest tensorized common-functional obstruction does not by
itself close the family at every size. For the local functional
\(\chi(H)=p+r\), the row cut has a triangular block of size four,
but the column compression has size at most one: the coefficient of
\(p\) in the column cut is \(N'\), of rank one. That size is attained.
The local rank-product contribution is therefore four, below the
factor eight required by \(d_P^3\).

More generally, the local tensors satisfy \(E=P+Q\). A normalized
functional has \(\chi(P)+\chi(Q)=1\). If both values are nonzero,
the rank-one coefficient in each opposite cut bounds both scalar
compressions by one. The two binary extreme choices give at most
four when triangular certificates are allowed. Thus polynomially
many copies of a fixed core cannot make this particular local
certificate into an all-size \(d_P^{3/2}\) obstruction.

Shared-parameter coefficient collisions and the assembled interactions
between marked-site submatrices remain available. The finite binary
scan alone does not decide those effects over extension fields or
across tensor powers. The separate compression certificate now
settles this specific 240-label span; larger two-even-site
constructions remain outside that theorem.
