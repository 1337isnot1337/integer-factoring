# Eight-site barrier from bounded Jordan modules

## Result and scope

The two-even-site construction with eight physical sites has

\[
d=\binom82 2^6=1792
\]

logical labels. Its full coefficient span has a common functional that gives a
row rank certificate of dimension **769,060** and a column scalar compression
of dimension **2,106**. Consequently, every decomposition of a tensor power of
its crossing tensor into elements of the full tensor coefficient span satisfies

\[
S\ge (769060\cdot2106)^{n/2}
  =1619640360^{n/2}
  >1792^{99n/70}
  >1792^{\sqrt2 n}.
\]

These statements hold over every field of characteristic two. They allow
arbitrary coefficient-space recombinations at every tensor power. They concern
this fixed seed family; they do not give a lower bound for general integer
factoring algorithms.

The comparison is exact:

\[
1619640360^{35}>1792^{99},\qquad
99^2=9801>9800=2\cdot70^2.
\]

The computed cost exponent is approximately 1.41538; the proof uses only the
integer comparisons above. No full eight-site matrix is constructed. The largest
matrix used by the checker has dimension 5,120.

## Coefficient functional and local factors

The [all-size core factorization](two-even-site-core-factorization.md) gives the
normalized polynomial family

\[
G_P(t,z)=\sum_{i=0}^{P-1}A_i t^i+
             z\sum_{i=0}^{P-3}B_i t^i,
\qquad A_0=E.
\]

Here \(u=q+r\), \(t=(p+q)/u\), and \(z=qr/u^2\); the original family is
homogeneous of degree \(P\). These variables describe polynomial coefficient
identities. They do not require division when a coefficient matrix is evaluated.

For \(P=8\), take the functional whose normalized scalar polynomial is

\[
\chi_8(G_8)=(1+t)^6=1+t^2+t^4+t^6.
\]

Thus the functional takes values

\[
(1,0,1,0,1,0,1,0;0,0,0,0,0,0)
\]

on the fourteen displayed directions. Its original homogeneous scalar is
\((q+r)^2(p+r)^6\), and \(\chi_8(E)=1\). The compression identities below imply
that this functional respects any relations among the displayed directions; an
independence assertion about all fourteen matrices is unnecessary.

At a free site, the normalized row and column factors are

\[
H_r(t)=I_4+tF,\qquad H_c(t)=I_4+tN,
\]

where \(F+I_4\) and \(N\) each have rank one and square zero. The constant
similarities checked in the
[seven-site package](two-even-site-P7-jordan-barrier.md) put these in the forms

\[
H_r(t)\sim (1+t)\oplus(1+t)
             \oplus\big((1+t)I_2+tJ\big),
\]

\[
H_c(t)\sim1\oplus1\oplus(I_2+tJ),\qquad
J=\begin{pmatrix}0&1\\0&0\end{pmatrix}.
\]

In particular, the entire row factor is upper triangular with all four diagonal
entries equal to \(1+t\). For \(m\) free sites, the summands with \(j\) Jordan
factors occur with multiplicity

\[
\binom mj 2^{m-j}.
\]

Polynomial coefficients with equal powers are combined before any compression
is tested.

## The three retained syndrome classes

Pair coordinates split both cuts into syndrome blocks. At eight sites these
syndromes have weights zero, two, or four. Different syndrome blocks are disjoint,
so their rank certificates can be added. Within a block, the proof retains only
the subspace described below.

| Syndrome class | Number of blocks | Retained family | Row dimension | Column dimension |
| --- | ---: | --- | ---: | ---: |
| Weight four | 70 | Entire four-site core with four free sites | 573,440 | 1,120 |
| Weight two | 28 | Five-site principal restriction with three free sites | 182,732 | 980 |
| Weight zero | 1 | Five-site principal restriction with three free sites | 12,888 | 6 |
| Total | 99 | Direct sum of the retained blocks | **769,060** | **2,106** |

### Weight four: entire blocks

A syndrome supported on four sites forces the exceptional edges to lie on those
sites. Its entire block is the 96-dimensional four-site core tensored with four
free local factors. The independently checked core row certificate is upper
triangular on 32 coordinates, all with diagonal functional \((1+t)^2\).
Tensoring the four full row factors gives

\[
32\cdot4^4=8192
\]

row coordinates with diagonal \((1+t)^6\) per block.

For the column cut, the Jordan decomposition reduces the calculation to modules
of dimensions \(96\cdot2^j\), for \(0\le j\le4\). The stored witness has
compression dimension 16 in the \(j=4\) module. The proof discards the other
modules. This gives 16 column coordinates per physical weight-four block.

### Weight two: principal restrictions

For a syndrome \(S=\{a,b\}\), choose three other physical sites as anchors. Keep
only pairs of exceptional edges

\[
\{a,c\},\{b,c\},\qquad c\text{ one of the three anchors},
\]

in either order. Together with the local single occupations, these coordinates
form the complete weight-two block of a five-site core, tensored with the three
remaining free factors. This is a **principal restriction** of the full
weight-two block. It is not asserted to equal the whole block.

The certificate chooses the three anchors separately for the row and column
cuts. The relative syndrome inside the ordered five-site core determines the
small core matrices. Different anchor choices are compatible because every
retained certificate has the same functional \(\chi_8\).

For the row cut, use the verified seven-site certificate for this five-site core
with two free factors. Each of its scalar compressions has functional
\((1+t)^5\). Tensor its two maps with \(I_4\) on the third free factor. The
result is upper triangular, with all diagonal entries equal to
\((1+t)^5(1+t)=(1+t)^6\). Its dimension is four times the old dimension.
This gives the new row bound without searching larger row matrices.

For the column cut, the checker constructs fresh eight-site modules from the
five-site core with \(j=0,1,2,3\) Jordan factors. Their dimensions are at most
\(384\cdot8=3072\). The nonzero stored column compressions are:

| Relative syndrome | Compression for \(j=2\) | Compression for \(j=3\) | Weighted column dimension |
| ---: | ---: | ---: | ---: |
| 3 | 0 | 8 | 8 |
| 5 | 4 | 20 | 44 |
| 6 | 0 | 8 | 8 |
| 9 | 4 | 22 | 46 |
| 10 | 4 | 20 | 44 |
| 12 | 0 | 8 | 8 |
| 17 | 0 | 8 | 8 |
| 18 | 4 | 22 | 46 |
| 20 | 4 | 22 | 46 |
| 24 | 0 | 8 | 8 |

The multiplicities for \(j=2,3\) are six and one. The stored anchor list selects
the best available relative syndrome for each physical syndrome and cut. Its
weighted sums are the weight-two totals in the first table.

### Weight zero: principal restriction

Retain the coordinates whose equal exceptional edges lie inside physical sites
\(0,1,2,3,4\). This gives the five-site zero-syndrome core of dimension 640,
with three free factors. As above, lifting the seven-site row certificate gives
\(4\cdot3222=12888\) row coordinates. The column certificate retains six
coordinates from the \(j=3\) module of dimension \(640\cdot8=5120\).

## Independent witness audit

The new checker reconstructs the four- and five-site matrices from the published
core sources. By default it regenerates the literal five-site basis. An optional
basis file is accepted only after its full coefficient digest matches the
published five-site certificate.

The discovery routine proposes coordinate subspaces using tagged elimination.
The audit uses the separate low-pivot dual solver from the
[common-functional compression package](common-functional-rank-compressions.md).
For each retained scalar module it reconstructs a left map \(W\), then checks
all entries of

\[
WM_iV=\chi(M_i)I.
\]

The audit checks 22,307,228 scalar compression entries. It independently
rechecks the seven-site row module witnesses used in the local lift. It also reconstructs the four-site row algebra and checks its 32-coordinate
triangular character block. The local constant similarities are checked directly.

The physical embedding audit compares **594 complete eight-site polynomial
entries** with literal contractions: 420 entries in entire weight-four blocks,
168 entries in weight-two principal restrictions, and six entries in the
zero-syndrome principal restrictions. These controls accompany the all-size
factorization proof; they are not a substitute for its universal argument.

The row result is an upper-triangular rank certificate with common diagonal. It
is not claimed to be a scalar compression of dimension 769,060. The column result
is a scalar compression of dimension 2,106.

## Tensor powers and mixed products

For every element \(X\) of this coefficient span, the row compression is upper
triangular with all diagonal entries \(\chi_8(X)\), and the column compression is
\(\chi_8(X)I_{2106}\). Tensoring these maps preserves both statements, including
for arbitrary linear combinations of tensor products.

In any decomposition of \(E^{\otimes n}\), the sum of the values under
\(\chi_8^{\otimes n}\) is one. At least one summand has nonzero value. Its row
and column ranks are at least \(769060^n\) and \(2106^n\), respectively. Its
single contribution already gives the stated lower bound on \(S\).

The same argument combines this seed with the previously certified label sizes
24, 26, 80, 240, and 672. If a mixed tensor product has total label count

\[
D=24^{n_{24}}26^{n_{26}}80^{n_{80}}240^{n_{240}}
  672^{n_{672}}1792^{n_{1792}},
\]

then every decomposition inside the corresponding tensor coefficient span obeys

\[
S>D^{99/70}>D^{\sqrt2}
\]

whenever at least one exponent is positive. The earlier seeds have stronger
certified exponents, so the eight-site comparison supplies the shared exponent.

## Reproduce

```bash
python tools/check_two_even_site_P8_compressions.py \
  --output /tmp/two-even-site-P8-compression.json
```

To reproduce the discovery choices as well as the independent audit:

```bash
python tools/check_two_even_site_P8_compressions.py \
  --discover \
  --output /tmp/two-even-site-P8-compression.json
```

An existing certified literal five-site basis can be passed with
`--basis5 /path/to/basis.json`. This changes neither the matrices nor the report.
The new package consists of this note, the checker, the exact compression
certificate, and its audit report. Its dependencies are the already published
four-, five-, and seven-site packages.
