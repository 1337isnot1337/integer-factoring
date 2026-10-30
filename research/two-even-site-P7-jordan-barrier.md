# A 672-label all-field certificate from small Jordan blocks

The seven-site two-even-site seed is now excluded over every
characteristic-two field and at every tensor power of its full
coefficient space. The proof uses selected syndrome restrictions and
small Jordan modules. It never constructs a full seven-site matrix.

The final row rank certificate has size **108,315** and the column
compression has size **1,150**. Their product is

\[
108315\cdot1150=124562250.
\]

For every crossing decomposition in the \(n\)-fold coefficient span,

\[
\boxed{
S\ge124562250^{n/2}>672^{10n/7}>672^{\sqrt2 n}.
}
\tag{1}
\]

The exponent comparison is exact:
\(124562250^7>672^{20}\) and \(10^2>2\cdot7^2\).
This extends the fixed-seed all-field frontier from 240 to 672 labels.
It does not improve the factoring exponent or rule out changed codes.

## 1. Full coefficient space and the common functional

The code has
\(\binom72 2^5=672\) logical labels. Its inputs and upper-site outputs
are those in the
[all-size core factorization](two-even-site-core-factorization.md).
That theorem gives the twelve formal normalized directions

\[
1,t,\ldots,t^6,z,zt,\ldots,zt^4,
\qquad t=(p+q)/(q+r),\quad z=qr/(q+r)^2.
\]

They contain the complete scalar polynomial coefficient span. Define
the common functional by the coefficients of \((1+t)^5\):

\[
\chi=(1,1,0,0,1,1,0,0,0,0,0,0).
\tag{2}
\]

On the original homogeneous polynomial family,

\[
\chi(G_7)=(q+r)^2(p+r)^5.
\]

The crossing coefficient has functional value one. The matrix
certificates below also ensure that the functional respects every
linear relation among the formal coefficient directions.

## 2. Exact free-site Jordan decompositions

In the two cuts, a normalized free single-mode site has matrices

\[
I_4+tF,\qquad I_4+tN',
\]

where \(F+I_4\) and \(N'\) are rank-one square-zero matrices.
The constant basis columns are respectively

\[
(e_0,e_3,e_1+e_2,e_1),\qquad
(e_1,e_2,e_0+e_3,e_0).
\]

In these bases, the row factor is the direct sum of two copies of
the scalar \(1+t\) and the two-dimensional block
\((1+t)I_2+tJ\). The column factor is the direct sum of two scalar
ones and \(I_2+tJ\), where

\[
J=\begin{pmatrix}0&1\\0&0\end{pmatrix}.
\]

For \(m\) free sites, the summand with \(j\) Jordan factors therefore
has multiplicity

\[
\binom mj2^{m-j}
\tag{3}
\]

and Jordan dimension \(2^j\). The row cut also has the scalar factor
\((1+t)^{m-j}\). The column scalar factor is one. These are constant
direct-sum decompositions, so the integer multiplicities count actual
copies of identical modules. They are not sampling weights.

All factors use the same parameter \(t\). The module builder multiplies
the polynomial matrices, merges equal degrees, and then constructs
the constant rank certificates. This retains the relevant coefficient
collisions between factors.

## 3. Which parts of the full matrices are used

The complete seven-site cuts have 57 syndrome blocks. Their sizes
are 21,504 for syndrome zero, 10,240 for each of the 21 weight-two
syndromes, and 6,144 for each of the 35 weight-four syndromes.

The construction uses the following pieces.

| Physical syndrome | Number used | Exact object retained | Dimension before Jordan decomposition |
|---|---:|---|---:|
| Weight four | 35 | Entire block: four-site core times three free sites | 6,144 |
| Weight two | 21 | Principal restriction to three chosen shared anchors: five-site core times two free sites | 6,144 |
| Zero | 1 | Principal restriction to even edges within five chosen sites: five-site core times two free sites | 10,240 |

The first row is an equality for the whole syndrome block. The other
two rows are principal submatrices; omitted coordinates can only
increase the full ranks.

For a weight-two syndrome \(S=\{a,b\}\), a row edge pair is
\((\{a,c\},\{b,c\})\), in either order. Choose three allowed anchors
\(c\) and retain only those coordinates in both the rows and columns
of the selected cut. Let \(U\) be \(S\) together with those anchors.
The resulting principal restriction is exactly the weight-two block
of the five-site core on \(U\), tensored with the free factors outside
\(U\).

The row and column cuts may choose different anchor sets. Each choice
is a constant restriction of its own flattening, and every choice
carries the same coefficient functional (2). No equality between the
two selected submatrices is required. The certificate lists both
anchor choices for each physical syndrome and verifies their relative
core-syndrome indices.

For syndrome zero, the retained five sites are \(\{0,1,2,3,4\}\).
The retained edge pairs are the equal pairs whose even edge lies
inside that set. This is also a principal restriction.

Different physical syndromes are genuine direct summands of the
original coefficient matrices. Their certificates can therefore be
combined by direct sum. Within a weight-two syndrome, only one
chosen principal restriction per cut is counted.

## 4. Constant scalar maps and a stronger triangular row block

The small modules use the
[direct subspace criterion](common-functional-rank-compressions.md).
For each module, selected coordinate vectors span \(V\) with

\[
V\cap\sum_i(M_i-\chi_iI)V=\{0\}.
\]

The independent checker reconstructs the left maps and tests
\(WM_iV=\chi_iI\) entry by entry. Combining these scalar maps gives
the intermediate totals below.

| Stage | Row scalar compression | Column scalar compression | Rank product |
|---|---:|---:|---:|
| Weight-four blocks only | 50,785 | 560 | 28,439,600 |
| Add two-anchor weight-two restrictions | 72,483 | 990 | 71,758,170 |
| Use three anchors and add the weight-zero restriction | 87,420 | 1,150 | 100,533,000 |

The first two stages do not establish the seed exclusion. The third
already passes the exact \(99/70\) sufficient exponent. These partial
certificates are retained so the contribution of each extension is
explicit and reproducible.

There is a stronger existing row certificate for the weight-four
core. Its adapted row representation contains a principal upper
triangular block of size 32 with diagonal polynomial \((1+t)^2\).
Each free row factor is upper triangular with diagonal \(1+t\) at
all four positions. Therefore the three free sites give a row
triangle of size \(32\cdot4^3=2048\) for each weight-four syndrome,
with diagonal \((1+t)^5\).

Across 35 syndromes this contributes 71,680 row dimensions, replacing
the scalar subtotal 50,785. The resulting final contributions are:

| Source | Row rank certificate | Column scalar compression |
|---|---:|---:|
| Weight-four blocks | 71,680 | 560 |
| Three-anchor weight-two restrictions | 33,413 | 584 |
| Weight-zero principal restriction | 3,222 | 6 |
| **Total** | **108,315** | **1,150** |

The row map is now a direct sum of scalar identities and upper
triangular matrices with the same functional on the diagonal. It is
not asserted to be a scalar compression of size 108,315.
For every coefficient-space tensor \(T\), \(\chi(T)\ne0\) makes this
row matrix nonsingular and makes the column compression nonsingular.
Thus \(r(T)\ge108315\) and \(c(T)\ge1150\).

## 5. All fields, all powers, and mixed products

Every basis change, restriction, scalar compression, and triangle
above is defined over \(\mathbb F_2\). These finite matrix identities
persist over every characteristic-two field.

Tensor products of the row maps remain upper triangular, with
\(\chi^{\otimes n}(T)\) on every diagonal position. Tensor products
of the column maps remain scalar identity compressions. In a
decomposition of the crossing, the functional values sum to one,
so at least one piece has nonzero value. Its contribution alone
gives (1), even when the pieces are arbitrary elements of the whole
tensor-product coefficient span.

Together with the 24-, 26-, 80-, and 240-label certificates, this
gives a mixed-product result. Let \(D\) be the product of the label
dimensions in any nonempty mixture of the five certified spaces.
Every crossing decomposition in the full mixed coefficient span
satisfies

\[
S>D^{10/7}>D^{\sqrt2}.
\]

This includes independent local parameter curves and arbitrary
coefficient recombinations. It does not apply automatically to
changed output sections, extra directions, or the eight-site seed.
It is not a lower bound for integer factoring.

## 6. Verification and reproduction

The largest matrix constructed by this checker has dimension 2,560.
The five-site core can be regenerated from literal contractions or
reused by exact digest. The four-site core is the previously audited
24-label basis. The checker verifies:

- the two constant local Jordan similarities;
- all module dimensions and integer multiplicities;
- all physical anchor restrictions and their relative core indices;
- 342 full seven-site polynomial entries, checked directly against
  the original local Pfaffian kernel;
- 23,413,104 scalar-compression entries;
- a separate reconstruction of the four-site row algebra and its
  32-dimensional triangular block;
- the exact integer exponent margin.

The literal seven-site controls include both cuts and all three
kinds of retained blocks. They compare full homogeneous polynomial
entries, not values at a few field points.

Files:

- [check_two_even_site_jordan_compressions.py](../tools/check_two_even_site_jordan_compressions.py);
- [two-even-site-P7-jordan-compression-certificate.json](../evidence/two-even-site-P7-jordan-compression-certificate.json);
- [two-even-site-P7-jordan-compression.json](../evidence/two-even-site-P7-jordan-compression.json).

Independent replay, including regeneration of the five-site literal core:

```sh
python tools/check_two_even_site_jordan_compressions.py \
  --output /tmp/two-even-site-P7-jordan-compression.json
```

If the already checked five-site basis is available, pass
`--basis5 /tmp/omega_two_even_P5/basis.json`. Its coefficient digest is
verified against the earlier 80-label certificate before reuse.
