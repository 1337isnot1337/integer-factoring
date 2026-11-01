# Mixed products of the certified coefficient spaces

## Corollary

Choose any nonempty finite list of local coefficient spaces from the table below.
Different tensor positions may use different listed spaces. Over any field of
characteristic two, let \(D\) be the product of their label counts and let
\(\mathcal V\) be their tensor-product coefficient space. Every decomposition
of the resulting crossing tensor into elements of \(\mathcal V\) satisfies

\[
\boxed{\quad S>D^{17/12}>D^{\sqrt2}.\quad}
\]

The decomposition may use arbitrary linear combinations across the entire tensor
coefficient space. In particular, the pieces need not be elementary tensors or
products of pieces from local decompositions.

| Local space | Labels \(d\) | Row rank base \(r\) | Column rank base \(c\) | Product \(rc\) |
| --- | ---: | ---: | ---: | ---: |
| Symplectic 24-label span | 24 | 384 | 43 | 16,512 |
| Symplectic 26-label span | 26 | 436 | 43 | 18,748 |
| Two-even-site, five sites | 80 | 2,072 | 286 | 592,592 |
| Two-even-site, six sites | 240 | 15,688 | 718 | 11,263,984 |
| Two-even-site, seven sites | 672 | 108,315 | 1,150 | 124,562,250 |
| Two-even-site, eight sites, strengthened row | 1,792 | 1,225,260 | 2,106 | 2,580,397,560 |
| Two-even-site, nine sites | 4,608 | 7,522,480 | 3,332 | 25,064,903,360 |
| Cyclic cubic outputs \(C=\lambda S\) | 7 | 14 | 24 | 336 |
| Jordan cubic outputs \(C=\lambda J\) | 7 | 14 | 28 | 392 |
| Cubic outputs with one fixed image \(C=uv^T\) | 7 | 9 | 28 | 252 |
| Mixtures of coordinate-image cubic outputs | 7 | 14 | 28 | 392 |

The cyclic column base 24 comes from a seven-point evaluation-grid argument.
It is not asserted to be a 24-dimensional scalar compression. The other bases
come from scalar compressions or simultaneous upper-triangular compressions with
a common diagonal functional.

For every table entry the exact inequality

\[
(rc)^6>d^{17}
\]

holds. For the cubic entries these are the three comparisons
\(252^6>7^{17}\), \(336^6>7^{17}\), and \(392^6>7^{17}\).
The final comparison uses \(17^2=289>288=2\cdot12^2\).

## Scope of the local choices

Each factor is one of the fixed coefficient spaces proved in the linked local
packages:

- [24/26-label character certificates](symplectic-seed-character-barrier.md);
- [five-site scalar compressions](common-functional-rank-compressions.md),
  [six-site certificate](two-even-site-P6-character-barrier.md), and
  [seven-site certificate](two-even-site-P7-jordan-barrier.md);
- [strengthened eight-site and nine-site row corners](two-even-site-row-corner-barrier.md);
- [cyclic cubic outputs](cubic-cyclic-all-power-barrier.md) and
  [Jordan cubic outputs](cubic-jordan-all-power-barrier.md);
- [fixed-image rank-one outputs](cubic-fixed-image-rank-one-barrier.md) and
  [coordinate-image mixtures](cubic-coordinate-image-mixtures.md).

For a rank-one factor, choose its image vector \(u\) once. Every covector \(v\)
and every coefficient recombination allowed by that theorem may vary within
that factor. Different tensor positions may choose different fixed vectors
\(u\). A union over distinct general image vectors inside one local factor is
outside this corollary.

The coordinate-image theorem separately permits all coordinate vectors \(e_i\)
to vary within its one local space. Its scope concerns mixtures of maps that
individually have one coordinate image. It does not include the cross terms of
an arbitrary complementary map with several nonzero output rows.

Likewise, a tensor factor can be cyclic or Jordan, but an arbitrary union of the
cyclic and Jordan spaces within one factor is not added by this corollary. These
are tensor-product closure statements. Enlarging a local coefficient space
requires its own compatible certificate.

## Proof with the cyclic factors included

For every noncyclic factor \(j\), let \(\lambda_j\) be its certified linear
functional, normalized by \(\lambda_j(E_j)=1\). Its row and column compressed
matrices are upper triangular, with \(\lambda_j\) on every diagonal entry.
Scalar compressions are a special case.

For a cyclic factor, write \(e_j\) for the row functional and \(v_j\) for the
second column diagonal form from its Fourier proof. Its row compression is
\(e_jI_{14}\). Its column restriction is simultaneously upper triangular,
with four copies of

\[
e_j+xv_j
\]

for each seventh root of unity \(x\). Extend the ground field to contain these
roots if necessary; all matrix ranks remain unchanged under the extension.

Suppose there are \(m\) cyclic factors. Define the product functional

\[
\Lambda=\bigotimes_{j\text{ noncyclic}}\lambda_j
 \ \otimes\!\bigotimes_{j\text{ cyclic}}e_j.
\]

Tensoring the fixed row maps gives an upper-triangular compression of dimension

\[
R=14^m\prod_{j\text{ noncyclic}}r_j,
\]

with every diagonal entry equal to \(\Lambda(T)\), even for an arbitrary
\(T\in\mathcal V\). Consequently \(\Lambda(T)\ne0\) implies
\(\operatorname{rank}R(T)\ge R\).

The tensor product of the column maps has diagonal entries given by

\[
F_T(x_1,\ldots,x_m)
 =\left(\bigotimes_{j\text{ noncyclic}}\lambda_j\right)
   \otimes\left(\bigotimes_{i=1}^m(e_i+x_iv_i)\right)(T),
\]

on the product of the seven-point root sets. Each value is repeated
\(4^m\prod_{j\text{ noncyclic}}c_j\) times. This polynomial has degree at most
one in each variable; its constant coefficient is \(\Lambda(T)\).

If \(\Lambda(T)\ne0\), the polynomial is nonzero. The product-grid lemma in
the cyclic proof says that it is nonzero on at least \(6^m\) grid points.
The column rank is therefore at least

\[
C=24^m\prod_{j\text{ noncyclic}}c_j.
\]

This argument uses the complete polynomial associated with \(T\). It does not
assume that \(T\), or its diagonal polynomial, factors across tensor positions.
For \(m=0\), the same statement is the common-diagonal compression argument
without any grid variables.

Finally, in a decomposition \(\sum_aT_a=\bigotimes_jE_j\), the sum of
\(\Lambda(T_a)\) is one. At least one piece has nonzero functional value.
That piece alone contributes at least

\[
\sqrt{RC}
 =336^{m/2}\prod_{j\text{ noncyclic}}\sqrt{r_jc_j}
 >\prod_jd_j^{17/12}=D^{17/12}
\]

to the rank cost. This proves the corollary.

## Arithmetic and dependency check

The local matrix identities remain in their original proof packages. The small
[corollary checker](../tools/check_mixed_common_functional_barrier.py) reads their
verified dimensions, pins the source-report hashes, and checks every exact
\((rc)^6>d^{17}\) margin. Its
[report](../evidence/mixed-common-functional-barrier.json) records the table and
its scope.

```bash
python tools/check_mixed_common_functional_barrier.py \
  --output /tmp/mixed-common-functional-barrier.json
```

The ten-site Frobenius partial certificate is deliberately absent from the table:
its currently certified exponent does not reach \(\sqrt2\).
