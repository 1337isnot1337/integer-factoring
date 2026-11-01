# An all-field, all-power barrier for the 240-label seed

The six-site, 240-label two-even-site seed also admits a common-functional
rank certificate. The full ten-direction coefficient space fails the
finite-seed criterion over **every characteristic-two field and every
tensor power**, including arbitrary recombinations in the entire tensor
span.

The constant row and column compression dimensions are

\[
a=15688,\qquad b=718,
\qquad ab=11263984.
\tag{1}
\]

Their product is below \(240^3\), so this certificate does not claim the
stronger \(d^{3/2}\) margin established for the smaller seeds. It does
prove the sufficient bound

\[
\boxed{
S\ge11263984^{n/2}>240^{37n/25}>240^{\sqrt2 n}
}
\tag{2}
\]

for every crossing decomposition in the \(n\)-fold coefficient span.
The exponent \(37/25=1.48\) is separated from \(\sqrt2\) by exact
integer inequalities.

## 1. Seed and functional

Use the code and exact ten coefficient directions from the
[four-site-core factorization](two-even-site-core-factorization.md).
With \(u=q+r\), \(v=p+q\), \(s=qr\), the normalized polynomial is

\[
u^{-6}G_6=
\sum_{i=0}^{5}A_i t^i+z\sum_{i=0}^{3}B_i t^i,
\qquad t=v/u,\quad z=s/u^2.
\]

The coefficient \(A_0\) is the crossing. In the ordered basis
\((A_0,\ldots,A_5,B_0,\ldots,B_3)\), define

\[
\chi=(1,0,0,0,1,0,0,0,0,0).
\tag{3}
\]

Its value on the polynomial family is

\[
\chi(G_6)=u^6+u^2v^4
=u^2(u+v)^4
=(q+r)^2(p+r)^4.
\tag{4}
\]

In particular, \(\chi(A_0)=1\).

## 2. Direct constant compressions

Apply the
[common-functional subspace criterion](common-functional-rank-compressions.md)
separately to each of the 31 syndrome blocks on both cuts. For each
block, the witness is a coordinate subspace \(U\) satisfying

\[
U\cap\sum_i(M_i-\chi_iI)U=\{0\}.
\tag{5}
\]

The certificate lists only the selected coordinate indices. From
these indices, the independent checker constructs the actual left
covectors, producing constant matrices

\[
W_rM_r(T)V_r=\chi(T)I_{15688},\qquad
W_cM_c(T)V_c=\chi(T)I_{718}
\tag{6}
\]

for every tensor \(T\) in the full ten-direction span. These are
binary matrix identities. They persist after extension to every
field of characteristic two. Neither flattening algebra needs to
have \(\chi\) as an algebra character.

The row witnesses contain 15,688 coordinate vectors; their total
difference-image dimension is 37,472. The column witnesses contain
718 coordinate vectors, with total difference-image dimension 5,123.
The direct-sum conditions and all generator identities are checked
block by block.

## 3. Tensor powers and exact separation

Tensoring (6) gives scalar identity compressions of sizes
\(15688^n\) and \(718^n\). Therefore

\[
\chi^{\otimes n}(T)\ne0
\quad\Longrightarrow\quad
r(T)c(T)\ge11263984^n.
\]

If \(A_0^{\otimes n}=\sum_jT_j\), then
\(\sum_j\chi^{\otimes n}(T_j)=1\). At least one summand has nonzero
functional, so its cost alone proves the first inequality in (2).
This argument allows arbitrary pieces in the tensor-product span,
not just pure tensors or coefficients of a common parameter curve.

The exact remaining comparisons are

\[
11263984^{25}>240^{74},\qquad
37^2>2\cdot25^2.
\]

Taking positive real powers yields
\(\sqrt{11263984}>240^{37/25}\) and \(37/25>\sqrt2\).
The checker records the integer difference in the first comparison.
No floating-point threshold enters the proof.

## 4. Mixed products of the four closed spaces

The previous certificates and this one give:

| Label dimension | Row certificate size | Column certificate size | Rank-product lower bound |
|---:|---:|---:|---:|
| 24 | 384 | 43 | 16,512 |
| 26 | 436 | 43 | 18,748 |
| 80 | 2,072 | 286 | 592,592 |
| 240 | 15,688 | 718 | 11,263,984 |

The first two row certificates are triangular restrictions; the
others are scalar compressions. Their tensor products retain the
same rank consequence when the product functional is nonzero.

For any nonempty mixed product with \(n_d\) copies of each listed
seed and total label dimension \(D=\prod_d d^{n_d}\), every crossing
decomposition in the full mixed tensor-product span satisfies

\[
S\ge\prod_d (a_db_d)^{n_d/2}
>D^{37/25}>D^{\sqrt2}.
\tag{7}
\]

This includes independent local parameter curves and extension
fields. It does not apply automatically to altered input or output
maps, additional tensor directions, or the seven-site construction.
It is not a lower bound for integer factoring.

## 5. Faster witness discovery and independent verification

The new finder preserves the earlier discovery tool and its fixed
certificates. It maintains a Gaussian basis of \(Z(U)+U\), tagged
by the desired coordinate-map values on \(U\). Every proposed
difference image must map to zero. A relation with nonzero tag
rejects the candidate; a new independent basis vector extends the
map. The proposed coordinate vector must map to a new unit vector.
This tests (5) incrementally without rechecking every previous
coordinate vector after each proposal.

Only one shuffled coordinate-vector pass with seed zero was needed.
No dense coefficient algebra was generated. The resulting witness
indices are independent of the later verification algorithm.

The checker reconstructs the six-site matrices from exact small
cores. It checks ten independent coefficient coordinates and repeats
346 direct literal polynomial controls. It then uses a separate
low-pivot elimination to reconstruct the left covectors and checks
all entries of (6): **89,439,480 row-cut entries and 225,160 column-cut
entries**. This is a direct finite matrix certificate for the
all-field theorem.

Files:

- [find_tagged_common_character_compressions.py](../tools/find_tagged_common_character_compressions.py);
- [check_two_even_site_compressions.py](../tools/check_two_even_site_compressions.py);
- [two-even-site-P6-compression-certificate.json](../evidence/two-even-site-P6-compression-certificate.json);
- [two-even-site-P6-compression.json](../evidence/two-even-site-P6-compression.json).

Independent replay:

```sh
python tools/check_two_even_site_compressions.py \
  --output /tmp/two-even-site-P6-compression.json
```

Deterministic discovery:

```sh
python tools/find_tagged_common_character_compressions.py \
  --pairs 6 --character 1000100000 \
  --output /tmp/two-even-site-P6-compression-certificate.json
```
