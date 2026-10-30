# Small row corners close the nine-site coefficient span

## Result

A calculation on five-site matrices of dimension at most **640** gives stronger
row certificates for the two-even-site family. Combined with the published
bounded column witnesses, it closes the full coefficient span at nine sites,
with **4,608 labels**:

\[
R_9=7522480,\qquad C_9=3332,\qquad R_9C_9=25064903360.
\]

For every characteristic-two field, every tensor power, and every decomposition
inside the full tensor coefficient span,

\[
S\ge25064903360^{n/2}>4608^{17n/12}>4608^{\sqrt2 n}.
\]

Both comparisons are certified by exact integer arithmetic:

\[
25064903360^6>4608^{17},\qquad17^2=289>288=2\cdot12^2.
\]

The same row calculation strengthens the eight-site certificate to

\[
R_8=1225260,\qquad C_8=2106,\qquad R_8C_8=2580397560,
\]

which also exceeds the exponent \(17/12\). The earlier eight-site package stays
unchanged; this note gives a separate strengthening.

The proof below also yields explicit common-functional rank bounds for every
\(P\ge8\). Those bounds do not exclude the target exponent at all sizes. The
specific all-power exclusions established here are the eight- and nine-site
seeds. These are exclusions of fixed coefficient spaces, not lower bounds for
general factoring algorithms.

## A compression criterion from joint nilpotence

Let \(M_0=I_n,M_1,\ldots,M_s\) be matrices over a field. Fix scalars
\(\chi_0=1,\chi_1,\ldots,\chi_s\). Choose maps

\[
V:\mathbb F^k\longrightarrow\mathbb F^n,\qquad
W:\mathbb F^n\longrightarrow\mathbb F^k,\qquad WV=I_k.
\]

Write

\[
B_i=WM_iV,\qquad N_i=B_i-\chi_iI_k.
\]

Suppose the operators \(N_i\) are jointly nilpotent: every sufficiently long
product of them is zero. Then the compressed matrices \(B_i\) have a common
upper-triangular basis, with every diagonal entry of \(B_i\) equal to
\(\chi_i\).

To see this, form the descending image filtration

\[
F_0=\mathbb F^k,\qquad F_{j+1}=\sum_iN_iF_j.
\]

It is descending because \(F_1\subseteq F_0\), and this inclusion propagates
under every \(N_i\). Joint nilpotence gives \(F_\ell=0\) for some \(\ell\).
Choose a basis by extending from the deepest nonzero filtration space toward
\(F_0\). Each \(N_i\) sends a basis vector into earlier filtration spaces, so
it is strictly upper triangular in this basis.

This criterion does **not** require the image of \(V\) to be invariant under the
original matrices. It also does not assert that compression preserves products
of original matrices. Only the actual compressed matrices \(WM_iV\) occur in
the criterion.

A useful way to select the two maps is to choose an idempotent \(e\), take
\(V\) to be a basis of its image, and let \(W\) annihilate its kernel while
satisfying \(WV=I\). Then \(VW=e\). The present calculation uses small powers
of one coefficient matrix to construct \(e\), but the criterion itself applies
to any maps with \(WV=I\).

For tensor products, compression does have the required identity:

\[
(W_1\otimes W_2)(M_i\otimes M_j)(V_1\otimes V_2)
=(W_1M_iV_1)\otimes(W_2M_jV_2).
\]

Tensor products and linear combinations of upper-triangular matrices remain
upper triangular. If all individual diagonals are the selected functional,
then every tensor-space element has the corresponding tensor functional on its
entire diagonal. No invariance assumption is introduced by this step.

## Exact five-site row corners

Use the normalized coefficient family from the
[all-size core factorization](two-even-site-core-factorization.md):

\[
G_P(t,z)=\sum_{i=0}^{P-1}A_it^i+z\sum_{i=0}^{P-3}B_it^i,
\qquad A_0=E.
\]

For the five-site row cut, the desired functional has scalar polynomial

\[
\chi_5(G_5)=(1+t)^3=1+t+t^2+t^3.
\]

On each syndrome block, set \(e=A_1^{2^a}\), using the smallest tested power
that is idempotent. Construct \(V\) from its image and \(W\) from its kernel.
The exact compressed difference filtrations terminate as follows.

| Relative syndrome | Original size | Idempotent exponent | Corner dimension | Successive image dimensions |
| ---: | ---: | ---: | ---: | --- |
| 0 | 640 | 4 | 640 | 640, 354, 92, 9, 0 |
| 3 | 384 | 2 | 384 | 384, 236, 52, 8, 0 |
| 5 | 384 | 4 | 256 | 256, 168, 48, 8, 0 |
| 6 | 384 | 4 | 384 | 384, 260, 80, 12, 0 |
| 9 | 384 | 4 | 128 | 128, 88, 24, 0 |
| 10 | 384 | 4 | 256 | 256, 184, 64, 8, 0 |
| 12 | 384 | 4 | 384 | 384, 260, 80, 12, 0 |
| 17 | 384 | 4 | 0 | 0 |
| 18 | 384 | 4 | 128 | 128, 88, 24, 0 |
| 20 | 384 | 4 | 256 | 256, 168, 48, 8, 0 |
| 24 | 384 | 2 | 384 | 384, 236, 52, 8, 0 |

If a nonzero relative syndrome has endpoints \(a<b\) inside the ordered
five-site core, its corner dimension is

\[
128\bigl(4-(b-a)\bigr).
\]

These are exact finite identities, not dimensions inferred from sampled ranks.
The certificate stores a basis already adapted to the nilpotent image
filtration. The independent audit reconstructs the original matrices and the
left maps, then directly checks the diagonal and every entry below it in every
compressed generator. It also recomputes each full image-dimension sequence.
Thus the triangular conclusion receives a direct matrix audit in addition to
the filtration argument.

The zero-dimensional corner at syndrome 17 is deliberately retained in the
report. The proof handles the corresponding extreme physical syndrome with an
older scalar witness, as described next.

## Embedding the corners into larger seeds

For \(P\) physical sites the two-even-site code has

\[
d_P=\binom P2 2^{P-2}
\]

logical labels. Pair coordinates in each cut split into syndrome blocks of
weights zero, two, and four.

At a free site the entire row factor is upper triangular, after a constant
similarity, with four diagonal entries equal to \(1+t\). The column factor
has the principal scalar compression

\[
e_0^T(I_4+tN)e_0=1+t.
\]

Here \(N\) is the published rank-one square-zero local column matrix. Therefore:

- A five-site row corner of dimension \(r\), tensored with \(P-5\) full free
  row factors, gives dimension \(r4^{P-5}\) and functional \((1+t)^{P-2}\).
- A seven-site scalar row certificate, tensored with \(P-7\) full free row
  factors, gives its old dimension times \(4^{P-7}\) and the same functional.
- An eight-site scalar column certificate, tensored with \(P-8\) copies of
  the local coordinate \(e_0\), retains its old dimension and has the same
  functional.

The original homogeneous functional in every case is

\[
\chi_P(G_P)=(q+r)^2(p+r)^{P-2}.
\]

### Weight-four blocks

Each entire weight-four block is its four-site core tensored with \(P-4\)
free factors. The existing row core has a triangular certificate of dimension
32, giving \(32\cdot4^{P-4}\) per physical block. The eight-site column
certificate supplies dimension 16 per block, and the local scalar lift retains
that dimension for larger \(P\).

### Weight-two principal restrictions

For syndrome \(S=\{a,b\}\), choose three anchor sites. Retain only the edge
pairs \(\{a,c\},\{b,c\}\), in either order, for those anchors \(c\). This
is the full five-site weight-two core on the chosen sites with \(P-5\) free
factors, embedded as a **principal restriction** of the physical block.

The row and column cuts choose their anchors independently. Both compressions
carry \(\chi_P\), so this does not obstruct the common-functional argument.
The exact physical anchor choices are stored in the certificate.

For rows, choose anchors that minimize the relative distance between \(a\) and
\(b\). For all but six physical syndromes, three anchors can lie outside the
interval, giving a relative distance of one and corner dimension 384. Three
syndromes require distance two, giving 256; two require distance three, giving
128. The single extreme syndrome \(\{0,P-1\}\) has relative distance four
and zero corner. For that syndrome, use the independently audited seven-site
scalar row certificate of dimension 1,163, then lift it by the full free row
factors.

For columns, lift the eight-site witnesses. Adjacent physical endpoint pairs
and the extreme pair have dimension eight. The two pairs \(\{0,2\}\) and
\(\{1,3\}\) have dimension 44. Every other pair can choose anchors that
realize a stored relative syndrome of column dimension 46:

- If \(a\ge2\), place two anchors below \(a\) and one between the endpoints.
- If \(a=1\) and \(b\ge4\), place one anchor below \(a\) and two between.
- If \(a=0\) and \(3\le b\le P-2\), place two anchors between the endpoints
  and one above \(b\).

These cases exhaust the nonadjacent pairs other than the two listed exceptions
and the extreme pair.

### Weight-zero principal restriction

Retain equal exceptional edges inside physical sites \(0,1,2,3,4\). The
five-site zero-syndrome row corner has dimension 640. With the free row factors
it contributes \(640\cdot4^{P-5}\). The lifted column witness contributes six.

## An explicit bound for every \(P\ge8\)

Adding the disjoint syndrome contributions gives the row certificate

\[
R_P=\left(2048\binom P4+6144\binom P2-9077\right)4^{P-7}
\]

and the column scalar compression

\[
C_P=16\binom P4+46\binom P2-38P+2.
\]

The row formula follows by summing the three nonzero core dimensions above,
the single 1,163-dimensional fallback, and the zero-syndrome contribution. The
column formula counts \(P\) blocks of dimension eight, two of dimension 44,
\(\binom P2-P-2\) of dimension 46, and the weight-four and zero blocks.

These formulas give, for every characteristic-two field and every tensor power,

\[
S\ge (R_PC_P)^{n/2}
\]

inside the corresponding tensor coefficient space. This all-size statement is
about an explicit certificate. Its ability to exceed the target is checked
separately at each size.

| Sites | Labels | Row certificate | Column compression | Rank product | Consequence |
| ---: | ---: | ---: | ---: | ---: | --- |
| 8 | 1,792 | 1,225,260 | 2,106 | 2,580,397,560 | Excludes cost exponent below \(17/12\) |
| 9 | 4,608 | 7,522,480 | 3,332 | 25,064,903,360 | Excludes cost exponent below \(17/12\) |
| 10 | 11,520 | 44,638,912 | 5,052 | 225,515,783,424 | This certificate does not reach \(\sqrt2\) |

The ten-site shortfall is real for this certificate: its product satisfies

\[
225515783424^5<11520^{14},\qquad(7/5)^2<2.
\]

This does not show that a cheap decomposition exists at ten sites. It only
identifies where these particular restrictions stop excluding the target.
The formulas have \(R_P=\Theta(P^4 4^P)\) and \(C_P=\Theta(P^4)\), so the
cost exponent supplied by these witnesses tends to one as \(P\) grows.

## Tensor powers and mixed seeds

For a fixed seed, every retained row compression is upper triangular with the
same functional on its diagonal, and the column compression is scalar with that
functional. Direct sums over distinct syndrome blocks preserve these statements.
Tensor products preserve them as well.

In any decomposition of \(E^{\otimes n}\), the sum of the tensor-functional
values is one. At least one summand has nonzero value, hence has row rank at
least \(R_P^n\) and column rank at least \(C_P^n\). Its single contribution to
\(S\) gives the bound above.

The same argument combines the certified label sizes

\[
24,\ 26,\ 80,\ 240,\ 672,\ 1792,\ 4608.
\]

For every nontrivial mixed tensor product of these fixed coefficient spaces,
with total label count \(D\), arbitrary decompositions inside the mixed tensor
span satisfy

\[
S>D^{17/12}>D^{\sqrt2}.
\]

For the eight-site factor this uses the strengthened row certificate from this
note. The earlier label sizes already have stronger certified exponents.

## Audit and reproduction

The new audit reconstructs every five-site corner, its left map, and its
triangular compressed generators. It also independently rechecks the specific
seven-site fallback row witness and all eight-site column modules that the
larger-size proof uses. The largest newly discovered row matrix has dimension
640; the largest dependency matrix has dimension 5,120.

Literal polynomial controls check the physical core embeddings at both eight
and nine sites. There are 594 eight-site controls and 978 nine-site controls,
covering entire weight-four blocks and the selected principal restrictions in
the other two classes. The all-size embedding claim follows from the published
core factorization proof and the explicit anchor construction, rather than an
extrapolation from these controls.

The source files are:

- [check_two_even_site_row_corners.py](../tools/check_two_even_site_row_corners.py);
- [two-even-site-row-corner-certificate.json](../evidence/two-even-site-row-corner-certificate.json);
- [two-even-site-row-corner.json](../evidence/two-even-site-row-corner.json).

Run the audit with literal five-site regeneration:

```bash
python tools/check_two_even_site_row_corners.py \
  --output /tmp/two-even-site-row-corner.json
```

Add `--discover` to reproduce the stored adapted bases. A previously generated
five-site basis can be supplied with `--basis5`; its full coefficient digest must
match the published five-site certificate.
