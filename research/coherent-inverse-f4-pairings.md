# F4 inverse systems and bilinear dual physical codes

The derivative inverse-system construction now works with genuine F4
coefficients at every occupied height. A second construction replaces its
coordinate retraction by any nondegenerate pairing with another homogeneous,
deletion-closed space. In the self-paired case, constant Gram units identify
the two complete regular Smith profiles.

For an explicit two-parameter quartic family on eight modes, the inverse
system has Hilbert function `(1,8,22,8,1)` on `ab=1`, and
`(1,8,28,8,1)` elsewhere, over every characteristic-two field. Six explicit
two-by-two blocks explain the entire jump. All sixteen F4 members have
completed exact two-cut profiles. Their numerical pressure fits, and four
additional Gram-paired probes, stay above `1/sqrt(2)`. These finite probes
do not exclude other inverse systems or pairings. The complete factoring
score remains **1/2**.

## 1. Field-valued inverse systems

Use occupation words, heights `h(A)=sum_(i in A)(2i+1)`, and commuting bit
deletions `D_i` from the
[binary construction](coherent-hot-compound-retractions.md). For finitely
many forms homogeneous in height and occupation, take the span of all their
iterated deletions. The leading-word downset proof in that note uses only
field elimination, so applies unchanged over F4 or any other field.

The [F4 constructor](../tools/coherent_inverse_systems_f4.py) normalizes
each pivot by its field inverse and performs reduced elimination. Its
columns are the resulting basis, and its retraction rows are the leading
occupation coordinates. It explicitly checks `LP=I`, both deletion-closure
conditions, and the common height and occupation of each label.

F4 is encoded by `0,1,2,3`, with `2=omega` and
`omega^2=omega+1`. Pairings are ordinary bilinear pairings. No conjugation
is applied. The mathematical construction is field-independent; this
implementation and the retained finite enumeration use F4.

At most `2^m` independent derivative vectors enter the work queue on `m`
modes. There are at most `m*2^m` single-deletion insertion attempts after
the initial top forms. Elimination and the subsequent exact compound
evaluation retain their full finite costs. This is an effective local
construction and supplies no new uniform factoring runtime by itself.

## 2. Pairing independently chosen spaces

Let the columns of `P` and `Q` be bases of two spaces of the same dimension,
both homogeneous in height and occupation and closed under all deletions.
Suppose their ordinary pairing matrix

\[
G=Q^{\mathsf T}P
\]

is invertible. Then

\[
L=G^{-1}Q^{\mathsf T},\qquad LP=I.
\]

The row space of `L` is exactly the transpose of the space spanned by `Q`,
so it is deletion closed. Distinct homogeneous cells pair to zero. The
inverse pairing therefore produces a retraction row of the same height
and occupation as each embedding column. This proves every hypothesis of
the existing complete two-cut compound construction.

This argument works over any field for which the stated pairing is
nondegenerate. The [implementation](../tools/coherent_orthogonal_f4.py)
uses exact F4 elimination, rejects singular pairings, and rechecks the
physical-map hypotheses. Changing the basis of the second space leaves
`L` unchanged. Choosing another second *space* can change the crossing.

For `Q=P`, the matrix `G=P^T P` is symmetric. Write the common ambient
compound restriction before the harmless column-pair swap as

\[
C=(P^{\mathsf T}\otimes P^{\mathsf T})E_B(P\otimes P).
\]

The zero and infinity regular matrices are respectively

\[
(G^{-1}\otimes G^{-1})C,
\qquad
(G^{-1}\otimes I)C(G^{-\mathsf T}\otimes I).
\]

Both are obtained from `C` by constant invertible matrices. Hence their
complete Smith orders agree. The Gram matrix preserves both gradings,
so the same constant-unit argument is compatible with occupation blocks
and the inherited height normalization. The routine `orthogonal_profiles`
computes one complete cut and explicitly labels the second as derived by
this theorem. It does not claim to compute the second cut independently.

## 3. A quartic catalecticant with one exceptional locus

On eight modes let

\[
g=x_0x_7+x_1x_6+x_2x_5+x_3x_4,
\qquad
f_{a,b}=\gamma_2(g)+a x_0x_3x_5x_6+b x_1x_2x_4x_7.
\]

Here `gamma_2(g)` means the sum of the six products of distinct displayed
pair terms. This is a divided square, rather than the ordinary square in
characteristic two. Every monomial has occupation four and height 32.
The eight occupation words, in increasing order, are
`60,90,102,105,150,153,165,195`; the coefficients are
`1,1,1,a,b,1,1,1`.

The order-two catalecticant is indexed by two-element subsets and has entry
`coefficient_(A union B)(f)` for disjoint `A,B`, and zero otherwise.
It splits into blocks pairing derivative height `h` with `32-h`.
For `h=8,12,14,18,20,24`, the block is, up to row and column ordering,

\[
\begin{pmatrix}1&a\\b&1\end{pmatrix}\oplus I
\]

or its transpose. The other blocks have determinant one. The central
four-by-four block is `I+J`, whose square is `I` in characteristic two.
Consequently the full determinant is `(1+ab)^6`. If `ab=1`, each of the
six displayed blocks has rank exactly one less, since its diagonal is
nonzero. Otherwise every block is invertible. The middle derivative
dimension is therefore 22 or 28 exactly, over every characteristic-two
field, including nonfinite extensions.

The eight first derivatives are nonzero and have distinct heights, so are
independent. The third catalecticant is their transpose under occupation
pairing and also has rank eight. The top and vacuum dimensions are one.
This proves the stated complete Hilbert functions and total dimensions
40 and 46. The exceptional locus removes one quadratic state at each of
heights `8,12,14,18,20,24`. A smaller physical alphabet alone does not
guarantee better crossing pressure.

## 4. Completed finite probes

The [retained sixteen-member family](../evidence/coherent-inverse-f4-paired-family.json)
contains four original source-pinned shards, their exact input text, maps,
complete orders, native precision receipts, and numerical pressure fits.
Every declared member completed; no partial Smith profile is used. The
best fit in this specified family is the binary member `a=b=0`, with
46 letters and `rho` approximately `0.7155813177`.

| Retraction | Top coefficients `(a,b)` | Dual top | Letters | Numerical rho |
|---|---|---|---:|---:|
| Coordinate | `(0,0)` | Leading downset | 46 | 0.7155813177 |
| Coordinate | `(1,1)` | Leading downset | 40 | 0.7181143032 |
| Self-paired | `(omega,omega+1)` | Same space | 40 | 0.7181191868 |
| Self-paired | `(omega,0)` | Same space | 46 | 0.7155853113 |
| Paired | `(1,0)` | `(omega,0)` | 46 | 0.7155853116 |
| Paired | `(0,0)` | Coefficient of word 60 changed to omega | 46 | 0.7155852012 |

The four Gram probes have separate receipts
[0](../evidence/coherent-gram-f4-case-0.json),
[1](../evidence/coherent-gram-f4-case-1.json),
[2](../evidence/coherent-gram-f4-case-2.json), and
[3](../evidence/coherent-gram-f4-case-3.json).
These values are floating-point one-variable minimizations used for
discovery. They are not rational rate certificates or proofs of global
optimality. Neither this family nor the larger class has been exhausted.
The coefficient choices do not improve on the existing
[96-letter construction](coherent-q96-pressure.md).

## 5. Verification and replay

Run the small exact controls with:

```sh
make verify-coherent-inverse-f4 verify-coherent-gram-f4
```

The [inverse-system checker](../tools/check_coherent_inverse_systems_f4.py)
and [receipt](../evidence/coherent-inverse-f4-controls.json) retain:

- Three exact specializations to the existing binary constructor.
- Forty-eight derivative reconstructions from independently computed binary
  coefficient planes, plus presentation and Frobenius covariance controls.
- All 908 literal F4 polynomial entries of a six-letter nonbinary example,
  including 85 nonbinary entries, assembled from independent permutation
  determinants of the ambient one-particle matrix.
- All 72 Smith orders of those two cuts, replayed with the separate Python
  fraction-free oracle, and two binary-cut specialization controls.
- Symbolic height-block determinants over `F2[a,b]` and the constructor's
  Hilbert function at every one of the sixteen F4 parameter pairs.

The [Gram checker](../tools/check_coherent_gram_f4.py) and
[receipt](../evidence/coherent-gram-f4-controls.json) verify both products
of each inverse Gram matrix, independence from the second basis, singular
pair rejection, and all 522 orders of both cuts on three small self-paired
codes. At least one case has a non-diagonal homogeneous Gram block.
The one-cut derived routine agrees with those separately computed cuts.

The larger family can be recomputed with:

```sh
python3 tools/search_coherent_inverse_systems_f4.py --modes 8 --particles 4 --height 32 --candidate-file evidence/coherent-inverse-f4-paired-inputs.json --max-trials 16 --seconds 7200 --output /tmp/coherent-inverse-f4-replay.json
python3 tools/search_coherent_gram_f4.py --case 0 --output /tmp/coherent-gram-f4-case-0.json
```

Repeat the last command for cases 1, 2 and 3. Times and numerical optimizer
metadata need not be byte-identical. Compare the input forms, maps, complete
orders and completion flags. The quick targets rerun the small exact
controls, rather than repeating all sixteen large native evaluations.
All checks here are local author checks; no independent mathematical
review or stronger leaderboard score is claimed. No new Lean work is
required for this checkpoint.
