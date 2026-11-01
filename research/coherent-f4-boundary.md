# Coherent boundary retractions over the four-element field

## Construction and current scope

The height-16 boundary construction now has a complete exact evaluator
over `F_4`. It covers all 5,440 rank-three idempotents on the four
two-particle boundary states. Frobenius and transposition reduce this
declared class to 1,416 representatives without changing either complete
Smith profile. A deterministic driver checks these representatives and
retains a checkpoint after every completed candidate.

The catalogue's `complete` field determines whether the whole class has
finished. A partial file certifies its completed prefix only. No full
F4 optimality claim or stronger factoring score follows from an unfinished
run. The reference remains the explicit coherent crossing with
`rho < 1781/2501` from
[the binary boundary construction](coherent-hot-compound-retractions.md).

The retained search stopped after 153 representatives, covering 530
projectors. Its 153 distinct profile pairs consist of the reference and
152 profiles with strictly larger product pressure at every scale. The
separate [protected-submatrix theorem](coherent-boundary-protected-floor.md)
then proved `rho >= 17/24 > 1/sqrt(2)` for this entire family, over every
characteristic-two field and even allowing all five height-16 directions.
The catalogue was stopped at its saved prefix because finishing the F4
enumeration cannot produce the required exponent improvement. It remains
an exact partial catalogue, with no full F4 optimum asserted.

The new evaluator reproduces both complete binary reference profiles.
It also agrees with independently constructed literal Pfaffian tensors
for two six-letter examples having genuinely nonbinary coefficients.
An independent exact interval check retains the reference rate certificate.
The numerical pressure fit is used only to order discoveries; every
all-scale comparison retained by the driver is an integer polynomial
certificate.

## 1. Every rank-three boundary projector gives a strict crossing

Work over `F_4=F_2[omega]/(omega^2+omega+1)`. Retain the 28 coordinate
occupation words of height below 16 and the four-particle word 15.
The remaining boundary has the ordered basis

\[
(e_{24},e_{36},e_{66},e_{129}).
\]

All four words have height 16 and occupation two. Let `w` be a nonzero
column vector in `F_4^4`, and let `ell` be a row vector satisfying
`ell(w)=1`. Then

\[
\Pi=I-w\ell
\tag{1}
\]

is a rank-three idempotent, with kernel the line spanned by `w` and
image `ker ell`. Conversely every rank-three idempotent has this form.
Normalize the first nonzero coordinate of `w` to one; its covector
`ell` is then unique for the given projector.

Choose an index `p` with `ell_p!=0`. For each `i!=p` put

\[
P_i=e_i+\frac{\ell_i}{\ell_p}e_p,
\qquad
L_i=e_i^{\mathsf T}+w_i\ell.
\tag{2}
\]

Direct multiplication gives `LP=I_3` and `PL=Pi`. Append these three
columns and rows to the 29 retained coordinate states. This gives a
32-letter physical embedding and retraction with exact `LP=I_32`.
Every paired vector has one height and one occupation number.

Every nonidentity bit deletion lowers height. All derivatives of a
boundary vector land in the retained coordinate space of height below
16; the same is true of the retained four-particle word. The two map
spaces therefore satisfy the deletion-closure hypotheses of the existing
coherent construction. These arguments hold over every extension of
`F_4`, and all coefficient calculations here take place in the fixed
four-element field.

Apply the maps to the same reflection tensor as in the binary theorem.
The physical endpoint is the scale-one crossing because `LP=I`. The
two natural cut matrices are computed separately:

\[
M_0=(L\otimes L)E_B(P\otimes P),
\qquad
M_\infty=(L\otimes P^{\mathsf T})E_B(L^{\mathsf T}\otimes P).
\tag{3}
\]

The triangular deletion factors remain unimodular after compression.
Consequently (3) gives the complete natural Smith profiles, including
all occupation blocks. No profile from one cut substitutes for the other.
The opposite common scalar shifts cancel in the product pressure, exactly
as in the binary construction.

Changing the three-dimensional boundary basis changes (3) only by
constant invertible row and column transformations. Since the boundary
heights coincide, its physical type data also remain unchanged. Thus
enumerating projectors covers every retraction in this declared class,
up to this irrelevant basis choice.

## 2. The exact finite class and its two symmetries

There are

\[
\frac{4^4-1}{4-1}=85
\]

possible kernel lines. For each normalized `w`, exactly `4^3=64`
covectors satisfy `ell(w)=1`. Hence (1) gives exactly 5,440 projectors.

Two commuting involutions preserve both profiles:

1. **Frobenius.** Square every field coefficient. This is an automorphism
   fixing the binary kernel. A determinant coefficient is zero exactly
   when its image is zero, so all determinantal valuations are preserved.
2. **Transpose.** Replace `Pi` by `Pi^T`. Up to boundary basis changes,
   this interchanges `P` with `L^T`. The zero matrix in (3) is transposed.
   At infinity the two physical tensor factors are interchanged. The
   reflection one-particle matrix has form `[[D,J],[J,D]]`, so it is
   symmetric and invariant under that interchange. Its compound matrix
   has the same property. The new infinity matrix is therefore related
   to the old one by simultaneous row and column permutations.

For normalized key `(w,ell)`, transposition is implemented exactly as

\[
(w,\ell)\longmapsto
\left(a^{-1}\ell^{\mathsf T},\;a w^{\mathsf T}\right),
\tag{4}
\]

where `a` is the first nonzero coordinate of `ell`.

The fixed-projector counts of Frobenius, transpose, and their product
are respectively 120, 64, and 40. The first count is the rank-three
binary class, `15*8`. A transpose-fixed rank-one complement is determined
by a line with `sum_i w_i!=0`, giving `85-21=64` lines. For the
Frobenius-transpose product, the complement is Hermitian; the quantity
`sum_i w_i^3` is nonzero exactly when the support has odd size. The
number of such projective lines is

\[
\frac{\binom41 3+\binom43 3^3}{3}=40.
\]

Burnside's formula therefore gives

\[
\frac{5440+120+64+40}{4}=1416
\tag{5}
\]

orbits. The checker independently enumerates every projector, checks its
rank and idempotence, verifies both involutions and their commutation,
and reconstructs all three fixed-point counts. A nonbinary orbit is
also evaluated independently in all four orientations, comparing all
8,192 resulting Smith orders.

## 3. Complete exact field-valued compression

The native evaluator accepts sparse physical vectors with coefficients
`1,omega,1+omega`. It computes the binary Gaussian minors and sums them
with the actual products of their F4 map coefficients. It checks every
compressed block endpoint over the same field before calculating Smith
orders. A missing visible pivot causes a larger series window or an
explicit failure; an incomplete profile is not recorded as a candidate.

The normalized reflection minors use the existing height-based degree
bound. Their one-particle entries are

\[
K_{ij}=(1+t)t^{i\mathbin{\mathrm{xor}}j}
\quad\text{within a side},\qquad
K_{ij}=[i=j]\quad\text{between sides}.
\]

For the selected row and column occupations, the complete unnormalized
minor is a monomial times this normalized determinant, and the determinant
degree is bounded by that monomial exponent. The evaluator retains this
bound and the monomial separately. Exact minors are cached by their two
occupation masks. The cache changes only the runtime; it does not change
the arithmetic or identify distinct physical candidates.

The independent small controls use the six-letter physical space spanned
by the vacuum, four single modes, and

\[
e_6+\omega e_9
\quad\text{or}\quad
e_6+(1+\omega)e_9.
\]

Coordinate evaluation at word 6 is the final retraction row. The checker
constructs both full tensors from literal ambient Pfaffians and computes
their profiles by separate Python fraction-free F4 elimination. Across
these examples it checks 5,184 literal endpoint entries and 144 complete
Smith orders. It then compares all 4,096 orders of the two binary
32-letter controls with the existing binary evaluator.

## 4. Pressure certificates and checkpoint interpretation

Every candidate has the same physical height multiset. Write its two
counting polynomials as `Z(r)` and `I(r)`. If the difference between its
product and the reference product has an expansion

\[
Z(r)I(r)-Z_{\rm ref}(r)I_{\rm ref}(r)
=r^a(1-r)^b\sum_{j=0}^N c_jr^j(1-r)^{N-j},
\quad c_j\ge0,
\tag{6}
\]

with at least one positive coefficient, then the candidate has strictly
larger pressure at every `0<r<1`. The type-restoration terms are identical,
so this comparison remains valid for every common physical type.

The driver uses integer arithmetic to construct (6) and reconstruct its
power-basis coefficients. It retains the quotient polynomial and the
Bernstein degree. If its finite list of Bernstein degrees does not certify
the sign, it marks the comparison unresolved. No sign follows merely from
a sampled numerical minimum. Equal product polynomials are recorded
separately.

The file stores the full two-cut order histograms for every distinct
profile pair and associates every completed projector with its profile.
Histograms retain integer multiplicities; expanding them gives all 1,024
orders at each cut. A deterministic source hash, catalogue hash, current
input, and completed prefix support interrupted-run resumption. The
readback checks every saved profile hash and exact pressure comparison.

Run the construction and independent controls with

```sh
python3 tools/coherent_f4_boundary.py \
  --seconds 7200 --output evidence/coherent-f4-boundary-orbits.json
python3 tools/check_coherent_f4_boundary.py \
  --catalogue evidence/coherent-f4-boundary-orbits.json \
  --output evidence/coherent-f4-boundary-controls.json
```

To continue an interrupted catalogue, use the same first command with
`--resume`. The source hashes must agree. For a quick readback of a later
checkpoint, add `--catalogue-only` to the checker. Full class coverage is
asserted only when `complete` is true and all 1,416 representatives,
covering all 5,440 projectors, have finished.

This work uses ordinary algebra and exact computational controls, with
substantial OpenAI assistant support. No new Lean theorem or independent
formal verification is claimed. The established factoring score remains
`1/2` unless a separate complete strict-pressure certificate crosses the
required threshold.
