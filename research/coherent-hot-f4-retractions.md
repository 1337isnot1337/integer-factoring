# Coherent physical retractions over F4

## Construction and current scope

The homogeneous, deletion-closed physical compression in
[the binary construction](coherent-hot-compound-retractions.md)
works with coefficients in F4 as well. This note supplies an exact
weighted implementation and independent coefficient controls. It does
not claim an improved factoring exponent or a completed search over
all F4 physical maps.

The purpose is constructive: field-extension parameters may permit
cancellations absent from the binary boundary catalogue. The existing
binary optimum does not settle that question. Both complete cut
profiles must be computed; improving only one cut is insufficient.

## 1. Paired maps and their endpoint

Use the field F4 with omega squared equal to omega plus one. Codes
0, 1, 2, 3 represent 0, 1, omega, 1+omega. Pairings are bilinear,
not Hermitian: no Frobenius conjugation appears in a contraction.

Let P map a q-dimensional label space into the occupation basis.
Let L map back, with LP=I. Require each paired column of P and row
of L to have one common physical height and occupation. Require
both their spans to be invariant under deleting each occupied bit.
The selected compound entries are

```math
M_{ij}(t)=\sum_{R,C} r_{i,R}\,c_{j,C}\det B(t)_{R,C},
```

with the same binary one-particle matrix B as the
[occupation-block constructor](https://github.com/1337isnot1337/integer-factoring/blob/8d3a2da11eca5e1a7d064b645c698719ac141098/research/hot-code-charge-profile-constructor.md).
The row and column coefficients now lie in F4. At t=1 the
one-particle matrix swaps the two mode sets. Multilinearity therefore
reduces the endpoint to the two LP identities, giving the literal
identity on every paired occupation block.

The deletion-closure argument is unchanged: the elementary Gaussian
factors act triangularly inside each selected space. Their inverse
factors do as well. Extending the coefficient field does not change
the polynomial identities or their invertibility. Thus the two
regularized natural cuts are obtained by using these maps:

- Zero cut: rows L tensor L, columns P tensor P, with the paired swap.
- Infinity cut: rows L tensor P, columns P tensor L, with the paired swap.

This is the same paired convention as the binary constructor. The
implementation validates LP=I, homogeneity, and both deletion closures
before computing either cut. It returns exactly q squared Smith
orders per cut, or raises an error.

## 2. Exact coefficient and series arithmetic

Each determinant remains binary. The native implementation computes
that determinant with the existing exact hot-kernel degree bound,
then multiplies it by its F4 map coefficient and accumulates two
binary coefficient planes. It does not encode omega as a polynomial
in t. Such an encoding would change valuations and is not valid here.

For a minor with k occupied modes, the hot-kernel factorization writes
its determinant as t to the power e times a binary polynomial of
degree at most e. Here e is the sum of the row and column mode indices
plus k. Homogeneous sparse vectors have the same e for every summand.
Consequently the complete weighted entry has degree at most 2e;
field-linear summation cannot increase that bound.

The native F4 DVR kernel then computes the complete Smith profile.
Weights are zero, so every reported pivot order is an exact integer.
The endpoint proves full generic rank. If the series window does not
resolve all pivots, the wrapper retries with a larger window and
ultimately fails rather than returning a partial profile.

## 3. Independent executable controls

Run `python3 tools/check_coherent_hot_f4.py` from the repository root.
The [receipt](../evidence/coherent-hot-f4-controls.json) records:

- 2,412 complete F4 polynomial entries compared against literal
  weighted contractions of the independently generated binary ambient
  compound blocks, including 738 entries with a nonzero omega plane.
- Both complete binary-limit profiles compared to the prior engine.
- Both profiles invariant under Frobenius conjugation and paired
  inverse rescaling, for four additional profile comparisons.
- Rejection of invalid endpoint, homogeneity, and deletion-closure maps.

The tests retain source hashes. They validate an exact finite
constructor; they are not evidence of a better numerical exponent.

## 4. A concrete field-extension search space

For the 32-label boundary, retain the existing 28 lower-height words,
the height-16 four-particle word 15, and a three-dimensional retract
of the two-particle boundary words 24, 36, 66, 129.

Choose a nonzero kernel vector k and covector c over F4 with c k=1.
Then Q=I+k c is an idempotent of rank three in characteristic two.
Choose a pivot p with c_p nonzero. An explicit basis and dual are

```math
P_i=e_i+(c_i/c_p)e_p,\qquad L_i=e_i^{\mathsf T}+k_i c,
\qquad i\ne p.
```

Direct multiplication gives c P_i=0 and L_i P_j=delta_ij. Their
product is Q. Every deletion of a boundary vector lands in the fully
retained lower-height space, so both maps satisfy the constructor's
closure requirement.

Normalize the first nonzero coordinate of k to one. There are 85
projective kernel lines over F4, and 64 covectors satisfying c k=1
for each normalized k. These give 5,440 distinct rank-three
idempotents. The corresponding binary class has 120. A finite sample
of this larger class is not an exhaustive exclusion or a global
optimality theorem.

`tools/search_coherent_hot_f4.py` records both complete cut profiles
for each completed trial. Its pressure minimization is numerical and
used only to identify candidates for later exact certification. The
complete factoring exponent remains 1/2.

Substantial OpenAI assistant support was used for this construction,
implementation, tests, and text.
