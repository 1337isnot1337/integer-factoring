# Polynomial couplings between two copies of a crossing

## Construction and scope

Let `B(w)` be a strict scalar crossing on `q` physical letters, with
polynomial degree at most `d`. Let `F,G` be arbitrary constant tensors
on the same four physical legs. For a polynomial `h` satisfying `h(1)=0`
and a field scalar `c`, define a crossing on ordered pairs of letters by

\[
T(w)_{(x_1,x_2)(y_1,y_2)(u_1,u_2)(v_1,v_2)}
 =B(w)_{x_1y_1u_1v_1}B(w)_{x_2y_2u_2v_2}
  +c h(w)F_{x_1y_1u_1v_1}G_{x_2y_2u_2v_2}.             \tag{1}
\]

At `w=1`, the correction is zero and the product has the literal crossing
endpoint on `q^2` letters. Its complete degree is at most
`max(2d,deg h)`. This is an ordinary fixed finite scalar seed, so the
existing effective tensor-power transfer applies once its two cut
profiles and a favorable pressure inequality have been proved.

The coupling can correlate the two physical copies. It is not assigned
the product of the separate copies' ranks. Both complete matrices in
(1) are constructed and evaluated after adding the correction.

This note provides a checked constructor and a finite search route. It
does not provide an improved factoring exponent. The numerical reference
remains above the required threshold.

## Exact weighted matrices

Use the same four integer physical height functions on both copies at
each place. If the one-copy cleared matrix is `A_s(w)` and its scalar
clearing shift is `a_s`, then the uncoupled two-copy cleared matrix is
exactly

\[
A_s(w)\otimes A_s(w),                                   \tag{2}
\]

up to the fixed row and column ordering of the two physical copies.
Its scalar clearing shift is `2a_s`. Consequently its complete local
orders are the pairwise sums of the one-copy orders. This gives a
complete baseline control, including multiplicities.

For a supported entry of `F tensor G`, the new polynomial is multiplied
by the actual sum of the eight physical heights. At infinity the full
polynomial is reversed at degree `2d` before this monomial is inserted.
At zero its original coefficient order is used. These operations retain
all leading and trailing zero coefficients and the complete degree cost.

The [implementation](../tools/two_copy_crossing_couplings.py) forms the
two matrices literally. The current native interface supports source
alphabets with `q^4<=1024`. Its integer inputs produce integer local
orders. Failure to resolve all promised pivots is recorded as unresolved;
it is not interpreted as a favorable rank defect. The strict endpoint
proves the required generic full rank in both cuts.

## A finite source-module catalogue

The current experiment uses the degree-18, five-letter reference from
[F4 frame synthesis](f4-recentered-crossing-lattices.md), giving 25
physical letters and complete degree at most 36. The three constant
triple slabs from the
[saturated source module](source-slab-module-saturation.md) can each be
placed at any of the five first physical labels. They supply fifteen
binary tensors `F` or `G`.

The correction polynomials are

\[
h(w)=w^s(1+w)^r,\qquad s\ge0,\quad r\ge1,
\quad s+r\le36.                                        \tag{3}
\]

Every Lucas coefficient in (3) is retained. Powers of two are especially
sparse, but no sparse approximation replaces a different power.

The finite driver selects shifts and powers for which the added entries
do not fall below the smallest reference order at either place. This is
a search restriction, not a theorem that other corrections are useless.
It includes the maximal permitted endpoint multiplicity and specified
smaller powers of two at both ends of the allowed shift interval.

Unordered pairs `F,G` suffice for this experiment: exchanging the two
physical copies is a relabeling preserving their additive height lists.
The nonzero F4 coefficients have two Frobenius classes, represented by
`1` and `omega`. Since the base tensors and height monomials are binary,
the choices `omega` and `omega^2` have identical complete orders.

With additive paired energies, the uncoupled 25-letter control has the
same numerical ratio as the original five-letter reference,

\[
\rho\approx0.7263479477044684.
\]

The complete pairwise-sum profile identity is exact; this displayed
minimization is numerical. The half-rate threshold is
`1/sqrt(2)`, approximately `0.70710678`, so the control is not an
improvement. A successful coupled witness would still require an
independent exact profile and rational pressure certificate.

## Independent controls and evidence

The [checker](../tools/check_two_copy_crossing_couplings.py) constructs
all scalar coefficients separately by direct convolution. It checks
7,073 literal endpoint entries and 14,146 independently weighted matrix
entries on source alphabets two and three. It then verifies every order
of the selected complete profiles with 17 binary block-Toeplitz ranks
after restriction of scalars from F4. Six tensor-product profile controls,
33 literal binomial-power controls, and six invalid-polynomial controls
are also retained.

```sh
python3 tools/check_two_copy_crossing_couplings.py > /tmp/two-copy-controls.json
cmp /tmp/two-copy-controls.json evidence/two-copy-crossing-coupling-controls.json
```

The [control record](../evidence/two-copy-crossing-coupling-controls.json)
contains the complete small tensors' supports, height lists, profiles,
and constructor hashes. The
[search driver](../tools/search_two_copy_crossing_couplings.py) saves all
completed recipes, profile classes, numerical costs, unresolved queries,
and explicit completion flags. A saved finite prefix is not a completed
catalogue or a proof of a global minimum.

No new Lean formalization is used.
