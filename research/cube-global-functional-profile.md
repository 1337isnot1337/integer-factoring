# A functional profile excludes globally favorable products at every block count

For the fixed actual m=5, weight-one, two-word Gaussian coefficient family,
**no scalar reweighting of globally favorable common-orientation products
restores the endpoint, at any block count B>=1**. Coefficients may lie in
any characteristic-two extension field. This strengthens the
[single-block separator](cube-favorable-reweighting.md) by using its full
functional support profile, rather than inferring that global selection
must be favorable separately in each block.

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338808)
fixes the same base family throughout. This does not exclude larger native
word alphabets, new coefficient recombinations before orientation, other
orientation/width rules or new tensor constructions. It supplies no stronger
complete factoring exponent. The score remains **1/2**.

## 1. A product-support lemma

Let `A_1,...,A_m` be four-parent tensors in each of finitely many common
orientations h. Write `P_h(A_i)` for the reoriented tensor. Let L be a
linear functional on the fixed parent tensor space and suppose
`L(E)=1` for the endpoint E. Suppose that for every h there is an axis
`a_h != h` such that every i with `L(P_h(A_i)) != 0` satisfies

~~~
width_(a_h)(A_i) > width_h(A_i).
~~~

Then E^tensor B is not in the linear span of the product tensors

~~~
P_h(A_(i_1)) tensor ... tensor P_h(A_(i_B))
~~~

whose common hidden axis h maximizes the sum of block widths.

**Proof.** The product functional L^tensor B has value one on
E^tensor B. Its value on a displayed product is the product of the B
individual functional values. If that product is nonzero, every factor
has strictly greater width on the same retained axis a_h than on h.
Summing gives a strictly greater total on a_h. Such a product is excluded
by the hidden-maximum rule. Thus the functional vanishes on every allowed
generator, and hence on their linear span. This contradicts value one on
the endpoint. The argument works over any field where the hypotheses hold.
It permits arbitrary coefficients and arbitrarily many allowed generators.

The support premise is stronger than merely vanishing on all single-block
favorable orientations. That weaker assertion alone does not imply the
product-support lemma.

## 2. Acquire the complete support premise for the actual family

Use the 22 coefficient tensors A_i of `G_mask*G_(mask xor255)` from
the [actual-family global-pairing audit](two-axis-global-pairing-repair.md).
The three axes here have the audit's most-significant-bit-first numbering.
The fixed parent labels have two binary modes, so E=E4 on alphabet size four.
Set `L(T)=T[240]+T[255]`. Its endpoint value is one.

All 66 oriented functional values are checked, including the unfavorable
orientations. For each hidden h, precisely **one** original piece has a
nonzero value. The complete support table is:

| Hidden axis h | First retained axis a_h | Coefficient exponent | Exact original cut ranks | Widths |
| --- | --- | --- | --- | --- |
| 0 | 1 | (0,24,0) | (1,16,1) | (0,4,0) |
| 1 | 0 | (0,0,24) | (16,1,1) | (4,0,0) |
| 2 | 0 | (0,0,24) | (16,1,1) | (4,0,0) |

The nonzero functional value in every row is one. The retained-axis
width exceeds the hidden-axis width by four, so the lemma applies for
every B. Its separating functional has 2^B assignment entries; this is
a proof certificate, not a claim of constant-size explicit expansion.

The conclusion also holds if the hidden maximum is chosen by **exact
product cut ranks**, or their logarithmic widths, instead of the sums of
individual rounded widths. A product with nonzero functional consists
entirely of the unique listed extreme piece. Its cut rank is exactly
16^B on a_h and one on the other two axes. Its hidden axis therefore
cannot be a maximum under these rules either.

It is essential that the generators are literal products of these base
pieces in one common orientation. The proof does not cover choosing an
orientation after taking new linear combinations with new rank bounds.
Nor does it turn a fixed base-family obstruction into a theorem about
all word coefficients of the original Gaussian construction.

## 3. Exact controls and reproduction

The [checker](../tools/check_cube_global_functional_profile.py) first
reconstructs the original coefficients/ranks and all 66 functional values.
It then constructs every two-block truth vector on 65,536 assignments.
There are 1,452 piece-pair/orientation cases; 624 satisfy the global
maximum-width rule. Their independently eliminated binary span has
dimension 110 and excludes E^tensor2. The four-entry product functional
vanishes on all 624 generators and has value one on the endpoint.
Each of the three unrestricted fixed-orientation endpoint identities
is a positive control. The
[receipt](../evidence/cube-global-functional-profile.json) retains the
full one-block support data and all derived counts/source hashes.

Run `make verify-cube-global-functional-profile`. The all-B theorem is the
ordinary support proof above, not an extrapolation from B=1,2. These
controls do not constitute independent mathematical acceptance or new Lean.
