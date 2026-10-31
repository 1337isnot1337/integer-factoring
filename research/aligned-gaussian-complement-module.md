# The complete polynomial module of the aligned Gaussian complements

## Result and scope

The [actual Gaussian complement construction](aligned-gaussian-complement-perturbations.md)
provides 65 endpoint-zero slabs at the retained five-letter reference.
Their rational-function span has dimension **37**. The complete
polynomial vectors in that span form a saturated module with the
following degree-reduced basis.

| Generator degree | 0 | 3 | 4 | 5 | 8 | 9 | 10 | 11 | 13 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Number of generators | 1 | 8 | 7 | 3 | 5 | 6 | 2 | 2 | 3 |

After restoring a factor `1+z` and allowing the five first-label
placements, the complete endpoint-zero polynomial family of degree
at most 18 has dimension **2,120** over every characteristic-two field.
Every direction is retained literally. This gives a larger search
space than separately dividing the 65 slabs by their scalar greatest
common divisors.

The source is aligned with the existing reference by a full coefficient
identity, `G_code(w)=w^16 T_5(w^2)`, proved in the preceding construction.
It is distinct from the earlier homogeneous weight-16 source and from
the mixed-weight neutral source. No improvement in the factoring
exponent is asserted here.

## 1. The rational and polynomial spaces

Let `S_alpha(z)` be the 125-entry complementary slab indexed by a
four-particle occupation outside the five selected code states. Its
physical coordinates are `y+5u+25v`. Define

\[
 V=\operatorname{span}_{\mathbf F_2(z)}\{S_\alpha(z)\},
 \qquad M=V\cap\mathbf F_2[z]^{125}.                    \tag{1}
\]

Exact rational row reduction gives `dim V=37`. The retained
transformation matrix reconstructs all 4,625 scalar entries of the
37-row rational echelon basis from the original slabs.

The complete bounded-degree calculation uses the pivot coordinates
of that rational basis. A polynomial vector of degree at most `d`
has polynomial pivot coordinates of degree at most `d`. Clearing the
denominator of each remaining coordinate produces exact binary
conditions: the remainder must vanish and the quotient must have
degree at most `d`. This gives every bounded polynomial vector in
(1), without imposing a denominator search limit.

Extracting successive generators from these complete spaces gives
the 37 degrees in the table. Their leading coefficient vectors are
independent, so every nonzero polynomial combination satisfies

\[
 \deg\Bigl(\sum_i a_i(z)G_i(z)\Bigr)
   =\max_{a_i\ne0}\bigl(\deg a_i+\deg G_i\bigr).        \tag{2}
\]

## 2. A finite certificate of complete saturation

The [bounded saturation lemma](source-slab-module-saturation.md)
reduces all-degree completeness to degree one below the largest
generator degree. Here that degree is 12.

At degree 12 there are `37*13=481` pivot-coordinate variables and
242 independent binary equations. Thus the complete polynomial
space has dimension 239. The displayed generators and all fitting
monomial multiples span exactly those 239 dimensions. The saturation
lemma therefore proves that these generators form a basis of the
entire module `M`, including every larger degree.

Together with (2), this gives the exact Hilbert function

\[
 \dim_{\mathbf F_2}M_{\le d}
    =\sum_{i=1}^{37}\max(0,d-\deg G_i+1).               \tag{3}
\]

All equations have binary coefficients. Their ranks and the
independence in (2) persist over every characteristic-two extension
field; the same saturation proof applies there.

## 3. The complete endpoint-zero family

The original slabs vanish at one. Polynomial saturation can introduce
vectors with nonzero endpoint, so the endpoint condition must be
restored explicitly. Its complete submodule is

\[
 M_0=(1+z)M.                                           \tag{4}
\]

One containment is immediate. Conversely, if `H` belongs to `M` and
`H(1)=0`, divide every coordinate by `1+z`. The quotient is polynomial
and stays in the rational span `V`, hence belongs to `M`.

For every physical first label `x`, insert one vector of `M_0` into
the slice with that label. The five placements have disjoint
coefficient supports. Therefore the complete number of directions
through tensor degree `D` is

\[
 5\sum_{i=1}^{37}\max(0,D-\deg G_i).                    \tag{5}
\]

At `D=18`, the sum is `5*424=2120`. The exported primitives are
`(1+z)G_i`, together with every fitting monomial shift. The checker
also verifies independence by eliminating all 2,120 literal dense
polynomial tensors, separately from the degree argument.

The 37 generator evaluations at `z=1` are independent. Nevertheless,
none of the five endpoint slabs `[v=x][y=u]` belongs to their span.
Thus the module is a complementary direction space. It does not
contain a substitute strict crossing merely by changing its endpoint
coefficients. A valid candidate has the form

\[
 T_5(z)+\Delta(z),\qquad
 \Delta(1)=0,\quad \deg\Delta\le18,                    \tag{6}
\]

and both complete cut profiles of (6) must be evaluated anew.

## 4. Reproduction and retained certificates

```sh
python3 tools/aligned_gaussian_module_directions.py evidence/aligned-gaussian-complement-directions-q5.json --output evidence/aligned-gaussian-module-directions-q5.json
python3 tools/check_aligned_gaussian_module_directions.py
```

The [constructor](../tools/aligned_gaussian_module_directions.py),
[independent checker](../tools/check_aligned_gaussian_module_directions.py),
[complete module record](../evidence/aligned-gaussian-module-directions-q5.json),
and [verification record](../evidence/aligned-gaussian-module-checks.json)
retain the frozen source hash, the reference digest, all rational
transformations, the complete polynomial basis, and every degree-18
direction descriptor. The checker replays the complete construction,
reconstructs the rational rows, tests every original and new generator
for rational-span membership, certifies degree-12 saturation, checks
the Hilbert function at additional degrees, verifies all endpoint
factors, and eliminates the full dense direction list.

The profile searches using this family are separate finite experiments.
The certificate above does not depend on a favorable search outcome or
on any floating-point optimization.

This construction and verification used substantial OpenAI assistant support.
