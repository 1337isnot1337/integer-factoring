# Fourth roots of saturated source slabs

## Exact result

The ten minimal generators in the
[saturated source module](source-slab-module-saturation.md) have only
exponents divisible by four. Their entrywise fourth roots give a new
polynomial module with a degree-reduced basis of degrees

\[
 0,0,0,2,2,2,3,4,4,4.
\]

Multiplying these generators by `1+w` gives ten endpoint-zero primitives
of degrees `1,1,1,3,3,3,4,5,5,5`. The complete degree-at-most-18 family
obtained by their polynomial combinations and free first physical label
has dimension **795** over every characteristic-two field. Its previous
uncompressed counterpart had dimension 480.

This creates a new family of candidate perturbations. It does not assert
membership in the original source operator image, a favorable cut
profile, or a stronger factoring exponent.

## Literal exponent compression

Write the previously certified generators as `G_i(w)`. Inspection of
every scalar coefficient proves that

\[
 G_i(w)=H_i(w^4)=H_i(w)^4,
 \qquad H_i\in\mathbf F_2[w]^{125}.
\]

The last equality is the Frobenius identity. The constructor obtains
`H_i` by dividing each nonzero exponent by four, and the checker inflates
all 1,250 scalar entries back to the original generators. This is an
identity of complete polynomial slabs, rather than an inference from
selected entries or sampled evaluations.

For every polynomial slab `H`, the vector `(1+w)H(w)` vanishes entrywise
at `w=1`. It can therefore be inserted at any first logical label and
added to any strict crossing on the same five-letter physical space.
The endpoint remains unchanged for every field coefficient. Its two
complete cut profiles must be recomputed after the addition.

The root operation is applied before restoring the single endpoint
factor. Taking an exponent root of `(1+w)G` itself would be invalid:
that polynomial generally has exponents not divisible by four.

## Complete polynomial module

Let

\[
 M'=\operatorname{span}_{\mathbf F_2(w)}\{H_0,\ldots,H_9\}
       \cap\mathbf F_2[w]^{125}.
\]

The exact rational row reduction has rank ten. A degree-reduced basis
of `M'` has the ten degrees displayed above. The leading coefficient
vectors are independent, so for every nonnegative integer `d`,

\[
 \dim_{\mathbf F_2}M'_{\le d}
   =\sum_{i=0}^{9}\max(0,d-d_i+1).
\]

Completeness follows from the bounded saturation lemma in the preceding
note. The largest basis degree is four. The complete polynomial space
of degree at most three has 40 pivot-coordinate variables and 21
independent equations, hence dimension 19. The displayed generators
and their fitting monomial multiples span exactly those 19 dimensions.
This excludes polynomial torsion at every degree and proves that the
displayed basis generates the entire intersection `M'`.

The endpoint-zero module is `(1+w)M'`. Free insertion at the five first
logical labels gives the following exact dimensions.

| Maximum tensor degree | 1 | 3 | 4 | 5 | 9 | 18 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Independent endpoint-zero directions | 15 | 60 | 95 | 145 | 345 | 795 |

All equations use binary coefficients. Their ranks, the independent
leading vectors, and the bounded saturation argument persist under
extension to any characteristic-two field.

## Reproduction

```sh
python3 tools/source_slab_root_directions.py evidence/source-slab-module-directions-h16.json --output evidence/source-slab-fourth-root-directions-h16.json
python3 tools/check_source_slab_root_directions.py
```

The [constructor](../tools/source_slab_root_directions.py),
[checker](../tools/check_source_slab_root_directions.py),
[complete primitive record](../evidence/source-slab-fourth-root-directions-h16.json),
and [verification record](../evidence/source-slab-fourth-root-checks.json)
retain the source hash, exponent roots, rational transformations,
polynomial basis, saturation certificate, and all 795 direction
descriptors. The checker independently reconstructs the roots and
rational rows, checks the degree-three polynomial space, verifies every
endpoint coefficient, and eliminates the complete literal dense set
of 795 directions. Two malformed root operations are rejected.

This construction and verification used substantial OpenAI assistant support.
