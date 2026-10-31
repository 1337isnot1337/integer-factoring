# Matching-value curves that escape weaker certificates

The [binary site-symmetry theorem](binary-site-symmetry-barrier.md)
rules out post-half improvement by both curves below. Their role is to
separate weaker algebraic certificates from the universal degree-four
certificate. Neither finite experiment gives a stronger factoring bound.

## Explicit construction

Over `F_p`, choose a complete rank-one row matching or column matching
at every non-endpoint projective point, and the strict crossing `E`
at point one. Each matching is an outer product of two invertible
two-by-two matrices in its low-rank cut, so its cut ranks are `(1,4)`
or `(4,1)`. The constructor uses a fixed random seed only to select
these invertible local matrices; all subsequent arithmetic is exact.

The degree-`p` interpolation basis is

```math
b_a(z)=1-(z-a)^{p-1}\quad(a\in F_p),\qquad
b_\infty(z)=z^p-z.
```

Finite point evaluations of `b_a` are Kronecker deltas, its leading
degree-`p` coefficient is zero, and `b_infinity` vanishes at every
finite point with leading coefficient one. Thus the tensor polynomial
has exactly the prescribed projective values. The full coefficients
and values are saved, rather than merely their pointwise ranks.

## Degree five: no common projected cofactor over any extension

The [degree-five seed](../evidence/degree-five-matching-candidate.json)
over `F_5` has row matchings at `0,2,3` and column matchings at
`4,infinity`. Its entry span has dimension six. The two degree-20
determinants are coprime, including at infinity, and both possible
constant partial-transpose intertwiner systems have full constraint
rank four.

To search for a common projected three-by-three minor, write it as
a weighted sum of the 16 cofactors in each cut. Each four-by-four
weight matrix must have rank at most one. Equality of the two
polynomials, together with value one at the endpoint, imposes 17
independent linear conditions on the 32 weights. Exact elimination
leaves 15 affine variables. The two weight matrices give 72 quadratic
two-by-two minor equations `G_i=0`.

The saved certificate specifies 1,452 multiplier monomials and verifies
the identity

```math
\sum_{i=1}^{72}Q_iG_i=1\quad\text{in }F_5[z_1,\ldots,z_{15}],
\qquad \deg Q_i\le2.
```

This identity proves nonexistence over every extension of `F_5`,
including the algebraic closure. It is stronger than enumerating
projection vectors over the base field. Singular 4.4.1 generated the
certificate with `liftstd`; a Python standard-library checker expands
the polynomial identity independently.

Yet the unrestricted cofactor spans have dimensions 10 and 12 and
intersection dimension six. A saved common polynomial with value one
at one has coefficients

```text
0 0 0 3 2 2 3 3 1 2 1 3 1 0 0 0.
```

Its two weight matrices have ranks three and two. Frobenius powers of
this polynomial, although not single projected minors, belong to
coefficient spaces of suitable Weyl quotients. The checker verifies
the polynomial identities at powers `1,5,25` and exact rectangular
multiplicities at powers through `625`. This is a useful intermediate
certificate; the universal degree-four theorem subsumes it here.

## Degree eleven: even the cofactor spans are disjoint

The [degree-eleven seed](../evidence/degree-eleven-matching-candidate.json)
over `F_11` has row matchings at `0,2,3,4,5,6`, column matchings at
`7,8,9,10,infinity`, and entry span dimension 12. Its cut determinants
are coprime. The combined cofactor coefficient matrix has shape
`34 by 32` and full column rank 32. The certificate stores 32 selected
rows and the inverse of that square submatrix. Both inverse products,
2,048 entries in total, are checked directly.

Thus there is no nonzero common cofactor polynomial over any field
extension. Nevertheless the universal degree-four polynomial from
the site-symmetry theorem still exists: it multiplies matrix entries
by cofactors. Disjoint cofactor spans are therefore insufficient as
a criterion for a viable scalar seed.

## Finite actual jets and reproduction

The complete cut ranks of all local jets through three tensor copies
are saved for both curves. Exact Hermite reconstruction checks all
`16+256+4096=4368` crossing entries per curve.

| Field and degree | One copy: max rank product | Two copies | Three copies | Schedule coverage |
| --- | ---: | ---: | ---: | --- |
| `F_5`, degree 5 | 12 | 225 | 3584 | Every minimal Hermite allocation: 210, 1365, 4845 |
| `F_11`, degree 11 | 12 | 165 | 2688 | One prefix-minimax DP schedule at each length |

The degree-eleven dynamic program minimizes a prefix rank bound.
It is not an exhaustive optimization after zero interpolation weights
are removed. These values have no asserted uniform rank upper bound.

Run the fast saved-certificate checks with

```sh
make verify-matching-curve-controls
```

To regenerate the Nullstellensatz certificate, with Singular installed:

```sh
python3 tools/check_degree_five_matching_candidate.py --write-singular /tmp/matching-five.sing
Singular -q /tmp/matching-five.sing > /tmp/matching-five-cas.log
python3 tools/check_degree_five_matching_candidate.py \
  --import-certificate /tmp/matching-five-cas.log \
  --candidate-output /tmp/matching-five-replay.json
cmp evidence/degree-five-matching-candidate.json /tmp/matching-five-replay.json
```

The generation and independent certificate check were both run, and the
reproduced candidate was byte-identical. Singular is optional for normal
verification of the stored polynomial identity. Recompute the finite
jet experiments with `make verify-matching-curve-jets`; the exact data
are [degree five](../evidence/degree-five-matching-jets.json) and
[degree eleven](../evidence/degree-eleven-matching-jets.json).

Implementation and verification: Codex, with substantial AI assistance.
The representation-theoretic inputs and their sources are identified
in the site-symmetry proof; no independent review is claimed.
