# Exact certificates for overlapping three-label subalphabets

Work over `F2` with five external labels. Let `E,R,C` denote the
crossing and the two other standard pairings. For a three-element
label set `A`, put

```math
H_A=(E+R+C)\,[x,y,u,v\in A].
```

This note closes all scalar entry spaces

```math
\operatorname{span}_{\mathbf F_2}
 (E,R,C,H_{A_1},\ldots,H_{A_m}),\qquad m\in\{2,3\},
```

where the sets `A_i` are distinct, for the **observable individual-word
rank spectral method**. Every observable transition table with entries
in one of these spaces has spectral base at least `13`, regardless
of its finite auxiliary dimension. Since `13^2>5^3` and
`sqrt(2)<3/2`, this exceeds the score-improvement threshold
`5^(sqrt(2))`.

The conclusion is specific to `F2` and to this spectral method.
It does not exclude direct rank savings after grouping words, and
it does not assert the same certificate over extension fields.

## 1. Six label orbits suffice

Replace each three-element set by its two-element complement.
A collection of sets becomes a simple graph on five labels, with
two or three edges. Relabelling the external alphabet is exactly
graph isomorphism, and reordering the `H_A` basis tensors changes
only the coefficient-coordinate order.

With two edges, the edges meet or are disjoint. With three edges,
the graph is a triangle, a three-edge star, a four-vertex path, or
a three-vertex path plus a disjoint edge. These are all cases:
an acyclic three-edge graph has either one connected component
with edges, giving one of the two three-edge trees, or components
with two and one edges. Three disjoint edges need six vertices.

Here are representatives, orbit sizes, and verified lower bases.
The notation `012` denotes the set `{0,1,2}`.

| Three-label sets, in coefficient order | Complement graph | Labelled families | Lower base |
| --- | --- | ---: | ---: |
| `012,013` | Adjacent edges | 30 | 15 |
| `012,234` | Disjoint edges | 15 | 16 |
| `012,013,014` | Triangle | 10 | 15 |
| `012,013,023` | Star | 20 | 13 |
| `012,013,024` | Four-vertex path | 60 | 14 |
| `012,234,013` | Three-vertex path plus edge | 30 | 14 |

The totals are `45=binom(10,2)` and `120=binom(10,3)`.
The checker independently enumerates all these families and all
120 label permutations to verify completeness.

## 2. Compact integer certificates

Write the coefficient vector as `(a,b,c,h,j)` or `(a,b,c,h,j,k)`
in the order `(E,R,C,H_(A_1),...)`. Each expression below describes
a diagonal linear detector: a term `w ell` means `w` copies of
the linear functional `ell` on the diagonal. Its rank on a tensor
is the sum of the weights whose functionals do not vanish.

| Three-label sets | Weighted functionals |
| --- | --- |
| `012,013` | `3a, 2(a+j), 2(a+h), 4(a+c), 4(a+b)` |
| `012,234` | `2a, 3(a+j), 3(a+h), 4(a+c), 4(a+b)` |
| `012,013,014` | `2a, (a+k), 2(a+j), 2(a+h), 4(a+c), 3(a+b), (a+b+k)` |
| `012,013,023` | `3a, 2(a+j), 4(a+c), 3(a+b), (a+b+k)` |
| `012,013,024` | `2a, (a+k), 3(a+j), 4(a+c), 3(a+b), (a+b+k)` |
| `012,234,013` | `a, (a+k), 3(a+j), 2(a+h), 4(a+c), 2(a+b), (a+b+k)` |

All functionals take value one at the crossing. The sums of the
weights are therefore the lower bases in the preceding table.
For each of the 31 or 63 nonzero coefficient vectors in its space,
the certificate satisfies the exact integer inequality

```math
\left(\sum_{\ell(z)\ne0} w_\ell\right)^2
 \le r(t(z))c(t(z)).
```

The [finite functional certificate theorem](finite-functional-rank-detector-certificates.md)
then gives the stated lower spectral base for every finite
observable table in that scalar entry space. The all-word endpoint
condition is essential, exactly as in that theorem.

## 3. Verification and its scope

`check_overlapping_subalphabet_detectors.py` is self-contained and
uses bitwise Gaussian elimination over `F2`. It constructs all tensor
entries directly, verifies the basis rank in each orbit, computes both
25-by-25 flattening ranks for all **314 nonzero tensor cases**, and
checks each squared integer certificate. It also checks the six-orbit
classification against all 165 labelled families.

The accompanying evidence stores every coefficient vector, both cut
ranks, and the corresponding detector rank. Verification uses no linear
programming or floating-point inequalities. Numerical optimization was
used only to propose the certificates before rounding to the displayed
small integer weights.

Repeated sets add no new scalar direction. Their cases reduce to two
distinct sets or to the previously certified one-set space. Thus allowing
repeated `A_i` cannot escape these conclusions for two or three terms.
