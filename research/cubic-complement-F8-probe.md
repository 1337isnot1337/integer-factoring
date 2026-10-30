# A bounded non-diagonal cubic probe over F8

The order-seven cyclic matrices have all their eigenvalues in `F8`.
This makes that field a useful separate test of the five-direction
[complementary-output construction](cubic-complementary-output-maps.md).
The scalar and diagonal cases are already subject to the
[syndrome sumset bound](syndrome-sumset-rank-cost.md); the new probe
therefore concentrates on non-diagonal complementary images.

The complete bounded experiment finds no strict seed. Every one of its
428 coefficient spaces has an affordable span missing the crossing.
The result concerns these maps over `F8`, not arbitrary output maps,
larger fields, tensor powers, or new recombinations across different maps.

## Field and output maps

Use

\[
\mathbb F_8=\mathbb F_2[z]/(z^3+z+1),
\]

with the binary integer codes `0,...,7` for polynomial coefficients.
The element `z`, encoded by 2, has multiplicative order seven.
For each recorded matrix `C`, the output map is

\[
\Pi_i=[e_iD_{i+1}]+\sum_j C_{ij}[f_jD_{j+1}],\qquad i,j\pmod7.
\]

This gives `Pi I=I_7` and `Pi J=C`. The crossing endpoint is checked
exactly. The five raw coefficient coordinates are
`(r^3,p r^2,p^2 r,q r^2,q^3)`.

## Exact bounded scope

| Complementary image family | Maps |
| --- | ---: |
| Scaled order-seven spectral projectors | 378 |
| Cyclic shift plus a scalar identity | 8 |
| Nilpotent Jordan shift plus a scalar identity | 8 |
| Seeded matrices of rank at most one | 16 |
| Seeded matrices of rank at most two | 16 |
| Zero and identity controls | 2 |
| **Total** | **428** |

For each nonempty proper subset `S` of the seven frequencies, the
spectral family uses

\[
C_{ij}=\lambda\sum_{k\in S}z^{-k(j-i)},
\qquad \lambda\in\{1,z,z^2\}.
\]

There are `126*3=378` such maps. The unscaled matrix is the cyclic
spectral projector onto `S`: the usual sum of seven roots of unity
works in characteristic two because the field element represented by
the integer seven is one. The shift and Jordan families each use all
eight scalars. The rank-one and rank-two families use the recorded seed
`20261030`; they are not exhaustive over matrices of those ranks.

## Full coefficient enumeration

For each map, the native evaluator enumerates all

\[
1+8+8^2+8^3+8^4=4681
\]

projective vectors in the five-coordinate parametrization. The total
is **2,003,468** vectors. The actual tensor space may have smaller
dimension for a control map; the formal parametrization is still
exhaustive. Zero tensors from its coefficient kernel are retained with
zero cost, so they cannot invalidate the span test.

The evaluator accepts every vector with `r(T)c(T)<256`. This is a
deliberately generous individual-piece gate: `16^70>7^99` and
`99^2>2*70^2` imply `16>7^sqrt(2)`. Every piece in a strict seed
decomposition would have to pass it.

Across all maps there are 4,376 affordable projective coefficient
vectors. Their spans have dimension three in 419 cases and dimension
four in nine cases. The crossing coordinate is outside the affordable
span in all 428 cases. This excludes every `F8` decomposition within
each recorded coefficient space.

## Verification

The C++ evaluator stores field rows in three bit planes. It computes
the row rank exactly, then stops the column elimination only after a
proven lower bound already exceeds the cost gate. Every reported
affordable vector therefore has both exact ranks.

The Python audit uses the opposite pivot order and recomputes every
affordable vector plus rejected controls: **5,658** vectors in total.
A separate checker recomputes all 428 affordable-span exclusions,
directly expands **72,030** physical polynomial coefficient entries
across six control maps, and reproduces six complete native screens
byte for byte. The source compiles with C++20 and
`-Wall -Wextra -Werror`.

```sh
c++ -std=c++20 -O3 -Wall -Wextra -Werror \
  tools/screen_cubic_gf8.cpp -o /tmp/cubic-gf8-screen
python tools/explore_cubic_complement_gf8.py \
  --native /tmp/cubic-gf8-screen \
  --output /tmp/cubic-complement-F8.json
python tools/check_cubic_complement_gf8.py \
  --input evidence/cubic-complement-F8.json \
  --native /tmp/cubic-gf8-screen \
  --output /tmp/cubic-complement-F8-audit.json
```

The explorer also supports `--replay evidence/cubic-complement-F8.json`
to regenerate every tensor basis and repeat every independent rank
control without rerunning the native enumeration.

- [Native `F8` rank evaluator](../tools/screen_cubic_gf8.cpp).
- [Construction, enumeration, and rank replay](../tools/explore_cubic_complement_gf8.py).
- [Physical expansion and archive checker](../tools/check_cubic_complement_gf8.py).
- [All maps and complete native screens](../evidence/cubic-complement-F8.json).
- [Independent audit certificate](../evidence/cubic-complement-F8-audit.json).
