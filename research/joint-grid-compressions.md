# Joint pair compression and a new physical grid

A crossing can be compressed on a pair of legs jointly and then re-encoded
as a different Cartesian alphabet. The joint maps need not factor into
separate maps on the old physical legs. This enlarges the existing
[physical split construction](physical-split-crossings.md) and permits
rectangular maps beyond the invertible
[joint conjugation family](joint-height-conjugation-floor.md).

Two complete coordinate families of the four-letter source have exact
obstructions. All 1,232 declared vacuum-containing `2 by 3` grids have
balanced natural pressure `rho >= 3/4`. All 1,152 declared full `4 by 4`
grids have no better pressure than the original source at any scale and
any two supporting weights. These conclusions extend to every
characteristic-two field for the same coordinate maps. They do not cover
general field-valued joint projections, different source tensors, changed
height conventions, or joint spaces omitting the vacuum. The factoring
score remains **1/2**.

## 1. The endpoint needs a joint split identity

Write the original crossing as the operator

\[
A(z)_{(X,Y),(V,U)}=T(z;X,Y,U,V),\qquad A(1)=I.
\]

Choose any constant joint embedding `P` and joint retraction `L` with
`LP=I_k`, where `k=ab`. Label the new basis by `(x,y)` in a Cartesian
`a by b` grid and put

\[
T'(z;x,y,u,v)=(LA(z)P)_{(x,y),(v,u)}.              \tag{1}
\]

Then `T'(1;x,y,u,v)=[x=v][y=u]` exactly. The maps can mix the original
two legs and need not commute with swap. The `(V,U)` order in (1) is
essential. In the old `(U,V)` order it includes an explicit swap before
the embedding. Polynomiality and effective construction are preserved
for fixed finite maps. A new local tensor must still pay its complete
construction cost.

Suppose a label `(x,y)` has total height `h_X(x)+h_Y(y)` and both maps
preserve this total. The zero regularization compresses by `L` and `P`
at the inherited heights. Its complete Smith profile must be computed
for a rectangular compression. The opposite cut is the actual reshuffle
of the new tensor in (1); joint matrix multiplication is not a row or
column equivalence after reshuffling. We retain that cut in full.

## 2. The declared coordinate domain is finite and complete

Use the literal degree-twelve four-letter reflection source with heights
`h=(0,3,1,4)`. The [source](../tools/check_triangular_auxiliary_row_fiber.py)
already supplies its full polynomial tensor. A coordinate embedding is
an injection

\[
f:(x,y)\longmapsto(X,Y),\qquad
h_X(x)+h_Y(y)=h(X)+h(Y).                            \tag{2}
\]

The retraction extracts the selected coordinates. We require nonnegative
sorted height lists beginning with one zero, and preserve the original
vacuum. Sorting loses no pressure possibilities because it merely
relabels symbols. Equation (2) bounds each individual height by eight.
Repeated positive heights are allowed.

The ambient height-fiber sizes at heights `0,...,8` are
`1,2,1,2,4,2,1,2,1`. If the grid needs `n_j` states in fiber `j`, there
are exactly `prod_j d_j!/(d_j-n_j)!` injections, provided `n_j<=d_j`.
Enumerating the sorted height lists and these ordered injections gives:

| Grid | Valid height pairs | Joint coordinate maps | Profile and height classes |
|---|---:|---:|---:|
| 2 by 3 | 58 | 1,232 | 214 |
| 4 by 4 | 3 | 1,152 | 23 |

The [constructor](../tools/joint_grid_compressions.py) checks every
endpoint entry before computing both complete profiles. It uses exactly
the original regularization: multiply the zero row by
`z^(h_X(x)+h_Y(y)-4)`, and reverse at infinity with degree
`h_X(v)+h_Y(u)+4`. All negative powers cancel entrywise. The opposite
matrix has row `(x,u)` and column `(v,y)`.

The original [joint-height source](joint-height-conjugation-floor.md)
is credited for the literal seed and regularization. The rectangular
split and new-grid enumeration are the additional construction here.
No reduced Gaussian block is substituted for the literal opposite cut.

## 3. Every declared 2-by-3 grid has balanced rho at least 3/4

For one profile class define `Z_0(x)=sum_i x^nu_i`,
`Z_infinity(x)=sum_i x^mu_i`, `F_X(x)=sum_a x^h_X(a)`, and similarly
`F_Y`. At `rho=3/4`, the balanced natural criterion is equivalent,
after writing `x=exp(-tau/3)`, to

\[
D(x)=\bigl(Z_0(x^3)Z_\infty(x^3)\bigr)^2
          -\bigl(F_X(x^4)F_Y(x^4)\bigr)^3\ge0.      \tag{3}
\]

All 214 integer polynomials in (3) have exact positive certificates.
Remove their leading monomial and expand at degree `N` as

\[
x^{-v}D(x)=\sum_{i=0}^N\beta_i x^i(1-x)^{N-i}.
\]

For a monomial coefficient `c_j`, its contribution to `beta_i` is
`c_j binom(N-j,i-j)`. If required, degree elevation replaces adjacent
coefficients by their sums. Every final integer `beta_i` is nonnegative,
and both end coefficients are positive. The largest certificate degree
is 216. Thus (3) is strictly positive for `0<x<1`.

The Gibbs entropy expression on the physical side increases with `rho`,
so every balanced pressure root is at least `3/4`. This is an all-scale
certificate for the complete declared coordinate family. The numerical
discovery minimum, about `0.7551305316`, is not used in the proof.
Arbitrary asymmetric supporting weights are outside this 2-by-3 result.

## 4. Every full 4-by-4 grid is dominated by the original source

Let `b=(0,1,3,4,5,6,7,8,8,9,10,11,12,13,15,16)`.
For every one of the 23 profile classes, the zero profile is exactly
`b`. If the opposite profile is `c` in increasing order, then

\[
\sum_{i=1}^j c_i\le\sum_{i=1}^j b_i
\qquad(1\le j\le16).                              \tag{4}
\]

The standard positive-part integral argument in the
[joint-height note](joint-height-conjugation-floor.md) turns (4) into
`sum exp(-t c_i) >= sum exp(-t b_i)` for every `t>0`.
Moreover all sixteen ambient coordinates are retained, so

\[
F_X(x)F_Y(x)=\left(\sum_{a=0}^3x^{h(a)}\right)^2.
\]

The physical partition product is therefore unchanged even for the two
nonstandard factorizations into height lists. Each weighted cut
contribution is no smaller at every scale. This proves no improvement
over the original source for every choice of the two positive supporting
weights and optimized physical types. It is a relative dominance result,
not a classification of all joint linear maps.

## 5. Independent profile replay and reproduction

The discovery receipts are
[2 by 3](../evidence/joint-grid-2x3.json) and
[4 by 4](../evidence/joint-grid-4x4.json). They retain every injection,
complete profiles, height lists, and source hashes. Numerical fits are
discovery diagnostics, not parts of either proof.

The [checker](../tools/check_joint_grid_compressions.py) reconstructs each
tensor from the literal source coefficients independently of the search
constructor's restored matrices. It repeats both weighted profiles using
the separate native engine: **51,648 orders** and **339,264 strict endpoint
entries**. The native calculation uses the proved full rank, integer
weights and a 512-coefficient window; every promised pivot must resolve.
The Python fraction-free discovery and native replay agree exactly.
Per-height injection counts prove completeness by the factorial formula,
and duplicate maps are rejected. The
[control receipt](../evidence/joint-grid-compression-controls.json) also
records every Bernstein certificate and all source/evidence hashes.

```sh
make verify-joint-grid-compressions
python3 tools/search_joint_grid_compressions.py --x 2 --y 3 --max-cases 2000 --seconds 900 --output /tmp/joint-grid-2x3.json
python3 tools/search_joint_grid_compressions.py --x 4 --y 4 --max-cases 2000 --seconds 900 --output /tmp/joint-grid-4x4.json
```

All tensor coefficients in these coordinate families are binary, so rank
and profile conclusions survive every characteristic-two extension.
This does not enumerate field-valued joint embeddings and retractions.
Checks are local author checks. No independent review, new Lean theorem
or improved factoring exponent is claimed.
