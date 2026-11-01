# An all-scale comparison for the six four-column orthogonal codes

## Result and scope

All six normalized F4 orthogonal scalings of the four-column rank-two
paired differential code have a strictly larger natural two-cut profile
product than the earlier coordinate-retraction code at height twenty,
at every finite positive pressure scale. This is an exact polynomial
comparison for the complete six-map class and this height truncation.

The result removes a specific candidate family from the natural-pressure
search. It does not exclude orthogonal physical maps in general, other
height cutoffs, field extensions, gauges or larger coefficient systems.
It gives no new factoring exponent.

The [existing protected-core theorem](coherent-protected-height25.md)
already places this class above the post-half threshold. The present
comparison is stronger within the finite class: all six choices are
dominated at every scale by one explicit earlier coordinate code.

## 1. The complete four-column class

Use physical pairs `(0,7),(1,6),(2,5),(3,4)`, each of height sixteen.
The four projective coefficient columns are `(1,t)` for t in F4.
Scale them by lambda, with the first scale normalized to one.

The [prior subset-Gram criterion](orthogonal-differential-pair-codes.md)
accepts eighteen nonzero scalings before this normalization. A common
nonzero scale leaves the kernel unchanged, so the six normalized cases
are

```text
(1,1,2,2), (1,1,3,3), (1,2,1,2),
(1,2,2,1), (1,3,1,3), (1,3,3,1).
```

Here 2 and 3 denote omega and omega plus one. For each scaling, form the
common differential kernel and its inverse-Gram retraction. The
[complete-profile driver](../tools/tensor_orthogonal_pair_profiles.py)
constructs the actual embedding and retraction, checks their validity,
and computes both cut profiles independently.

The retained run completes all six scalings at height cutoffs sixteen
and twenty, on twelve modes. The respective dimensions are 31 and 48.
It resolves 39,180 Smith orders in total. The two natural cuts agree in
every case, as predicted by the prior orthogonal change-of-basis theorem.

The numerical pressure minima across the six scalings are approximately
0.7122177405 at height sixteen and 0.7115540628 at height twenty. These
numbers are discovery diagnostics; the next section gives the exact
all-scale statement at height twenty.

## 2. A matched coordinate control

The control uses the unscaled coefficient matrix

```math
C=\begin{pmatrix}1&1&1&1\\0&1&\omega&1+\omega\end{pmatrix}
```

and the coordinate retraction of its differential kernel. This is the
height-twenty, two-equation control in
[the earlier F4 interior-code experiments](f4-differential-pair-codes.md).

All seven maps have exactly the same multiset of 48 physical heights.
This follows both from the uniform rank-two matroid Hilbert formula and
from the constructed bases. The comparator checks the actual multisets.

Let C_0(r), C_infinity(r) be the two control counting polynomials and
let Z_lambda(r) be the common counting polynomial of the orthogonal
map at scaling lambda. These polynomials contain every order:

```math
C_0(r)=\sum_{i=1}^{48^2}r^{a_i},\qquad
C_\infty(r)=\sum_{i=1}^{48^2}r^{b_i},\qquad
Z_\lambda(r)=\sum_{i=1}^{48^2}r^{c_{\lambda,i}}.
```

The control's two profiles were recomputed from its explicit maps with
the existing native engine. All 4,608 orders agree exactly with the pinned
earlier control. This is a fresh run of that engine, not a second profile
implementation. The orthogonal profiles come from the completed native
run described above.

## 3. Six exact positive differences

For each of the six scalings, exact integer polynomial arithmetic gives

```math
Z_\lambda(r)^2-C_0(r)C_\infty(r)
=r^{44}(1-r)Q_\lambda(r),                            \tag{1}
```

where Q_lambda has degree 162 and every one of its 163 coefficients is
strictly positive. The complete coefficient lists are retained in the
[certificate](../evidence/tensor-orthogonal-pressure-dominance.json).
The checker reconstructs (1) coefficient by coefficient.

For every `0<r<1`, (1) is strictly positive. Because all physical heights
agree, the physical entropy and mean-height costs agree for matched
types. At every fixed finite scale, replacing the orthogonal map by
the coordinate control strictly decreases the channel part of the
natural balanced two-cut pressure. Thus none of the six orthogonal maps
can improve that natural-pressure bound at any scale.

This conclusion uses the product of both complete channel polynomials.
It does not infer domination from one cut, determinant sums or fitted
numerical minima. Additional nonnatural weights are outside the claim.

## 4. Replay scope

Run the quick exact coefficient replay with

```sh
python3 tools/check_tensor_orthogonal_pressure_dominance.py
```

This reads the complete saved orthogonal profiles and the saved coordinate
control, checks their physical heights, and reconstructs all six positive
identities. It does not rebuild the Smith profiles by default.

Use `--recompute-control` to reconstruct the coordinate profiles, or
`--recompute-all` to reconstruct those profiles and all six height-twenty
orthogonal profile pairs before checking the identities. The certificate
records how many Smith orders were actually reconstructed in that run.
The initial retained certificate used `--recompute-control` and records
4,608 reconstructed orders. The
[complete orthogonal source receipt](../evidence/tensor-orthogonal-pair-profiles-small.json)
records the original 39,180 orders and the input maps.

The construction and comparison build on revision
`941cc1d8f56e6e80ec45938f1d660532de2f5046`. The previously proved F4
differential construction, orthogonal criterion and native profile
engine remain credited to their linked notes. Substantial OpenAI
assistant assistance was used for the search, exact comparison and text.
