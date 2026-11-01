# Every piece count in W18 has balanced parameter above three quarters

## Result and scope

Let `W` be the specified 18-dimensional F4 subspace of the ordered
two-copy coefficient envelope on nine physical letters, defined in the
[construction note](fused-f4-envelope.md). Every decomposition of the
crossing `E` with all pieces in `W`, after discarding zero pieces, satisfies

```math
 S=\sum_i\sqrt{r_i c_i}\ge\sqrt{162}+2\sqrt{56}>27>9^{\sqrt2}.
```

The [balanced finite-block parameter](finite-block-crossing-amplification.md)
is `rho_star=log(S)/log(81)`, so it is **strictly greater than `3/4`**
for every such decomposition, regardless of piece count. The same
piece-count refinement strengthens the existing full 33-dimensional
**binary** certificate. The full 33-dimensional **F4** space, asymmetric
criteria, higher powers and physical type restrictions remain open.
The complete factoring score remains **1/2**.

## Complete finite coverage

Frozen `e0f8d423` completed all **17,179,869,184** affine-one positions
in 512 contiguous chunks. Exactly 39 points have rank product at most
501; their minimum is exactly **162**. Native time was 15,975.44 seconds.
The [affine-one packet](../evidence/fused-f4-one-census.json.xz) retains
every summary, retained point, hash, registration and as-run source.

The [companion affine-zero census](../evidence/fused-f4-zero-census.json.xz)
already excludes every nonzero product at most 55, giving floor **56**,
without claiming this floor is sharp. It covers the same number of
coordinate positions, with 5,726,623,061 normalized native tests, two
other scalar representatives skipped per orbit, and zero counted once.

The retained functional has `f(E)=1`. Affine one gives one representative
of every projective point with nonzero `f`; affine zero supplies the
remaining `(4^17-1)/3` representatives. Together the native threshold
tests cover exactly **22,906,492,245 projective points**, equal to
`(4^18-1)/3`. Coordinate visits, native tests and scalar skips are separate
counters. This is the proper space `W`, not the full F4 envelope.

## Argument for every number of pieces

At least one nonzero piece of a decomposition has nonzero `f`. Scaling
it to `f=1` preserves both ranks, so its product is at least 162. Every
other nonzero piece has product at least 56, by the two completed scans.

For one piece the cost is 81. For two pieces, normalize one having
nonzero functional value. If its product exceeds 501, the total cost is
at least `sqrt(502)+sqrt(56)>29`. Otherwise its representative belongs
to the 39-point catalogue. Every nonzero F4 scalar is tested against
its crossing complement. All **117** resulting two-piece costs have
integer lower bound at least **84**, obtained by summing the integer
square-root floors of the two products. Both cases exceed 28.

For three or more pieces, the nonzero-functional piece costs at least
`sqrt(162)`, and two other pieces cost at least `sqrt(56)` each. This
sum lies strictly between 27 and 28:

```math
 \sqrt{162}>25/2,\quad \sqrt{56}>29/4,
 \qquad \sqrt{162}<13,\quad \sqrt{56}<15/2.
```

All four inequalities are exact integer-square comparisons. Thus the
same displayed lower bound holds in every case. Finally `sqrt(2)<3/2`
gives `9^sqrt(2)<27`. This excludes the stated balanced rank-sum target;
it does not price all asymmetric rank criteria.

## Independent checks and reproduction

The [combined checker](../tools/check_fused_f4_decomposition_obstruction.py)
checks all coverage intervals, source/input hashes, retained Gray
coordinates, span membership and functional values. It independently
computes every retained rank and all three scalar-complement ranks by
full dense F4 elimination: **312 full matrix-rank calculations**. It
also reruns the separate exact affine-zero normalization-count checker.
The [receipt](../evidence/fused-f4-decomposition-obstruction-controls.json)
retains every rank and complement.

The quick gate does not independently repeat the billions of native
threshold rejections. Their frozen scanner, inputs, registration and
complete chunk ledgers remain available for a full independent replay.
The proof combines that finite catalogue with the piece-count argument
above; it is not a new Lean formalization or independent mathematical
acceptance.

```sh
make verify-fused-f4-obstruction
```

## Binary full-space refinement

The existing [33-dimensional binary certificate](fused-ordered-scalar-envelope.md)
has the same affine-one minimum 162 and nonzero floor 56. Its listed
complements have product at least 3,969, hence each listed two-piece
cost is at least `sqrt(162)+63>28`. The outside-catalogue and three-piece
arguments above apply unchanged. Therefore its previous `sqrt(504)`
lower bound strengthens to `sqrt(162)+2sqrt(56)>27`, and its balanced
finite-block parameter also exceeds `3/4`.

The combined gate replays the old binary audit provenance and all 84
retained full ranks before applying this refinement. It does not repeat
the separate exhaustive binary scan or infer a full-space F4 theorem.
Original PR59 provenance and source history are preserved.
