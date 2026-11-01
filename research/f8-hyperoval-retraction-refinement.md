# The same F8 code with a better retraction

## Exact result

Keeping the same 96-letter size, physical heights and hyperoval image,
a different left inverse gives the strict bound

```math
\rho < \frac{355002893599}{500000000000}=0.710005787198.
```

This improves the [earlier F8 certificate](f8-hyperoval-crossing-bound.md)
`0.710005835`. Both complete natural profiles improve at every positive
scale. The result is still above `1/sqrt(2)`, so the complete factoring
exponent remains `1/2`. This is not a global optimum over retractions or
arbitrary alphabet sizes, and it is not a practical factoring speedup.

## 1. Change one coordinate readout

Retain the same 93 coordinate directions and the same rank-three F8
hyperoval image on the ten height-25 triple words from the earlier note.
Let `eta` be the F8 element encoded by seven, namely `1+x+x^2` in
`F2[x]/(x^3+x+1)`.

Replace the finite hyperoval columns by

```math
(1,a+\eta,(a+\eta)^2)^{\mathsf T}\quad(a=0,1,\ldots,7),
```

and retain the two final columns `(0,1,0)` and `(0,0,1)`.
This is an invertible change of basis on the same three-dimensional image.
Read coordinate columns **7, 8, 9**, rather than 0, 8, 9.
They correspond to occupation words **265, 517, 1027** and form the
identity matrix in the new basis. Thus `LP=I` remains exact.
All changed vectors have height 25 and occupation three. Their deletions
land in the fixed lower core, so both map spaces remain deletion closed.

The physical projector changes by only rank one. With `e_j` denoting
coordinates on the ten-word boundary, put

```math
v=\sum_{j=0}^{7}e_j,
\qquad
\varphi=e_7^*+e_0^*+\eta e_8^*+\eta^2e_9^*.
```

Then `Q_new=Q_old+v*phi`, and `phi` annihilates the entire retained image.
The checker verifies all 100 boundary entries of this identity.
The lower 93 coordinate directions are unchanged. The improvement
therefore comes from the retraction, not from a larger code or alphabet.

## 2. What the completed 120-case screen establishes

All reported orders use the deflated variable `t=w^2`, as in the predecessor.
Every three-column subset of the hyperoval is invertible. Normalizing
those columns gives 120 coordinate retractions of the fixed image.
The [completed census](../evidence/f8-hyperoval-charge3-retraction-screen.json)
computes both **charge-three blocks only**, through order 63. Each block
has 932 channels. The unresolved tails remain explicit.

There are six low-profile classes, of sizes `32,52,8,8,16,4`.
The four jointly best truncated profiles use
`{6,7,8}`, `{6,7,9}`, `{6,8,9}`, and `{7,8,9}`.
Relative to the earlier `{0,8,9}` choice, they remove one zero-cut channel
at order 61, one zero-cut channel at order 63, and two infinity-cut
channels at order 63. These channels move beyond the cutoff.

This is a complete low-order census, **not** a complete-profile census
or proof of a global pressure optimum. The
[source manifest](../evidence/f8-hyperoval-retraction-census-manifest.json)
retains the semantic data hash and source fingerprints. Timing fields
are excluded from that hash. Two selected candidates were then rebuilt
at full precision; their charge-three results agree with the census.

## 3. Complete constructions and entry controls

Both `{6,7,8}` and `{7,8,9}` were reconstructed with cutoff 127.
Every one of their 18,432 orders resolved, with no tails remaining.
The complete records are
[678](../evidence/f8-hyperoval-retraction-678.json) and
[789](../evidence/f8-hyperoval-retraction-789.json).
The numerical screening values were approximately 0.7100057889695 and
0.7100057871974; only the exact certificate below is used for the result.

For each candidate, a separate subset-determinant calculation checks
166 complete polynomial entries, using 342 independent original-B
minors. The best candidate includes 17 nonbinary entries, 15 with a
nonzero third field plane. These are selected-entry controls, not a
second full matrix reconstruction. See the
[789 controls](../evidence/f8-hyperoval-retraction-789-entry-controls.json)
and [678 controls](../evidence/f8-hyperoval-retraction-678-entry-controls.json).

Compared with the earlier F8 construction, the best zero-cut profile
first improves at prefix 2,883 and gains 856 in total order. The infinity
profile first improves at prefix 3,318 and gains 626. Every prefix sum
is at least the earlier value. The best profile also dominates the
other fully reconstructed candidate in both cuts, with total gains
400 and 158.

The hinge-sum argument from the earlier note therefore proves strict
improvement of both natural partition sums at every positive scale.
This does not compare arbitrary physical gauges or certify that the
other 118 complete profiles are worse.

## 4. Exact pressure and uniform physical type

Set

```math
A=355002893599,\quad B=500000000000,\quad
z=1392830724074/1392830724075.
```

Using the complete orders, the checker proves
`4A log Q-B(log Z_0+log Z_infinity)>0`, with the same definitions as
in the earlier F8 certificate. The 192-bit lower gap numerator is

```text
17467037659011421086900957553589821363482329483680215394580
```

in units of `2^-192`. The 224-bit replay also passes.
A nearby test of `0.710005787197`, at its explicitly adjusted base,
has strictly negative gap. This is a parameter-specific negative control,
not a global lower bound on the construction.

The [pressure receipt](../evidence/f8-hyperoval-retraction-refinement-pressure.json)
includes an explicit rational physical type with denominator `2^36`
and minimum count 2,813,778. Its separate entropy inequality is strictly
positive. Thus the bound has the same fixed-seed uniform word construction
as the prior certificate, rather than relying only on a real-valued type.

## Reproduction and remaining scope

Run `make verify-f8-retraction-refinement` for source checks, independent
entries, census validation and exact pressure. Run
`make verify-f8-retraction-refinement-full` to reconstruct the best full
profile without overwriting its canonical record.
The full low-order census can be rerun with
`python3 tools/search_f8_hyperoval_retractions.py`.

Only coordinate retractions were exhaustively screened at low order.
The complete affine space of left inverses of this fixed rank-three image
has 21 field parameters and was not exhaustively searched.
The [protected-core barriers](jensen-tail-protected-core-floor.md) still
prevent this boundary-only family from crossing the factoring threshold.
No new Lean formalization or independent mathematical acceptance is claimed.
