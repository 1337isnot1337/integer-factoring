# F4 extension of the ordered fused coefficient envelope

This constructs the full **33-dimensional F4 scalar extension** of the
literal two-copy coefficient envelope, with exact ranks at both cuts.
It extends the binary representation from
[the symmetric envelope](fused-scalar-projection-envelope.md) and
[PR59](https://github.com/1337isnot1337/integer-factoring/pull/59) at
`4c7805d6279f4f3742a90fb0b0138c337dda031e`. The latter's binary obstruction
does not imply a field-extension obstruction. No improved factoring
exponent is claimed here; the complete score remains **1/2**.

The completed controls and newly registered larger scans have different
coverage. The large scans below are **not complete at registration**.

The **affine-zero half is now complete**: all `4^17=17,179,869,184`
coordinate positions, in 512 contiguous chunks, with no nonzero vector
of rank product at most 55. This gives rank product **at least 56**
on the specified 17-dimensional functional kernel inside `W`. It is
not a claimed sharp minimum. The affine-one half remains incomplete.

The run performed **5,726,623,061 exact native threshold tests**, skipped
**11,453,246,122** other scalar representatives and counted zero once.
The [complete packet](../evidence/fused-f4-zero-census.json.xz) retains
all chunk summaries, hashes, registration and as-run sources. The
[checker](../tools/check_fused_f4_zero_census.py) validates the entire
ledger and independently reconstructs each chunk's zero and scalar-skip
counts using [binary linear systems](../tools/fused_f4_coverage.py).
Sixty-four brute-force small cases control that separate count algorithm.
It does not independently repeat all native rank-product tests.
Run `make verify-fused-f4-zero-census`; its
[receipt](../evidence/fused-f4-zero-census-controls.json) records the
actual checks. No full-envelope, all-field or factoring improvement
follows from completing this half alone.

## Representation and genuine decomposition target

Use the 33 retained binary basis tensors, in their original order. An F4
point is `a+omega*b`, represented by two 33-bit coordinate planes, with
`omega^2=omega+1`. Constant field extension preserves linear independence.
Both 81-by-81 flattenings retain the exact 25-block charge splitting.
The [implementation](../tools/fused_f4_envelope.py) checks every basis
matrix for off-block entries before adding block ranks.

For actual pieces summing to the crossing, the
[finite-block amplification criterion](finite-block-crossing-amplification.md)
charges `S=sum_i sqrt(r_i*c_i)`. On this nine-letter alphabet, an
improvement needs `S<9^sqrt(2)`. A single low-product vector is not a
decomposition and cannot establish that target. The implementation
explicitly checks the coordinate sum of proposed decompositions.

The inherited F16 encoding is `F4[u]/(u^2+u+omega)`, with F4 embedded as
digits `0,1,2,3`. The F4 constant-coefficient projection has binary
digits `chi(z),chi(omega*z)`. It is not the ordinary field trace.
All four such projections fixing one give genuine nine-piece
decompositions of the crossing; their complete ranks are retained as
controls and their cost is much too large.

## Completed finite field controls

The [packet](../evidence/fused-f4-finite-controls.json.xz) retains:

- All 484 distinct projective points `a+omega*b` from zero, the crossing
  and the 21 cheap affine-one binary vectors in frozen PR59. Of these,
  462 are nonbinary. The minimum product remains 162; the nonbinary
  minimum is 342. All three scalar choices and their crossing complements
  are retained. The cheapest nontrivial two-piece cost is approximately
  84.7279, far above the target.
- All 61,440 points in 15 distinct affine F4 cosets of the extension of
  the six-dimensional symmetric cheap binary span, based at those 21
  binary vectors. They give 696 rank pairs and minimum product 162.
  This is a proper finite family, not the whole field extension.

Run `make verify-fused-f4-finite-controls`. The
[checker](../tools/check_fused_f4_finite_controls.py) reconstructs both
coverage ledgers, all line/complement ranks and every coset rank-class
representative with full dense elimination: 1,392 independent full
matrix ranks. It does not recompute every native coset rank.

The separate [backend checker](../tools/check_fused_f4_envelope.py)
checks all 256 two-by-two F4 matrices, 162 independent full 81-by-81
ranks, scalar invariance, all four literal projection decompositions,
the registered subspace and three complete small Gray-code scans.
The [receipt](../evidence/fused-f4-envelope-controls.json) records the
actual checks. Run `make verify-fused-f4-envelope`.

## Registered larger exact scans

Let `W` be the F4 span of the crossing and all 21 retained cheap
affine-one binary vectors. It has dimension **18**. This definition,
the exact basis, input bytes, source hashes, seed and enumeration mask
are frozen in the [manifest](../evidence/fused-f4-subspace-manifest.json).
The checker proves equality with the declared span in both directions.
Since all these vectors are binary, their binary independence also
proves their F4 independence by separating the two coefficient planes.

The retained binary functional `f` has mask 1572866 and `f(E)=1`.
Each affine half `f=1` or `f=0` of `W` has `4^17=17,179,869,184` points.
Register two threshold scans:

| Affine value | Rank-product threshold | Coverage |
| --- | ---: | --- |
| one | 501 | Every point; already normalized by its nonzero functional value |
| zero | 55 | Every coordinate position; only the unique scalar-normalized representative of a nonzero F4 orbit requires a rank computation |

In affine zero the other two scalar multiples have identical cut ranks.
They are counted separately as skips; zero itself is counted once.
The final expected skip count is `2*(4^17-1)/3`. This distinction prevents
counting scalar-normalized rank queries as if all positions were evaluated.

The [scanner](../tools/scan_fused_f4_subspace.cpp) traverses a binary
Gray code on 34 bits, XOR a registered mask. Its 17 low and 17 high
bits are exactly the free F4 coordinates. The directions are independent,
so every affine point occurs once. For each non-skipped point it adds
exact block ranks. Their partial product is nondecreasing, so a rejection
above threshold is sound. Retained vectors have both full ranks.

The [runner](../tools/run_fused_f4_subspace.py) checkpoints contiguous
chunks with hashes, validates each retained point against its Gray
position, and independently checks every new rank-pair representative
by full dense elimination. It retains as-run source texts and distinguishes
incomplete from complete summaries. Resume rejects changed sources,
binary, input or registration. The canonical binary input is
[retained](../evidence/fused-f4-subspace-input.bin) and independently
reconstructed by the focused gate.

The registration alone is no exclusion theorem. Even full completion
would concern this specified 18-dimensional F4 space, rather than the
whole 33-dimensional envelope or other coefficient fields. Any final
decomposition obstruction needs the actual resulting catalogues and
a separate argument covering arbitrary numbers of pieces.

Substantial OpenAI assistant support was used for the implementation,
controls and text. No new Lean formalization is included.
