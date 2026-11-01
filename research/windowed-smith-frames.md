# Exact finite windows for adaptive binary Smith frames

The complete-entry affine search can rebuild local frames without keeping
unbounded intermediate polynomials. The new backend returns exactly the
coefficients of the original regular frames needed by every constructed
equation. Its stopping condition checks both pivot resolution and coefficient
precision; resolving all pivots alone is insufficient.

This enables the adaptive changed-energy search. It changes neither the
feasible coefficient fiber in a specified frame nor the factoring score,
which remains **1/2**.

## Precision guarantee

The existing polynomial algorithm selects the least valuation pivot and
performs regular row and column operations using its unit part. Suppose
the same operations are performed modulo `w^P`. A pivot of order `a`
determines its unit and shifted quotient modulo `w^(P-a)`. Every affected
entry of the active source matrix has order at least `a`, so discarding
these undetermined coefficients changes that matrix only modulo `w^P`.
Previously cleared rows and columns remain zero in the required positions.
Thus its pivot path agrees with the original algorithm whenever the
selected pivots resolve below `P`.

Let `M` be the largest resolved pivot order. The accumulated regular
left and right frames agree with the original frames modulo `w^(P-M)`:
all operations have nonnegative polynomial orders, so prior errors cannot
move to a lower degree, and each new quotient loses at most `M` digits.
The argument is an induction on the actual pivot path, not an assumption
about undiscovered orders.

With row allowance `b`, the affine solver requests coefficients only
below `M+b`. Therefore the explicit condition

```math
P>2M+b
```

guarantees every requested frame coefficient is identical to the original
polynomial algorithm. Since every gauged candidate matrix is polynomial,
these low frame coefficients determine all requested coefficient equations.

The [backend](../tools/windowed_binary_smith.py) doubles `P` until all
pivots resolve and this inequality holds. A nonzero determinant has degree
at most the sum `D` of the row degree maxima, hence `M<=D`. The finite
ceiling `P=2D+b+1` is sufficient for every full-rank source. An unresolved
or singular matrix raises an error instead of supplying missing coefficients.
The returned truncated frames are also checked directly against the source.

The [windowed constructor](../tools/windowed_full_alphabet_lattices.py)
inherits the original solving, sampling and literal-witness methods.
The historical complete-affine source remains byte-for-byte unchanged,
preserving its pinned search receipts. Call
`WindowedFullAlphabetTwoCutLattice` explicitly to select the new backend.

## Exact and measured controls

Run `make verify-windowed-smith-frames`. The
[checker](../tools/check_windowed_binary_smith.py) compares 1,044 frame
and diagonal coefficients on 36 complete matrices, including dense
polynomial-unit examples with nonzero Smith orders. An independent native
profile engine checks all order lists. Six complete affine systems have
identical equations and matching feasibility results. A high-order example
forces precision retries, including a stage where all pivots resolve but
the coefficient adequacy condition still fails. A singular source is rejected.
The [receipt](../evidence/windowed-binary-smith.json) retains every input.

On the actual dense five-letter adaptive source, the backend built all
**58,175** coefficient equations with 26,250 binary coefficient variables
in about **7.66 seconds** on the current host. Both places needed precision
128, returning frame coefficients through precisions 67 and 65.
The literal reference satisfies its row bounds and both complete profiles
were reconstructed. The new negative target remains infeasible in this
particular frame; this is not general target infeasibility.

The earlier untruncated construction had not completed when it was stopped
after more than fifteen minutes. The [packet](../evidence/windowed-dense-q5.json.xz)
retains the measurements, exact source hashes, source registration and
stopping scope. This is an observed construction improvement, not a
universal runtime theorem or a completed dense old/new output comparison.
The exact precision argument and the small full-system controls justify
coefficient equivalence independently of that unfinished baseline.

Substantial OpenAI assistant support was used. No Lean formalization or
improved uniform factoring exponent is claimed.
