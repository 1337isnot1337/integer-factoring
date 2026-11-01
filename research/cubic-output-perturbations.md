# A gate for changing a few cubic output rows

Allowing output monomials from other syndromes escapes the exact block
formula for [syndrome-preserving cubic outputs](cubic-radical-output-maps.md).
However, changing only a few logical output rows leaves a large unchanged
subtensor. This gives a separate all-field gate.

## Untouched-output theorem

Start with the odd-dimensional cubic input family and any permitted
syndrome-preserving output map from the linked theorem. Change the output
functionals on at most `s` logical labels, in any way that preserves the
crossing endpoint. Keep the inputs fixed. Put `m=d-s`.

For every decomposition of the resulting crossing into coefficient-space
pieces, over every field of characteristic two,

\[
S\ge\frac{m(m-1)}{\sqrt2}.
\]

**Proof.** Restrict all four logical indices to the `m` unchanged output
labels. Both cut ranks can only decrease under this restriction. The
restricted endpoint remains the crossing on those labels. Every
unordered pair of retained labels has the original two-by-two row block
and scalar column block, so the same pair charge gives
`sqrt(r c) >= sqrt(2)` times the number of active retained pairs.
Each of the `binom(m,2)` pairs must occur in at least one decomposition
piece. Summing gives the bound. The original input vectors can still
involve all `d` physical pairs; restriction concerns logical indices,
and their unchanged block formulas remain valid. In particular, `m`
need not be odd. \(\square\)

For `d=7,s=1`, the lower bound is `30/sqrt(2)`, already above
`7^sqrt(2)`. The checker records the exact integer comparison

\[
450^{35}>7^{99},\qquad 99^2>2\cdot70^2.
\]

Thus **every single-output-row change is excluded over all fields of
characteristic two**, including changes outside the old syndrome
family. A construction must change at least two output rows to escape
this particular seven-label gate. The hypothesis concerns the number
of changed rows, not the rank of a perturbation spread over many rows.

## Complete single-monomial experiment

The exploratory search fixes the cyclic coordinate outputs
`Pi_i=[e_i D_(i+1)]`. It adds one triple-single monomial to one output
row. The triple can use any three distinct pairs and any of their
eight low/high orientations. Such a correction annihilates every
input vector, so `Pi I=Id` is unchanged. It also retains the two zero
endpoint coefficients; the literal polynomial checker verifies all
endpoint entries before computing ranks.

At seven labels there are exactly

\[
7\binom73\,8=1960
\]

single corrections. All are included. Every corrected polynomial
tensor differs from the baseline. They produce 57 distinct complete
binary rank tables, all with coefficient dimension four and integer
shortest-path cost lower bound 38. No correction improves on the
baseline decomposition cost `38.7490157328`. The all-field theorem
above now explains why a strict seed is impossible even for much more
general changes to a single output row.

The experiment is exhaustive only for the stated single-monomial
corrections. It does not enumerate simultaneous changes to two or more
output rows. The source retains every correction, both cut ranks for
the best and baseline spaces, and a digest of each complete rank table.

```sh
python tools/explore_cubic_off_syndrome_outputs.py \
  --output /tmp/cubic-off-syndrome-output-corrections.json
```

- [Literal tensor and complete binary search](../tools/explore_cubic_off_syndrome_outputs.py).
- [Complete correction records](../evidence/cubic-off-syndrome-output-corrections.json).
