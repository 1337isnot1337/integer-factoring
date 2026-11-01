# Low-order coherent profiles with certified tail bounds

## Result and purpose

A complete Smith profile is not required to obtain a rigorous upper bound
on its channel partition function. The new
[partial constructor](../tools/coherent_hot_f4_partial.py) computes every
natural coherent F4 Smith order at most an integer cutoff K. It separately
reports the number of remaining orders. Those orders are not claimed known.

If the full cut has n channels, the retained orders are nu_1,...,nu_r,
and `0<t<1`, then

\[
\sum_{i=1}^r t^{\nu_i}
\ \le\ Z(t)\ \le\
\sum_{i=1}^r t^{\nu_i}+(n-r)t^{K+1}.                 \tag{1}
\]

The upper expression can replace Z in a sufficient pressure certificate.
A strict inequality proved with this overestimate remains valid for the
actual crossing. The lower expression is a lower bound at that scale;
a numerical optimization of it is not a certified global pressure floor.
No improved factoring exponent is asserted here.

The original full-profile constructor and its historical source pins are
unchanged. This is a separate optional engine.

## 1. Discard entries that cannot affect the requested orders

Use the natural coherent compound matrix from
[coherent-hot-compound-retractions.md](coherent-hot-compound-retractions.md),
with homogeneous height- and occupation-preserving F4 maps. For a row and
column in one common occupation block, let H_R,H_C be their pair heights.
The corresponding entry has the form

\[
t^e P(t),\qquad e=(H_R+H_C)/2,\qquad \deg P\le e.     \tag{2}
\]

The exponent e is a nonnegative integer. The known binary minor formula
and the field-linear compression prove (2), including possible cancellation
that only raises the valuation.

If `e>K`, the entry is zero modulo `t^(K+1)` and need not be constructed.
If every entry in a charge block has this property, the entire block is
skipped and all its channels are placed in the unknown tail. Otherwise,
the engine constructs the relevant binary minors exactly and accumulates
both F4 coefficient planes. It uses enough precision to include degree K.
For every constructed minor, the degree bound in (2) prevents a hidden
higher coefficient from affecting the calculation.

Both actual cuts are treated separately. The outer constructor checks
LP=I, both deletion closures and the two homogeneous grades before using
the compound reduction. A skipped block is not called an evaluated
endpoint matrix.

## 2. DVR elimination at a finite cutoff

The [native engine](../tools/coherent_hot_f4_partial.cpp) reduces the matrix
modulo `t^(K+1)` and repeatedly chooses an entry of least finite valuation.
It normalizes its unit factor and clears its row and column by regular
invertible operations. All divisions remove only the chosen monomial
valuation and invert a unit. Each update is reduced modulo `t^(K+1)`.

These are valid DVR Smith operations. If a residual entry remains, its
least valuation is the next Smith order and is at most K. If the residual
matrix is zero modulo `t^(K+1)`, every remaining Smith order is greater
than K. No missing pivot is interpreted as a zero polynomial or as a
rank-deficient endpoint.

The full crossing has generic rank n because its endpoint is invertible.
Therefore there are exactly `n-r` remaining finite orders, each at least
`K+1`. Since powers of t decrease with the exponent, (1) follows.

The public result fields are `retained_orders_t`, `tail_count` and
`cutoff`. In particular, the engine does not label repeated cutoff values
as a complete profile.

## 3. Use in pressure searches and certificates

For a fixed rational scale, (1) is an exact rational interval for the
channel sum. It applies independently to both natural cuts, so their
upper sums can be used in the existing two-cut entropy inequality. The
physical alphabet, height distribution and endpoint tensor stay unchanged.

Increasing K can tighten the estimate. Search code can first use a small
cutoff, then increase it for candidates whose upper pressure estimate is
competitive. Any final claim must still certify the pressure sign and
physical type with outward arithmetic or exact integer comparisons.
The partial profile alone is not a new bound.

This procedure can avoid many high-order determinants and elimination
steps. It gives no universal speedup: the benefit depends on the retained
rank and the distribution of heights. Large unskipped blocks still require
substantial memory and time.

## 4. Exact controls

The [checker](../tools/check_coherent_hot_f4_partial.py) constructs complete
reference profiles for two ten-letter F4 differential kernels and a
43-letter orthogonal paired kernel. At cutoffs `0,7,20,40`, it compares
every retained order and tail count in both cuts with the complete result.
All 24 cut comparisons agree.

For each comparison, exact rational evaluation at `t=1/2` and `t=9/10`
checks both sides of (1), giving 48 rational partition comparisons. Six
invalid cutoff requests are rejected. The
[receipt](../evidence/coherent-hot-f4-partial-controls.json) records the
cases and the sources of both algorithms.
