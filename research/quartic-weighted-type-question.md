# A bounded weighted-type diagnostic for the fixed quartic interaction

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338662)
changes the physical restriction of one fixed tensor from the
[completed quartic panel](quartic-interaction-panel.md). Its complete
factoring score remains **1/2**. This document registers and closes a search, not
a pressure inequality or an improved bound.

The support `S149={X0,Y0,U0,V1}` flips the occupation parity of both
mode copies while preserving total parity. The panel's common equal-bit
frequency cannot distinguish the four physical letters fully. A
correlated four-letter type and independently weighted physical gauges
might suppress the new low-order sector while using the lifted tail.
That is the specific mechanism being tested.

The all-block Gaussian profile floor does not apply to this non-Gaussian
tensor. The six-case panel floor only covers its stated common frequency
and unweighted scale; it does not settle these types and gauges. The
reflection-coordinate floor has a different seed and inherited gauges.
The independently permuted binary-class theorem also needs a premise
check because the interaction couples the classes at a common site.
These scope differences are not a claim to escape every invariant.

## Fixed inputs, budget and stopping condition

Run [the driver](../tools/search_quartic_weighted_type.py) from a clean
frozen checkout, with an output path outside that checkout. It reads the
exact 256 entries of only `S149:x4+y4` and its Gaussian control from the
[panel receipt](../evidence/quartic-interaction-panel.json). No fields,
supports, degrees or tensor entries vary. The native F4 oracle also works
on these binary coefficients and independently checks the endpoint.

The recorded budget is sixteen starts per tensor and cut orientation,
four groups and sixty-four starts in total, with seed
`0x51554152544943`, at most 250 optimizer iterations per start and 50,000
completed native queries in total. Heights lie in `[-6,6]` and log(s)
in `[0.05,4]`. There are four natural occupation starts, four starts
with unequal copy occupations, and eight prescribed seeded
perturbations. Four independent 25-coordinate finite-difference
gradient controls are charged to the same query budget.

For the native raw log partition product `F`, let `w_a` be the summed
height of letter a on the eight end/leg combinations, with `w_0=0`.
At the critical rate `c=4/sqrt(2)`, the minimized margin is

~~~
F - c*log(sum_a exp(-w_a/c)).
~~~

The Gibbs probabilities `pi_a=exp(-w_a/c)/sum_b exp(-w_b/c)` optimize
the common four-letter type at this rate. This eliminates the type
variables from the numerical search, without assuming their two bit
frequencies are equal or independent.

The driver records exact input tensors, source commit/tree and hashes,
environment versions, every successful query's parameters, probabilities,
complete spectra and gradient, all native rejections, and per-start
termination reasons. A compressed query log has both a compressed and
uncompressed SHA-256 hash. Its polynomial arithmetic is exact, but
weights, pivot comparisons, pressures and optimization use floating
point. A finite negative is not a global all-gauge theorem.

Stop immediately on a margin below `-1e-7`. Such a numerical witness
requires exact rational weighted-minor replay, a positive rational
physical type, certified logarithmic inequalities, and the complete
uniform construction/transfer before a stronger bound can be claimed.
Otherwise stop after this one budgeted run and preserve the scoped
numerical negative. Do not add random starts to polish the decimal.
No Lean formalization or subagents are used.

## Recorded result: no strict witness, panel closed

The clean frozen source was `3ddcf36f9ab335893f0662010e5b5783340dc37f`,
tree `e119e09624810d814a9f91588f715a53c3805fee`. The
[run receipt](../evidence/quartic-weighted-run.json),
[all 5,995 query records](../evidence/quartic-weighted-run.queries.jsonl.gz)
and [complete stdout](../evidence/quartic-weighted-run.stdout.txt.xz)
retain the full experiment. All 64 starts were attempted, with 5,995
completed queries, zero native rejections and no margin below `-1e-7`.
The four finite-difference controls have maximum coordinate error below
`1.7e-9`. Wall time was 1.90 seconds on this machine.

All sixteen starts finished normally in each group except the Gaussian
control at cut one: eight finished normally, seven reached the iteration
limit and one had abnormal line-search termination. Every query remains
in the archive; numerical convergence is not asserted for those starts.
The minimum recorded margin was about `5.78570e-5`. Its heights all
equal the upper bound six and its Gibbs probability of letter zero is
about `0.9999998721`, with entropy about `2.30e-6` nats. This small
absolute margin arises near a pure-type boundary, not from a certified
near-critical pressure ratio.

Run `make verify-quartic-weighted-run` to check the archive's hashes,
source hashes, query counts, bounds, all recorded partition arithmetic
and the per-start accounting. This is an integrity check of numerical
evidence, not an independent exact weighted-Smith replay or a proof of
global optimality. The registered experiment is closed without expanding
the start/height budget. The complete factoring score remains **1/2**.
