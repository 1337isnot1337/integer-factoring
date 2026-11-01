# A bounded weighted-type diagnostic for the fixed quartic interaction

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338662)
changes the physical restriction of one fixed tensor from the
[completed quartic panel](quartic-interaction-panel.md). Its complete
factoring score remains **1/2**. This document registers a search, not
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
