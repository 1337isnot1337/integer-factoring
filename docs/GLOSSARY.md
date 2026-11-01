# Glossary

These terms introduce the half-rate proof and the research notes linked below.
For the current proof status, use the [README result table](../README.md#results-and-evidence).

## Input bit length

For an integer `N >= 2`, its input bit length is the number `n` of binary
digits needed to write it, so `n = floor(log2 N) + 1`.
Runtime bounds are functions of `n`, not of the numerical value `N`.
See the [problem statement](../README.md) and the
[half-rate theorem](../research/half-rate-factoring.md#theorem-and-evidence-level).

## Logarithmic-saving exponent alpha

In a bound `T(n) = exp(O(n / (log n)^alpha))`, `alpha` measures the
logarithmic saving inside the exponent; a larger `alpha` gives a smaller
asymptotic upper bound of this form.
The half-rate result has `alpha = 1/2`, meaning a square-root logarithm in
the denominator, not a halving of practical runtime.
See the [half-rate theorem](../research/half-rate-factoring.md#theorem-and-evidence-level)
and the [general rate formula](../research/pieces-to-factoring.md#1-the-precise-hypothesis).

## Crossing piece

A crossing piece is one four-label array in a sum that reconstructs the
crossing constraint `[x = v][y = u]` (each bracket is 1 when its equality
holds and 0 otherwise), possibly multiplied by a common nonzero scale.
Each piece must also have two explicit factorizations, corresponding to
the row and column cuts; the normalized modular construction calls its
pieces `A_s(x,y,u,v)` and has scale one.
See the [crossing-family requirements](../research/pieces-to-factoring.md#1-the-precise-hypothesis)
and [modular crossing identity](../research/modular-crossing.md#3-exact-crossing-and-its-two-factorizations).

## Rank width

For a factorization with `R` intermediate channels, its binary width is
`ceil(log2(max(1,R)))`: this counts bits needed to label channels, rather
than the channels themselves.
The transfer bounds the sum `c_i^row + c_i^col` of the two factor widths
for each piece; the modular note instead writes channel counts `R_s,C_s`
and takes their rounded logarithms, and these counts can exceed the actual
matrix ranks.
See the [width-sum hypothesis](../research/pieces-to-factoring.md#1-the-precise-hypothesis)
and [modular width estimate, equation (14)](../research/modular-crossing.md#5-what-the-state-counts-prove-about-contraction).

## Uniform construction

A uniform construction uses fixed programs to generate the required
objects for every allowed input or label width, with the stated computation
and encoding costs included.
It cannot replace those programs with a separate unexplained choice for
each size, or assume free access to tables that have not been written.
See the [effective-family hypothesis](../research/pieces-to-factoring.md#1-the-precise-hypothesis)
and the [fixed evaluator and its input costs](../research/half-rate-factoring.md#1-the-prior-evaluator-with-its-actual-input-charged).

## Finite certificate

A finite certificate records exact, reproducible evidence for specified
finite cases, such as the crossing identities and ranks checked at selected
widths.
It supports those cases and does not by itself prove a uniform construction
or an improved asymptotic factoring bound.
See the [modular exact checks](../research/modular-crossing.md#8-exact-checks-and-the-remaining-target)
and the [README evidence distinctions](../README.md#results-and-evidence).

## Conditional transfer

A conditional transfer proves that an algorithmic conclusion follows if
an explicitly stated family of constructions satisfies all its hypotheses.
Here an effective crossing family with fixed rational `rho` in `[1/2,1)`
gives `alpha = log2(1/rho)`, but the implication does not supply that family:
the crossing identities, factorization widths, coefficient bounds, and
uniform program costs still have to be established.
See the [transfer hypothesis and theorem](../research/pieces-to-factoring.md#1-the-precise-hypothesis)
and its [conditional targets](../research/pieces-to-factoring.md#6-exact-targets-for-a-new-crossing-construction).

Notation across these notes: `lambda` is a local binary label width, distinct
from the factoring input length `n`; the modular note has `m = 2^lambda`
labels, partition weight `h`, and fixed prime characteristic `p`.
Its index `s` labels pieces, while the general transfer uses `i`;
`rho` controls the width bound and `alpha` describes the resulting runtime.
See the [modular label convention](../research/modular-crossing.md#3-exact-crossing-and-its-two-factorizations)
and the [transfer notation](../research/pieces-to-factoring.md#1-the-precise-hypothesis).
