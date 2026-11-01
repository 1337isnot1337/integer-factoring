# Arbitrary Pascal exponent permutations with independent physical gauges

## Constructive search interface

The [permutation constructor](permuted-pascal-zero-spectra.md) and the
[weighted Clifford evaluator](weighted-pascal-flag-crossings.md) can be used
together outside the explicit Boolean-isotonic infinity cone. The
[search driver](../tools/search_weighted_pascal_flags.py) now accepts
`--images`, the complete exponent permutation, instead of a coordinate
weight list. Every one of the eight physical mode-height vectors is free.
An optional invertible constant flag shear is also supported. Each
objective evaluation reconstructs both complete Clifford profiles and
charges both independent matching-pair physical energies.

This is a constructive search tool and independent finite validation.
There is no new contraction certificate or factoring exponent. The
[all-scale cone comparison](permuted-pascal-cone-dominance.md) still applies
inside its stated cone; the search here does not assume that cone or its
explicit infinity profile.

## Why the existing evaluator applies

For every exponent permutation the actual matrices satisfy
`AB=(1+u)^(2L)P` and `BPA=(1+u)^(2L)I`. These are exactly the identities
used to construct the two Clifford matrices and the polynomial inverse.
The existing uniform determinant-one kernel, physical matching endpoint,
and Clifford identities therefore remain available. The evaluator reads
the actual polynomial entries and Laurent shifts. It never substitutes
the exponent list for a cut profile. Its generic Clifford rank is full
by the same inverse identity, independently of physical gauges.

Both matching-pair energies are retained even when they differ. Arbitrary
real search heights are exploratory; rationalizing them and independently
checking the resulting actual profiles and a strict pressure inequality
would be necessary for a bound claim.

## Independent original-polynomial controls

The [checker](../tools/check_permuted_pascal_free_gauges.py) constructs
seven small actual kernels, including nonmonotone exponent permutations,
independent perturbations of the mode-height lists, and two
constant-flag shears. Pure Python polynomial Smith elimination on the
original Fock matrices agrees with the reconstruction from both complete
Clifford profiles: **2,144 full Fock orders** in total. The optional native
comparison also agrees at every integer Clifford order. These controls
are exact and include negative physical heights.

```sh
python3 tools/check_permuted_pascal_free_gauges.py --native \
  --output /tmp/permuted-pascal-free-gauges.json
```

The [receipt](../evidence/permuted-pascal-free-gauges.json) retains the
physical heights, complete Clifford lists, original-Fock comparison
counts, and source hashes. Omit `--native` to run the independent
standard-library checks without compiling the native engine.

An example exploratory search is:

```sh
python3 tools/search_weighted_pascal_flags.py --images 0 1 3 2 \
  --restarts 8 --iterations 400 --seconds 120 --seed 932119 \
  --output /tmp/permuted-pascal-search.json
```

The initial eight-start control reduced the numerical pressure from
approximately `0.84145` at inherited heights to `0.80683`. This remains
above the half-score threshold, and it is neither a certified optimum
nor a new factoring bound. The next constructive question is whether
other permutations, flags, or physical alphabet restrictions improve
both complete cuts after all costs are charged.

This work reuses and credits the preceding permutation, weighted Pascal,
Clifford, and original-Fock constructors. Substantial OpenAI assistant
support was used in the implementation, computation, and text. No new
Lean formalization is claimed.
