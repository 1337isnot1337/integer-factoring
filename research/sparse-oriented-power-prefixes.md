# Larger power-prefix gains from sparse oriented residuals

## Result

The [sparse selector](../tools/sparse_oriented_power_search.py) combines
incremental oriented coprime atoms with exact signed-graph reduction and
integer unit elimination. On two certified inputs it exposes useful
relations substantially earlier in the same prime-power stream than the
existing exact-key collision selector:

| Input size | Cutoff | Refined queries | Original queries to its first factor | Query ratio |
| --- | ---: | ---: | ---: | ---: |
| 64 bits | 127 | 2048 | 3983 | 1.94 |
| 80 bits | 257 | 10240 | 38818 | 3.79 |

At each refined stopping prefix, the original selector returns no factor.
All queried powers, residual ideal identities, and integer presentation
maps replay independently. The new search neither factors residual norms
nor tests their primality.

These are query and relation-supply gains. The current implementation is
slower in wall time than the complete original traversal. Its repeated
large-product updates and unit-elimination passes dominate the 80-bit
run. No stronger uniform factoring exponent follows from these examples.

## 1. Exact composition of the components

The initial prime-form base, starting exponents, round-robin ordering,
leading-norm decomposition, and canonical residual keys match
`factor_with_cofactor_collisions`. Every norm query is retained, including
smooth norms and repeated residual keys. A new residual key is inserted
once into the [incremental oriented basis](incremental-oriented-coprime-basis.md).

At each declared checkpoint, the selector exports sparse signed atom
words and substitutes them into all supplied witness rows. It then uses
the [signed residual graph theorem](signed-residual-graph-elimination.md) to eliminate
tree and unit-anchor rows. Larger residual words pass through the same
integer substitution. They are not discarded. Remaining unit residual
pivots are eliminated with explicit integer presentation maps.

The oriented refinement preserves the evaluation of every old torsion
class. The graph and unit steps preserve the full integer presentation.
Only residual columns absent from every relation are removed as free
direct summands. Finally, a maximal rationally independent subset of
genuine rows is passed to the existing relation-halving extractor. Equal
rational row spans have equal saturations, so this last step preserves
the evaluated torsion image. It does not assert equality of the integral
row lattices.

The extractor verifies its supplied relations using generic quadratic
ideal HNF composition. The search additionally checks every partial
leading-norm identity. Its saved certificate includes all original
witnesses, sparse atom words, transported rows, both complete integer
presentation maps, and the selected extraction rows.

The bounded policy can skip extraction if a compressed core exceeds its
declared column cap. It still retains the rows and records that skip.
Neither reported run skips a checkpoint: their column cap is 256.

## 2. Certified inputs and reproducible policies

The input generator uses the existing public deterministic seed
`integer-factoring-fresh-prime-cofactor-20261101`. It constructs two
half-size primes with recursive Lucas certificates. The factor-search
function receives only their product and its bounded policy. Input
construction and its certificate work are separate from the factoring
algorithm.

| Bits | Input N | Certified factors |
| --- | --- | --- |
| 64 | 14132880935215926553 | 3332910829, 4240401757 |
| 80 | 1035490702726305096363691 | 1017821676341, 1017359648351 |

The 64-bit policy uses cutoff 127, a maximum of 8192 queries, and
checkpoints every 256 queries. The 80-bit policy uses cutoff 257, a maximum
of 32768 queries, and checkpoints every 1024 queries. Both use a 256-column
extraction cap. The first experiment ran on the 64-bit input; the larger
80-bit run followed it. These finite exploratory controls are not claimed
as a blinded benchmark or a uniform success theorem.

## 3. The 64-bit prefix

After 2048 queries, the original method has 1806 residual keys and six
independent fixed-base rows. Its 242 candidate collisions include 236
dependent rows, and exact extraction finds no factor.

The same supplied keys refine to 1854 coprime atoms with 2664 nonzero
signed exponent entries. Of the 2048 transported rows, 1955 have residual
absolute weight at most two. All 93 larger words are transported too.
Graph and unit elimination yield 140 genuine rows on the 17 original
base generators; their rank is 17. Extraction returns factor 3332910829
from the ambiguous form `(3786656293,907490928,3786656293)`.

The independent checker verifies all 2048 original powers, all 1806
literal residual ideal products, 2048 input-relation identities for each
presentation map, and the final factor. It separately runs the original
selector at exactly the same query cap. It then reproduces the complete
original traversal, which first factors the input after 3983 queries.

## 4. The 80-bit prefix

After 10240 queries, the old method retains 9776 residual keys and only
one independent base row. It has 464 collision candidates, 463 of which
are dependent. Exact extraction finds no factor at this cap.

The oriented basis has 10079 atoms and 15045 nonzero signed exponent
entries. The 10240 transported rows include 809 words of residual absolute
weight greater than two. The complete graph substitution retains them.
After 764 additional unit pivots, the core has 30 rows on 22 generators:
21 base generators and one residual atom. Its rank is 22.

The retained atom has norm 24691 and canonical root 5505. Its residual
coefficients in the final rows are zero or minus two, so the unit pass
cannot remove it. Keeping this double-root information preserves the
complete presented torsion. This does not claim that the residual
generator is necessary to obtain the factor by every valid elimination.

Exact extraction returns 1017359648351 from the ambiguous form
`(508795331173,231013995,508795331173)`. The independent replay checks all
10240 literal powers, 9776 literal residual ideal products, and both
full presentation maps. The complete original traversal first obtains
the same factor after 38818 queries.

## 5. Complete measured costs and remaining work

All generated powers and all checkpoints from the start are charged.
The original comparison uses the same cutoff and exact candidate
verification. Its complete traversal is also rerun independently to
reproduce the first-factor query count and mathematical certificate.

| Input | Refined full search | Original complete search | Refined query gain |
| --- | ---: | ---: | ---: |
| 64 bits | 1.232 seconds | 0.262 seconds | 1.94 times fewer |
| 80 bits | 67.006 seconds | 2.168 seconds | 3.79 times fewer |

These are retained single-run wall times on a shared machine, not stable
cross-machine performance claims. Independent replay and final JSON
compression occur after the reported refined-search interval. During the
search, however, every graph and unit presentation is checked at every
checkpoint. The timing does not omit those certificate costs.

On the 80-bit run, incremental atom insertion takes 18.630 seconds,
graph reduction 1.607 seconds, unit elimination 45.516 seconds, sparse
transport 0.267 seconds, and extraction 0.144 seconds. The remaining time
includes factor-base construction, power generation and witness checks.
The atom tree's final product is 271098 bits. Its gcd and multiplication
operand-bit totals are both above 2.6 billion; small gcd-decision counts
alone do not describe this cost.

This profile identifies concrete next steps: avoid multiplying two large
ancestor products on every small leaf change, and avoid rescanning and
renumbering every sparse word for each unit pivot. Sparse indices and
presentation-certificate size must remain charged. None of these
engineering targets supplies the missing uniform relation-search bound.

## 6. Artifacts and replay

- [Search implementation](../tools/sparse_oriented_power_search.py).
- [Independent power, ideal and presentation checker](../tools/check_sparse_oriented_power_search.py).
- [64-bit complete compressed certificate](../evidence/sparse-oriented-power-64.json.xz)
  and [readable summary](../evidence/sparse-oriented-power-64-summary.json).
- [80-bit complete compressed certificate](../evidence/sparse-oriented-power-80.json.xz)
  and [readable summary](../evidence/sparse-oriented-power-80-summary.json).

The complete artifacts include the certified inputs and original
full-traversal comparison records. The summaries include source hashes,
all structural counts, exact factor outputs, and the measured costs.

```sh
python3 tools/check_sparse_oriented_power_search.py \
  --replay evidence/sparse-oriented-power-64.json.xz
python3 tools/check_sparse_oriented_power_search.py \
  --replay evidence/sparse-oriented-power-80.json.xz
```

To rerun the larger search and write a new complete receipt:

```sh
python3 tools/check_sparse_oriented_power_search.py \
  --bits 80 --cutoff 257 --max-queries 32768 --checkpoint-interval 1024 \
  --output /tmp/sparse-oriented-power-80-rerun.json.xz
```

This implementation and analysis used substantial OpenAI assistant support.
