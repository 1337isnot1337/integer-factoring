# Fresh certified inputs for prime-cofactor relation halving

Date: 2026-11-01. Status: finite exact computations; no stronger uniform
factoring bound is claimed.

## Fresh 160-bit comparison

Both the fixed-base and selective-promotion implementations factor

```text
1179574670119658454958945004503742038025773331723
  = 1131805077794960945876933 * 1042206554168995786239631.
```

The input was generated from the recorded deterministic seed
`integer-factoring-fresh-prime-cofactor-20261101`, with the two highest
bits of each prime candidate set. The generator chooses the next odd
certified prime. Its recursive Lucas certificates and rejected-candidate
record were committed before these search results. This is a fresh
generated control, not a benchmark input chosen after seeing a successful
factorization.

The input-generation helper uses its finite control factorizer to build
prime certificates. The relation-halving search receives the integer
and uses zero calls to that factoring helper. Input generation and the
search therefore have separate computational roles.

Both runs use cutoff `Y=16381`, residual bound `Y^2`, a power window of
400,000, and an initial base of 911 primes. Selective promotion admits
new residual primes up to `16Y=262096` into subsequent trial divisions.
The input, initial base, residual bound, and traversal rule are shared.

| Measurement | Fixed base | Selective promotion |
|---|---:|---:|
| Norm queries | 7,792,462 | 5,453,561 |
| Admissible partial identities | 20,870 | 23,462 |
| New residual generators | 19,959 | 22,551 |
| Fixed-base relation rank | 911 | 911 |
| Augmented relation rank | 20,870 | 23,462 |
| Maximum definition depth | 129 | 152 |
| Maximum elimination coefficient bits | 1,019 | 1,013 |
| Maximum saturation product bits | 11,225 | 56,550 |
| Recorded search seconds | 221.26 | 356.71 |

Selective promotion uses about 30.01% fewer norm queries, but takes
about 1.61 times as long in this shared-machine run. It creates more
dependencies, more trial divisions, and substantially larger saturation
products. These measurements support a query-count improvement on this
input, not a wall-clock improvement or a uniform asymptotic claim.

## Exact certificates and replay

The [prime-cofactor saturation theorem](prime-cofactor-saturation.md)
supplies the algebraic justification. Rough residuals below `Y^2` are
certified prime by their lack of factors up to `Y`; the search makes no
residual primality-test calls. The retained dependency graph is
eliminated over the integers, preserving the exact relation identities
and saturation exponents.

Independent replays check all 20,870 fixed-base partial identities and
all 23,462 selective identities, their weighted dominance bounds,
19,959 and 22,551 residual-prime certificates, and 911 integer graph
eliminations for each run. The selective replay additionally checks
8,715 definitions with multiple noninitial dependencies and 11,147
promoted-prime occurrences. Both extracted divisors divide the original
integer exactly.

The full evidence is compressed deterministically with gzip timestamp
zero. Compression changes no mathematical data. The driver reads both
plain JSON and `.json.gz` replay files. Source hashes inside the run
records identify the search implementation before gzip transport support
was added; that update changes only file input and output.

- [Certified input](../evidence/prime-cofactor-fresh-160-input.json)
- [Fixed-base certificate](../evidence/prime-cofactor-fresh-160-fixed.json.gz)
- [Selective-promotion certificate](../evidence/prime-cofactor-fresh-160-selective.json.gz)
- [Comparison and uncompressed hashes](../evidence/prime-cofactor-fresh-comparison.json)
- [Generation, search, and replay driver](../tools/check_fresh_prime_cofactor_control.py)

```sh
python3 tools/check_fresh_prime_cofactor_control.py \
  --replay evidence/prime-cofactor-fresh-160-fixed.json.gz
python3 tools/check_fresh_prime_cofactor_control.py \
  --replay evidence/prime-cofactor-fresh-160-selective.json.gz
```

The separately committed 176-bit input is a further test target. A
completed result is added only after its search and exact replay finish.
The established complete factoring score remains `alpha=1/2`.
