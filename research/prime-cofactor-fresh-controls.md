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

The 160-bit evidence uses gzip with timestamp zero; the larger 176-bit
certificate uses XZ compression. Compression changes no mathematical data.
The driver reads plain JSON, `.json.gz`, and `.json.xz` replay files.
Source hashes inside the run records identify the search implementation
before these transport updates, which change only file input and output.

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

## Fresh 176-bit completion

Selective promotion also factors the separately committed 176-bit input:

```text
74210458568770352552128983761565394602512667014226413
  = 291627194339293422246369983 * 254470296355250926433520211.
```

This run uses `Y=32749`, residual bound `Y^2`, promotion cutoff `16Y`,
and a power window of 800,000. It completes after 23,998,215 norm queries
and 1,423.76 recorded search seconds. Its 1,776-prime initial base reaches
full rank. The augmented system has 51,902 generators and the same rank.
The largest definition depth is 148, the largest elimination coefficient
has 1,143 bits, and the largest saturation product has 115,022 bits.

An independent replay verifies all 51,902 partial identities and their
dominance bounds, 50,126 residual-prime certificates, and 1,776 integer
graph eliminations. It also checks 20,022 definitions with multiple
noninitial dependencies and 25,708 promoted-prime occurrences. The
37 recursive input-prime certificates and 100 rejected input candidates
are verified again. This is one completed policy run, not a matched
176-bit comparison with fixed-base traversal.

The [full 176-bit certificate](../evidence/prime-cofactor-fresh-176-selective.json.xz)
retains all witnesses. Its uncompressed SHA-256 is
`158148dddac3263de8a9e116a730fd1b169ff4608ea696676a10f01dc1a9aacb`.

```sh
python3 tools/check_fresh_prime_cofactor_control.py --bits 176 \
  --replay evidence/prime-cofactor-fresh-176-selective.json.xz
```

A fresh certified 192-bit input was also committed before its search.
Further completed results are added only after exact replay. The
established complete factoring score remains `alpha=1/2`.
