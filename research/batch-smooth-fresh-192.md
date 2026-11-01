# A certified 192-bit batched selective-promotion run

## Result

The predeclared batched selective-promotion search factors

$$N=5190550906325328607235154905347093439262037985110177641581$$

as

$$N=73990305065426046823116548401
       \cdot70151770583126743163701981181.$$

It returns the second prime after **48097746 generated leading norms**,
of which **46616135** are scanned. The remaining **1481611** generated
norms are charged work beyond the earlier stopping witness in their
batches. Search time is 4100.936 wall seconds and 3715.588 CPU seconds.

The complete retained certificate passes a separate saved-artifact
replay. This is a larger finite factoring control. It does not improve
the uniform logarithmic saving exponent of `1/2`, and there is no
matched unbatched 192-bit traversal from which to claim a speedup.

## Declared policy and charged work

The [certified input](../evidence/prime-cofactor-fresh-192-input.json)
and [bounded policy](../evidence/batch-smooth-fresh-192-policy.json)
were retained before the run. The input's 40 recursive prime
certificates and 138 rejected preceding candidates replay exactly.
Input-fixture generation is separate from the factoring algorithm;
the search receives only N.

| Policy parameter | Value |
| --- | ---: |
| Initial smoothness cutoff | 65521 |
| Per-generator window | 1600000 |
| Maximum generated norms | 120000000 |
| Maximum generators | 250000 |
| Promotion cutoff | 16 times the initial cutoff |
| Maximum smooth-filter batch size | 32 |
| Address-space bound | 3584 MiB |

The [batch filter](batch-smooth-filtering.md) uses exact product and
remainder operations with the classical smooth-part method. Within
each phase its smoothness product is fixed, candidates keep their
original order, and the first admissible witness is selected. The
geometric batch schedule and all generated-ahead candidates remain
charged. A state-changing accepted witness ends the frozen phase.

The complete run obeys

$$48097746=46616135+1481611.$$

Its query count includes the rightmost term. Discarding those generated
norms would understate the performed work. The retained source hashes
match all nine guarded source files after independent replay.

## Exact relation and ideal certificates

The initial base has 3299 primes. The successful result retains
98975 generators and 98968 independent augmented identities. Its
fixed-base relation rank is **3292**, so the factor appears before
full fixed-base rank. No full-rank claim is needed or made.

The [segmented compressed certificate](../evidence/batch-smooth-fresh-192-selective.json.xz.manifest.json)
and [readable summary](../evidence/batch-smooth-fresh-192-summary.json)
retain these exact checks:

| Independent saved-certificate check | Count |
| --- | ---: |
| Partial ideal identities | 98968 |
| Weighted dominance inequalities | 98968 |
| Prime residual certificates | 95676 |
| Integer DAG eliminations | 3292 |
| Negative residual orientations | 48512 |
| Identities with multiple noninitial dependencies | 35303 |
| Occurrences of promoted primes | 44715 |

The maximum definition depth is 148. The largest fixed relation entry
and elimination coefficient have 1124 bits. The largest smoothness
product has 217022 bits. These sizes and all exact integer elimination
work remain part of the computation.

The search invokes neither a factoring oracle nor a residual primality
test. Its retained prime residuals have the certificates required by
the roughness and bounded-norm rules of the existing selector.
The input's supplied prime factors are used only by fixture verification,
not by the search.

The search interval includes candidate generation, filtering, promotions,
relation construction, and extraction. Independent certificate replay
and serialization happen afterwards. Replay verifies every retained
identity and the exact final factor; it does not repeat the entire
48-million-candidate search. The source-guarded traversal records the
charged candidate counters.

The complete XZ file has 22214068 bytes. Its compressed SHA-256 is
`52354c1ea02b6e90f1ee39b7244cebe29302e0ad9df83a51b3277bcdd728b80e`.
Streaming decompression independently reproduces the uncompressed hash
`4c0b059fe944d8e8bc020a8b4b0c84894749fc0e7c3866ea10f1ba293fb76f0f`.

The archive is retained as three contiguous binary segments, each at most
8 MiB. The manifest records every segment's offset, length and SHA-256,
plus the complete archive checksum. The
[assembler](../tools/assemble_factoring_certificate.py) checks all of
these before making a completed archive available. Reconstructing into
a fresh temporary directory reproduces the exact byte-identical XZ file.

## Replay and bounded follow-on

```sh
python3 tools/assemble_factoring_certificate.py \
  evidence/batch-smooth-fresh-192-selective.json.xz.manifest.json
python3 tools/check_fresh_prime_cofactor_control.py --bits 192 \
  --replay evidence/batch-smooth-fresh-192-selective.json.xz
```

The precommitted [208-bit policy](../evidence/batch-smooth-fresh-208-policy.json)
allows its follow-on only after the 192-bit search returns a proper
factor and this saved certificate passes independent replay. Both
conditions have passed, and the follow-on has started under its
unchanged finite limits. That is an active attempt, not a completed
208-bit factorization result.

This implementation and analysis used substantial OpenAI assistant support.
