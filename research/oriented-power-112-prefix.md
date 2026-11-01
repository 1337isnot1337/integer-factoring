# A completely replayed 112-bit oriented prefix

## Result

The fixed [optimized selector](../tools/optimized_sparse_oriented_power_search.py)
factors the certified input

$$N=3271277229293916941510651942419709
 =54488953565822267\cdot60035603828255527$$

after 73728 tested leading norms. The complete original exact-key selector
needs 96752. Both complete searches and the integer presentation
certificates are retained and independently reproduced.

The oriented search uses fewer queries but substantially more time on
this fixture. It takes 573.002 wall seconds and 259.123 CPU seconds;
the original selector takes 79.989 wall seconds and 35.585 CPU seconds.
This is a larger finite relation-supply control. It does not establish
a uniform saving or improve the current score of `1/2`.

## Fixed search and final presentation

The deterministic fresh-input generator supplies 30 recursive prime
certificates and eight rejected preceding candidates. Fixture construction
is outside the factoring algorithm; the search receives only N.
The policy has cutoff 4093, at most 131072 tested norms, checkpoints
every 8192 queries, a 512-column extraction cap, and a 1536-MiB
address-space bound. No checkpoint is skipped. The ninth checkpoint
is the first successful one.

The frozen implementation uses small-leaf product-tree quotients,
sparse expanded unit maps, and exact rational row selection. Source
hashes captured at search time agree with replay. Every checkpoint
reduces and verifies the complete signed graph, transports the larger
residual words, verifies integer unit substitutions, and tests both
relation systems. All repeated work remains charged.

| Quantity at the successful checkpoint | Value |
| --- | ---: |
| Tested leading norms | 73728 |
| Distinct residual ideal keys | 72687 |
| Pairwise coprime atoms | 75434 |
| Nonzero signed atom usages | 98572 |
| Graph-eligible rows | 71144 |
| Larger residual rows transported | 2584 |
| Additional integer unit pivots | 2505 |
| Free component generators removed from the torsion view | 2914 |
| Final core rows | 619 |
| Final core generators | 289 |
| Independent core rank | 289 |

The core consists of the original 289 base generators. All 2584 larger
residual rows pass through exact integer substitutions. Free direct
summands are removed only from the torsion view. The extractor receives
independent original rows; rational rank testing does not divide class
relations or replace them by unsupported rational combinations.

The resulting ambiguous form is
`(54488953565822267,0,60035603828255527)`. Extraction uses two saturation
steps and three kernel passes. Four tops are evaluated: two identities,
one rejected top, and the useful involution. No class-number determinant
or supplied annihilator is used.

## Complete independent comparison

The [compressed receipt](../evidence/optimized-oriented-power-112.json.xz)
contains every queried witness, every signed atom word, the complete
graph and unit maps, the input certificates, and the original selector's
successful traversal. The [summary](../evidence/optimized-oriented-power-112-summary.json)
records the measurements, checks, and source hashes.

Independent replay checks all 73728 powers using generic HNF composition
and all 72687 literal integral ideal identities. The graph certificate
checks 75723 old-generator identities and all 73728 input-relation
identities. The unit certificate independently verifies that same
original presentation. Exact rational row selection and the final
extraction are rerun.

At the identical cap, the original selector has rank 184, 1041 candidate
relations, 857 dependent relations, and no factor. Its complete traversal
finds the same proper factor at query 96752, with rank 277 and 95372
residual keys. Independent execution reproduces that entire result.

| Method | Tested norms to factor | Wall time | CPU time |
| --- | ---: | ---: | ---: |
| Oriented prefix | 73728 | 573.002 s | 259.123 s |
| Original exact-key traversal | 96752 | 79.989 s | 35.585 s |

The query ratio is about 1.31. Other research processes were active,
so CPU and elapsed times are separate. These are complete single runs.
The oriented run spends 250.276 wall seconds in extraction, including
rational row selection, and 165.979 seconds maintaining the global
atom tree. Its final product has 2743294 bits. Those costs remain in
this historical receipt; later component improvements are separate.

Post-search independent replay and serialization are outside the search
interval, as in earlier controls. Initial powers, next-power updates,
factor-base construction, rejected candidates, and every checkpoint
remain in the total. The query counter counts tested norms, using the
same convention as the original selector.

## Replay

```sh
python3 tools/check_optimized_oriented_fresh_run.py \
  --replay evidence/optimized-oriented-power-112.json.xz
```

To repeat the bounded search:

```sh
python3 tools/check_optimized_oriented_fresh_run.py \
  --bits 112 --cutoff 4093 --max-queries 131072 \
  --checkpoint-interval 8192 --max-extraction-columns 512 \
  --memory-limit-mib 1536 --output /tmp/oriented-112-rerun.json.xz
```

A stronger uniform theorem still needs a bound on how soon every input
supplies useful residual-ideal relations. This finite prefix does not
provide that missing bound.

This implementation and analysis used substantial OpenAI assistant support.
