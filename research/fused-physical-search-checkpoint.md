# Completed fused physical-basis search

The declared two-hour search on the three- and four-letter fused seeds
is complete. It did not reach the factoring threshold. Its best
recorded balanced pressures are:

| Physical alphabet | Fusion parameter m | theta | Best recorded rho |
| --- | --- | --- | --- |
| 3 | 2 | 2 | 0.8915621555195623 |
| 4 | 3 | 4 | 0.9647087871380579 |

These values concern the searched fused seeds, physical bases, gauges,
common paired types, and place allocations. They are not global minima
and do not improve the repository's best factoring exponent.

## 1. Search and time accounting

The source is
[`search_fused_physical_bases.py`](../tools/search_fused_physical_bases.py).
It uses the full exact fused coefficient matrix, applies independent
constant physical basis changes on the two paired legs, checks the
literal crossing endpoint, and computes complete weighted local Smith
profiles. The degree and common physical type means are included in
every reported pressure.

The proposal families include identity, ordinary and quantum Pascal,
inverse Pascal, ordinary and geometric evaluation bases, and elementary
row mutations. Three retained place partitions across the two alphabets supply the
starting allocations. Full physical-letter heights, common paired
types, and allocations are trained within the declared numerical
search bounds. Completed optimizer calls do not imply global
optimization.

The environment interrupted the original worker. Its atomic prefix was
preserved, and
[`continue_interrupted_worker.py`](../tools/continue_interrupted_worker.py)
ran only the unspent declared budget, with a separately recorded random
seed. The harness continued proposal numbering and carried the retained
fused winners. It did not claim to recover the lost in-memory random
generator state.

| Segment | Saved wall seconds | Attempted trials | Completed trials |
| --- | ---: | ---: | ---: |
| Preserved prefix | 2678.168740536 | 563 | 563 |
| Continuation | 4522.019976077 | 486 | 485 |
| Combined | 7200.188716613 | 1049 | 1048 |

The final trial reached the declared wall limit. The count excludes
unsaved work after the original prefix, and the continuation does not
double-count the prefix time. The combined report retains every saved
trial number and outcome.

## 2. Archived inputs and exact winner replay

The [compact summary](../evidence/fused-physical-basis-search-summary.json)
contains all trial outcomes, exact physical basis matrices for both
winners, all winning numerical parameters, complete local records,
source hashes, and explicit interruption lineage. The complete reports
and original baseline are retained in deterministic gzip archives:

- [Prefix](../evidence/fused-physical-basis-prefix.json.gz)
- [Continuation](../evidence/fused-physical-basis-continuation.json.gz)
- [Baseline](../evidence/fused-physical-basis-baseline.json.gz)

The checkpoint tool
[`checkpoint_fused_physical_search.py`](../tools/checkpoint_fused_physical_search.py)
checks the archive lineage, contiguous proposal numbering, all frozen
source hashes, and complete budget accounting. It then rebuilds both
winning physical tensors and checks all 337 endpoint entries.

For each winning local profile, it interprets the saved heights and
outer scale as exact decimal rationals. It normalizes the group
allocations exactly and derives each local scale from the common outer
scale and exact per-root allocation. It also checks every saved
allocation and scale relation before this normalization. The independent fraction-free
bordered-minor oracle computes every Smith order over F16. The resulting
11 complete profiles contain 148 channels. Exact normalized rational
paired types supply their type shifts. The final logarithmic pressure
is recomputed numerically and agrees with the recorded values within
2e-8. This is an exact local-order replay with a floating final cost,
not an interval threshold certificate.

The [replay record](../evidence/fused-physical-basis-search-replay.json)
contains every exact rational order and type mean. The full native
search records remain available for later proposal selection, without
requiring the search to be repeated.

```sh
python tools/checkpoint_fused_physical_search.py
```

The default command reads the three committed gzip archives directly.
Explicit `--prefix`, `--continuation`, and `--baseline` arguments also
accept their uncompressed forms. The historical absolute baseline path
is retained as provenance and is not required for replay.
