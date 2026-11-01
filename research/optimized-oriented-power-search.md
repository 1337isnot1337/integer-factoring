# Exact component acceleration and a certified 96-bit prefix

## Result

The [optimized selector](../tools/optimized_sparse_oriented_power_search.py)
uses the unchanged power-search body with two exact replacements:
[small-leaf quotient updates](quotient-update-oriented-basis.md) in the
atom tree, and [sparse unit elimination](sparse-unit-presentation.md).
Every witness, atom word, relation, presentation certificate, stopping
query, and extracted factor agrees with the retained 64/80-bit baseline.
Only runtime measurements and tree implementation counters change.

The combined replacements reduce the retained 80-bit search time from
67.006 seconds to 9.694 seconds, with all checkpoint verification still
included. A separate certified 96-bit input factors after 20480 tested
norms, while the complete original exact-key selector needs 35106.
The optimized method still takes more CPU time than that original
selector on this fixture. These results improve finite relation supply
and the cost of the new implementation; they do not improve the uniform
factoring score beyond `1/2`.

## 1. Frozen search body and exact matched receipts

The baseline search and all of its source guards remain unchanged.
The adapter creates an isolated function namespace with exactly two
dependency bindings replaced. Its Python code object and keyword
defaults are identical to the reference function. The original module's
global namespace is not modified. This keeps the witness order,
checkpoint schedule, stopping rule, and verification calls identical.

The [matched checker](../tools/check_optimized_oriented_power_search.py)
compares the complete resulting record with the baseline. It permits
differences only in the declared timing fields, atom-tree implementation
counters, and component labels. It then independently replays the powers,
literal ideal products, integer presentation maps, and factor extraction.

The compact [matched receipt](../evidence/optimized-oriented-power-search.json)
refers to the existing complete compressed certificates by SHA-256 and
stores only those allowed metadata replacements. Reconstruction yields
the complete optimized record. No mathematical witness is omitted or
replaced by an unchecked digest.

### Measured matched costs

| Fixture | Reference component implementation | Optimized implementation | Measurement |
| --- | ---: | ---: | --- |
| 64 bits, 2048 queries | 3.349 s | 1.945 s | Medians of three paired repeats, alternating order |
| 80 bits, 10240 queries | 67.006 s | 9.694 s | Retained complete reference run and identical-result optimized run |

The paired 64-bit CPU medians are 1.573 and 0.940 seconds. The optimized
80-bit run uses 8.565 CPU seconds. The shared machine was running other
research work, so wall and CPU times differ. The earlier 64-bit baseline
artifact records a different single-run wall time; the paired experiment
is retained separately rather than replacing that historical measurement.

All presentation checks remain inside the search interval at every
checkpoint. Independent post-search replay and certificate serialization
are outside that interval, as in the baseline. These timings compare the
two implementations of the oriented method. The older exact-key selector
remains faster on the retained inputs despite using more tested norms.

## 2. A larger certified input

The 96-bit input is

$$N=61655890914495222928111725881
 =221484168306413\cdot278376063562237.$$

The existing deterministic input generator produces recursive prime
certificates; 22 certificates and 29 rejected preceding candidates replay.
The search receives only N. Its bounded policy uses cutoff 1021, at most
65536 norm queries, checkpoints every 2048 queries, and a 256-column
extraction cap. No checkpoint is skipped.

At its first successful checkpoint, the search has:

| Quantity | Value |
| --- | ---: |
| Tested leading norms | 20480 |
| Distinct residual keys | 19705 |
| Coprime atoms | 20402 |
| Nonzero signed atom usages | 27584 |
| Graph-eligible rows | 19600 |
| Larger residual rows transported | 880 |
| Additional integer unit pivots | 849 |
| Final core rows | 321 |
| Final core generators | 89 |
| Independent core rank | 89 |

The core retains 88 base generators and one residual generator. All 880
larger residual words pass through the exact substitution. The remaining
generator is retained in the presentation; no rational elimination is
used to remove it. The extractor obtains the ambiguous form
`(221484168306413,0,278376063562237)` and returns its proper leading factor.

The independent replay checks all 20480 powers using generic HNF
composition, all 19705 literal residual ideal factorizations, and both
complete presentation certificates. The graph certificate verifies
20490 old-generator identities and 20480 input-relation identities.
The unit certificate independently preserves that same original system.

At the identical query cap, the original exact-key selector has 19705
residual keys, 775 candidate relations, 737 dependent relations, rank 38,
and no factor. Its separately retained full traversal first factors the
input after 35106 queries with rank 86. Independent execution reproduces
that entire result exactly.

## 3. Costs and scope

The query ratio is `35106/20480`, about 1.71. The retained full search
costs are:

| Method | Queries to factor | Wall time | CPU time |
| --- | ---: | ---: | ---: |
| Optimized oriented prefix | 20480 | 43.028 s | 18.982 s |
| Original exact-key selector | 35106 | 10.816 s | 4.818 s |

Thus fewer queries do not yet imply a faster complete search. The
optimized run's recorded wall-time components are 10.672 seconds for
incremental atom insertion, 10.387 for graph reduction, 10.726 for unit
elimination, 3.289 for extraction, and 2.025 for sparse transport. The
remaining time includes factor-base construction, power generation,
witness checks, and process scheduling effects. The atom tree's final
product has 633029 bits; all gcds, exact divisions, multiplications,
certificate operations, and repeated checkpoints remain charged.

The query tables count tested leading norms, matching the original
selector's counter. Construction of the initial powers and next-power
updates is included in total runtime. It is not relabeled as additional
tested norms. No previously rejected candidate or earlier checkpoint is
removed from the recorded work.

These controls establish a new, reproducible finite supply of useful
relations and cheaper processing of that supply. They do not bound
uniformly how long every input needs to reach a useful torsion class.
That missing search bound remains necessary for a stronger logarithmic
saving exponent.

## 4. Replay

The [larger-run checker](../tools/check_optimized_oriented_fresh_run.py)
retains the full source-hashed adapter, input certificates, every queried
witness, and both complete presentation maps in the
[96-bit compressed receipt](../evidence/optimized-oriented-power-96.json.xz).
The [readable summary](../evidence/optimized-oriented-power-96-summary.json)
contains all measured costs, counts, comparison results, and source hashes.

```sh
python3 tools/check_optimized_oriented_power_search.py \
  --replay evidence/optimized-oriented-power-search.json
python3 tools/check_optimized_oriented_fresh_run.py \
  --replay evidence/optimized-oriented-power-96.json.xz
```

To repeat the complete larger search under its original bounded policy:

```sh
python3 tools/check_optimized_oriented_fresh_run.py \
  --bits 96 --cutoff 1021 --max-queries 65536 --checkpoint-interval 2048 \
  --max-extraction-columns 256 --output /tmp/optimized-oriented-96-rerun.json.xz
```

This implementation and analysis used substantial OpenAI assistant support.
