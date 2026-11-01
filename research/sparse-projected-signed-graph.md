# Sparse intermediate graph projections

## Result

The [separate graph implementation](../tools/sparse_projected_signed_graph.py)
keeps the complete four-map integer presentation certificate and stores
its projected row view sparsely. Every generator definition, relation
source, inverse identity, and output row agrees with the frozen graph
implementation. Expanding the sparse view reproduces its dense view
exactly.

This removes an intermediate rectangular allocation that stopped a
bounded 128-bit run. It changes presentation processing and proof
storage, with no new relation-supply or factoring-exponent claim.

## 1. The observed memory boundary

The first [declared 128-bit policy](../evidence/prepared-oriented-fresh-128-policy.json)
uses the prepared atom basis, small internal integer pivots, compact
unit traces, and a 2048-MiB address-space bound. Its last completed
checkpoint has 114688 tested norms, 120653 atoms, and refined rank 85
in a 511-column core. The old exact-key rank is 68. No factor is found.

At the next 131072-query checkpoint, the graph verifier raises
`MemoryError` while constructing a second complete dense projected
matrix for equality comparison. The
[failure receipt](../evidence/prepared-oriented-fresh-128-bounded-failure.json)
retains the exact traceback and links the
[unchanged progress log](../evidence/prepared-oriented-fresh-128-progress.jsonl)
by hash. No complete or partial raw witness record was saved.

The last checkpoint's 108.282 seconds are only a lower bound on that
failed attempt's work. Its final timing was not recorded by the old
driver. The attempted 131072 queries are not replaced by the smaller
last-completed count, and this attempt is not presented as a success.

## 2. The projected view is redundant storage

Let the graph output have m rows and q retained nonfree columns. The
existing certificate already stores each output row as an integer
column dictionary. If pi is the order-preserving map from retained
old columns to projected columns, the sparse view stores

$$\{\pi(j):a_{ij}\mid a_{ij}\ne0\}$$

for row i. This is exactly the dense row with its zero entries omitted.
No generator, relation, or coefficient is changed. Free residual
columns still occur in no output relation, and their direct summands
are removed only from the torsion view.

The verifier checks each projected dictionary against the complete
output row. It still replays all four integer-map identities. The
compact unit pass consumes the complete sparse output rows and
generator indices, so it does not need the intermediate dense view.
Its final extraction matrix remains the same ordinary dense matrix.

For s nonzero projected entries, the view uses m sparse row records
and s entries instead of m q cells. Indices, integer coefficient
lengths, and all other certificate maps remain charged. Verification
creates at most one projected-row dictionary at a time, rather than
a second m-by-q array.

The implementation also discovers unused columns from the support
set, reuses supplied read-only row containers, and releases completed
construction arrays before verification. Already normalized integer
dictionaries can be read without another copy. JSON dictionaries are
still normalized and checked. These changes do not alter any map.

## 3. Complete controls

The [component checker](../tools/check_sparse_projected_signed_graph.py)
and [receipt](../evidence/sparse-projected-signed-graph.json) retain:

- 2570 complete reference-certificate equalities, including all 1728
  three-vertex signed graphs and 420 mixed larger-row systems;
- 3596 independent finite-character counts and two ideal-factor controls;
- JSON round trips and unchanged supplied input rows;
- nine invalid-input rejections, four corrupted-map rejections, and
  four additional corrupted-projection rejections;
- complete equality with the retained 64-bit and 80-bit graph certificates.

A separate 1024-row supplied presentation has 1027 intermediate
columns. Its 1051648 dense projection cells become 3950 sparse entries.
All graph maps and the subsequent unit output agree exactly. The unit
verifier replays 1024 isomorphisms and 138783 row substitutions while
checking the complete original graph again.

## 4. Integration and a declared retry

The [new adapter](../tools/projected_oriented_power_search.py) uses the
same search code object and four exact dependency replacements.
Its [matched checker](../tools/check_projected_oriented_power_search.py)
expands the graph view only on bounded comparison fixtures, requires
equality of every mathematical field, then independently replays all
retained powers, ideals, maps, unit steps, and factor extraction.

The [complete matched receipts](../evidence/projected-oriented-power-search.json)
retain the same 2048-query and 10240-query successes on the 64-bit
and 80-bit fixtures. Recorded search wall times are 0.400 and 3.387
seconds. Checker CPU measurements also include in-memory JSON
canonicalization of the returned record; they are not identical in
scope to the selector's internal wall timer. These single runs are
controls for exact integration, not a uniform runtime assertion.

The [first sparse 128-bit retry policy](../evidence/projected-oriented-fresh-128-policy.json)
keeps the same input, cutoff, query cap, checkpoint interval, extraction
cap, and address-space bound. Its 29 source hashes were fixed before
the first retry query, after complete component and selector controls.
The run restarts from the first query. Work from the failed attempt
remains separately retained; no saved prefix is assumed to exist.

This retry passes the formerly failed 131072-query checkpoint and
completes 245760 queries, reaching refined rank 418 versus old rank
220 in 511 columns. It then reaches the same address-space guard at
the attempted 262144-query checkpoint, during sparse graph construction.
The [second failure receipt](../evidence/projected-oriented-fresh-128-bounded-failure.json)
and [progress log](../evidence/projected-oriented-fresh-128-progress.jsonl)
retain that separate attempt. The new driver records 475.496 wall
seconds, 468.490 CPU seconds, and maximum RSS 1992.055 MiB. These
outer-attempt timings include all attempted work and failure cleanup.
Again, no raw partial witness record is available and no factor is found.

A [third declared policy](../evidence/projected-oriented-fresh-128-3072-policy.json)
uses the identical source files, input, cutoff, checkpoint interval,
and extraction cap. It raises the query bound to 327680 and the
address-space guard to 3072 MiB. The third run restarts from query one;
both failed attempts and their costs remain separately retained. It is
an active bounded attempt, and any resulting saved certificate still
requires the separate independent replay checker before a success claim.

```sh
python3 tools/check_sparse_projected_signed_graph.py
python3 tools/check_projected_oriented_power_search.py \
  --replay evidence/projected-oriented-power-search.json
```

This implementation and analysis used substantial OpenAI assistant support.
