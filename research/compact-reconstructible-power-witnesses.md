# Compact reconstructible witnesses for exact power prefixes

## Result and scope

The exact oriented selector can store each power witness as a triple:

\[
  (\text{base word},\ \text{residual key index},\ \text{orientation}).
\]

The input integer, factor-base policy and query index determine every omitted
field. The independent verifier reconstructs the prime, exponent, target form
and full partial norm split. It then checks the retained triple against those
values. No factorization or residual primality query is introduced.

The new producer also clears the five final-checkpoint record references before
constructing the next checkpoint. This removes the simultaneous retention of
those preceding and succeeding records. The live atom state, old-key selector,
current witnesses and all temporary construction work remain charged.

These are storage changes to a finite exact search. They do not bound relation
supply, improve the uniform factoring exponent, or show that query savings imply
a reduction in total running time. The earlier bounded 128-bit failures remain
separate evidence; this implementation does not change their recorded outcomes.

## Why the omitted fields are determined

Let the sorted initial prime base have size \(d\), let \(p_j\) be its entry
in column \(j\), and let \(e_j\) be the least exponent for which
\(p_j^{e_j}\) exceeds the fixed reduced-form height bound. At zero-based query
\(q\), write

\[
  q=ud+j,\qquad 0\le j<d.
\]

The selected prime is \(p_j\), the exponent is \(e_j+u\), and the target is
the reduced ideal-class power \(\mathfrak p_j^{e_j+u}\). The algorithm obtains
successive targets by prime-form composition. The verifier computes every
target independently by binary powering with generic ideal-HNF composition.

The target, fixed known primes and their orientations determine the exact
smooth part, its signed base word, and the oriented residual ideal. Residual
keys receive their indices at their first occurrence. Thus the entire key
sequence is determined by the same query prefix.

After subtracting the decoded known-prime word from the queried prime power,
the producer retains the resulting base word \(w_q\), the residual index
\(i_q\), and its sign \(\varepsilon_q\). A missing residual uses a null
index. If the final oriented atom decomposition of that key is \(v_{i_q}\),
the exact transported relation is

\[
  w_q-\varepsilon_q v_{i_q}=0.
\]

This expression is exactly the row produced from the older full witness.
Its derivation uses the whole integer word, including signs and multiplicities.
No division, rational replacement or parity-only approximation is applied to
the emitted relation.

A direct factor found in the current target still has a separate terminal
witness containing its prime, exponent and target. Such a target is counted in
the query total and is not also inserted as a partial-norm witness.

## Query and checkpoint equivalence

The reference search body is pinned at SHA-256
`57a5178c9afb628ef102776e4353c68763366e80cb104a70287cb306c0000c07` in
`tools/sparse_oriented_power_search.py`.

The syntax checker compares the entire function after three declared edits:

1. Replace the full witness dictionary with the triple above.
2. Use the corresponding triple-decoding transport helper.
3. Clear the preceding final-checkpoint record references immediately before
   exporting the next atom basis.

Only the function name, docstrings and added descriptive result metadata are
ignored in this comparison. The query loop, earliest direct-factor test, old
exact-key comparison, target advancement, checkpoint condition, graph call,
unit call, row selection, extraction and stopping conditions remain identical.

The component bindings are the already checked prepared bucket basis,
small-magnitude primitive row selector, sparse intermediate graph projection,
and compact unit trace. Their proofs and controls are in
[the prepared-basis note](batched-coprime-bucket-basis.md),
[the row-selection note](small-pivot-integer-row-selection.md),
[the sparse-projection note](sparse-projected-signed-graph.md), and
[the compact unit-trace note](unit-presentation-traces.md).

No later query or checkpoint reads the omitted witness fields. At each
checkpoint, both representations therefore produce the same atom inputs,
transported integer rows, graph maps, unit pivots, selected original rows and
factor. Clearing an obsolete record changes object lifetime, not these values.
The deallocation and replacement work occurs inside the charged transport
phase. This is not a resumable-checkpoint guarantee: an unsuccessful allocation
can still leave no complete raw witness artifact.

## Independent replay

`tools/compact_witness_oriented_verifier.py` retains the complete mathematical
checks from the earlier verifier, with compact-witness decoding substituted at
the input boundary. It checks:

- the recursively certified input primes and deterministic input-generation
  prefix;
- every literal generic-HNF power, every reconstructed partial split, the
  complete residual-key order, and total trial divisions;
- each supplied residual norm as a literal product of oriented atom ideals;
- the complete transported rows and all four integer graph-presentation maps;
- every exact local unit substitution, with the graph certificate independently
  checked again;
- the selected original rows using an independent rational rank routine;
- the final relation-halving extraction and the original exact-key traversal at
  the same query cap.

When a complete original traversal is supplied, replay also reproduces its
entire saved mathematical result. Reconstructed compact triples are reused
after their individual checks; full omitted split dictionaries are not retained
for every query. Transported rows are compared one row at a time instead of
building a second serialized copy of the whole matrix.

The fresh-run driver uses the same declared-policy checks and bounded failure
logger as the preceding adapter. Its separate source manifest includes the new
producer, verifier, controls and drivers. It does not change any source guarded
by an already active search.

## Retained finite evidence

The complete matched 64-bit and 80-bit runs retain the same stopping queries,
factors and mathematical payloads as their original reference certificates.
Only the declared metadata, exact component proofs and witness representation
are replaced. Both new certificates undergo full independent replay.

| Input bits | Queries | Factor | Full witness JSON bytes | Compact witness JSON bytes | Full logical cells | Compact logical cells |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 2,048 | 3,332,910,829 | 637,814 | 54,441 | 101,759 | 21,675 |
| 80 | 10,240 | 1,017,359,648,351 | 3,390,155 | 264,539 | 487,396 | 99,986 |

The byte counts use canonical JSON. Logical cells count scalar values,
container nodes and dictionary keys. These exact representation counts are not
process RSS measurements. In particular, the table does not predict the peak
memory of atom refinement, graph construction or certificate replay.

The additional controls compare complete search results on 25 bounded policies
for six certified inputs from 32 to 42 bits. They exercise all three stopping
methods: direct target factor, oriented presentation factor, and query-cap
exhaustion. Six corrupted encodings, signs, indices or base words are rejected.
The source-backed receipt records the exact count of reconstructed powers and
the search CPU and elapsed times. Timing excludes JSON conversion and replay;
other research processes were active, so no stable speedup is claimed here.

The evidence is
`evidence/compact-witness-oriented-power-search.json`. Its patches reference
the complete retained 64-bit and 80-bit archives by SHA-256, preserving every
mathematical field rather than supplying only summary rows.

## Reproduction

```sh
python3 tools/check_compact_witness_oriented_power_search.py \
  --replay evidence/compact-witness-oriented-power-search.json
```

To generate a new receipt and repeat the complete matched searches:

```sh
python3 tools/check_compact_witness_oriented_power_search.py \
  --output /tmp/compact-witness-oriented-power-search.json
```

Fresh declared policies use `tools/compact_witness_oriented_fresh_run.py` and
`tools/check_compact_witness_oriented_fresh_run.py`. Input generation is
separate from the factoring algorithm. A saved search result becomes verified
evidence only after the latter tool completes its independent replay.
