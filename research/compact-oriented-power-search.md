# Compact proofs in the exact oriented power selector

## Result

The [compact selector](../tools/compact_oriented_power_search.py) combines
three exact component replacements in the frozen oriented power search:
small-leaf quotient updates, primitive integer rank tests, and local unit
substitution traces. The retained 80-bit run has the same 10240 queried
witnesses, the same signed atom words, the same graph presentation, the
same original independent relation rows, and the same factor as the
complete reference run. Its final unit proof takes 21639 JSON bytes,
compared with 3113278 bytes for the expanded presentation maps.

This is an implementation and proof-storage improvement. It does not
bound the supply of useful relations uniformly and does not improve the
current logarithmic saving exponent of `1/2`.

## 1. Three exact replacements

The adapter constructs an isolated Python function namespace. The search
code object, its defaults, and its complete control flow stay identical
to [the retained implementation](../tools/sparse_oriented_power_search.py).
Exactly three dependency bindings change:

| Component | Exact replacement | Reason it preserves the result |
| --- | --- | --- |
| Atom product tree | [Small-leaf quotient updates](quotient-update-oriented-basis.md) | Every ancestor product changes by exact division by its old leaf and multiplication by its new leaf. |
| Independent row selection | [Primitive temporary residuals](primitive-integer-row-selection.md) | Integer residuals are nonzero rational multiples of Gaussian residuals; the emitted relations are untouched original rows. |
| Unit presentation certificate | [Local unit substitution trace](unit-presentation-traces.md) | Every recorded unit pivot gives an independently checked integer presentation isomorphism. |

The original module globals are unchanged. The binding checker requires
the exact search code object and the exact three replacement names. No
candidate, checkpoint, or verification call is removed from the search.

The graph reducer already verifies its complete integer presentation
certificate at each checkpoint. The compact unit adapter checks that
verification flag and avoids immediately repeating that same graph
verification inside the unit component. It still verifies every unit
substitution locally. Independent replay of a saved result verifies the
graph again, then every unit substitution, and finally the factor
extraction.

The primitive row basis divides only an internal rank-test residual by
its content. It returns the same original rows as rational Gaussian
selection. In particular, a relation with a common factor is never
replaced by an unproved divided relation.

## 2. A matched complete certificate

The [receipt](../evidence/compact-oriented-power-search.json) refers to
the immutable [80-bit reference](../evidence/sparse-oriented-power-80.json.xz)
by file hash. It stores the new unit trace and a small declared metadata
patch. The checker reconstructs the complete compact result from that
reference before independent replay.

The only permitted differences are elapsed times, tree implementation
counters, component labels, and the unit proof representation. The final
unit presentation must agree exactly in all of these fields: generator
indices, output rows, torsion rows, torsion generator indices, free
columns, base size, and operation counts. Every other mathematical field
of the complete search record must remain equal. Thus the receipt does
not replace any power, ideal identity, relation, or final presentation
by an unchecked hash.

For

$$N=1035490702726305096363691
  =1017821676341\cdot1017359648351,$$

the run uses cutoff 257, a 32768-query bound, checkpoints every 1024
queries, and a 256-column extraction cap. It returns the factor
`1017359648351` after 10240 tested norms. The final core has 30 rows and
22 generators, with independent rank 22. All 809 larger residual words
remain in the presentation; 764 unit pivots are replayed exactly.

The independent checker replays all 10240 powers and all 9776 literal
residual ideal identities. It checks the complete original graph and
all 764 unit isomorphisms, including 9969 row substitutions. The final
ambiguous form is
`(508795331173,231013995,508795331173)`.

The old exact-key selector still has no factor at that same query cap.
Its retained full traversal needs 38818 tested norms and is reproduced
independently. This preserves the previously established query benefit;
it does not create a new query improvement over the oriented baseline.

## 3. Measured costs and limitations

The complete compact search takes 11.675 wall seconds and 5.132 CPU
seconds in the retained run. The earlier quotient-tree and expanded-unit
implementation takes 9.694 wall seconds and 8.565 CPU seconds on the same
mathematical search. Other research processes were active during both
runs. The lower CPU cost is useful evidence for these exact replacements;
these two single runs do not establish a stable wall-time speedup.

The compact run's recorded wall components are 2.925 seconds for atom
insertion, 0.625 for sparse transport, 4.846 for graph reduction, 1.389
for compact unit elimination, and 0.226 for extraction. Construction of
the factor base, power updates, all tested candidates, and every repeated
checkpoint remain inside the total search interval. Independent replay
and final serialization are outside it, as in the reference protocol.

The proof reduction concerns the incremental unit certificate. The
complete original graph certificate and all queried witnesses remain
available and checked. It is not a claim that the whole search proof has
only 21639 bytes. The [unit trace component study](unit-presentation-traces.md)
separately measures final-pass time and compressed certificate size.

The search counter counts tested leading norms. Initial power
construction and next-power updates are charged in runtime but are not
renamed as additional tested norms. No stopping-query credit comes from
discarding generated or checked work.

## 4. Replay

```sh
python3 tools/check_compact_oriented_power_search.py \
  --replay evidence/compact-oriented-power-search.json
```

To repeat the complete bounded search and then verify it:

```sh
python3 tools/check_compact_oriented_power_search.py \
  --bits 80 --output /tmp/compact-oriented-rerun.json
```

The checker requires every declared source hash and the retained full
reference hash to agree before using the compact receipt.

This implementation and analysis used substantial OpenAI assistant support.
