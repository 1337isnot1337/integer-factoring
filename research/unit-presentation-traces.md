# Compact traces for integer presentation isomorphisms

## Result and scope

The [compact unit pass](../tools/unit_presentation_trace.py) preserves the
same full integer presentation as the
[expanded-map implementation](sparse-unit-presentation.md). It retains
the supplied graph certificate and records each subsequent unit pivot by
three integers: its stable row index, column index, and sign. A separate
verifier reconstructs every elementary substitution from those records.

On the retained 80-bit example, the additional unit certificate falls
from 3,113,278 to 21,638 bytes in canonical JSON. Compressed sizes fall
from 251,508 to 7,096 bytes. The common graph certificate remains
necessary and is excluded from both sides of this incremental comparison.
No information about the evaluated torsion is discarded.

This changes the certificate representation and implementation cost.
The final rows, retained generators, pivot choices, and factor-search
queries are unchanged. No better uniform factoring exponent follows.

## 1. A unit pivot is an integer isomorphism

Suppose a presented abelian group has generators `x_0,...,x_(d-1)` and
one relation

\[
 s x_j+\sum_{i\ne j}u_i x_i=0,\qquad s\in\{-1,1\}.
\]

Set `e=-s sum_(i!=j) u_i x_i`. Replace `x_j` by `e` in every other
relation, then remove that generator and its defining relation. This
gives an isomorphic integer presentation.

Indeed, the forward map fixes every surviving generator and sends
`x_j` to `e`. The inverse map includes the surviving generators in
the original quotient. One composition is literally the identity.
For the other, the only generator difference is

\[
 x_j-e=s\left(s x_j+\sum_{i\ne j}u_i x_i\right),
\]

which is an integer multiple of the pivot relation. Substitution sends
the pivot to zero and sends every remaining relation to its recorded
new row. Any additional row that becomes zero may be removed. These
arguments prove both well-definedness and the inverse identities over
the integers; rational rank alone is not used.

Composing the elementary maps proves an isomorphism for any valid
finite pivot sequence. Their evaluation into the class group commutes
with the maps because every substitution uses a genuine relation.
The coefficient must be exactly `1` or `-1`. In particular, an isolated
double-root relation remains part of the presentation.

## 2. What the certificate stores

The original signed-graph certificate is retained separately. It
already supplies and verifies all four maps from the complete supplied
relation matrix to the graph presentation. The compact extension stores:

- a SHA-256 binding to the graph's base size, surviving original
  generator indices, and complete output relation matrix;
- one triple `(row, column, sign)` for each unit pivot;
- the final surviving row and generator indices and the final rows;
- the list of free columns and the resulting torsion projection.

The fingerprint identifies the input to the local trace; the full
graph verification establishes its relationship to the original
relations. Correctness follows from exact replay on the supplied
matrix and does not require a collision-resistance assumption.

The row and column identifiers remain stable until the end. Therefore
a verifier can obtain a pivot's entire expression from the current
row without storing the expression a second time. The final row and
column renumbering is monotone and occurs once.

All residual generators absent from every final relation span free
direct summands. The complete presentation still includes them. The
separate `torsion_rows` projection drops those summands only when
passing the evaluated torsion image to the factor extractor.

The composite expanded maps remain mathematically defined. They can
be materialized by replaying the substitutions with the original map
bookkeeping if requested. The compact verifier does not need to write
them all to certify the composition of the elementary isomorphisms.

## 3. Construction, independent verification, and cost

The constructor uses sparse row incidence and a lazy heap. It follows
the same deterministic support/column/row pivot rule as the expanded
implementation. It no longer substitutes into every original-generator
word or expands its original-row source expression.

The verifier uses neither that incidence index nor that heap. For
each recorded pivot, it checks that the row is live, the residual
column has not been removed, and the coefficient is the recorded
unit sign. It then scans the live rows and performs the exact integer
row subtraction using a copy of the pivot. It checks the complete
final presentation and verifies that no unit residual pivot remains.
It accepts any valid complete pivot sequence; correctness does not
depend on reproducing the producer's heap implementation.

If there are `k` pivots, `r` initial rows and `d` initial generators,
the trace triples use `O(k log(r+d+2))` bits, besides the final output
and the shared graph certificate. The verifier's row scans cost
`O(kr)` incidence tests, and each actual row subtraction charges the
pivot support and all integer operand sizes. Sparse fill, coefficient
growth, the producer's stale heap entries, fingerprint construction,
and final serialization remain charged. No linear-time claim is made.

The default verifier checks the complete original graph as well as the
unit trace. In a search that has just verified that same graph, the
unit component may set `verify_graph=False` to avoid repeating that
completed check. The returned record explicitly states whether the
graph was checked in this particular call. A standalone replay uses
the default full check.

## 4. Exact controls and retained measurements

The [checker](../tools/check_unit_presentation_trace.py) includes 432
small presentations and JSON round trips, 581 local unit isomorphisms,
138 independent finite-character counts, and eight rejected corrupted
traces. Every final presentation and retained-generator selection
agrees with the expanded reference. Every input graph is unchanged.

The two complete retained examples use all original supplied rows and
verify their complete graph maps before checking the trace. The 80-bit
example performs 764 unit pivots and 9,969 actual row substitutions.
Its final 22 nonfree generators still include the retained residual
column of that unit-reduced presentation; no claim is made
that every different valid factor-extraction strategy needs it.

| Fixture | Expanded JSON | Compact JSON | Expanded xz | Compact xz |
| --- | ---: | ---: | ---: | ---: |
| 64 bits | 281,705 B | 11,208 B | 33,156 B | 2,292 B |
| 80 bits | 3,113,278 B | 21,638 B | 251,508 B | 7,096 B |

These are the additional unit certificates only. Both methods share
the complete original witnesses, oriented atom words, and graph maps.

Three matched runs measure the unit component after graph verification.
Both methods perform their own required unit checks inside the timed
interval. Method order alternates. On the 80-bit fixture, median CPU
time is 0.951199 seconds for expanded maps and 0.192337 seconds for
the compact trace. Median wall times are 2.140370 and 0.432677 seconds
on the shared machine. Standalone graph replay is performed separately;
it is required once and is not claimed to be free.

The [receipt](../evidence/unit-presentation-trace.json) contains both
complete compact traces, all source and input hashes, every timing
sample, and the exact byte counts.

```sh
python3 tools/check_unit_presentation_trace.py
python3 tools/check_unit_presentation_trace.py --benchmark-repeats 3 \
  --output evidence/unit-presentation-trace.json
```

Whole-search integration and larger input experiments have separate
receipts. This note establishes the certificate composition and its
measured component costs, without assuming a bound on relation supply.
