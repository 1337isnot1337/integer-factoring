# Direct unit traces for supplied integer presentations

## Result

A sparse relation matrix can be reduced by exact unit substitutions without
first constructing expanded graph maps. The trace retains the original row
and column identifiers of every pivot. A separate verifier replays each local
integer isomorphism using a column-incidence index.

This preserves the **entire presented abelian group**. After the substitutions,
columns absent from every remaining relation are recorded as free direct
summands. Only the torsion-extraction view omits those free residual columns.
Relations with coefficient two, or any other nonunit coefficient, are retained.

The supplied relation matrix remains an input to this theorem. There is no
bound on how quickly a factoring algorithm supplies it, and no improvement to
the uniform factoring exponent is claimed. The direct reduction need not have
the same rows or pivot order as a graph-first reduction.

## Exact local substitution

Let a current presentation have generators \(x_1,\ldots,x_g\) and integer
relation rows. Suppose row \(r\) has coefficient
\(\varepsilon\in\{-1,1\}\) in an eligible residual column \(j\). Its
relation is

\[
 \varepsilon x_j+\sum_{k\ne j}r_kx_k=0,
 \qquad
 x_j=-\varepsilon\sum_{k\ne j}r_kx_k.
\]

Remove generator \(x_j\) and its defining row. For every other row \(v\),
replace it by

\[
  v-v_j\varepsilon r.
\]

The new coefficient in column \(j\) is zero. This is an integer presentation
isomorphism: one map includes all retained generators, and the reverse map
sends the removed generator to its displayed integral expression. The pivot
relation verifies the two compositions. The change preserves every remaining
relation and has an inverse obtained by adjoining the generator and its
defining unit relation.

No division by two or by a nonunit occurs. A relation \(2x=0\) remains visible
and continues to present its two-torsion. The pair \(2x=0,3x=0\) also remains
unchanged by the unit-only procedure; this method does not claim to produce a
complete Smith or Hermite normal form.

Composing the certified local isomorphisms proves the result without expanding
their composite generator maps. In a class-group application, the original
relation rows are already verified identities. Each local substitution then
commutes with evaluation. Retained generators still refer to their original
supplied ideals by unchanged column indices.

## Producer and certificate

`tools/direct_unit_presentation.py` applies the existing compact unit producer
to the original matrix, viewed as the identity presentation. Its deterministic
pivot heap and incidence update code are unchanged. Only the fingerprint
serializer is replaced with a streaming implementation of the same byte
sequence, avoiding a second complete serialized matrix.

The trace records each pivot as its original row index, original column index
and sign. It also records surviving original generator indices, final rows,
unused residual columns, and the torsion-extraction projection. No graph
certificate, expanded old-generator word, or expanded old-row combination is
required.

The input fingerprint binds the retained artifact to its supplied matrix. It
does not replace an algebraic check: the verifier starts from the actual
matrix and checks every substitution before comparing the complete output.

## Independent incidence verifier

For each residual column \(j\), maintain the set \(I_j\) of active rows
whose coefficient there is nonzero. The invariant is exact:

\[
  i\in I_j\quad\Longleftrightarrow\quad
  i\text{ is active and its coefficient in column }j\text{ is nonzero}.
\]

Initialization enumerates the original row supports. At each recorded pivot,
the verifier checks that its row is active, its column is an unremoved residual
generator, and its coefficient equals the recorded sign. The affected rows are
exactly \(I_j\) excluding the pivot row.

The verifier is written separately from the producer's incremental updates.
For each affected row it removes **all** old residual-column memberships,
computes the complete new row over the union of the old and pivot supports,
and inserts **all** new memberships. Every untouched row has coefficient zero
in the pivot column and remains unchanged. Removing the pivot row's memberships
and checking that its column has no remaining incidence completes the step.
This proves the invariant inductively.

At the end, the verifier independently rebuilds the entire final incidence
index and checks equality. It checks that no eligible unit coefficient remains,
renumbers only once, reconstructs the complete output rows, and checks every
retained original generator index and free-column projection. It does not
consult the producer's heap, version numbers or pivot-choice logic.

Stable row slots preserve original indices throughout. After checking a pivot,
the verifier releases its obsolete current row. Original input rows remain
available and unchanged.

## Cost statement

The old full-scan verifier tests every active row at every pivot. The new
verifier visits only rows containing the pivot column. Let \(A_s\) be that
affected set at step \(s\), let \(r_s\) be its pivot row, and write

\[
 U_{s,i}=|\mathrm{supp}(v_i)\cup\mathrm{supp}(r_s)|
 \quad(i\in A_s).
\]

Apart from input normalization, fingerprint encoding and final output, its
sparse arithmetic work is bounded by the sums of
\(|\mathrm{supp}(r_s)|\), \(|A_s|\log(|A_s|+1)\), and
\(\sum_{i\in A_s}U_{s,i}\). The logarithmic term charges the independent
sorting of affected row indices. Each coefficient update uses exact integer
multiplication and addition; actual coefficient bit lengths must also be
charged. The receipt counts substitutions, support visits and maximum replayed
row support.

This is a bound in terms of actual sparse updates, not a universal near-linear
bound. Fill, coefficient growth, the producer's stale heap entries, and the
final dense extraction matrix can all be large. The unchanged producer retains
its row slots and charges their stored contents; no bound based only on active
rows is asserted. The independent verifier stores the original input, current
rows, their incidence sets, the trace, and final output. It does not store
expanded composite presentation maps.

## Finite controls

The checker validates 607 bounded presentations against the unchanged producer
and its full-scan verifier. It checks 476 unit isomorphisms, 775 row
substitutions, 607 JSON round trips, 1,938 independent character counts, eight
corrupted-trace rejections and six invalid-input rejections. Character counts
compare the full original and reduced groups, and separately restore the
explicit free-factor contribution after projection.

Controls include negative unit pivots, empty rows, empty presentations, free
columns, coefficient-two torsion, and the nonunit pair \(2x,3x\). They also
check input immutability and byte-for-byte equality of the streaming and
reference fingerprints.

The retained complete 64-bit and 80-bit supplied prefixes give:

| Input bits | Original rows | Original generators | Unit pivots | Final rows | Nonfree columns | Rank | Extracted factor |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 2,048 | 1,871 | 1,794 | 138 | 17 | 17 | 3,332,910,829 |
| 80 | 10,240 | 10,100 | 9,749 | 28 | 22 | 22 | 1,017,359,648,351 |

These presentations differ from the graph-first presentations while preserving
the complete original groups. Both retain a proper factor under independently
checked class-group extraction. The 80-bit replay performs 11,008 row
substitutions and 94,581 support visits, instead of scanning all 10,240 row
slots at every one of 9,749 pivots.

The receipt, `evidence/direct-unit-presentation.json`, contains full traces and
references original supplied matrices by archive SHA-256. Its component times
include the incidence verifier. Input loading, final rank selection and factor
extraction are separate; these are not complete factoring times.

## Complete power-prefix adapter

`tools/direct_unit_oriented_power_search.py` keeps the entire compact-witness
query function unchanged. Two dependency bindings replace the graph step with
a dimension-checked identity descriptor and the later unit pass with the
direct reduction. Every norm test, atom update, checkpoint, incidence replay,
rank selection and factor extraction remains charged to the search.

The literal-power verifier also remains in place. A private constant-space
range supplies the retained-prefix length expected by its older graph
interface; this range is not saved or interpreted as a graph map. The saved
descriptor contains only the initial generator count, total generator count
and input row count. Those dimensions are checked against the actual atom
basis and supplied rows. The verifier then calls the independent direct unit
replay on that entire original matrix.

The complete matched 64-bit and 80-bit searches still stop at 2,048 and 10,240
queries with the factors in the table. Their full candidate witnesses, atom
words, original relation matrices and old exact-key comparisons equal the
retained references. The final presentations are allowed to differ and are
checked by their own integer-isomorphism traces. This is stronger than merely
comparing ranks or extracted factors.

The full search receipt is
`evidence/direct-unit-oriented-power-search.json.xz`. It includes both complete
new result records and independent verification. Search CPU times exclude JSON
conversion and saved-artifact replay, and include every checkpoint's incidence
verifier. No stable end-to-end speedup is inferred from one concurrent run.

Separate declared policies use `tools/direct_unit_oriented_fresh_run.py` and
`tools/check_direct_unit_oriented_fresh_run.py`. They preserve the existing
input and source-hash guards and bounded-failure reporting. The new code does
not modify any source used by an already active search.

## Reproduction

```sh
python3 tools/check_direct_unit_presentation.py \
  --replay evidence/direct-unit-presentation.json
python3 tools/check_direct_unit_oriented_power_search.py \
  --replay evidence/direct-unit-oriented-power-search.json.xz
```

The local isomorphism proof builds on
[compact unit-presentation traces](unit-presentation-traces.md). Direct
application and incidence verification remove the need to first construct the
[expanded graph presentation maps](signed-residual-graph-elimination.md).
