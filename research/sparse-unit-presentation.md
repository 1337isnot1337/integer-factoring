# Sparse unit elimination with complete integer certificates

## Result and scope

The [incidence implementation](../tools/sparse_unit_presentation.py) returns
exactly the same integer presentation, pivot sequence, and certificates as
the [reference unit pass](../tools/signed_residual_unit_elimination.py).
It keeps stable identifiers until the final output and updates only the
rows and generator words that contain the current pivot. On three matched
80-bit final-pass runs, the median falls from 13.116 seconds to 0.853
seconds, including the complete certificate verification.

This removes repeated scans and renumberings from the implementation.
It does not change the relation stream, supply new relations, or establish
a better uniform factoring exponent. Coefficient growth, sparse fill,
expanded certificates, and all verification work remain charged.

## 1. The exact substitution

Use the notation of the
[signed residual graph theorem](signed-residual-graph-elimination.md).
An output row with residual coefficient `s` equal to `1` or `-1` has the form

\[
 s y_j+u=0,\qquad y_j=-s u.
\]

Let its source expression in the original rows be `b`. Substituting this
identity into a remaining row with coefficient `a` at column `j` removes
that coefficient and adds `-a s u` to the other coefficients. Its source
certificate changes by `-a s b`. An old-generator word with coefficient
`a` at that column receives the same word substitution, and its definition
certificate changes by `a s b`.

These are the identical integer operations used by the reference pass.
The pivot relation becomes zero under substitution. Any other row that
becomes zero may also be discarded. Removing their indices from the
old-row expressions preserves the input-relation identity. Every surviving
output generator remains one of the original generators.

Consequently the four complete maps still satisfy

\[
 RW=AS,\qquad SV=BR,\qquad I-WV=DR,\qquad VW=I.
\]

They prove an isomorphism of the full integer presentations, commuting
with their evaluation into the class group. A residual coefficient `2`
is never treated as a unit. The 80-bit retained example still has its
nonfree residual generator and all double-root information.

## 2. Stable identifiers preserve the deterministic pivot rule

The reference chooses the lexicographically least triple
`(row support size, residual column, row index)` among unit coefficients.
Its deletions renumber surviving rows and columns monotonically. Therefore
the same comparisons give the same choice when the original identifiers
are kept throughout the pass.

For each residual column, the new pass keeps two incidence sets: the
output rows containing it and the old-generator words containing it.
Only those vectors receive a substitution. Their incidence memberships
are changed precisely when a coefficient appears or cancels. Source
certificates are updated along with the corresponding vector.

A heap contains the unit candidates. Each row has a version counter.
After a row changes, its new candidates enter the heap with the new
version; older entries are ignored when removed. A pivot row is marked
inactive and removed from its incidence sets. Thus the least valid heap
entry is the same pivot the full reference scan would find.

Induction over pivots proves equality of all live vectors and source
certificates after interpreting their stable identifiers. The final pass
renumbers survivors once and updates component roots and free columns.
The complete serialized mathematical output agrees with the reference,
including its list of pivot sources and all four presentation maps.

## 3. Work and space accounting

Let a pivot have `e` other nonzero entries and a source certificate with
`s` nonzero entries. If it occurs in `a` other rows and `g` generator words,
its substitutions use `O((a+g)(e+s))` coefficient operations, with their
actual integer operand sizes included. Updating incidence memberships
requires the same touched supports. Removing the pivot also charges its
own support. Unaffected vectors receive no substitution or renumbering.

Every changed row is scanned to insert its new unit candidates. If `H`
heap entries are inserted in total, heap maintenance takes
`O(H log(1+H))` identifier comparisons and `O(H)` worst-case stored entries,
including stale entries not yet removed. The initial and final passes
charge every stored vector and certificate entry. Sparse fill can be
large, and an expanded source expression can be much larger than its
resulting relation. No essentially-linear bound is asserted.

The original unit-elimination bit-size analysis continues to apply,
because the arithmetic operations and pivot sequence are unchanged.
The implementation is a reduction of unnecessary data movement, with
the integer certificate output retained in full.

## 4. Independent controls

The [checker](../tools/check_sparse_unit_presentation.py) compares both
passes on 606 presentations, including signed pivots, extra rows,
cancellation, isolated free columns, empty inputs, and JSON-loaded maps.
It verifies complete integer maps and independently compares finite
character counts on small presentations. Every input certificate is
checked to remain unchanged.

It then replays both complete retained power-prefix certificates at
64 and 80 bits. Equality includes every output row, generator word,
old-row expression, source certificate, component record, pivot step,
and torsion projection. The ordinary checker always performs these
replays; matched timing runs are optional.

```sh
python3 tools/check_sparse_unit_presentation.py
python3 tools/check_sparse_unit_presentation.py --benchmark-repeats 3 \
  --output evidence/sparse-unit-presentation.json
```

The small controls produce 1,212 complete reference equalities, replay
651 unit pivots, and check 142 independent finite character counts.
The retained final-pass measurements are:

| Input fixture | Unit pivots | Reference median | Sparse median | Ratio |
| --- | ---: | ---: | ---: | ---: |
| 64 bits, 2,048 rows | 81 | 0.210785 s | 0.045798 s | 4.60 |
| 80 bits, 10,240 rows | 764 | 13.115864 s | 0.853242 s | 15.37 |

Each median uses three runs with alternating method order. Both methods
copy the same saved input and verify all four integer presentation maps
inside the timed interval. These are shared-machine measurements, not
cross-machine performance guarantees. The output comparison occurs
outside the timed interval and includes every certificate entry.

The [receipt](../evidence/sparse-unit-presentation.json) retains source
hashes, fixture hashes, full equality digests, verification counts, and
the measured samples. Timings concern the final presentation pass only;
the full search is a separate measurement.
