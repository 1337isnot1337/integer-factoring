# The twelve-site cut trace form has a cut-specific radical

## Result

The twelve-site coefficient table and the diagonal trace forms determine
the trace form of every product on the cut module. Evaluating the eight
primitive central idempotents gives their ranks modulo five. The new
structural fact is that the cut trace form is not the regular trace form:

```text
row cut trace form rank       9
column cut trace form rank    9
cut trace radical dimension   1
cut trace radical             span(C_9)
regular trace form rank       8
regular trace radical         two-dimensional
```

The unique cut-trace-radical direction is coefficient vector `C_9`. It has
regular trace `4`, so it is not in the regular trace radical. The cut
module therefore has a trace-invisible direction that the regular
representation does not see.

The coefficient label is `C_9`, not `C_8`: the table basis is
`0,2,3,4,5,6,7,8,9,10`, so basis position `8` is coefficient index `9`.
An earlier version of this note and of the checker used `C_8`; that was an
off-by-one indexing error. The peer `R,T` extension table independently
records extension by `C_9` as dimension `6`, while extension by the actual
`C_8` is dimension `7`.

The idempotent trace residues in the two orientations are

```text
row:    0, 3, 1, 2, 2, 1, 2, 3
column: 0, 1, 0, 1, 1, 3, 1, 4
```

for regular ranks `1,1,2,1,2,1,1,1`. These are the ranks modulo five.
The first scalar idempotent has trace zero in both cuts; each dual block
has even trace, as required. The residues are checksums for the exact
rank computation, not the ranks themselves.

## Why this matters

- The cut-module trace radical is not the algebra's regular trace radical.
  An argument that identifies the two, or that uses only the regular
  idempotent certificate to infer cut ranks, is incomplete.
- Coefficient direction `C_9` is invisible to the cut trace form. It lies
  outside the five-dimensional coefficient subalgebra `<C_0,T>` and
  extends it only to dimension `6`. It is therefore not one of the
  minimal completing pairs `(2,6)`, `(3,5)`, `(3,6)` found by the peer
  `R,T` scan: the cut trace radical is a non-completing coefficient
  direction, not a route to the full table by itself.
- The exact row and column ranks, the two-cut pairing, and the structured
  noncentral decomposition remain open. A direct `F5` rank computation,
  a characteristic-zero model, or an independent rank certificate is
  still required; the naive large-prime lift fails
  ([obstruction](bounded-spin-characteristic-five-lift-obstruction.md)).

No factoring exponent changes. The complete bound remains `1/2`.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_cut_trace_form.py \
  --output evidence/bounded-spin-twelve-cut-trace-form.json
make verify-bounded-spin-twelve-cut-trace-form
```

The checker rebuilds the row and column diagonal trace forms, multiplies
through the exact table, verifies the rank-nine cut trace form, identifies
the radical with `C_9`, checks its extension against the peer `R,T`
coefficient certificate, and evaluates the eight idempotents. It does not
build a cut operator or claim an exact rank.
