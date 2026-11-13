# The twelve-site cut trace form has a large radical

## Result

The twelve-site coefficient table and the diagonal trace forms determine
the trace form of every product on the cut module. Evaluating the eight
primitive central idempotents gives the values of the cut trace functionals
on them modulo five. The corrected structural facts are

```text
row cut trace form rank       3
column cut trace form rank    3
cut trace radical dimension   7
regular trace form rank       8
regular trace radical         two-dimensional
```

The row and column forms share the same radical basis. Coefficient
direction `C_9` lies in the radical, and the first primitive central
idempotent is exactly `4 C_9`; nevertheless the radical is seven- rather
than one-dimensional. `C_9` has regular trace `4`, so it is not in the
regular trace radical. The cut module therefore has a large
trace-invisible subspace that the regular representation does not see.

The idempotent trace-functional residues in the two orientations are

```text
row:    0, 1, 0, 1, 0, 0, 0, 4
column: 0, 4, 0, 1, 0, 0, 0, 1
```

for regular ranks `1,1,2,1,2,1,1,1`. Each list sums to one modulo five,
matching the normalized trace of the coefficient unit. These are the
values of the cut trace functionals modulo five, not module dimensions;
for an idempotent projection, however, the trace equals the rank, so a
nonzero functional value or a nonzero diagonal matrix entry gives a rank
lower bound. The residues are checksums for the exact rank computation,
not the exact ranks themselves.

## Correction and normalization

An earlier version of this note reported rank `9`, a one-dimensional
radical spanned by `C_9`, and residues `0,3,1,2,2,1,2,3`. That computation
evaluated the trace forms at coordinate `k` instead of at the coefficient
label `basis[k]`; since the table basis is `0,2,3,4,5,6,7,8,9,10`, the
omitted labels `1` and `11` shifted every entry after the first. The
earlier `C_8` and `C_9` labels were both consequences of that indexing
error. The present checker maps the basis explicitly.

The table normalization also has to be fixed. The raw row trace of the
coefficient unit is `4`, while the raw column trace is `1`; the row trace
vector is therefore divided by `4` (equivalently multiplied by `4`) and
the column vector is unchanged. After normalization both unit traces are
one and both idempotent residue lists sum to one. The trace-form rank and
radical are unchanged by this nonzero scaling.

The diagonal trace note's ten-site control uses the same basis mapping and
already fixes the orientation: the row form at `(a,b,b,a)`, not the
transposed pairing `(a,b,a,b)`, recovers the published primary dimensions
modulo five. The column form here is retained as a separate checksum and
is normalized independently.

## Why this matters

- The cut trace radical is seven-dimensional, not the regular trace
  radical. An argument that uses only a single invisible direction, or
  that identifies the cut radical with the regular radical, is incomplete.
- The first central idempotent is `4 C_9`; its zero residue is explained
  by the fact that `C_9` lies in the cut trace radical. A direct nonzero
  diagonal entry shows that the idempotent has nonzero rank; since its
  trace is zero modulo five, its rank is a positive multiple of five.
- `C_9` lies outside the five-dimensional coefficient subalgebra
  `<C_0,T>` and extends it only to dimension `6`. It is not one of the
  minimal completing pairs `(2,6)`, `(3,5)`, `(3,6)` found by the peer
  `R,T` scan. The cut trace radical is therefore not a route to the full
  coefficient table by itself.
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

The checker rebuilds the row and column diagonal trace forms, maps the
table basis into the coefficient labels, normalizes the coefficient unit,
verifies the rank-three trace form and its seven-dimensional radical,
checks that `C_9` lies in the radical and that idempotent `0` is `4 C_9`,
checks the `C_9` extension against the peer `R,T` certificate, and
evaluates the eight idempotents. It does not build a cut operator or claim
an exact rank.
