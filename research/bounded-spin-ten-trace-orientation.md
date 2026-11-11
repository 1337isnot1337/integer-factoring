# Ten-site trace control: row ranks, column functionals, and path swap

## Result

The ten-site operator model distinguishes the row-cut ranks, the column
trace functional, and an explicit path-pair swap.  Let `e_i` be the seven
ten-site primitive central idempotents, `E_i` their row-cut matrices,
`C_i` their column-cut matrices, and `P` the permutation of ambient path
pairs that swaps the two paths.  The exact traces are

```text
i                  0   1    2    3   4    5    6
rank(E_i)         34   1  225  300  35  406  155
trace(E_i) mod 5   4   1    0    0   0    1    0
(a,b,b,a) mod 5    4   1    0    0   0    1    0
trace(C_i) mod 5   4   1    0    0   0    1    0
trace(P E_i) mod 5 4   1    0    0   0    4    0
(a,b,a,b) mod 5    4   1    0    0   0    4    0
```

Thus the `(a,b,b,a)` pairing evaluates the ordinary row trace and is a
row rank residue.  The `(a,b,a,b)` pairing evaluates `trace(P E_i)`.

The path-pair swap `P` is **not** the action of the coefficient unit
`C_0`.  The coefficient identity `e_i C_0=e_i` therefore says nothing
about `trace(P E_i)`, and `(a,b,a,b)` is not an `R`-character of `e_i`.
This is consistent with the
[column-trace retraction](bounded-spin-twelve-column-trace-correction.md).

The ordinary column idempotent traces are `[4,1,0,0,0,1,0]`, equal to the
row traces here, while the ordinary column ranks are

```text
931, 1156, 356, 581, 771, 776, 715.
```

The coefficient idempotents are row-cut projections; they are not
column-cut projections.  The column trace functional is not the column
rank, and the relation between them remains open.

## Consequence for the twelve-site cut trace

The corrected twelve-site **row** residues are
`0,1,0,1,0,0,0,4`, with the same provenance as the ten-site row trace.
The corrected **column** residues are `0,4,0,1,0,0,0,1`; they are values
of the column trace functional, not column ranks.  Any claim that uses
them as column ranks is withdrawn pending an operator-level twelve-site
or transfer-matrix column computation.

This leaves the corrected row residues and the separator's `T=1`
two-component split.  It removes the previously claimed column splits
and the uncertified `T=0` and `T=4` splits.  The exact `T`-primary module
dimensions are `1485,3281,3155`; they are upper bounds for sums of the
individual central-idempotent ranks, not the ranks themselves.

## Reproduction

```sh
python3 tools/check_bounded_spin_ten_trace_orientation.py
make verify-bounded-spin-ten-trace-orientation
```

The checker rebuilds the ten-site coefficient tensor and both cut
operators, forms the seven idempotent matrices, computes ordinary traces
and ranks, evaluates both diagonal pairings, and compares the ordinary
pairing with `trace(P e)`.  No twelve-site operator or score change is
claimed.
