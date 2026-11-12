# Retraction: the column trace is not an R-twisted trace

## Result

The earlier twelve-site `R`-primary split is **withdrawn**.  The two cut
trace orientations in the [cut trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
are not `tr(e_i)` and `tr(e_i R)`.

In the coefficient algebra, `R=C_0` is the unit and every primitive central
idempotent satisfies

```text
e_i R = e_i.
```

The cut action is an associative algebra representation, so the same
identity holds for the operators on the cut module:

```text
action(e_i) action(R) = action(e_i).
```

Therefore `tr(action(e_i) action(R)) = tr(action(e_i))`.  The row and
column trace residues differ, so the column form cannot be the
replica-swap-twisted trace.  The claimed decomposition into `R=+` and
`R=-` multiplicities was an algebraic misinterpretation and is not a
structural fact.

The column form is a second diagonal trace functional, not the
column-cut **rank** and not `tr(e_i R)`.  The exact
ten-site column-cut rank pairs are
`[931,1156,356,581,771,776,715]` in
[the idempotent-S control](../tools/check_bounded_spin_idempotent_s.py);
their residues do not match the column trace residues either.  Also, the
trace functional evaluated at `C_0` is `4` in the row orientation and
`1` in the column orientation, while `C_0` is not the cut-module
identity; these traces are therefore not module dimensions.  The
relation between the trace functionals and the column-cut ranks remains
open.

## Exact `T`-primary dimensions

The coarse `(T,R)` certificate records exact joint generalized dimensions
`1320,165,2015,1266,670,2485`.  Summing the two `R` signs gives the exact
`T`-primary module dimensions

```text
T=0: 1485    T=1: 3281    T=4: 3155,    sum 7921.
```

They are dimensions of the cut module, and their residues are
`0,1,0 mod 5`.  They are not the trace-functional values `2,3,4` in the
`T`-power table; the row trace functional has total `4 mod 5` at `C_0`,
while the module dimension is `7921 = 1 mod 5`.  The exact dimensions are
upper bounds for the sum of the central idempotent ranks inside each
`T` block, not the individual ranks.  Reproduce with
`make verify-bounded-spin-twelve-T-primary-exact`.

## What is retained

- The corrected exact twelve-site cut trace residues are recorded in the
  [rank-3 cut trace certificate](bounded-spin-twelve-cut-trace-form.md):
  both orientations have rank `3` and a seven-dimensional radical, and
  `C_9`, not `C_8`, is the radical direction carried by the first
  idempotent. The earlier `C_8` label was an indexing artifact; only the
  `R`-twist interpretation of the second orientation is withdrawn.
- The separator's `T`-group split and the `T,S` character data are
  unaffected.
- No `R` character, R-pure refinement, or R-mixed obstruction is claimed.

## Reproduction

The retraction is checked by the symbolic identity `e_i*C_0=e_i` in the
exact coefficient table and by the ten-site column-cut rank control:

```sh
make verify-bounded-spin-twelve-RT-coefficient
make verify-bounded-spin-twelve-cut-trace-form
python3 tools/check_bounded_spin_idempotent_s.py
```

No full cut operator or score change is involved.
