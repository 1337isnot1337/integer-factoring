# The twelve-site separator splits all T groups modulo five

## Result

Use the twelve-site coefficient algebra and the separator combination
`S=C_3+C_5` from the [separator direction](bounded-spin-twelve-RT-coefficient.md#the-separator-direction-at-the-algebra-level).
The joint `T,S` characters of the eight primitive central idempotents are

```text
idempotent   0     1     2     3     4     5     6     7
(T,S)       (0,2) (1,0) (1,2) (4,3) (4,0) (1,3) (0,4) (1,1)
```

The [corrected cut trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
records the idempotent trace-functional residues modulo five:

```text
row:    0, 1, 0, 1, 0, 0, 0, 4
column: 0, 4, 0, 1, 0, 0, 0, 1.
```

Here the column residues are a second diagonal trace orientation, not a
separate column-cut rank and not `tr(e_i R)`; see the
[column-trace retraction](bounded-spin-twelve-column-trace-correction.md).
Grouping by the `T` eigenvalue gives

```text
T=0:  idempotents 0,6         row residues 0,0       column residues 0,0
T=1:  idempotents 1,2,5,7     row residues 1,0,0,4   column residues 4,0,0,1
T=4:  idempotents 3,4         row residues 1,0       column residues 1,0
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. The nonzero residues already split `T=1` (`S=0` and `S=1`) and
give one nonzero component in each of `T=0` and `T=4`. The remaining
primaries with zero residues are resolved by explicit diagonal matrix
entries next.

## Direct diagonal witnesses

For an idempotent action, a nonzero diagonal matrix entry proves that the
operator is nonzero. The checker evaluates three such entries directly
from the path-pair polynomial forms:

```text
idempotent  path pair   row entry   column entry   consequence
0           (38,82)     4           1              rank >= 5
4           (0,0)       2           2              rank >= 5
6           (0,0)       4           4              rank >= 5
```

Because each idempotent action is a projection over `F5`, its trace equals
its rank. Idempotents `0`, `4`, and `6` have zero trace residue in both
orientations, so a nonzero entry forces their ranks to be positive
multiples of five, hence at least five. The nonzero residues of
idempotents `1` and `7` (`T=1`) and idempotent `3` (`T=4`) already give
nonzero components there.

Therefore every `T` group splits into two nonzero separator components in
both orientations:

```text
T=0:  S=2 (idempotent 0) and S=4 (idempotent 6)
T=1:  S=0 (idempotent 1) and S=1 (idempotent 7)
T=4:  S=3 (idempotent 3) and S=0 (idempotent 4)
```

This is a modulo-five splitting certificate. It does not give the exact
dimensions of the joint `T,S` components, their exact ranks, a two-cut
rank product, or a factoring exponent. The rank lower bounds are five for
idempotents `0,4,6` and one for idempotents `1,3,7`; idempotents `2` and
`5` are nonzero but are not needed for the splits.

## Relation to the withdrawn R interpretation

The attempted `R`-primary split is [withdrawn](bounded-spin-twelve-column-trace-correction.md).
In the coefficient algebra `R=C_0` is the unit, so `e_i R=e_i`; the second
trace orientation is therefore not a replica-swap-twisted trace `tr(e_i R)`
and does not give `R=+` and `R=-` multiplicities. No `R` character,
`R`-pure refinement, or `R`-mixed obstruction is claimed here.

## What changed

The first separator checkpoint established only character separation at
the coefficient-algebra level. A later revision restricted the cut-module
conclusion to `T=1` because the trace residues of idempotents `0`, `4`,
and `6` vanish. The present certificate removes that restriction with
explicit nonzero diagonal entries: the trace residue alone is not the
only possible witness. All three `T` groups are now certified to split in
both orientations. Exact component dimensions and the two-cut pairing
remain open.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
make verify-bounded-spin-twelve-separator-split
```

The checker reads the RT coefficient receipt and the corrected cut
trace-form receipt, checks their common basis, verifies the three `T`
groups, and evaluates explicit nonzero diagonal entries of idempotents
`0`, `4`, and `6`. It records all three certified `T`-group splits in
both orientations and the resulting rank lower bounds. No full cut
operator or exact rank computation is performed.
