# The twelve-site separator splits the T=1 group modulo five

## Result

Use the twelve-site coefficient algebra and the separator combination
`S=C_3+C_5` from the [separator direction](bounded-spin-twelve-RT-coefficient.md#the-separator-direction-at-the-algebra-level).
The joint `T,S` characters of the eight primitive central idempotents are

```text
idempotent   0     1     2     3     4     5     6     7
(T,S)       (0,2) (1,0) (1,2) (4,3) (4,0) (1,3) (0,4) (1,1)
```

The [corrected cut trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
records the idempotent rank residues modulo five:

```text
row:    0, 1, 0, 1, 0, 0, 0, 4
column: 0, 4, 0, 1, 0, 0, 0, 1.
```

Here the row residues are the ordinary row-cut trace and the column
residues are a second diagonal trace orientation, not a separate
column-cut rank. Neither orientation is `tr(e_i R)`; see the
[column-trace retraction](bounded-spin-twelve-column-trace-correction.md).
Grouping by the `T` eigenvalue gives

```text
T=0:  idempotents 0,6         row residues 0,0       column residues 0,0
T=1:  idempotents 1,2,5,7     row residues 1,0,0,4   column residues 4,0,0,1
T=4:  idempotents 3,4         row residues 1,0       column residues 1,0
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. Therefore `T=1` splits into two row components and two column
components, namely the primaries with `S=0` and `S=1`. The trace residues
alone leave the second `T=4` primary (`S=0`) and the second `T=0` primary
(`S=4`) ambiguous; the direct matrix entries below certify both. Thus
each `T` group has two nonzero certified central components. This is a
nonzero-component certificate, not a claim that the two components
exhaust the full `T`-primary block.

## The first idempotent is nonzero

Idempotent `0` has zero trace residue in both orientations, but it is not
forced to zero. Its coefficient vector is supported only at coefficient
`9`:

```text
e_0 = 4 C_9.
```

At the twelve-site path indices `(38,82)`, the diagonal entries are

```text
row    (a,b,b,a):  [z^9] P = 1,   e_0 entry = 4,
column (a,b,a,b):  [z^9] P = 4,   e_0 entry = 1.
```

Both are nonzero, so the idempotent operator is nonzero in both
orientations. Since its trace residue is zero, its rank is a positive
multiple of five. This certifies one `T=0` component, with separator
character `S=2`; the other `T=0` idempotent `6` (`S=4`) still has zero
residues and may be zero or a positive multiple of five. Hence `T=0` is
not split by the present certificate.

## Direct witnesses for idempotents 4 and 6

The remaining two primaries also have nonzero row-cut entries. At the
diagonal path pair with path index `0` (the return path
`(1,1,1,-1,1,-1,1,-1,1,-1,-1,-1)`), the row-cut diagonal entries are

```text
idempotent 4 (T=4,S=0):  2,
idempotent 6 (T=0,S=4):  4.
```

Together with `e_0=4C_9` at path indices `(38,82)`, this gives two
nonzero central components in every `T` group:

```text
T=0:  e_0 (S=2), e_6 (S=4)
T=1:  e_1 (S=0), e_7 (S=1)
T=4:  e_3 (S=3), e_4 (S=0).
```

The `T=1` witnesses can also be read from the nonzero row residues; the
new entries are for idempotents `4` and `6`. This is a modulo-five
nonzero-component certificate. It does not give the exact dimensions of
the joint `T,S` components, their ranks, whether two components exhaust
a block, a two-cut rank product, or a factoring exponent. It uses the
corrected cut-trace certificate and the coefficient-algebra joint `T,S`
characters.

## Relation to the withdrawn R interpretation

The attempted `R`-primary split is [withdrawn](bounded-spin-twelve-column-trace-correction.md).
In the coefficient algebra `R=C_0` is the unit, so `e_i R=e_i`; the second
trace orientation is therefore not a replica-swap-twisted trace `tr(e_i R)`
and does not give `R=+` and `R=-` multiplicities. No `R` character,
`R`-pure refinement, or `R`-mixed obstruction is claimed here.

## What changed

The previous separator checkpoint established the character separation
at the coefficient-algebra level. This note transfers the nonsingular
part of that separation to the cut-module partition. Direct diagonal
witnesses now certify two nonzero central components in every `T` group,
including the previously ambiguous primaries `4` and `6`. The unresolved
targets are the exact ranks and the two-cut pairing; nonzero components
need not exhaust their `T`-primary blocks. The second trace orientation
is not `tr(e_i R)`, so no `R` refinement is available from these residues.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
python3 tools/check_bounded_spin_twelve_separator_witnesses.py
make verify-bounded-spin-twelve-separator-split
make verify-bounded-spin-twelve-separator-witnesses
```

The first checker reads the RT coefficient receipt and the corrected cut
trace-form receipt, checks their common basis, verifies the three `T`
groups, and evaluates the two explicit diagonal entries that make
idempotent `0` nonzero. The second evaluates the diagonal entries for
idempotents `4` and `6` at path index `0`. They record nonzero certified
components in all three `T` groups; exact ranks and exhaustion of the
blocks remain open. No full cut operator or exact rank computation is
performed.
