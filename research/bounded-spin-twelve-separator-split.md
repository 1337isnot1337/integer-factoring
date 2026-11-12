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

Here the row residues are ordinary ranks and the column residues are the
replica-swap-twisted traces, not a separate column-cut rank. Grouping by
the `T` eigenvalue gives

```text
T=0:  idempotents 0,6         row residues 0,0       column residues 0,0
T=1:  idempotents 1,2,5,7     row residues 1,0,0,4   column residues 4,0,0,1
T=4:  idempotents 3,4         row residues 1,0       column residues 1,0
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. Therefore `T=1` splits into two row components and two column
components, namely the primaries with `S=0` and `S=1`. The `T=4` group has
only one nonzero residue in each orientation (`S=3`), so its split is
**not** certified. The `T=0` group has only one certified nonzero
primary, as explained next, so it is also **not** certified to split.

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

This is a modulo-five splitting certificate. It does not give the exact
dimensions of the joint `T,S` components, their ranks, a two-cut rank
product, or a factoring exponent. It uses the corrected cut-trace
certificate and the coefficient-algebra joint `T,S` characters.

## Relation to the R refinement

The column trace is the replica-swap-twisted trace, so the corrected row
and column residues also give the `R=+` and `R=-` multiplicities modulo
five. The corrected result is much weaker than the earlier
pre-correction claim: no nontrivial primary is certified mixed. Primaries
`1` and `7` are pure `R=-`, primary `3` is pure `R=+`, and the other five
primaries are undetermined because both residues vanish. See
[the R-primary note](bounded-spin-twelve-R-primary-split.md).

## What changed

The previous separator checkpoint established the character separation
at the coefficient-algebra level. This note transfers only the `T=1`
nonsingular part of that separation to the cut-module partition. The
direct diagonal witness makes idempotent `0` nonzero but does not split
`T=0`; `T=4` also remains uncertified. The unresolved targets are now the
exact ranks of idempotents `0,6` (both multiples of five) and the dual
idempotent `4`.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
make verify-bounded-spin-twelve-separator-split
```

The checker reads the RT coefficient receipt and the corrected cut
trace-form receipt, checks their common basis, verifies the three `T`
groups, and evaluates the two explicit diagonal entries that make
idempotent `0` nonzero. It records the one certified `T=1` split and the
unresolved `T=0` and `T=4` groups. No full cut operator or exact rank
computation is performed.
