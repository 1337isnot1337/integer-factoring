# The twelve-site separator splits two coarse blocks modulo five

## Result

Use the twelve-site coefficient algebra and the separator combination
`S=C_3+C_5` from the [separator direction](bounded-spin-twelve-RT-coefficient.md#the-separator-direction-at-the-algebra-level).
The joint `T,S` characters of the eight primitive central idempotents are

```text
idempotent   0     1     2     3     4     5     6     7
(T,S)       (0,2) (1,0) (1,2) (4,3) (4,0) (1,3) (0,4) (1,1)
```

The [cut trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
records the idempotent rank residues modulo five:

```text
row:    0, 3, 1, 2, 2, 1, 2, 3
column: 0, 1, 0, 1, 1, 3, 1, 4.
```

Grouping by the `T` eigenvalue gives

```text
T=0:  idempotents 0,6         row residues 0,2       column residues 0,1
T=1:  idempotents 1,2,5,7     row residues 3,1,1,3   column residues 1,0,3,4
T=4:  idempotents 3,4         row residues 2,2       column residues 1,4
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. Therefore `T=1` splits into at least two row components and at
least three column components, and `T=4` splits into two components in
each orientation. The `T=0` group is **not certified**: idempotent `0`
has residue zero in both orientations, so its rank could be zero or a
positive multiple of five.

This is a modulo-five splitting certificate. It does not give the exact
dimensions of the joint `T,S` components, their ranks, a two-cut rank
product, or a factoring exponent. It uses the cut-trace certificate's
identification of the residues as ranks modulo five; the ten-site control
for that identification is checked in
[the diagonal trace note](bounded-spin-diagonal-trace.md).

## What changed

The previous separator checkpoint established the character separation
at the coefficient-algebra level. This note transfers the two nonsingular
parts of that separation to the cut-module coarse partition: the `T=1`
and `T=4` blocks are genuinely split by the separator, not merely
distinguished as algebra characters. The unresolved `T=0` block is now a
specific target: decide whether idempotent `0` has zero rank or a positive
rank divisible by five. Either answer is compatible with the present
residue certificate.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
make verify-bounded-spin-twelve-separator-split
```

The checker reads the RT coefficient receipt and the cut trace-form
receipt, checks their common basis, verifies the three `T` groups, and
records the four certified splits and the two unresolved `T=0`
orientations. No cut operator or exact rank computation is performed.
