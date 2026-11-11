# The twelve-site separator splits the T-eigenspace groups modulo five

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

Grouping by the `T` eigenvalue alone gives

```text
T=0:  idempotents 0,6         row residues 0,2       column residues 0,1
T=1:  idempotents 1,2,5,7     row residues 3,1,1,3   column residues 1,0,3,4
T=4:  idempotents 3,4         row residues 2,2       column residues 1,4
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. Therefore `T=1` splits into at least two row components and at
least three column components, and `T=4` splits into two components in
each orientation.

The `T=0` group has idempotent `0` residue zero in both orientations, so
the residues alone do not decide whether its rank is zero. The coefficient
vector of idempotent `0` is supported only at coefficient `9`:

```text
e_0 = 4 C_9.
```

At the twelve-site path indices `(38,82)`, the diagonal entries are

```text
row    (a,b,b,a):  [z^9] P = 1,   e_0 entry = 4,
column (a,b,a,b):  [z^9] P = 4,   e_0 entry = 1.
```

Both are nonzero, so the idempotent operator is nonzero in both cuts. Its
trace residue is zero, hence its rank is a positive multiple of five; in
particular `T=0` also splits into the two distinct `S` characters `2` and
`4`.

Thus all three `T`-eigenspace groups split in both orientations modulo
five. These groups are unions of the six joint `(R,T)` coarse blocks in
[the coarse spectrum](bounded-spin-coarse-twelve.md); the `R` labels of
the eight primaries are not determined by the coefficient table, because
`R=C_0` is the coefficient unit. Consequently this note does **not**
certify a split of any individual joint coarse block; that requires the
`R` characters of the same primaries, or a joint `R,T,S` primary
computation.

This is still only a splitting certificate. It does not give the exact
dimensions of the joint `T,S` components, their ranks, a two-cut rank
product, or a factoring exponent. It uses the cut-trace certificate's
identification of the residues as ranks modulo five; the ten-site control
for that identification is checked in
[the diagonal trace note](bounded-spin-diagonal-trace.md).

## What changed

The previous separator checkpoint established the character separation
at the coefficient-algebra level. This note transfers that separation to
the `T`-eigenspace partition of the cut module. All three `T` groups are
genuinely split by the separator, not merely distinguished as algebra
characters. The former `T=0` obstruction is removed by an explicit
nonzero diagonal entry of idempotent `0`; its rank is a positive multiple
of five, not forced to zero.

The next target is the `R` character of each primary, which would turn
these `T`-group splits into splits of the six joint `(R,T)` coarse blocks.
Exact component dimensions and the two-cut pairing remain open.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
make verify-bounded-spin-twelve-separator-split
```

The checker reads the RT coefficient receipt and the cut trace-form
receipt, checks their common basis, verifies the three `T` groups, and
evaluates the two explicit diagonal entries that make idempotent `0`
nonzero. It records all six `T`-group splits and leaves the `R` refinement
open. No full cut operator or exact rank computation is performed.
