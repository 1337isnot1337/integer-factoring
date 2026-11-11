# The twelve-site separator splits one coarse block modulo five

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

Grouping by the `T` eigenvalue gives

```text
T=0:  idempotents 0,6         row residues 0,0       column residues 0,0
T=1:  idempotents 1,2,5,7     row residues 1,0,0,4   column residues 4,0,0,1
T=4:  idempotents 3,4         row residues 1,0       column residues 1,0
```

A nonzero rank residue certifies that the corresponding primary image is
nonzero. Therefore `T=1` splits into two row components and two column
components. The `T=4` group has only one nonzero residue in each
orientation, so its split is **not** certified. The `T=0` group has zero
residues on both idempotents in both orientations, so it is completely
unresolved: each rank is either zero or a positive multiple of five.

This is a modulo-five splitting certificate. It does not give the exact
dimensions of the joint `T,S` components, their ranks, a two-cut rank
product, or a factoring exponent. It uses the cut-trace certificate's
identification of the residues as ranks modulo five after the basis and
unit normalization.

## Correction

An earlier version of this note certified splits of both `T=1` and `T=4`,
and treated the `T=0` ambiguity as confined to idempotent `0`. Those
conclusions used the pre-correction cut-trace residues, which came from
evaluating the trace forms at coordinate `k` instead of at coefficient
label `basis[k]`. The corrected residues remove the `T=4` split and put
both `T=0` idempotents in the ambiguous class. The separator still gives
two certified `T=1` components, so the coefficient-level separation is
not lost; its cut-module consequences are just weaker than previously
stated.

## What changed

The previous separator checkpoint established the character separation
at the coefficient-algebra level. This note transfers the one nonsingular
part of that separation to the cut-module coarse partition: the `T=1`
block is genuinely split by the separator, not merely distinguished as an
algebra character. The unresolved `T=0` and `T=4` groups are now specific
targets: decide the exact ranks of idempotents `0,6` (both multiples of
five) and idempotents `3,4` (the dual block `4` is still completely
ambiguous).

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_separator_split.py
make verify-bounded-spin-twelve-separator-split
```

The checker reads the RT coefficient receipt and the corrected cut
trace-form receipt, checks their common basis, verifies the three `T`
groups, and records the one certified split and the unresolved groups. No
cut operator or exact rank computation is performed.
