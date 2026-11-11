# The twelve-site primaries are not R-pure modulo five

## Result

The two cut trace orientations differ by the replica swap `R`.  On the
twelve-site cut module the row trace is the ordinary idempotent trace and
the column trace is the `R`-twisted trace:

```text
row residue    = tr(e_i)     = r_i^+ + r_i^-   (mod 5)
column residue = tr(e_i R)   = r_i^+ - r_i^-   (mod 5),
```

where `r_i^+` and `r_i^-` are the dimensions of the `R=+1` and `R=-1`
parts of the image of the primitive central idempotent `e_i`.  Since
`R^2=1` and `e_i` commutes with `R`, the two multiplicities are recovered
modulo five by

```text
r_i^+ = (row_i + column_i)/2,   r_i^- = (row_i - column_i)/2.
```

The [cut trace certificate](../evidence/bounded-spin-twelve-cut-trace-form.json)
gives the residues; the resulting split is

```text
i  T  rank  tr(eR)  r+  r-  status
0  0    0      0     0   0  undetermined
1  1    3      1     2   1  mixed
2  1    1      0     3   3  mixed
3  4    2      1     4   3  mixed
4  4    2      1     4   3  mixed
5  1    1      3     2   4  mixed
6  0    2      1     4   3  mixed
7  1    3      4     1   2  mixed
```

Thus every nontrivial primitive central idempotent has **both** `R` signs
with nonzero multiplicity modulo five.  Idempotent `0` is the only
undetermined case: its rank is a positive multiple of five, and both of
its `R` residues vanish, so its `R` split is a pair of multiples of five
(or a pure multiple of five in one sign).

The sum of the eight idempotents is the coefficient unit `C_0`, not the
cut-module identity; consequently these residues are ranks of the
primary **images**, not the dimensions of the six joint `(R,T)` coarse
blocks.  In particular the group sums below are not coarse dimensions:

```text
T=0: rank residue 2, R+ residue 4, R- residue 3
T=1: rank residue 3, R+ residue 3, R- residue 0
T=4: rank residue 4, R+ residue 3, R- residue 1
```

## Control and scope

At ten sites the same calculation is an exact control.  The row residues
equal the seven published primary dimensions modulo five, and the column
residues equal `+d_i` for `34,1,225,300,35,155` and `-d_i` for `406`.
The signs are the known joint `(T,R)` coarse labels: the two `T=2`
primaries `300,35` lie in the `R=+` coarse block `335`, while `406` is
the `T=0,R=-` primary.  This fixes the orientation of the two forms and
recovers the known pure `R` signs of the ten-site primaries.

The column form is not the **column-cut rank**.  The exact ten-site
column-cut rank pairs are `[931,1156,356,581,771,776,715]` in
[the idempotent-S control](../tools/check_bounded_spin_idempotent_s.py),
whose residues do not match the column trace residues above.  The two
forms are the row and column diagonal orientations of the same row-cut
matrix; the column cut is a different operator.

At twelve sites the result is only a modulo-five congruence.  It does
not give exact dimensions, exact row or column ranks, a rank product, or
a two-cut pairing.  It does show that the central primaries cannot
refine the joint `(R,T)` coarse blocks into pure `R` pieces: seven of the
eight primaries are mixed, and the remaining one is undetermined.  The
noncentral separator action, or a direct exact rank computation, is
still required.

No factoring exponent changes; the complete bound remains `1/2`.

## What changed

The [separator note](bounded-spin-twelve-separator-split.md) left the
`R` character of each primary open.  This note computes that character
modulo five and answers it negatively: the primaries are mixed, not
pure.  This is a structural obstruction to using the central idempotents
alone to split the six `(R,T)` coarse blocks.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_R_primary_split.py \
  --output evidence/bounded-spin-twelve-R-primary-split.json
make verify-bounded-spin-twelve-R-primary-split
```

The checker recomputes the ten-site control, reads the twelve-site cut
trace and separator receipts, solves the two R-sign residues for each
primary, and records the group sums and the mixed/undetermined status.
