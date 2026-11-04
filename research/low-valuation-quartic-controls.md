# Lower-valuation quartic controls at q8

## Result and scope

The [colour-asymmetric order-twelve floor](colour-asymmetric-order12-floor.md)
assumes distinct-colour quartic gates of valuation at least twelve.  The
remaining local escape is to lower that valuation.  This note records a
finite q8 control for the natural determinant-shaped family.  It does not
prove a lower-valuation theorem and does not recompute a pressure
certificate.

Keep the spatial-mode support rule and alpha matching, and replace the
three outer exponents of the quartic coefficient

```math
(12,20),\ (16,16),\ (20,12)                            \tag{1}
```

by

```math
(a,32-a),\ (16,16),\ (32-a,a),\qquad a=6,9,10,11.      \tag{2}
```

The middle exponent, degree `32`, endpoint cancellation and support rule
are unchanged.  For each value of `a`, both cuts are rebuilt with all
eight charge blocks and all occupations.  The exact low prefixes are:

| coefficient | low Smith prefix |
| --- | --- |
| `a=12` canonical | `0^1,6^14,12^35,15^84,16^70,17^14,18^208` |
| `a=6` | `0^1,6^70,12^283,15^84,16^98,18^184` |
| `a=9` | `0^1,6^14,9^56,12^35,15^248,16^98,18^86` |
| `a=10` | `0^1,6^14,10^56,12^35,16^346,18^26` |
| `a=11` | `0^1,6^14,11^56,12^35,16^98,17^248,18^26` |

The order-twelve increment is `35` for `a=9,10,11` and rises to `283`
for `a=6`; it is never below the canonical `35`.  Lowering the valuation
moves mass into orders `6`, `9`, `10` or `11` and changes the returning
layers.  It does not lower the protected order-twelve block.

## Why this is a relevant control

The Gaussian `P=Q` order-twelve block is independent of the quartic
coefficient.  A distinct-colour quartic cannot act on that block, which
uses only two colours.  Lowering the quartic valuation therefore adds
lower-order off-diagonal and returning channels; it does not supply a
term that can cancel the diagonal Gaussian rank-four block.  The finite
prefixes are consistent with that mechanism: the canonical order-twelve
count is unchanged for `a=9,10,11`, and the very low `a=6` case adds
further order-six and order-twelve channels.

The control does not cover arbitrary lower-valuation coefficients,
repeated-colour quartic gates, a different leg grouping or the weighted
transfer.  It does not show that every lower-valuation gate is harmless;
it shows that this natural determinant-shaped lowering does not improve
the fixed q8 profile.

## What is and is not proved

- This is a finite q8 control for four coefficient choices, not a
  theorem about all lower-valuation gates.
- The pressure certificate is not recomputed: additional low orders
  would only be left in the conservative tail, and the result is not a
  new crossing or non-crossing claim beyond the tabulated prefixes.
- The next changed assumptions are repeated-colour gates, a different
  leg grouping, a changed pair kernel, or the weighted/unequal-scale
  transfer.  No complete-bound score changes; it remains `1/2`.

## Reproduction

```sh
make verify-low-valuation-quartic-profiles
```

The target runs the checker and `make verify-docs`.  The
[checker](../tools/check_low_valuation_quartic_profiles.py) asserts both
cuts and all tabulated prefixes; the
[receipt](../evidence/low-valuation-quartic-profiles.json) records the
ten axis-coefficient cases.
