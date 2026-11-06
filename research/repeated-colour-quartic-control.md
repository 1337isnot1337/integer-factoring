# Repeated-colour quartic control at q8

## Result and scope

The [colour-asymmetric order-twelve floor](colour-asymmetric-order12-floor.md)
uses distinct-colour quartet gates.  The remaining local gate that can act
on the two-colour diagonal `P=Q` block is a four-vertex support with two
colours repeated twice, `a,a,b,b`.  This note records a finite q8 control
for the strongest cancellation form of such a gate.  It is not a theorem
about all repeated-colour gates and does not recompute pressure.

Replace the pair-product leading term on a repeated two-colour support by

```math
H=\text{pair-product}\cdot(1+x).                       \tag{1}
```

This has leading order twelve, endpoint zero, and shifts the leading
pair-product contribution to order fourteen.  Both cuts are rebuilt with
all eight charge blocks and all occupations.  The exact low prefixes are:

| control | low Smith prefix |
| --- | --- |
| canonical distinct gates | `0^1,6^14,12^35,15^84,16^70,17^14,18^208` |
| repeated gate plus canonical distinct gates | `0^1,6^14,12^35,14^56,16^98,17^10,18^208` |
| total-four cancellation, distinct gates disabled | `0^1,6^14,12^91,18^308` |

The order-twelve increment stays `35` with the canonical distinct gates.
The repeated gate moves fifty-six channels to order fourteen and changes
the returning layers, but it does not lower the protected diagonal block.
With the distinct gates disabled, the same cancellation leaves the
Gaussian order-twelve channel `91`; the leftover channels are one-particle
and off-diagonal, not the protected four-occupation block.

## Mechanism and limits

A repeated `a,a,b,b` gate can act on the diagonal `P=Q` block, so it is the
one local support pattern not covered by the distinct-colour protection.
The finite control shows that cancelling the diagonal pair product is not
enough to lower the prefix: the canonical distinct gates already remove the
off-diagonal/one-particle excess, and the protected rank-four block is not
the whole order-twelve channel.  The total-four cancellation alone leaves
the Gaussian `91`, which is higher than the canonical `35`.

Repeated-colour gates preserve separate colour parity.  The q>=34
parity-preserving sector obstruction therefore still gives pressure at
least `3/4`; a q8 repeated-colour gain would be finite only and would not
improve the asymptotic complete bound.  This control closes the last
local gate pattern at q8 but does not supply a uniform repeated-gate
theorem or a different leg grouping, kernel, or weighted transfer.

## Reproduction

```sh
make verify-repeated-colour-quartic-control
```

The target runs the checker and `make verify-docs`.  The
[checker](../tools/check_repeated_colour_quartic_control.py) asserts both
cuts and all three tabulated prefixes; the
[receipt](../evidence/repeated-colour-quartic-control.json) records the
six axis-control cases.  No complete-bound score changes; it remains
`1/2`.
