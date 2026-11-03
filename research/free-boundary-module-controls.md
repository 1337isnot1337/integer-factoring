# Free boundary vectors do not lower the q=2 module controls

## Result and scope

The [cancellation module automaton](cancellation-module-smith-automata.md)
uses a fixed terminal vector by default.  Releasing both boundary
vectors (`free_boundaries=True`) does not expose a subthreshold `q=2`
certificate in the recorded controls.

With free boundaries, the best sampled values are:

| family | best `B` | scale |
| --- | --- | --- |
| elementary scalar | `256/27` | `3` |
| ordered pairing | `16` | `1` |
| mixed nilpotent loop | `16` | `1` |
| explicit module crossing | `6153085/408942` | `3/2` |

The explicit crossing value is `15.046...`; the fixed-boundary value was
`1625/108=15.046...`.  All four remain above the `q=2` module target
`2^(2*sqrt2)=7.10299...`.  The search samples precisions `1..6`, the two
cuts, both sides, the three recorded heights, and scales
`1,3/2,2,3,5`.

This is a finite control result, not a theorem.  It shows that the
fixed-terminal masking obstruction is not the only reason the q=2
module lane stalls: releasing the boundary vectors alone does not
produce a crossing.  A new local alphabet, a different channel algebra,
or a genuine gauge/recombination is still required.

## Verification

`python3 tools/check_free_boundary_module_controls.py` reproduces the
table and the exact fractions.  The [saved record](../evidence/free-boundary-module-controls.json)
contains the same data.
