# Short-flag rank-one-row-cut units for the seven-letter quadratic inverse

**Current status: historical structural screen, not an active score route.**
The later [box-minor theorem](bare-seven-to-twelve-letters.md) closes the
bare criterion through twelve physical letters, and the
[adjugate constraint](adjugate-seven-letter-pressure-floor.md) excludes the
asymmetric profile.  The unit screen below records a structural fact about
the quadratic-inverse flag; it does not reopen the bare `q=7` route.

## Result and scope

The [quadratic-inverse screen](inverse-degree-crossing-compatibility.md)
starts from a `q=7` transpose-like unit whose kernel flag has maximum
level 16.  Its best degree is 24 and its best first-jet multiplicity is
20, so the first-jet pressure excludes that family.

The unit's row cut has rank one.  Writing the row-cut constant as
`a b^T` gives a complete rank-one parameterization of nearby units:
`U = row_cut^{-1}(a b^T)`.  The base unit has

```text
a = 0x1010101010101,   b = 0x1030303030303.
```

Flipping bit 17 or 25 of `b` preserves invertibility and the rank-one
row cut, and lowers the nilpotency index of `I-U` from 16 to 8.  The
single-flip screen covers 98 flips; 42 are invertible and nilpotent, and
the minimum index is 8.  No single flip reaches index 7.

For the bit-17 witness the kernel flag has maximum level 8.  The
quadratic constructor has degree `2*(8-1)-k = 14-k` when `k` levels
receive weight one:

| one-weight levels `k` | degree | best first-jet multiplicity |
| --- | --- | --- |
| 1 | 13 | 12 |
| 2 | 12 | 18 |

The degree-13, multiplicity-12 point is new: the first-jet pressure
`(49/s^2)(1+12/s)` at `s=3/2` is below the seven-letter threshold, so
this family is not excluded by the first-jet bound.  It is not a
crossing, however: the `D=1` physical matrix has `K(0)=U` invertible,
so its zero cut has all 49 Smith orders zero.  The target profile's
positive zero orders require auxiliary dimension greater than one.

Thus the short-flag unit is an infinity-side building block, not a
complete local crossing.  A `D>1` construction would be needed to
realize the historical profile, but the bare criterion is already
closed; this screen is retained as a structural record.  Further score
work must target a different transfer criterion or gauge, with its own
analysis.

## Verification

`python3 tools/check_short_flag_seven_letter_units.py` checks the base
rank-one factorization, the 98 single flips, invertibility, nilpotency,
the minimum index, the two candidate degrees and first-jet ranks, and
the all-zero `D=1` zero cut.  The [saved record](../evidence/short-flag-seven-letter-units.json)
reconstructs the same data.
