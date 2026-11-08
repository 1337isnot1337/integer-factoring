# A height-gauge q115 map numerically crosses the two-cut target

## Question

The q96 F8 hyperoval map has a mode `11` that no height-25 boundary
triple uses. Its standard height is `23`. Can the height of this unused
mode be changed as a gauge, and can the induced low block cross
`1/sqrt(2)` without changing the boundary?

## Construction

Start from the complete q96 profile. Assign the mode weights

```math
W=(1,3,5,7,9,11,13,15,17,19,21,11),
```

so mode `11` has height `11` instead of `23`; mode `5` also has height
`11`. The ten boundary triples do not use mode `11`, so they remain at
height `25` and the boundary image is unchanged. Under `W`, nineteen
words of height below `25` are absent from the q96 support. Add each as
a single-word coordinate to both the embedding and retraction:

```text
2050,2051,2052,2053,2054,2055,2056,2057,2058,2059,2060,2061,
2064,2065,2066,2067,2080,2081,2112.
```

The resulting q115 map has all `110` words of weighted height below
`25`. It satisfies `LP=I`, deletion closure, and height homogeneity for
`W`. The standard-height validation also passes; the added words are
not free in the standard grading, but they are free in the new grading.

## Result

At cutoff `47`, each cut retains `935` exact orders and leaves `12290`
unresolved tails. The known-order fit is a lower bound; assigning every
unresolved tail order `48` gives an upper bound. The numerical fits are

```math
\rho_{\rm known}=0.6999902708909971,
\qquad
\rho_{\rm upper}=0.7067324702960336
< \frac1{\sqrt2}=0.7071067811865475.
```

The upper margin is `3.7431e-4`. Thus the complete natural two-cut
pressure of this q115 map is numerically below the target. For
comparison, the q96 map at the same weights but without the nineteen
added words has complete numerical pressure `0.7064473232`; the q115
map is the maximal low-core version and is the better transfer
candidate. The standard q96 pressure is `0.7100058347`.

## Why this escapes the floor

The common-core/Jensen floors assume the standard mode heights and the
literal height-below-24 core. The new grading changes the low block and
gives mode `5` and mode `11` the same height. The q115 map retains the
full weighted low core, so this is not a truncation: it is a different
protected structure with a different height function. The floor
theorems do not apply as stated. The q94 deletion and the q96 low-block
pilots remain above target; the gauge does not delete letters.

## Height sensitivity and the duplicate

The crossing is sensitive to the mode-11 height. Replaying the same
q115 orders with a different `w11` gives the following tail-augmented
upper pressures:

| `w11` | weighting | upper `rho` | below target? |
|---:|---|---:|:---:|
| 11 | duplicate of mode 5 | `0.7067324703` | yes |
| 12 | even height | `0.7082128683` | no |
| 13 | next odd distinct | `0.7091841819` | no |
| 14 | even height | `0.7098097697` | no |
| 15 | odd distinct | `0.7102083474` | no |

The q96 map without the nineteen added words has complete pressure
`0.7064473232` at `w11=11`; it is not low-core maximal. The q115 map
adds the missing weighted-low words, which raises the pressure to the
quoted candidate, but makes the low block maximal under `W`.

Only the duplicate `w11=11` crosses. The next odd distinct choice
`w11=13` is above target. Thus the candidate requires two physical
letters of height `11` and no letter of height `23`; this is not the
standard one-letter-per-odd-height height polynomial `Q_H` used in the
protected-core floor. The floor does not exclude it, but duplicate
heights must be justified in the physical two-cut map and in the
uniform transfer. If the transfer excludes repeated heights, this
candidate is blocked.

## Open obligations

This is a **numerical candidate crossing of the natural two-cut
pressure**, not a score. The remaining checks are:

- establish whether repeated mode heights are admissible; the only crossing has duplicate height 11, while the next odd distinct height 13 fails;
- resolve the complete q115 profile if the duplicate is accepted (the
  cutoff-47 upper already lies below target, so this is confirmation);
- prove that the height assignment `W` is an admissible physical gauge,
  including the duplicate height `11` and the uniform transfer with all
  costs;
- check the transfer from the two-cut product to the ordered factoring
  algorithm; no complete exponent change is claimed.

The q115 map has `115` letters rather than `96`, so the cost comparison
is not automatic.

## Reproduce

```sh
make verify-f8-q115-height-gauge-crossing
```

The checker rebuilds the q115 map from the q96 receipt, validates
`LP=I`, deletion closure, the weighted low core, and height homogeneity,
then replays the saved cutoff-47 orders and both fits. The receipt is
[`evidence/f8-q115-height-gauge-crossing.json`](../evidence/f8-q115-height-gauge-crossing.json).
