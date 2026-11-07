# A second mode deletion does not rescue the q94 low-block drop

## Question

The q96 F8 map has twelve modes. Dropping mode `11` gives the q94 map;
its complete pressure `0.7100458192` is above q96 and does not cross.
Does deleting a second mode after mode `11` remove enough load-bearing
high orders to lower the pressure below the threshold? This tests
whether the deletion descent continues past the best single-mode drop.

## Screen

For each `m=0,...,10`, delete every q96 letter whose embedding or
retraction support uses mode `11` or mode `m`. All eleven resulting maps
validate `LP=I` and deletion closure. Resolve each at cutoff 47; the
known-order fit is a lower bound on the complete pressure.

| pair | q | known-order lower `rho` | tail-augmented upper `rho` |
| --- | ---: | ---: | ---: |
| `(11,10)` | 90 | `0.7100573927808442` | `0.7107399247277675` |
| `(11,9)` | 88 | `0.7102142859093288` | `0.7107984024550345` |
| `(11,8)` | 84 | `0.7104980448425247` | `0.7109342548869693` |
| `(11,7)` | 80 | `0.7110373272685115` | `0.7112905154287581` |
| `(11,6)` | 76 | `0.7119376483780790` | `0.7120365482422370` |
| `(11,5)` | 72 | `0.7134341326069014` | `0.7134536635684312` |
| `(11,4)` | 68 | `0.7159842656391489` | `0.7159856479478495` |
| `(11,3)` | 62 | `0.7207315352001564` | `0.7207315541031940` |
| `(11,2)` | 58 | `0.7313380106127316` | `0.7313380106249538` |
| `(11,1)` | 52 | `0.7656071830015672` | `0.7656071830015672` |
| `(11,0)` | 48 | `0.8190614097976898` | `0.8256215650891732` |

Every lower bound is above `1/sqrt(2)`. The best pair is `(11,10)` at
`0.7100574`, which is already above the **completed** q94 pressure
`0.7100458`. Adding the second deletion therefore makes the pressure
worse even at the lower-bound level; no pair can cross.

## Interpretation

The low-block deletion direction is monotone in the wrong direction:
the first deletion (mode 11) was the mildest, and every second deletion
is worse. This agrees with the q45 minimal-core negative and with the
load-bearing role of the removed height-23 words in the completed q94
profile. Further mode deletions are not a crossing route. A useful
low-block change must substitute coordinates or change the ambient
structure, not delete another mode.

## Scope

The screen covers the eleven pairs containing mode 11 only. Pairs not
containing mode 11, triples, arbitrary coordinate substitutions, and
complete profiles are outside it. The fits are numerical lower bounds,
not exact certificates; every one is above the threshold, so the
corresponding complete maps cannot cross. The complete factoring
exponent remains `1/2`.

## Reproduce

```sh
make verify-f8-two-mode-drop-screen
```

The checker rebuilds the eleven principal submaps, validates the maps,
and replays the saved cutoff-47 order hashes and fits. It does not
regenerate the profiles. The receipt is
[`evidence/f8-two-mode-drop-screen.json`](../evidence/f8-two-mode-drop-screen.json).
