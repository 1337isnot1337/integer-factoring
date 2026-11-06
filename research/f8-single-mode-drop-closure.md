# Every single-mode drop of the q96 F8 map is excluded

## Construction

The q96 F8 hyperoval map uses twelve hot modes and 96 paired letters.
For each mode `m=0,...,11`, delete every letter whose embedding or
retraction support contains a word using bit `m`. The surviving map is
a principal submap: `LP=I` restricts, and the deletion-closure check is
revalidated on the kept vectors. All twelve resulting maps are valid.
Mode `11` drops exactly the two letters at indices `90,91`, reproducing
the q94 pilot.

## Cutoff-47 lower bounds

Each map was resolved at cutoff 47. The known-order fit is a **lower**
bound on the complete pressure: unresolved orders only add positive
terms to `Z`, so the full pressure cannot be smaller. The table is
ordered by this lower bound.

| dropped mode | q | known-order lower `rho` | tail-augmented upper `rho` |
| ---: | ---: | ---: | ---: |
| 11 | 94 | `0.7099289418715207` | `0.7107270882301813` |
| 10 | 92 | `0.7100127987350093` | `0.7107413220852488` |
| 9 | 90 | `0.7101737073892955` | `0.7107987943356104` |
| 8 | 86 | `0.7104656188170568` | `0.7109333918728320` |
| 7 | 82 | `0.7110242451287223` | `0.7112918662843903` |
| 6 | 78 | `0.7119268764551838` | `0.7120332352446683` |
| 5 | 74 | `0.7134302643763426` | `0.7134512883792019` |
| 4 | 70 | `0.7159835082890482` | `0.7159849957056230` |
| 3 | 64 | `0.7207314770574176` | `0.7207314975692547` |
| 2 | 60 | `0.7313380098989488` | `0.7313380099123261` |
| 1 | 54 | `0.7656071830015263` | `0.7656071830015263` |
| 0 | 49 | `0.8182248402897341` | `0.8254155465585895` |

Every lower bound is above `1/sqrt(2)`. The best single-mode drop is
mode `11` at `0.7099289419`, still `2.82e-3` above the threshold. Hence
**no single-mode principal-submap drop of the q96 map can cross**. The
lower bounds get worse as more modes are removed; the q49 mode-0 drop
is far above.

This sharpens the earlier q94 pilot: mode `11` is not merely one tested
drop, it is the best of all twelve, and its lower bound excludes it.
It also illustrates why cutoff-31 screens are not decisive: mode `7`
looked best at cutoff 31 (`0.7147345`) but its cutoff-47 lower bound is
`0.7110242`, worse than several other drops.

## What remains

The census covers single-mode deletions only. Multi-mode deletions,
coordinate substitutions, a different rank-three image outside the
height-25 theorem, and a different ambient tensor remain outside it.
The q96 core words are load-bearing: deleting low letters raises the
pressure, and the q45 minimal core is much worse. The next low-block
search should **change** the surviving coordinates or the ambient
structure, not merely delete letters from the q96 map.

The complete factoring exponent remains `1/2`; no score, exact pressure
certificate, or uniform transfer is claimed.

## Reproduction

```sh
make verify-f8-single-mode-drop-census
```

The target rebuilds all twelve principal submaps, validates `LP=I` and
deletion closure, and replays the saved cutoff-47 order hashes and fits.
It does not regenerate the profiles. The source-pinned receipt contains
the orders, heights, tails, hashes and fits for every mode.
