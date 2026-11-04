# A minimal deletion-closed F8 core does not lower the pressure

## Question

The q96 F8 hyperoval map retains 81 height-below-24 core words, but its
ten height-25 boundary triples only delete through a much smaller set.
Can the unused core words be removed, reducing the physical capacity and
the numerical natural pressure?

## Exact quotient

The ten boundary words are

```math
52,76,82,97,138,145,262,265,517,1027,
```

with mode sets

```math
(5,9,11),(5,7,13),(3,9,13),(1,11,13),(3,7,15),
(1,9,15),(3,5,17),(1,7,17),(1,5,19),(1,3,21).
```

Their deletion closure has 52 words: the empty word, the 41 lower
subsets, and the ten triples. The q96 map has single-word coordinates for
the 42 lower words but not for the ten triples; the triples occur only
inside the three non-single hyperoval embeddings. Keeping those three
embeddings and their retractions, and the 42 lower singles, gives a
deletion-closed q45 map on 11 modes with maximum height 25 and height sum
696. The construction is the smallest deletion-closed coordinate set
containing the fixed hyperoval boundary image: any such set must contain
the 42 lower subsets, and the three non-single embeddings are needed for
the rank-three image.

The map passes `LP=I`, deletion closure, and the exact construction
checks. Its complete natural profile at cutoff 127 has `2025` orders per
cut and no unresolved tails. The numerical natural pressure is

```math
\rho \approx 0.7156911629839379,
```

which is above the q96 value `0.7100057871974` and the q99 height-39
value `0.7099930431440509`, and above `1/sqrt(2)` by about `0.008584`.
Removing the 52 q96 core words that are not needed by the boundary image
therefore makes the pressure worse, not better. The unused core words are
load-bearing for the pressure balance.

## Why the low core is saturated

The q96 map already contains all eleven height-24 words, and the only
height-25 words are the ten triples, which are the fixed hyperoval
boundary image and lie in the span of the three non-single embeddings.
There is therefore no missing low word below height 26 that could be
added to the core. A four-mode height-26 word such as `46=(3,5,7,11)`
does couple to the boundary triples, but its cutoff-47 fit is
`0.7098643`, slightly above the q99 value `0.7098582` at the same
cutoff. The remaining freedom is high-word extension, which is the
plateau already seen at q99 and q36.

## Matched-cutoff low-block diagnostics

Cutoff-47 fits are comparable only at the same cutoff. The complete q99
profile is `0.7099930`, but its cutoff-47 fit is `0.7098582`; comparing a
partial fit to a completed profile is biased. At cutoff 47:

| map | change from q96 | fit |
| --- | --- | ---: |
| q96 | none | `0.7098760` |
| q99 | add `160,192,224` | `0.7098582` |
| q94 | drop mode 11 | `0.7099289` |
| q95 | drop word `31` | `0.7098923` |
| q97 | add word `46` | `0.7098643` |
| q97 | q94 plus `160,192,224` | `0.7099123` |

The q99 extension is the best at the matched cutoff. Dropping mode 11,
dropping the redundant height-25 word, adding a boundary-coupled
height-26 word, and composing the mode drop with the q99 extension are
all worse. The q94 lower-bound argument
([mode-11 pilot](../research/f8-hyperoval-q94-mode11-pilot.md)) is
consistent with this: its complete pressure is at least the cutoff-47
fit, which is above the q99 full value. The q45 quotient shows the
opposite failure mode: deleting all unused core words is much worse.
The low core is saturated
([admissible high words](f8-admissible-high-words.md)), and the fixed
height-25 boundary and height-24 core leave only high-word extensions,
which stay on the plateau. This closes the tested low-block variants,
not the F8 family as a mathematical impossibility.

## Scope

The profile is a complete finite computation for the specified q45 map.
The pressure is a numerical fit, not an exact rational certificate; the
complete factoring exponent remains `1/2`. The map has the fixed F8
hyperoval image and the common balanced natural cut scale. No uniform
physical type, transfer, global optimality, or score change is claimed.

Reproduce with

```sh
make verify-f8-minimal-core
```

The [receipt](../evidence/f8-hyperoval-minimal-core.json) contains the
closure words, source indices, complete orders, source pins, and
numerical pressure. The fast checker replays the map and recorded orders;
it does not regenerate the profile. The adjacent q99 result is in
[the height-39 extension note](f8-hyperoval-height39-extension.md).
