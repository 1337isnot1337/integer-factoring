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

## Adjacent 11-mode diagnostic

Dropping mode 11 from the q96 map while retaining all remaining core
words gives a valid q94 map on 11 modes. Its cutoff-47 fit is `0.709929`,
which is *above* the q99 value `0.709858` at the same cutoff; the earlier
comparison to the full q99 profile was cutoff-biased. The q45 quotient
shows that core size alone is not the mechanism: a smaller
deletion-closed core can lose more pressure than the removed capacity
saves.

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
