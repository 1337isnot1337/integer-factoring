# All disjoint f3 placements span the five-label graded crossing

## Exact finite result

The [graded seam construction](invertible-fusion-graded-seams.md) used
free prefixes followed by disjoint `f_3` triples. Its unfiltered span
contains the crossing at two labels but misses it at five and thirteen
labels. The five-label failure is not caused by the prefix placement.

Enumerate every set of disjoint three-strand `f_3` placements on the
`2n` central sites, including triples that cross the middle seam. For
each placement take every boundary-degree coefficient tensor, both cut
orientations, and the two identity pairings. Exact Gaussian elimination
over `F16` gives:

| half-length `n` | labels | disjoint placements | nonzero pieces | all-placement span | contains crossing | prefix span | prefix contains |
| ---: | ---: | ---: | ---: | ---: | :---: | ---: | :---: |
| 4 | 2 | 13 | 68 | 13 | yes | 13 | yes |
| 6 | 5 | 60 | 674 | 475 | **yes** | 80 | no |

At five labels the enlarged span has dimension `475` rather than `80`
and contains the crossing. An explicit recombination uses `442`
nonzero coefficient pieces; the checker reconstructs it and verifies
the identity entry by entry.

This is a finite exact statement over the fixed field. It does not
assert a uniform family, a rank-sum optimum, or a stronger factoring
score. The `442`-piece recombination is far above the five-label
finite-block budget `5^(1/sqrt(2)) < 4`; it is a spanning certificate,
not a cost certificate.

## What changed in the scope

The parent note correctly reported that the **prefix** span misses the
five- and thirteen-label crossings. The present census shows that at
five labels the failure is removed by non-prefix disjoint placements.
It does not address thirteen labels, where the full placement census
has `277` disjoint placements and was not run here.

The remaining constructive target is therefore sharper:

1. run the all-placement census at thirteen labels;
2. if the span still misses, add genuinely fused projectors rather than
   more disjoint `f_3` placements;
3. if it contains the crossing, find a grouped representation with
   rank sum below the finite-block threshold, not merely a large
   spanning recombination.

The last step is the only one that could affect the complete bound.

## Reproduction

```sh
python3 tools/check_all_disjoint_fusion_placements.py \
  --output evidence/all-disjoint-fusion-placements.json
make verify-all-disjoint-fusion-placements
```

The receipt records the placement counts, exact span dimensions,
membership decisions, the explicit recombination, and source hashes.
No new Lean formalization or independent mathematical acceptance is
claimed.
