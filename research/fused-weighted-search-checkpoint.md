# Completed fused gauge search and physical-basis controls

## Scope

The completed search below explores complete physical-letter gauges and
paired types for four fused seeds. It improves their retained numerical
mixed-place bounds. It does not prove a global optimum, a new uniform
crossing certificate, or an improved integer-factoring exponent.

The search keeps the degree charge, all physical type shifts, and the
entropy denominator. Its local orders are computed with exact finite-field
arithmetic. Floating optimization selects the continuous parameters.
An independent replay then replaces each retained decimal parameter by
the corresponding exact rational number and checks all winning Smith
orders again. The logarithmic cost recomputation remains floating point.

## 1. A completed thirty-minute run

The declared wall-time budget was 1,800 seconds; the recorded elapsed
time was 1,800.7386 seconds. The process attempted 1,295 trials and
completed 1,294. The final trial was interrupted at the budget. Every
reported per-alphabet winner came from a completed trial.

The run used the three retained best unweighted place partitions for
each seed, optimized common paired physical types, and alternated linear
and unrestricted physical-letter heights. It also varied allocations and
included inactive places in some trials. The objective is the fixed
allocation criterion from the
[mixed-place theorem](mixed-place-hermite-crossings.md), with the
[weighted local oracle](fused-weighted-place-profiles.md).

| Alphabet | Seed `(m, theta)` | Completed linear trials | Completed full-height trials | Best retained numerical `rho` | Winning trial |
|---:|:---:|---:|---:|---:|---:|
| 3 | `(2, 2)` | 162 | 162 | 0.948414858180 | 300 |
| 4 | `(3, 4)` | 162 | 162 | 0.964784599892 | 437 |
| 5 | `(4, 4)` | 162 | 161 | 0.971769760852 | 38 |
| 6 | `(5, 4)` | 162 | 161 | 0.976475588576 | 639 |

All four winning trials used unrestricted physical-letter heights.
They are exploratory upper-bound candidates for the stated criterion.
Neither unsuccessful optimizer termination nor exhausted restart coverage
is evidence of a lower bound on the possible optimum.

The compact record retains every trial outcome, elapsed time, number of
oracle queries, optimizer-stage status, best cost, and a hash of its best
parameter vector. It retains the complete parameters and local profile
records for each per-alphabet winner. The redundant full run is identified
by its SHA-256 digest and byte length. The compact record is 363,948 bytes;
it avoids storing repeated nonwinning profiles in the repository.

The replay checks 29 complete local profiles containing 697 channels.
Every exact rational order agrees with its recorded floating order to
the stated tolerance, and every recomputed full cost agrees with the
retained objective. It also checks disjoint place supports and excludes
the crossing endpoint. Rational paired types are normalized exactly
before the physical mean shifts are recomputed.

The sources and retained records are:

- [`search_fused_mixed_weighted.py`](../tools/search_fused_mixed_weighted.py)
- [`compact_fused_weighted_search.py`](../tools/compact_fused_weighted_search.py)
- [`check_fused_weighted_winners.py`](../tools/check_fused_weighted_winners.py)
- [`fused-mixed-weighted-summary.json`](../evidence/fused-mixed-weighted-summary.json)
- [`fused-mixed-weighted-winner-replay.json`](../evidence/fused-mixed-weighted-winner-replay.json)

Replay the retained winners with:

```sh
python3 tools/check_fused_weighted_winners.py \
  evidence/fused-mixed-weighted-summary.json \
  --output evidence/fused-mixed-weighted-winner-replay.json
```

## 2. A physical-basis oracle for further searches

The full physical-basis adapter acts on the four tensor legs by

\[
 A_x\,B_y\,B^{-T}_u\,A^{-T}_v.
\]

It preserves the strict crossing endpoint. In the source matrix
convention it conjugates `F` by `B tensor A`. Each transformed seed is
checked at the endpoint before any local profiles are evaluated.
The two cost cuts are then rebuilt in full. They need not retain the
charge-block structure of the original fused matrix.

The fast path restricts scalars from `F16` to `F4`. With the same physical
weights on the two coordinates of each field element, every complete
Smith order occurs twice. The log partition sum therefore subtracts
`log(2)`. Its physical-height subgradient is obtained from the selected
rows and columns of the complete pivot chain. Its local-scale
subgradient uses the increments of the base determinantal valuations.

The independent exact path uses the
[fraction-free local Smith oracle](compressed-fusion-root-resonance.md)
on the original `F16` matrix. It accepts rational heights and scales.
The present adapter supports rational `F16` places and infinity. It
does not silently approximate higher-degree closed places.

The structured basis constructors include ordinary Pascal matrices,
quantum Pascal matrices, and Vandermonde evaluation matrices. Their
inverses are ordinary exact finite-field inverses; no semisimplicity
assumption is used.

The frozen precontrols cover three pairs of physical bases for each of
the three-letter and four-letter seeds. Across both cuts and the local
places zero, two, and infinity, 36 complete native profiles agree with
independent exact rational profiles. Ninety physical-height or
local-scale gradient coordinates agree with centered differences, and
1,011 endpoint entries are checked. These controls validate the oracle
used in the ongoing bounded search; they do not report that search as
completed.

The sources and exact-control evidence are:

- [`fused_physical_basis_oracle.py`](../tools/fused_physical_basis_oracle.py)
- [`check_fused_physical_bases.py`](../tools/check_fused_physical_bases.py)
- [`fused-physical-bases.json`](../evidence/fused-physical-bases.json)

```sh
python3 tools/check_fused_physical_bases.py \
  --output evidence/fused-physical-bases.json
```

The continued numerical search must retain one common physical basis
and one common pair of physical types for all reconstruction pieces.
The separate
[adaptive gauge-catalogue theorem](weighted-adaptive-gauge-catalogues.md)
allows the local rank-bound gauge to change between jets, while retaining
that common tensor and type restriction.
