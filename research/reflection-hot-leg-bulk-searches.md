# Completed separate-leg bulk searches

## Result and scope

The [separate-leg bulk theorem](reflection-hot-leg-bulk-gauges.md)
permits independent nonnegative cold-bit heights on the four physical
legs at each cut, subject to its stated ordering cone. Two finite
searches applied these gauges to 25 exterior polynomial defects on
the first logical mode of the reflection-compressed Gaussian.
Neither search found a critical pressure below its Gaussian baseline.

These records use the ordering cone proved in that theorem. They do
not search every ordering admitted by the later
[all-order bulk construction](reflection-hot-all-bulk-gauges.md).
The Smith profiles and their complete polynomial sources are exact;
the height and scale optimizations are numerical exploration.

## 1. The finite family and optimization box

The driver is
[`search_reflection_hot_leg_bulk_gauges.py`](../tools/search_reflection_hot_leg_bulk_gauges.py).
At each interval length it tests the same 25 coefficient rows:
the unperturbed Gaussian, the fixed hot-pair canceler, the seven
exceptional rows from the complete natural-gauge screen, and 16
distinct pseudorandom rows over \(\mathbf F_4\), with seed 133748.
Each row has all eight even exterior coefficients available on the
four selected physical variables. The recorded row and resulting
polynomial factor specify the tensor exactly.

Each candidate receives four optimizer starts. The natural start has
all extra heights zero. The other starts set the common boundary
bulk and its ordered excess to \((1,0),(0,1),(1,1)\), respectively.
The nominal function-evaluation budget is 100 per start; the record
retains each optimizer's actual evaluation count and termination
message.

At each boundary, the smaller leg bulk lies in \([0,8]\), and the
larger leg has an additional height in \([0,16]\). Eight nonnegative
hot-state slacks lie in \([0,24]\). Their conversion to relative hot
heights is the exact compensation parameterization in the theorem.
The scale \(\ell=\log s\) lies between one quarter and three times
the optimized natural Gaussian scale for that interval. Both
spectator partitions and both physical paired-type partitions are
included in every query.

## 2. Complete results

| Interval \(M\) | Candidates | Starts | Actual pressure queries | Best pressure | Best non-Gaussian pressure |
|---:|---:|---:|---:|---:|---:|
| 64 | 25 | 100 | 8,009 | 0.7075072113149573 | 0.7455867247441637 |
| 128 | 25 | 100 | 8,072 | 0.7072351451779789 | 0.7313702628621523 |

In both rows, the best candidate is the unperturbed Gaussian at its
natural heights. Both values exceed the target \(1/\sqrt2\).
The saved records are
[`reflection-hot-leg-bulk-search-M64.json`](../evidence/reflection-hot-leg-bulk-search-M64.json)
and
[`reflection-hot-leg-bulk-search-M128.json`](../evidence/reflection-hot-leg-bulk-search-M128.json).
They contain every candidate, each optimizer termination, complete
active orders, spectator orders, physical energies, and source-width
metadata for the best point reached on that candidate.

The \(M=64\) calculation completed before the shared execution
environment was interrupted. The \(M=128\) record had completed
nine candidates at that interruption. After recovery, the saved
prefix was checked against the deterministic coefficient sequence,
and the remaining 16 candidates were completed with the unchanged
candidate routine and the same four starts. The `resume` field
records that boundary explicitly. Its total execution time excludes
the period when the environment was unavailable.

The source hashes in each record identify the version that executed
that record. In particular, the \(M=64\) record predates the later
backward-compatible branch hooks and the explicit interval-scaled
source-precision request. Its original hashes are retained. The
published driver uses a complete source request of
\(\max(512,8M)\); this is 512 at \(M=64\), matching that earlier
calculation, and 1,024 at \(M=128\). Independent extended original-cut
checks of the shared evaluator are recorded in the theorem package.

These finite boxes and local optimizer starts do not establish a
global obstruction for separate-leg heights, other coefficient rows,
or larger intervals. They also do not constitute an all-order bulk
search.

## 3. Replay

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_reflection_hot_leg_bulk_gauges.py --interval 64 --starts 4 --maxfun 100 --random 16 --seed 133748 --output /tmp/reflection-hot-leg-bulk-search-M64-replay.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_reflection_hot_leg_bulk_gauges.py --interval 128 --starts 4 --maxfun 100 --random 16 --seed 133748 --output /tmp/reflection-hot-leg-bulk-search-M128-replay.json
```

The completed evidence is retained separately from the structural
verification records. A rerun may differ in floating-point optimizer
termination while using the same exact candidate tensors and
complete Smith-profile evaluator.

A separate [saved-point replay](reflection-hot-all-bulk-searches.md)
reconstructs all 50 candidate best points in these two records, along
with the later all-order records. It reproduces their complete spectra
and pressure costs using the published evaluator.
