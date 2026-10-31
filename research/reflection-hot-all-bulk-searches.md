# Completed searches across every bulk-order cone

## Result and scope

The [all-order bulk construction](reflection-hot-all-bulk-gauges.md)
admits every ordering of nonnegative physical-leg bulk heights. Three
completed finite searches used this enlarged family, including cold
spectators and both physical paired-type costs. None found a pressure
below the Gaussian baseline.

The searches enumerate all 16 pairs of cut branches for each tested
coefficient row. The continuous heights and scale are optimized only
from the recorded starts and within stated boxes. Thus the results
do not prove a global obstruction for arbitrary heights or defects.
The coefficient sources and Smith profiles are exact; the pressure
optimization is numerical.

## 1. A compensation parameterization in every cone

At each of the two boundaries of a cut, choose which of its two
physical legs has the larger bulk. Four independent choices across
both cuts give \(2^4=16\) cones. If \(l\ge0\) is the smaller
bulk and \(e>0\) the ordered excess, the two leg bulks are
\(l\) and \(l+e\). For each leg \(a\), introduce a nonnegative
hot-state slack \(s_a\) and use relative hot height

\[
h_a(1)=s_a-b_a-l,\qquad h_a(0)=0.
\tag{1}
\]

Every removed hot-cold internal Gaussian factor then satisfies the
proved regularity inequality. Its added order is at least
\(b_a+l+h_a(1)=s_a\). Internal factors involving two selected hot
variables remain in the exact active tensor.

Each numerical optimization has 17 variables: the logarithm of the
scale, four smaller bulks, four excesses, and eight slacks. The two
local cross pivots and the simultaneous physical-leg swaps select
the applicable proved regular basis. A direct finite-difference
check covers all 17 derivatives in all 16 cones: 272 comparisons,
with maximum discrepancy below \(2.4\cdot10^{-10}\).

## 2. Fixed-height boxes at M=64 and M=128

The [first driver](../tools/search_reflection_hot_all_bulk_gauges.py)
tests the Gaussian, the fixed hot-pair canceler, the seven exceptional
natural-profile rows, and four distinct random F4 coefficient rows
with seed 133749. The same 13 tensors are used at both interval
lengths. Each permits every even exterior coefficient on the four
selected physical variables.

The smaller bulks lie in \([0,8]\), excesses in
\([10^{-6},16]\), and slacks in \([0,24]\). The positive excess
floor keeps each optimizer inside its designated ordering cone.
The exact natural-height point is evaluated separately. The scale
\(\ell=\log s\) lies between one quarter and three times the
optimized natural Gaussian scale. Each cone has two starts, at
\((l,e)=(0,10^{-6})\) and \((1,1)\), with zero slacks.
The nominal optimizer evaluation budget is 70 per start.

| Interval | Coefficient rows | Cones | Starts | Pressure queries | Best pressure | Best non-Gaussian pressure |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 13 | 208 | 416 | 25,844 | 0.7075072113149573 | 0.7476474310261014 |
| 128 | 13 | 208 | 416 | 25,917 | 0.7072351451779789 | 0.7314765705605997 |

The query counts include one separate natural-height query per
candidate. Every planned cone and start completed. An optimizer
termination can report its budget or line-search limit; the records
retain every such message and actual evaluation count. Completion
of the finite schedule does not imply convergence to a global
minimum within each box.

The complete records are
[`reflection-hot-all-bulk-search-M64.json`](../evidence/reflection-hot-all-bulk-search-M64.json)
and
[`reflection-hot-all-bulk-search-M128.json`](../evidence/reflection-hot-all-bulk-search-M128.json).
For every cone they retain the best physical heights actually
evaluated, the selected local branches, active and spectator orders,
physical energies, optimizer results, complete-source metadata, and
executed source hashes.

## 3. Heights proportional to the interval

The [second driver](../tools/search_reflection_hot_scaled_bulk_gauges.py)
tests the Gaussian and the first two exceptional rows at \(M=128\).
It uses smaller bulks in \([0,M/2]\), excesses in
\([10^{-6},2M]\), and slacks in \([0,2M]\). Three starts per cone
use \((l,e)=(0,10^{-6}),(M/16,M/8),(M/4,M/2)\).
The scale window is unchanged, and the nominal evaluation budget is
100 per start.

All **three candidates, 48 cones, and 144 starts** completed, using
**11,954 pressure queries**. The Gaussian remains best at
0.7072351451779789. The best non-Gaussian value is
0.7325370173031095. This explores heights on the scale of the hot
physical energies while including their full cold-bit cost.

The canonical
[`complete record`](../evidence/reflection-hot-scaled-bulk-search-M128.json)
comes from a fresh run after an initial numerical proposal was
rejected. In that first attempt, the optimizer returned a slack of
\(-2.117582368135751\cdot10^{-22}\), outside its stated nonnegative
box. The unchanged strict regularity gate rejected the resulting
negative margin. That point was not scored. The original temporary
prefix output was not retained.

The corrected search projects optimizer proposals to the stated box
before evaluating them. The saved best vectors and physical heights
are those projected values actually evaluated. In the complete run,
seven proposals required projection, each by at most
\(2.117582368135751\cdot10^{-22}\). The offending slack becomes
exactly zero, satisfying the original regularity condition. No
evaluator inequality or theorem tolerance was changed. The
[`diagnostic`](../evidence/reflection-hot-scaled-bulk-rejected-query.json)
retains the rejected vector and its exact floating regularity margin.

## 4. Replay and limits

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_reflection_hot_all_bulk_gauges.py --interval 64 --random 4 --maxfun 70 --output /tmp/reflection-hot-all-bulk-M64-replay.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_reflection_hot_all_bulk_gauges.py --interval 128 --random 4 --maxfun 70 --output /tmp/reflection-hot-all-bulk-M128-replay.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_reflection_hot_scaled_bulk_gauges.py --interval 128 --maxfun 100 --output /tmp/reflection-hot-scaled-bulk-M128-replay.json
```

The frozen structural checker already replays literal original
weighted cuts in every branch pair. These new files record searches
using that evaluator. Their negative findings are limited to the
specified tensors, starts, scale windows, and height boxes. No
improved factoring exponent or all-height no-gain theorem follows.

## 5. Reconstructing the saved best points

The [saved-point checker](../tools/check_reflection_hot_bulk_search_records.py)
reconstructs every candidate best point and every cone best point in
these three records, together with the candidate best points in the
separate-leg records. It rechecks the complete schedule counts and
all 272 chain-rule comparisons. All **543 saved point entries**
reproduce their pressures, active and spectator spectra, and physical
energies. Some candidate best points also appear as cone best points.
The [replay record](../evidence/reflection-hot-bulk-saved-point-replay.json)
has zero maximum pressure discrepancy on this platform and retains
the hashes of all five input records. This replay does not rerun the
optimizers or certify their global optimality.

```bash
OPENBLAS_NUM_THREADS=1 make verify-bulk-search-records
```
