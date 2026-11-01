# Bounded ternary physical-basis search

## Outcome

The completed search retained the explicit ternary family's numerical value
`rho = 0.8203391769145085`. It did not find a better candidate. The rigorous
result remains the separately certified
[`rho < 821/1000` witness](ternary-basis-mixed-crossing.md), which is above
the threshold needed to improve the factoring exponent.

The search varied endpoint-preserving bases on both physical pairs,
independent paired types, all physical-letter heights at three places,
the row-place allowance, and the common positive pressure scale. It used
the row places `2,infinity` and the column place `0` over `F3`.
The finite search is not a proof of global optimality for these choices.

## Complete retained record

| Stage | Recorded loop time | Bases | Completed optimizations | Best numerical `rho` |
| --- | ---: | ---: | ---: | ---: |
| Initial search | 257.742832498 seconds | 48 | 384 | 0.8209477018467242 |
| Warm search from the exact family | 1540.254401347 seconds | 268 | 2144 | 0.8203391769145085 |
| Total | 1797.997233845 seconds | 316 | 2528 | 0.8203391769145085 |

The initial search was interrupted after a completed checkpoint to add the
explicit family's warm start. Its exact interruption time was not recorded;
the table gives the last saved elapsed time. The warm stage completed its
full 1540-second budget and finished the current optimization. Every
recorded basis has eight completed optimization records.

The [compact archive](../evidence/ternary-physical-basis-search.json)
preserves every retained trial, every basis, the best parameter vector at
each basis, iteration and termination flags, and both original progress-file
hashes. It also records snapshot modification times. Those timestamps are
file-write times, not reconstructed process-start times. The initial
program predates the optional warm branch; the stored source hashes
identify the durable reproducer rather than an unavailable historical
byte-for-byte snapshot of that earlier program.

## Mathematical and numerical checks

Each basis transformation uses `A,B,(B^-1)^T,(A^-1)^T` on `x,y,u,v`,
respectively, and checks all 81 endpoint entries. Complete polynomial
minor valuations are computed over `F3` by exact Laplace recursion.
Weighted local Smith orders then use successive minimum-valuation Schur
pivots. The floating weights choose a chamber; they do not certify a
real-valued optimum.

During the two stages, independent block-Toeplitz calculations checked
1008 weighted profiles. Finite differences checked 336 objective
gradients. These recorded controls cover negative heights, both cost
cuts, the two finite places, and homogeneous infinity orders.

The [replay checker](../tools/check_ternary_basis_search.py) verifies the
entire trial accounting and every basis inverse. It independently
reconstructs each stage's best transformed seed, recomputes its complete
minor valuations, and reproduces both retained pressure values. The
[replay evidence](../evidence/ternary-physical-basis-replay.json) records
the result and the archive hash. This replay does not rerun the entire
timed search or claim that every numerical optimizer reached a minimum.

## Reproduction

The durable implementations are
[`search_ternary_physical_bases.py`](../tools/search_ternary_physical_bases.py)
and [`ternary_mixed_pressure.py`](../tools/ternary_mixed_pressure.py).
The random seeds were `21831` and `29102`, respectively. The first command
below supplies a comparable bounded initial run; its stopping time will
not reproduce the deliberate interruption exactly.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/search_ternary_physical_bases.py \
  --seconds 258 --restarts 8 --seed 21831 \
  --output /tmp/ternary-basis-initial.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/search_ternary_physical_bases.py \
  --seconds 1540 --restarts 8 --warm --seed 29102 \
  --output /tmp/ternary-basis-warm.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/check_ternary_basis_search.py \
  --output /tmp/ternary-basis-replay.json
```

The search uses NumPy and SciPy. The proof and exact certificate for the
`821/1000` witness remain independent of these numerical optimizers.
