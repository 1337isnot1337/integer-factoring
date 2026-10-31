# Complete finite screens of symmetric polynomial tiles

## Result and scope

Four exact polynomial catalogues contain 3,830 distinct normalized
profiles. They were used in 32 searches over scale and aggregate charge,
with a separate filter choice for each tile's actual physical energies.
Every reported best point selected the Gaussian baseline in all tiles.
No strict crossing candidate was found.

The searches optimize the positive root of the complete paired-type
pressure. The [independent replay](../evidence/symmetric-filter-catalogue-check.json)
reconstructs every polynomial witness and re-evaluates every reported
root, physical entropy, and full critical pressure. It does not certify
global optimality over continuous parameters or exclude profiles outside
the finite catalogues.

## 1. Exact catalogues

The [binary](symmetric-multilevel-filter-synthesis.md) and
[prime-state](symmetric-prime-multilevel-filters.md) reductions permit
these finite coefficient searches:

| Characteristic | Modes per tile | Source ramifications | Exact degree vectors | Exact feasibility queries | Distinct normalized profiles |
|---:|---:|---|---:|---:|---:|
| 2 | 4 | 1, 3, 5 | 991 | 68,009 | 517 |
| 2 | 8 | 1, 3 | 10,314 | 1,612,776 | 2,802 |
| 3 | 2 | 1, 3 | 436 | 25,108 | 294 |
| 3 | 4 | 1 | 304 | 42,862 | 217 |

The binary catalogues include every occupation-weight support mask at
the stated degree levels, fractional degree offsets realized by odd
ramification, and a fixed small family of multilevel bounds. The qutrit
two-mode catalogue includes every support mask. Its four-mode catalogue
uses the complete support-mask list stored in its manifest: singleton,
prefix, suffix, digit-binomial, and deterministically selected masks.
It is not an enumeration of every four-mode support mask.

For every degree vector, three stated row priorities maximize beta-order
targets by exact feasibility queries. This produces a finite Pareto
sample, not all possible order profiles. Each retained profile stores
its source degree bounds, source order bounds, endpoint colors, source
ramification, and a hash of its literal polynomial lift. The normalized
profile is `(B_w/m,D_w/m)`, compared with exact rational arithmetic.

The complete records are:

- [binary four-mode tiles](../evidence/symmetric-tiled-filters-r4.json);
- [binary eight-mode tiles](../evidence/symmetric-tiled-filters-r8.json);
- [qutrit two-mode tiles](../evidence/symmetric-qutrit-tiled-filters-r2.json);
- [qutrit four-mode tiles](../evidence/symmetric-qutrit-tiled-filters-r4.json).

## 2. Actual physical mode partitions

Each catalogue is tested on four explicit disjoint partitions:

1. Consecutive physical energies within each tile.
2. Interleaved energies spanning the interval.
3. Reflected low/high energy pairs.
4. A deterministic shuffled partition.

The full list of actual positive odd energies is stored for every tile.
The generating polynomial for that list supplies its occupation-weight
probabilities. No equal-energy or continuum approximation is used.

Binary intervals `n=64,128` are tiled exactly. Qutrit intervals `n=27,81`
may leave one or three modes outside the chosen two- or four-mode blocks.
Every such mode is explicitly listed as an untouched spectator and
remains in the complete Gaussian and physical partitions.

Profiles from different source ramifications share a common odd
ramification. Their literal polynomial exponents are scaled by an
integer, as proved in the synthesis notes. The largest common
ramification in these records is 15.

## 3. Genuine positive roots

The floating objective solves

\[
G_{n,p}(J,t;\rho)+\sum_{\rm blocks}\log\mathcal R_{\rm block}(J,t)=0
\]

for its positive `rho`, then searches over `tau=nt` and `J/n`. The
baseline is always available with relative log cost exactly zero.
Profiles are selected independently in each tile because their complete
relative costs multiply.

The following values are the best recorded roots, repeated for every
tested block size and partition at the corresponding characteristic and
interval length:

| Characteristic | Interval length | Best recorded `rho` | Full pressure at `rho=1/sqrt(2)` | Physical entropy at the root |
|---:|---:|---:|---:|---:|
| 2 | 64 | 0.7075072433 | 0.0195069230 | 24.3623484 |
| 2 | 128 | 0.7072351467 | 0.0109546876 | 42.6737442 |
| 3 | 27 | 0.7091192809 | 0.0736842022 | 18.3317317 |
| 3 | 81 | 0.7074661596 | 0.0302753013 | 42.1323774 |

Every returned best point uses zero nonbaseline tiles. These roots
therefore reproduce finite Gaussian values and do not improve the
critical threshold. The positive physical entropies show that the
reported near-threshold values do not come from an empty-entropy limit.
The numerical values are search results rather than exact global minima.

## 4. Replay and reproducibility

```sh
python3 tools/search_symmetric_tiled_filters.py --block-size 4 --ramifications 1,3,5 --lengths 64,128 --output /tmp/symmetric-tiles-r4.json
python3 tools/search_symmetric_tiled_filters.py --block-size 8 --ramifications 1,3 --lengths 64,128 --iterations 35 --output /tmp/symmetric-tiles-r8.json
python3 tools/search_symmetric_prime_tiled_filters.py --block-size 2 --ramifications 1,3 --lengths 27,81 --output /tmp/qutrit-tiles-r2.json
python3 tools/search_symmetric_prime_tiled_filters.py --block-size 4 --ramifications 1 --lengths 27,81 --output /tmp/qutrit-tiles-r4.json
python3 tools/check_symmetric_filter_catalogues.py --output /tmp/symmetric-catalogue-replay.json
```

The [catalogue checker](../tools/check_symmetric_filter_catalogues.py)
reconstructed all 3,830 saved polynomial lifts, including their exact
orders, full degrees, endpoint colors, and polynomial hashes. It checked
all 32 complete root evaluations and the actual energy partitions.
The largest reported root residual is below `3e-10`; the replay allows
`1e-7` to accommodate floating solver and platform differences.

An earlier exploratory whole-symmetric 32-mode search was interrupted
after 210 of 290 planned degree vectors. It is not part of these
catalogues, counts, or conclusions. The exact natural-charge relative
gain in the separate four-mode theorem remains valid; the finite tiled
searches show that this gain did not improve the best recorded full
root in the tested families.
