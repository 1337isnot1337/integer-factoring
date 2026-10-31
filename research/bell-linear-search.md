# Exact linear-pencil calibration and larger finite searches

The [Bell construction](bell-diagonal-crossings.md) gives an exact crossing
for every ternary vector `v` by setting `s(z)=1+(z-1)v`. Its row eigenvalues
have only three forms over `F_3`: `1`, `z`, and `-(1+z)`. The column forms
are specified by the same three-value histogram of `Tv`.

This note records finite searches and no improved factoring exponent.
The objective is the largest rank product at tensor length **12**.
The retained witness is separately evaluated at lengths
`1,2,4,8,12,16,24,32`. These normalized finite ranks omit the logarithmic
slack in the asymptotic interface and are not leaderboard scores.

## Exact coefficient-rank formula

If `(a,b,c)` counts the three row forms, then the coefficient of degree
`j` has exact rank

```math
r_j=\sum_{n_0+n_1+n_2=t}
 \binom{t}{n_0,n_1,n_2}a^{n_0}b^{n_1}c^{n_2}
 \left[0\le j-n_1\le n_2\right]
 \left[\binom{n_2}{j-n_1}\not\equiv0\pmod3\right].
```

The formula follows by multiplying scalar eigenpolynomials and counting
their physical eigenwords. The sign from the third form is nonzero and
does not affect the rank. Crucially, only scalar polynomial coefficients
are reduced modulo three; the integer multinomial and alphabet
multiplicities remain dimensions. The column rank uses the histogram
of `Tv`. The literal-word controls in `make verify-bell` independently
check the formula for all 81 one-bit pencils through length four.

## Completed schedules

| Local dimension | Scope | Histogram pairs encountered | Best finite rho at t=12 | Same witness at t=32 |
| --- | --- | ---: | ---: | ---: |
| 2 | All `3^4=81` ternary vectors | 12 | 0.824557733 | 0.858063233 |
| 4 | All `2^16=65,536` binary vectors | 60 | 0.877280823 | 0.895793508 |
| 8 | Structured starts and 65,536 seeded local proposals over `F_3^64` | 4,441 | 0.910359791 | 0.919762089 |
| 16 | Structured starts and 32,768 seeded local proposals over `F_3^256` | 6,772 | 0.931327752 | 0.935942685 |

Only the first two schedules are exhaustive in their displayed coefficient
domains. The latter two are heuristic searches. The evaluation counters
include restart evaluations and duplicates, whereas the exhaustive input
counts above count the actual enumerated domain. No piece-count minimality
or larger-domain optimality is inferred.

There is a useful negative control: the dimension-two retained witness
has finite normalized rho `0.701838731` at length two, below `1/sqrt(2)`.
The complete-span theorem nevertheless excludes that asymptotic rate.
Its ranks at longer lengths visibly increase. A promising isolated rank
profile cannot stand in for a uniform family estimate.

## Saved evidence and replay

- [Dimension two](../evidence/bell-linear-d2.json).
- [Dimension four](../evidence/bell-linear-d4.json).
- [Dimension eight](../evidence/bell-linear-d8.json).
- [Dimension sixteen](../evidence/bell-linear-d16.json).

Each record retains the complete vector, its transformed vector, the
histograms, every exact integer rank, seed, schedule, and evaluation count.

```sh
python3 tools/check_bell_pencils.py
python3 tools/search_bell_pencils.py --bits 3 --steps 65536 --output /tmp/bell-d8.json
python3 tools/search_bell_pencils.py --bits 4 --steps 32768 --output /tmp/bell-d16.json
```

The checker reconstructs all 32 saved two-cut profiles and repeats both
exhaustive calibrations. The larger searches can be regenerated with the
displayed deterministic seeds; elapsed time is not an invariant field.

The next constructive step uses higher-degree polynomial vectors and an
all-length weighted rank bound, rather than continuing to extrapolate
finite linear-pencil profiles. Work used substantial OpenAI Codex assistance.
