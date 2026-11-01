# Completed four-letter Laurent-series family search

## Completed scope and result

The complete 65,536-word four-letter family

```math
N(w)=w^{-1}G_1+w^{-2}G_2,
\qquad
G_j=\begin{pmatrix}0&B_j\\0&0\end{pmatrix},
\qquad B_j\in\operatorname{Mat}_2(\mathbb F_4),
```

has been constructed and screened. All 32 retained words completed
their declared full-gauge refinement runs. The best recorded parameter
remains the undeformed source value

```text
rho = 0.7310809359734765.
```

The best nonzero word among the refined candidates has
`rho = 0.7546022282548617`. Both values exceed the post-half target
`1/sqrt(2)`. This is a completed finite search checkpoint, not a global
optimization theorem for its gauges or types.

The [full compressed archive](../evidence/reflection4-laurent-series-search.json.gz)
retains every one of the 65,536 screen outcomes, all source hashes,
the complete optimization outcomes, and the winner parameters. The
recorded elapsed time is 1,589.898 seconds. The declared time cap was
4,200 seconds; the finite family and all selected refinements finished
before that cap.

## 1. Exact family and numerical screen

The [triangular jet theorem](triangular-reflection-jet-families.md)
proves that all these frames are polynomial. The two matrices have
common upper-block support, so every product `G_i G_j` vanishes and
`N(w)^2=0`. The [complete series constructor](../tools/laurent_series_coupled_frames.py)
applies `I+N` on the two Bell-normalized primal legs and its transpose
on the two dual legs, then restores the original physical coordinates.
It checks the full endpoint and actual degree after all four actions.

The [search source](../tools/search_reflection4_laurent_series.py)
encodes the eight field coefficients as base-four digits of an integer
from zero through 65,535. The least significant digit comes first in
the order `(pole,row,column)`, with poles one and two, rows zero and
one, and columns two and three. Word zero is the undeformed seed.

Each word was evaluated at two declared natural-gauge scales, using
complete cut profiles and common Gibbs-optimized paired types. The
32 best screen values selected the refinement list. Each refinement
optimized all physical-height coordinates and the common scale, with
the balanced supporting angle and common types. It used both a warm
start and the original natural start, with the stored iteration and
time limits. These are numerical optimization runs; completion does
not imply that a global optimum was found.

The best nonzero refined word is 8,192. Its screen value was about
`0.8271990`, and refinement reduced it to about `0.7546022`. That
improves the sampled gauge for this word, but still loses to the
undeformed source and supplies no stronger factoring bound.

## 2. Independent replay

Run

```sh
python3 tools/check_laurent_series_search.py \
  --output evidence/laurent-series-search-replay.json
```

The [checker](../tools/check_laurent_series_search.py) verifies that all
65,536 indices occur exactly once, validates the archived source hashes,
and checks that the 32 completed refinements match the retained list.
It freshly reconstructs those words and additional fixed controls.

The independent coefficient calculation uses exponent dictionaries and
a separate field multiplier. Since all `G_i G_j` vanish, it factors
each single-leg action into `I+w^-1 G_1` and `I+w^-2 G_2`. It retains
all negative intermediate exponents through the complete four-leg
action. Only the final tensor is required to be polynomial. The original
Bell changes are then restored, and every resulting coefficient is
compared with the archived constructor.

Finally, every refined parameter vector is evaluated in a fresh
complete-profile oracle. The [replay record](../evidence/laurent-series-search-replay.json)
reports exact tensor comparisons separately from floating pressure
agreement. The archive supports the stated finite coverage and recorded
best values; it does not exclude better gauges for unrefined words.
