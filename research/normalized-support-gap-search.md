# Search the support gap with its angle scale normalized

**Registered heuristic, 2026-11-01.** The actual supporting value `R`
still has to be strictly below one, with an independent exact certificate,
before any improved recursion or factoring score is claimed. The complete
score remains **1/2**. This is a different numerical objective on the same
[twelve-letter literal graph family](q12-reflection-graph-search.md).
The previous registered raw-support walks continue unchanged at their
frozen source `10aed41`.

Let `a=cos(theta)`, `b=sin(theta)` and `y=log(tan(theta))`. Minimize

```
F = (R-1)/(a*b) = 2*cosh(y)*(R-1).
```

Since the multiplier is positive, `F<0` is exactly `R<1`. The derivative
is the existing support derivative multiplied by `2*cosh(y)`, plus
`2*sinh(y)*(R-1)` in the log-tangent coordinate. This is an optimizer
preconditioner, not a new mathematical bound or convexification.

Raw-support fits have repeatedly approached the edge of the declared
angle interval while remaining above one. The
[fixed-gauge axis theorem](polynomial-crossing-axis-obstruction.md)
does not prohibit moving gauges, nor does it prove this heuristic useful.
The normalization aims to compare positive gaps across angles without
giving a small positive gap at a small weight an automatic advantage.
Its performance and the resulting tensors must be assessed separately.

[Independent controls](../tools/check_normalized_support_gap.py) compare
54 real finite derivatives on an actual nonbinary twelve-letter tensor,
at three angles and across all eight physical-height groups, scale and
log tangent. The largest error was below `9e-10`. The actual rate and
the normalized objective have matching threshold signs in the declared
controls. The [receipt](../evidence/normalized-support-gap-controls.json)
is a numerical implementation check, not a rigorous pressure certificate.

## Registered run

The [manifest](../evidence/q12-normalized-gap-manifest.json) freezes its
source hashes and OS-random seeds before execution. All source hashes
are checked at startup. The same literal source, complete 144-channel
cuts, 88 free heights, physical costs and native rejection rules are used.

Three initial fixed-profile fits and free-height fits use the new
objective. From their common best point run two 1,800-second tracks:

- A single-coordinate graph walk with the supporting angle fixed at
  `pi/4`. Its fixed-angle starting point is actually refitted before
  acceptance comparisons. Its rate is `R=sqrt(2)*rho`.
- A multi-coordinate graph walk changing 2, 4 or 8 distinct field
  coefficients, with the free angle and normalized-gap objective.

The second track starts from the initial frozen point, not the first
track's discoveries. All height, scale, angle and per-trial bounds match
the predecessor manifest: heights `[-128,128]`, log scale
`[-9,log(4)]`, log tangent `[-2,2]`, 100 iterations and ten seconds per
trial, checked between queries. The first track fixes log tangent to
zero. The initial stage has 90 seconds total. Annealing temperatures
`(.001,.01,.05,.2)` now apply to objective differences, with 25-proposal
blocks and a restart from each track's best every 100 proposals.

Every trial preserves both its actual support and optimizer objective,
parameters, full best/last native details and rejected coordinates. The
smallest actual support query is recorded separately from the best
normalized objective. A negative floating objective is only a candidate;
fresh independent exact profiles and outward pressure remain necessary.
These different track modes and random histories do not constitute a
controlled speed or optimizer-effectiveness comparison.

```
make verify-normalized-support-gap
python3 tools/search_q12_reflection_graphs.py --manifest evidence/q12-normalized-gap-manifest.json --output /path/normalized-gap.json
```

## Completed normalized-gap tracks at `67a0666`

The [complete ledger](../evidence/q12-normalized-gap-complete.json.xz)
retains both 1,800-second tracks and all initial/fixed-angle starts.
The balanced track attempted 174 proposals and 11,166 complete queries;
the free-angle track attempted 173 proposals and 6,423 queries. Including
starts, the run made 18,246 successful complete queries with zero
rejections in 3,631.94 seconds.

No optimizer reported convergence: 347 attempts reached the time cap,
and four reached the iteration cap. Their bounded best/last points and
statuses are retained. Neither track improved its initial normalized
objective. Its retained best support is about `1.013947250012`, with
normalized gap about `0.0278945000234`, above zero. The separately
recorded minimum-support query must not be confused with this best-gap
point. These outcomes provide no exact pressure certificate, global
optimum or improved factoring score.

Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI assistant
support. The existing complete Smith, Gibbs/support gradient and literal
reflection/split implementations retain their attribution. No new Lean
work or independent mathematical acceptance is claimed.
