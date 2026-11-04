# Mixed-pair rank across a 40-bit cycle

**Status.** Finite diagnostic on one fixed 40-bit input. No uniform pair
supply, rank-window, or complete-bound theorem is proved. The complete saving
exponent remains **1/2**.

This extends the two-target certificate in
[mixed-pair rows on a 40-bit cycle](mixed-pair-gap-cycle.md), which in turn
continues [target orders and mixed rows](target-order-mixed-rows.md). The
construction and the cycle are due to that work; the only new object here is a
sweep over all outside targets in the same cycle, recording where the pair-row
rank succeeds and fails at window `3`, and a larger window-`4` check that
closes the failures.

## Construction and sweep

Use the same input and cutoff as the two-target certificate:

\[
 N=1000036000099=1000003\cdot 1000033,\qquad D=-N,\qquad Y=177.
\]

The reduced forms of the prime form `5` form a complete cycle of order
`h=198850`, and the `22` split primes through `177` lie in it. Take every
split prime `p` with `177<p<=2500` whose reduced form lies in the same cycle,
and form the `23`-column target set

\[
 \{5,7,19,\ldots,173\}\cup\{p\}.
\]

For each unordered pair of targets and each exponent `k=1,2,3`, test whether
the reduced leading coefficient of `(qr)^k` factors over the split primes
through `177`. When it does, the exact decoder gives the signed word `v` and
the relation row

\[
 k(e_q+e_r)-v.
\]

Pool the rows by window, compute the rank profile, and at the first full-rank
window select a maximal minor and run the two-part extractor. Every selected
row is replayed; every selected minor is checked to be a multiple of `h`.

## Results

The window-`3` receipt uses `max_target=2500`. There are `171` outside
targets in the cycle and none are skipped.

- Full rank first at window `2`: `38` targets.
- Full rank first at window `3`: `129` targets.
- No full rank at window `3`: `4` targets.
- The two-part extractor returns `1000003` for all `167` full-rank targets.

The only two targets in this sweep with index above
`D^{1/4}=1000.008...` reproduce the two-target certificate:

| target | order | index | first full window | selected minor | extractor |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1117 | 194 | 1025 | 2 | `4096h` | 1000003 |
| 1289 | 97 | 2050 | 3 | `1000h` | 1000003 |

The four rank failures are:

| target | order | index | target pair hits | rank profile `(window: rank/rows)` |
| ---: | ---: | ---: | ---: | :--- |
| 353 | 99425 | 2 | 0 | `1:0/0`, `2:20/20`, `3:22/44` |
| 937 | 198850 | 1 | 0 | `1:0/0`, `2:20/20`, `3:22/44` |
| 1511 | 99425 | 2 | 0 | `1:0/0`, `2:20/20`, `3:22/44` |
| 2339 | 19885 | 10 | 0 | `1:0/0`, `2:20/20`, `3:22/44` |

The four failures are **isolated target columns**: there is no pair hit
involving the target at any exponent `1..3`, so the target column is
identically zero and the base rows alone have rank `22`. They are not
connected-graph rank defects or determinant cancellations. The two full-rank
gap targets each have five target pair hits.

The window-`3` pattern is not permanent. A one-off window-`6` run on the same
`2500` range closes all four failures at window `4`. The committed window-`4`
receipt extends the range to `10000`: `608` outside targets, no skips, full
rank first at window `2` for `202`, at window `3` for `399`, and at window `4`
for `7`; there are no failures, and the two-part extractor returns `1000003`
for all `608`. The seven window-`4` targets are the four window-`3` failures
plus `3583`, `5857`, and `8527`. The isolated targets acquire pair hits at
window `4`.

The companion 40-bit note scans gap targets exhaustively through `100000` and
finds ten; the window-`4` sweep through `10000` sees only `1117` and `1289`,
each full at window `2`/`3` with seven target pair hits. This sweep is not an
exhaustive gap-target census.

The companion note now also proves an **isolated-column repair** lemma: if a
column is zero in every pair row, adjoining the diagonal row `M e_i` with a
common annihilator `M` raises the rank by one and multiplies the determinant
by `M`. A control on all split primes through `500` fills the zero column
`229` this way and extracts `1000003`. Thus the window-`3` isolated columns
have two finite repairs: a window-`4` pair hit (no annihilator factor) or a
common-annihilator row (cost about `sqrt(h)` group operations and a factor `M`
per isolated column). The second is proved but not free; the first is finite
here and not yet uniform.

## What the sweep says

- At the two-target certificate's window `3`, pair rows alone are not
  universal on this cycle: `4` of `171` targets have no pair hit and rank `22`
  rather than `23`. A uniform pair-only first-hit/rank theorem at window `3`
  cannot be asserted. The failure is a supply failure for the target pair, not
  a connected-graph rank defect, and all four close at window `4`.
- Every window-`3` failure has index at most `10`, far inside the small-index
  density range. The only two targets outside that range in the sweep, with
  index `1025` and `2050`, are full rank. This is a finite hint of
  **complementarity** at window `3`: density or order rows cover the
  small-index failures, while pair rows cover the medium/large-index gap side.
  The window-`4` receipt closes the same failures with pair rows alone, so the
  finite data does not force the hybrid. Neither observation is a uniform rank
  lower bound.
- The selected minor multiples for full-rank targets vary widely, including
  `4096`, `1000`, `2032`, `-8192`, and many other values. The determinant is
  always a multiple of the cycle order `h`, as required by the maximal-minor
  annihilator interface, but no uniform small cofactor is visible. The
  `4096h`/`1000h` factors in the two-target certificate are not representative
  of a uniform constant.
- The window-`3` failure at an index-`1` target shows that index alone does
  not guarantee a pair hit. The failure set is a property of the pair-row
  supply, not just the target order/index pair. At window `4` the same target
  acquires a hit.

The open obligation is unchanged in kind: prove a uniform supply of
admissible pair rows and a rank or lattice-index lower bound. The finite
window-`4` receipt makes a uniform window-`4` pair-hit/rank theorem the
sharper target; a lower-cost density or annihilator row for isolated targets
is the alternative. The common-annihilator repair is proved, but its
`sqrt(h)` cost and determinant factor are exactly what a uniform supply would
need to avoid. This sweep supplies neither; it records where the finite
construction stops.

## Reproduction

```sh
python3 tools/check_mixed_pair_rank_sweep.py \
  --max-target 2500 --max-window 3 \
  --output /tmp/integer-factoring-mixed-pair-rank-sweep.json
make verify-mixed-pair-rank-sweep

python3 tools/check_mixed_pair_rank_sweep.py \
  --max-target 10000 --max-window 4 \
  --output /tmp/integer-factoring-mixed-pair-rank-sweep-w4.json
make verify-mixed-pair-rank-sweep-w4
```

The receipts are `evidence/mixed-pair-rank-sweep.json` (window `3`) and
`evidence/mixed-pair-rank-sweep-w4.json` (window `4`).
