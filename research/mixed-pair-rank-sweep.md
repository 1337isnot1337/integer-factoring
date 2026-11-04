# Mixed-pair rank across a 40-bit cycle

**Status.** Finite diagnostic on one fixed 40-bit input. No uniform pair
supply, rank-window, or complete-bound theorem is proved. The complete saving
exponent remains **1/2**.

This extends the two-target certificate in
[mixed-pair rows on a 40-bit cycle](mixed-pair-gap-cycle.md), which in turn
continues [target orders and mixed rows](target-order-mixed-rows.md). The
construction and the cycle are due to that work; the only new object here is a
sweep over all outside targets in the same cycle, recording where the pair-row
rank actually succeeds and where it fails.

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

The committed receipt uses `max_target=2500`. There are `171` outside targets
in the cycle and none are skipped.

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

A one-off run with `--max-target 10000` (not part of the committed receipt)
found `608` outside targets, `7` rank failures, and the same two targets with
index above `D^{1/4}` in that range. The companion 40-bit note scans gap
targets exhaustively through `100000` and finds ten; this sweep is not an
exhaustive gap-target census. Every additional failure again had index `1`, `2`, or
`10`. The pattern is stable under the larger sweep, but it remains a finite
observation on this cycle.

## What the sweep says

- Pair rows alone are not universal at window `3` on this cycle: `4` of `171`
  targets have no pair hit and rank `22` rather than `23`. A uniform pair-only
  first-hit/rank theorem at this window cannot be asserted from the two-target
  certificate. The failure is a supply failure for the target pair, not a
  connected-graph rank defect.
- Every rank failure here has index at most `10`, far inside the small-index
  density range. The only two targets outside that range in this sweep,
  with index `1025` and `2050`, are full rank. This is a finite hint of **complementarity**: density
  or order rows cover the small-index side, while pair rows cover the
  medium/large-index gap side. It is not a proof of that division, and it does
  not supply a uniform rank lower bound.
- The selected minor multiples for full-rank targets vary widely, including
  `4096`, `1000`, `2032`, `-8192`, and many other values. The determinant is
  always a multiple of the cycle order `h`, as required by the maximal-minor
  annihilator interface, but no uniform small cofactor is visible. The
  `4096h`/`1000h` factors in the two-target certificate are not representative
  of a uniform constant.
- The rank failure at an index-`1` target shows that index alone does not
  guarantee a pair hit. The failure set is a property of the pair-row supply,
  not just the target order/index pair. In this sweep, however, every target
  without a pair hit is density-covered.

The open obligation is unchanged: prove a uniform supply of admissible pair
rows and a rank or lattice-index lower bound, or prove that the rank-deficient
targets are always covered by the density/order side. This sweep supplies
neither; it records where the finite construction stops.

## Reproduction

```sh
python3 tools/check_mixed_pair_rank_sweep.py \
  --max-target 2500 --max-window 3 \
  --output /tmp/integer-factoring-mixed-pair-rank-sweep.json
make verify-mixed-pair-rank-sweep
```

The receipt is `evidence/mixed-pair-rank-sweep.json`.
