# A single boundary e-gate spans the thirteen-label crossing

## Exact finite identity

The thirteen-label crossing does not require the disjoint `f_3`
placement family. Insert one pure `e` gate at the **central spatial
site** `site=n-1`, in the **first time layer** of the graded seam,
and take all boundary-degree slices. Over `F16` with `q=8`, their
scalar sum is exactly the crossing.

For a half-length `n`, let `T_b` be the degree-`b` coefficient tensor
of that single insertion. The checker finds

```math
\sum_b c_n T_b=\mathbf1_{P=S}\mathbf1_{Q=R},
```

with

| `n` | labels | nonzero slices | support | coefficient `c_n` |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 1 | 2 | 1 | 13 |
| 4 | 2 | 4 | 4 | 7 |
| 6 | 5 | 6 | 6 | 3 |
| 8 | 13 | 8 | 8 | 13 |
| 10 | 34 | 10 | 10 | 11 |

Each row is an exact coefficientwise reconstruction, not a rank or
cost estimate. At `n=8` the identity uses **eight pieces**, compared
with the `4086` nonzero pieces in the full disjoint-`f_3` census. It
uses no `f_3` placement at all.

The `n=10` row is a heavy spot-check (`tools/check_boundary_e_gate_length10.py`).
It spans in dimension `10` with coefficient `11` and has `b=0` cut
ranks `(1156,1)`, so the full-first-cut profile persists at the next
length. The full per-slice cut-rank table is not computed there because
it is expensive and does not change the rank-sum conclusion.

## The central site is isolated

At `n=8`, testing the same single-`e` construction at every spatial
site `0,...,13` gives:

- the central site `7` spans the crossing in dimension `8`;
- every other site spans only dimension `7` and misses the crossing.

Thus this is not a generic linear-algebra accident of the candidate
set. The crossing is carried by the central spatial channel, while
the earlier disjoint-`f_3` census spread thousands of pieces across
the whole seam and still missed it.

## Why this is not yet a cost certificate

The identity is a spanning certificate. It does not satisfy a
finite-block rank-sum budget. At `n=8` the two cut ranks of the
degree slices are

```text
b     cut1  cut2
0      169     1
1      169     8
2      169    43
3      169   155
4      168   165
5      169   169
6      169   167
7      161   167
```

The first cut is full (`d^2=169`) on every nonzero slice. The naive
rank-sum value is

```math
S=\sum_b\sqrt{r_bc_b}\approx964.33,
```

far above the one-label thresholds `13^(1/sqrt2)~=6.13` and
`13^sqrt2~=37.61`. No rank-sum or complete factoring improvement is
claimed.

The useful structural feature is the second-cut profile. The `b=0`
slice has second-cut rank one, and removing it drops the first-cut
rank of the remaining sum from `169` to `78`. This suggests a staged
decomposition: isolate low-second-cut-rank pieces and transport the
residual. The next target is a **fused transport or grouping identity
for this single central gate** with a controlled two-cut rank sum.

## Relation to the disjoint placement census

The [all-disjoint census](all-disjoint-fusion-placements.md) showed
that the entire disjoint central `f_3` family at thirteen labels has
span dimension `4014` and misses the crossing; its separating
functional is supported on the first boundary coordinate. The present
identity shows that the obstruction is removed by a different local
operator, not by adding more `f_3` placements. This is the concrete
escape from the boundary-supported dual certificate.

The identity is finite and field-specific. It does not assert a
uniform all-length formula for `c_n`, a uniform rank bound, or a
strict endpoint with constant identity term. Those remain separate
obligations.

## Reproduction

```sh
python3 tools/check_boundary_e_gate_crossing.py \
  --output evidence/boundary-e-gate-crossing.json
make verify-boundary-e-gate-crossing

# Heavy length-ten spot-check, not part of the default verify target.
python3 tools/check_boundary_e_gate_length10.py \
  --output evidence/boundary-e-gate-length10.json
```

The receipt records the exact coefficients, support, site census,
both cut ranks of every slice, and source hashes. No new Lean
formalization or independent mathematical acceptance is claimed.
