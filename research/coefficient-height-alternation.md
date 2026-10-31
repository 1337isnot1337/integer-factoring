# Exact coefficient walks after physical-height optimization

## Result and scope

Two complete 30-minute search budgets varied every degree-bounded F4
coefficient around previously optimized four-symbol references. Together
they completed 61,041 target-walk steps, evaluated 2,187 full profiles,
and recorded 1,432 distinct profile classes. Neither improved its reference
ratio. All eighty retained literal tensors and their full cost calculations
were independently replayed.

The runs had requested 50,000 steps each, but stopped at their explicit
time budgets before reaching that count. The records preserve
`completed_requested_trials: false`. These are finite searches at selected
frames and heights; they do not exclude the larger construction family.

## 1. Alternation between physical heights and complete coefficients

The earlier complete-entry affine search produced non-Gaussian F4 tensors
at a degenerate physical-energy choice. Optimizing all independent
physical-letter heights moved two of those tensors into substantially
better chambers. Their rounded integer-height references are
[`affine-degenerate-witness-0-integer-chamber.json`](../evidence/affine-degenerate-witness-0-integer-chamber.json)
and [`affine-degenerate-witness-1-integer-chamber.json`](../evidence/affine-degenerate-witness-1-integer-chamber.json).

Their binary coefficient planes have exactly the same two complete
profiles as the original F4 tensors at the recorded integer heights.
The [finite regular-frame checker](finite-precision-smith-frames.md)
verifies that equality independently. Those binary tensors therefore
supply exact polynomial frames for another full F4 coefficient search.

Both new affine models have degree 24 and 6,144 independent F4 variables,
with the strict crossing endpoint fixed. Common row floors eight below
the reference orders reduce them to 694 and 684 free F4 coordinates.
Every row target is allowed up to twelve above its reference order.
The coefficient tables are complete through those caps.

Each target walk lowers a selected group of constraints, then greedily
raises rows in a randomized order while maintaining exact feasibility.
Three F4 witnesses are sampled at each new feasible target vector.
Every reported cost uses both complete native Smith profiles and both
physical type-energy lists. Row-target bounds alone are never reported
as a realized pressure witness.

## 2. Completed time budgets

| Reference | Completed steps | Distinct target vectors | Full profile queries | Distinct profile classes | Best rho, including reference |
|---|---:|---:|---:|---:|---:|
| Integer chamber 0 | 27,104 | 331 | 993 | 759 | 0.7450485723725101 |
| Integer chamber 1 | 33,937 | 398 | 1,194 | 673 | 0.7449315952026363 |

Each displayed best ratio is the unchanged reference value. Both exceed
the critical threshold 1/sqrt(2), approximately 0.70710678. The records
therefore contain no factoring improvement.

The first run used 1,800.10 seconds and the second 1,800.06 seconds,
including construction and evaluation. Their elapsed times and unfinished
step requests are separate fields in the evidence; a completed time
budget is not presented as a completed 50,000-step catalogue.

## 3. Full replay and continuation

The replay checker verifies source hashes, the exact reference files,
the row windows, all saved strict endpoints, full untruncated row-frame
products, complete profiles, and recomputed roots. It also independently
reconstructs the two reference profiles through binary Smith calculations
and sixteen block-Toeplitz ranks.

```sh
make verify-coefficient-height-walks
```

The full records are
[`full-alphabet-alternating-walk-0.json`](../evidence/full-alphabet-alternating-walk-0.json)
and [`full-alphabet-alternating-walk-1.json`](../evidence/full-alphabet-alternating-walk-1.json),
with the replay summary in
[`full-alphabet-alternating-replays.json`](../evidence/full-alphabet-alternating-replays.json).

The reusable checker
[`check_full_alphabet_walk_records.py`](../tools/check_full_alphabet_walk_records.py)
also handles other alphabet sizes and explicitly timed runs. Larger Laurent
windows, full-floor constraint releases, other physical heights, and
nonbinary frames are distinct searches. Failure in these two recorded
chambers does not decide them.
