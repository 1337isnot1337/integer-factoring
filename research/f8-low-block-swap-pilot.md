# The word-31-for-46 swap pilot is dominated and cannot cross

**Result.** Swap the unique five-particle height-25 word `31` out of the q96
F8 hyperoval map and add the admissible height-26 word `46`. A matched
cutoff-47 profile resolves `935` of `9216` orders per cut, leaves `8281`
tails, and fits numerical pressure `0.7098808038921344` with
`scale_t=0.24576929194426034`. The [saved profile](../evidence/f8-low-block-swap-pilot.json)
is reproducible with `make build-f8-swap-pilot`; the [checker](../tools/check_f8_low_block_swap_pilot.py)
reconstructs the maps from the q96 receipt and reruns the pressure.

The fit is above `1/sqrt(2)` by about `0.002774` and is dominated by the
word-46-addition pilot `0.7098643133625308`. Because the map has a
height-zero letter and all orders are nonnegative, the truncated profile is
a [lower bound](truncated-profile-lower-bound.md) on the complete pressure.
The complete swap therefore cannot cross, and no full build is warranted.

This extends the matched cutoff-47 diagnostics of the
[low-block pilot census](f8-hyperoval-minimal-core.md#matched-cutoff-low-block-diagnostics)
by one multi-word vector: dropping the boundary word and adding a high word
does not improve on the addition alone. It remains a local vector change,
not the different rank-three image needed for a structural low-block
improvement. No exact profile, certificate, or factoring exponent changes;
the score remains `1/2`.
