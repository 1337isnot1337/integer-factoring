# Free physical heights and supporting angles on the retained F4 pencil

The completed pencil census used fixed physical heights. Its determinant
ceiling excludes a sufficient asymmetric escape at those heights; it does
not exclude other paired energies or free physical gauges. A new registered
experiment optimizes those heights together with the asymmetric supporting
angle. A supporting value below **1** would be a candidate for exact replay,
using the existing [supporting-line](asymmetric-support-line-crossings.md)
and [rectangular transfer](asymmetric-modular-transfer.md) proofs.

No such exact improving candidate is claimed at registration. The complete
factoring score remains **1/2**.

## Selection and search space

The [manifest](../evidence/two-place-free-gauge-asymmetric-manifest.json)
selects one representative from each of the 67 distinct complete natural-
profile classes in the [finished F4 census](two-place-f4-census.md). It
chooses the smallest numeric representative appearing among the saved
chunk representatives for that class. This is not a claim that it is the
smallest member of the entire class. Equality of profiles at the old
heights does not imply equivalence under new heights.

Each selected tensor is reconstructed from the explicitly credited PR57
pencil and its two F4 coefficient planes. The manifest fixes all selections,
literal-value hashes, source archive hash, natural profiles, seed, bounds,
seven starting angles and time budgets before the new search. Its 469
listed trials traverse all cases for one restart before beginning the next.
Each optimizer call has fifteen seconds and at most 1,000 iterations; the
global ceiling is two hours. Attempted trials, timed-out calls, optimizer
returns and global time-budget completion are separate counters.

All thirty-two nonreference physical heights are free, giving independent
potentials on the four legs at both places. Physical costs on both matching
pairs are recomputed at every query. The normalized height bound is 32
in absolute value; intrinsic scale ranges from `exp(-9)` to 4 and the
supporting angle ranges from `atan(exp(-8))` to `atan(exp(8))`. These
bounds leave the finite search incomplete outside the declared domain.

## Existing oracle and stable coordinates

The experiment reuses the existing complete
[asymmetric pressure oracle](../tools/asymmetric_crossing_pressure.py).
For unit weights `(a,b)`, it computes the positive paired Gibbs root

```math
a\log\sum_i e^{-\nu_{\infty,i}/a}
+b\log\sum_i e^{-\nu_{0,i}/b}
=R\left(\log\sum_xe^{-G_X(x)/R}+\log\sum_ye^{-G_Y(y)/R}\right).
```

The [coordinate wrapper](../tools/scaled_asymmetric_crossing_pressure.py)
uses actual heights `h=t*H`, `t=exp(log_t)`, and
`theta=atan(exp(log_tan))`. Its scale derivative is the dot product of
the native derivative with all actual height and scale variables; its
angle derivative includes `sin(theta)*cos(theta)`. These changes preserve
the objective and every complete profile. The predecessor sources remain
unchanged.

Run `make verify-scaled-asymmetric-pressure`. The
[checker](../tools/check_scaled_asymmetric_pressure.py) reconstructs all
67 selected literal tensors and hashes. Three binary/nonbinary cases have
150 complete orders independently checked by fraction-free elimination at
new integer gauges. It also verifies balanced normalization
`R=sqrt(2)*rho` and 204 finite-difference derivative comparisons. The
[receipt](../evidence/scaled-asymmetric-pressure-controls.json) retains
the controls and exact source hashes. Numerical derivatives are distinct
from exact profile checks and from a pressure certificate.

The [driver](../tools/search_two_place_free_asymmetric.py) retains exact
source commit/tree, as-run code, full best literal tensor, complete profiles,
both paired types and costs, individual optimizer outcomes and rejected
native queries. First-restart inputs are checked against every saved natural
profile. All native profile windows must resolve; failures are counted as
rejected queries. A numerical value below one still needs independent
untruncated profiles and an outward rational pressure inequality before
supporting a factoring improvement.

Substantial OpenAI assistant support was used. No new Lean formalization
or stronger factoring exponent accompanies this registration.
