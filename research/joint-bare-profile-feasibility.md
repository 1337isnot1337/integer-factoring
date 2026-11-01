# A combined compatibility gate before bare tensor construction

**Later theorem:** the [box-adjugate word inequalities](box-adjugate-word-caps.md)
exclude the integer point below in 52 oriented checks. They apply the
box theorem to the derived adjugate crossing, extending the families
used in this experiment. The original seven-family checks and exact rate
certificate remain valid historical calculations. This point is no longer
an open polynomial construction target; the bare thirteen-letter criterion
itself remains open.

## Registered question

Base: `d159502`, followed by the historical thirteen-letter target audit.
The complete score remains `1/2`. This is a necessary-profile diagnostic,
not a polynomial crossing construction or a stronger complete bound.

The old twelve-letter integer points omitted box budgets; the old
thirteen-letter box point omitted stronger jet constraints. Before spending
effort on realizing a staircase, ask whether a thirteen-letter profile can
pass **all** those currently available families while retaining formal bare
rate below `1/sqrt(2)`. The proved bare exclusion through twelve letters
applies there. Thirteen letters lie outside that exclusion, but the recorded
profile itself has now been rejected by newer constraints.

The next diagnostic fixes equal cut profiles, the old normalized thermal
scale `t=2000 log(414548/413547)`, and degree one after normalization.
Include box, cell/shared, hereditary adjugate, rectangular first-jet,
multichannel and higher-word conditions. The higher-word conditions imply
the weaker second-jet conditions; any proposed point must still pass their
literal independent replay. No physical gauges or changed endpoint conventions
are introduced.

For each smaller square alphabet, retain the ordered words that attain the
old profile's first channel costs, including multiplicities. Their affine
costs give a global **upper** bound on each minimum-word prefix. Substitution
therefore supplies a sufficient linear condition for that nonlinear cap.
This is an inner approximation of the necessary relaxation. A failed fit
cannot exclude another linearization, another profile, or an actual tensor.

Run one convex fixed-scale fit from the uniformly reduced old profile,
with at most 150 iterations and 120 seconds. Continue only if a candidate
has formal rate at least `1e-5` below `1/sqrt(2)` and all compiled constraints
pass. Such a candidate requires rounding, exact replay of every original
family, a rigorous rate certificate, and eventually an actual polynomial
realization before any factoring claim. Otherwise archive the bounded result
and change mechanism; do not enlarge this profile fit or sweep alphabet sizes.

No Lean or independent mathematical review is claimed. Substantial OpenAI
assistant support.

## Exact replacement target

The registered fit finished after 32 iterations and 32 function queries in
9.03 seconds; compiling the constraints took 12.14 seconds. The floating
objective was `0.7059116982`. This is a numerical feasible fit, not an
optimality certificate. The [complete fit receipt](../evidence/joint-q13-profile-fit.json)
retains all iterates, the full proposed profile, limits, environment and
source hash. Its 5,815 distinct linear rows include every compiled family;
168 additional rows impose sortedness. Sixty abstract controls verify
1,296 word-prefix majorant comparisons; these controls are not tensors.

Multiplying the candidate by `2000*(1-1e-10)` and rounding every order
down gives the following equal cut profiles, with `M=2000` and zero
physical energies. Superscripts denote multiplicities:

```
0,498,727,867,983,1039,1145,1289^6,1557^13,1766^13,
1900^13,1977^13,2080^6,2081^7,2091^26,2153^13,
2230^13,2306^13,2365^13,2369^12,3260.
```

The [full integer point](../evidence/q13-joint-word-target.json) passes
the original box and square-box adjugate constraints and all seven families
in the following exact replay. Both cut orientations coincide for this
equal-profile point.

| Family | Checks | Minimum slack |
| --- | ---: | ---: |
| Cell/shared prefix | 170 | 0 |
| Oriented box | 182 | 14 |
| Oriented hereditary square adjugate | 1816 | 1864 |
| Oriented inherited second jet | 1816 | 1864 |
| Oriented higher-word jet | 1816 | 968 |
| Rectangular first jet, one identical orientation | 8918 | 663 |
| Multichannel cofactor, one identical orientation | 8918 | 592 |

At the original rational scale `s=414548/413547`, the
[independent exact checker](../tools/check_q13_joint_word_target.py) computes

```
B = s^2000 (sum_i s^(-mu_i))^2
```

as a rational number. It rounds `B` upward to denominator `2^96`, then
uses the existing rational logarithm bounds to prove

```
10000 log(B) < 4*7061 log(13),
2*7061^2 < 10000^2.
```

Consequently `rho_formal < 7061/10000 < 1/sqrt(2)`. The
[receipt](../evidence/q13-joint-word-target-controls.json) retains the
actual dyadic upper bound, rational logarithm bounds and compiler hashes;
the numerical value about `0.7060680472` is display information only.
Run `make verify-q13-joint-word-target` for the exact replay and the short
majorant controls. This gate does not rerun the numerical optimizer.

For completeness, the majorant construction selects actual ordered words
of length at most `m-1`, together with arbitrarily many copies of the
formal cap cost `m*s_2`. Each selected word has a fixed nonnegative integer
count vector. The sum of any selected `k` costs bounds above the sum of the
`k` smallest costs. Thus the chosen source-minimizing prefixes are global
affine upper bounds even when the minimizing words change. Words with the
same count vector retain their distinct multiplicities. Substituting an
upper bound into the necessary word inequality supplies a sufficient
linear condition for that inequality; it does not characterize the entire
necessary relaxation.

At the experiment's source this replaced the old point with one passing
the **seven declared constraint families**. The subsequently composed
box-adjugate theorem supplies a further structural obstruction: the
`7 by 13` box requires `30588956<=30576000`, which is false. The focused
checker now replays that exclusion as well as all original successful
checks and the formal rate certificate. No profile data or discovery
history is removed. The seven-family relaxation remains open, but this
point cannot be a polynomial crossing. No further profile-only optimization
is run here. The complete factoring score remains `1/2`.
