# Full physical gauges for selected balanced Gaussian codes

## Result and scope

The finite searches in this note vary all 32 independent physical
heights and the common scale for nine exact five-state Gaussian tensors.
They do not improve the retained reference ratio
`0.7263479477044683`. The best nonreference value found is
`0.729549834747631` for an asymmetric physical code.

The better nonreference values occur close to the boundary of the type
simplex. In the best asymmetric state, the fifth Y/U symbol has mass
about `1.23e-12`. In the best symmetric state, the fifth symbol has mass
about `1e-15` on both sides. These are approaches to smaller physical
alphabets. They do not establish stronger interior five-symbol seeds.

All tensor coefficients and the complete Gaussian-minor reconstruction
are exact. Heights, weighted valuations, Gibbs types, scales, and the
reported ratios are floating-point search evidence. Neither a global
minimum nor a strict pressure-sign certificate is claimed.

## Exact source tensors and parameter space

The sources are the balanced four-particle tensors from
[the symmetric code screen](balanced-four-particle-codes.md) and
[the asymmetric code screen](asymmetric-four-particle-codes.md).
The first physical code is paired between X and V; the second is paired
between Y and U. Every chosen alphabet consists of five distinct
four-particle masks in eight modes. At the strict endpoint, the Gaussian
kernel is a block swap, so every selected tensor satisfies

\[
T_{xyuv}(1)=\delta_{xv}\delta_{yu}.
\]

The constructor retains every one of the 625 polynomial entries. It
compresses the even powers in the original Gaussian variable and removes
the common entry gcd, which evaluates to one at the endpoint. The saved
records contain the literal resulting coefficients, selected masks,
energies, gcd, and coefficient digest.

At each of the two local cuts, assign five heights to each physical leg.
The height of symbol zero is fixed to zero in each of the eight vectors.
This leaves 32 real height parameters. The thirty-third parameter is a
positive common valuation scale. X and V may have different local
heights, as may Y and U. The final physical energy for a symbol is the
sum of its four appearances over both cuts. Thus the optimization also
allows distinct X/V and Y/U Gibbs types.

For fixed weighted orders \(\sigma_\infty,\sigma_0\) and physical
energies \(g_X,g_Y\), the replay solves

\[
\log\!\sum_i e^{-\sigma_{\infty,i}}
+\log\!\sum_i e^{-\sigma_{0,i}}
=2\rho\left(
\log\!\sum_x e^{-g_X(x)/(2\rho)}
+\log\!\sum_y e^{-g_Y(y)/(2\rho)}\right).
\]

This is the paired Gibbs-type form used by the rectangular crossing
oracle. It charges the complete physical energies and both complete
matrix profiles.

## Nine-source finite search

The [driver](../tools/search_balanced_code_full_gauges.py) selects the
reference, the four best distinct nonreflection symmetric coefficient
tensors, and the four best distinct asymmetric coefficient tensors from
the earlier finite screens. It tries three starts for each source:
the recorded natural gauge, a two-zero energy target, and a uniform
physical-height start. Each SLSQP run has a 250-iteration limit, height
bounds `[-20,20]`, and scale bounds `[0.0001,2]`.

The [complete record](../evidence/balanced-code-full-gauge-search.json)
contains all 27 runs. The table reports the best successfully evaluated
complete profile encountered during each source's three runs.

| Family | X/V code | Y/U code | Earlier natural ratio | Best full-gauge ratio |
|---|---|---|---:|---:|
| Reference | 15,23,43,51,77 | same | 0.7263479477 | 0.7263479477 |
| Symmetric | 15,23,43,51,228 | same | 0.7482523811 | 0.7310818099 |
| Symmetric | 15,23,43,77,228 | same | 0.7540849221 | 0.7399943833 |
| Symmetric | 15,23,43,51,226 | same | 0.7624095996 | 0.7310940237 |
| Symmetric | 15,23,43,77,226 | same | 0.7654235726 | 0.7399750709 |
| Asymmetric | 15,23,43,51,77 | 15,23,43,51,228 | 0.7415036110 | 0.7311055681 |
| Asymmetric | 15,23,43,51,77 | 15,23,43,77,228 | 0.7458433484 | 0.7343789753 |
| Asymmetric | 15,23,43,51,77 | 15,23,43,51,226 | 0.7527790240 | 0.7296237140 |
| Asymmetric | 15,23,43,51,77 | 15,23,51,77,228 | 0.7541898896 | 0.7452958312 |

There were 14,650 optimizer-reported objective evaluations and 945
recorded failed oracle evaluations. A failed oracle call receives a
penalty and cannot become a retained witness. The driver stores the
best successful complete-profile evaluation separately from the final
optimizer iterate. It also records optimizer success flags and failure
messages. Several useful retained states came from runs that reached
the iteration limit; their validity as evaluated states does not imply
convergence to a local minimum.

The full nine-source schedule completed in about 45 seconds. Its
scope is exactly the saved source selection and finite start schedule.

## Refinement and type-boundary behavior

The [refinement wrapper](../tools/refine_balanced_code_gauges.py) exports
the literal best nonreference symmetric and asymmetric tensors, together
with their saved starting parameters. It then invokes the existing
paired-gradient rectangular optimizer for 32 restarts per source, a
maximum of 20,000 evaluations per restart, and height bounds `[-80,80]`.
The time budget was 1,200 seconds per source; both finite restart
schedules finished much earlier. The provenance records retain the
exact invocation, source hashes, seed hashes, starting hashes, and
result hashes.

| Refinement | Attempted profile queries | Failed queries | Best recorded ratio | Fifth-symbol type mass |
|---|---:|---:|---:|---|
| Asymmetric Y/U code 15,23,43,51,226 | 17,401 | 2,717 | 0.729549834747631 | Y/U: `1.23435e-12` |
| Symmetric code 15,23,43,51,228 | 12,256 | 727 | 0.7310816701980001 | X/V: `1.13424e-15`; Y/U: `1.37766e-15` |

The asymmetric X/V type remains genuinely supported on all five
symbols in the numerical record: its fifth mass is about `0.00804`.
The Y/U side collapses toward four symbols. The symmetric refinement
collapses toward four symbols on both sides. Therefore the numerical
improvement over those sources' natural gauges should be interpreted
as a boundary effect. No positive lower bound on the new-symbol type
mass was imposed in these finite runs.

The two complete provenance records are
[asymmetric](../evidence/balanced-asymmetric-gauge-refinement-provenance.json)
and [symmetric](../evidence/balanced-symmetric-gauge-refinement-provenance.json).
Each links its seed, start, and search record.

## Independent replay

The [checker](../tools/check_balanced_code_full_gauges.py) reconstructs
all 625 entries of each of the nine source tensors as literal strict
Gaussian minors. It independently assembles the two complete 25-by-25
matrices for each saved gauge and sends their row and column weights to
the general matrix oracle. It then solves the paired Gibbs equation
separately, without using the rectangular oracle's type routine.

The [retained replay](../evidence/balanced-code-full-gauge-replays.json)
passed all 5,625 literal minor identities, all 29 saved gauge states,
and all 58 complete weighted matrix profiles. The 29 states consist
of the 27 initial-run witnesses and the best state from each refinement.
Saved types, physical energies, weighted orders, and ratios agree
within the checker's stated floating-point tolerance.

This verifies the exact tensor construction and the numerical replay
of the retained states. It does not certify the floating-point
optimizer's global search, all intermediate trial profiles, or an exact
pressure inequality.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/search_balanced_code_full_gauges.py --output /tmp/balanced-code-full-gauge-search.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_balanced_code_full_gauges.py evidence/balanced-code-full-gauge-search.json --refinement evidence/balanced-asymmetric-gauge-refinement-provenance.json --refinement evidence/balanced-symmetric-gauge-refinement-provenance.json --output /tmp/balanced-code-full-gauge-replays.json
```

The first command repeats an exploratory floating-point schedule and
need not reproduce every optimizer iterate across numerical libraries.
The second replays the immutable saved witnesses and their source
provenance.
