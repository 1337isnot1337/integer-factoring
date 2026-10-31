# Ordinary-proof frontier session, 31 October 2026

## Outcome and the score target

The supported deterministic factoring exponent remains **1/2**:

\[
T(b)=\exp\!\left(O\!\left(\frac{b}{\sqrt{\log b}}\right)\right),
\]

where \(b\) is input bit length. This session has not constructed a
strict crossing witness below \(\rho=1/\sqrt2\). The existing
[quantitative transfer](pieces-to-factoring.md) would turn a complete
effective family at a fixed rational \(\rho\) into the exponent
\(\alpha=\log_2(1/\rho)\). That transfer requires the full crossing,
factorizations, coefficient bounds, and uniform construction costs.
A favorable local rank or a floating relative improvement alone does
not meet its hypotheses.

The work below consists of ordinary proofs, exact finite-field
constructions, outward-rational certificates, and explicitly scoped
numerical searches. No new Lean formalization accompanies this session.
The prior [constructive checkpoint](constructive-frontier-checkpoint.md)
remains available as the earlier research state.

## 1. Much larger exact finite certificates

The [accelerated interval certificate](accelerated-interval-certificates.md)
now reaches **\(2^{1636}-1\) modes**, with a rational upper bound
less than **\(2\cdot10^{-980}\) above \(1/\sqrt2\)**. It approaches
the threshold from above, so it does not improve the score.

The useful algorithmic advance is an exact finite moment identity
with a positive remainder. Binary exponentiation evaluates the long
geometric products, and a triangular difference table encloses their
logarithms. The certificate uses a number of integer operations
polynomial in the logarithm of the mode count. It can also certify
very long Gaussian spectator products around a fixed defect.
This does not provide a compact explicit listing of the enormous
physical alphabet.

The [all-scale balanced obstruction](all-scale-balanced-obstruction.md)
also removes the earlier scale restriction for its two stated repairs.
The [all-digit interval floor](truncated-symmetric-interval-frontier.md)
proves positive natural critical pressure for every digit count,
interval length, and positive scale. These are uniform results for
their declared families, not lower bounds for factoring itself.

The [arithmetic Cauchy-power theorem](arithmetic-cauchy-powers.md)
also determines every invariant factor of
\(H_{ij}=(T+d(i+j))^{-m}\), when \(T\) is a \(p\)-adic unit
and \(e=v_p(d)\ge1\):

\[
s_r=2re+v_p(r!)+v_p((m)_r),\qquad 0\le r<n,
\]

where \((m)_r=m(m+1)\cdots(m+r-1)\). The result holds for
every prime, positive integer power, and matrix size. It supplies exact
spectral data; a cheaper crossing or scalar decoder is still needed
to turn it into a factoring improvement.

## 2. Constructive cores around every prime power

The [prime-power selected-mode theorem](prime-power-joint-identity-filters.md)
uses weighted reflected Pascal identities to support arbitrary selected
identity indices at \(n=p^k\). Its zero core depends on the number of
selected directions, rather than the whole interval length.
Joint polynomial filters retain complete coefficient sources and
the full physical cost.

The new neighboring-length constructions are:

| Interval and selection | Actual active zero modes | Auxiliary Weyl matrix order | Main qualification |
|---|---:|---:|---|
| \(n=p^k\), any \(r\) selected indices | \(2r\) | \(p^r\) | Exact prime-power reflected basis |
| \(n=p^k-d\), \(r\) mixed initial/terminal indices | \(2r+d\) | \(p^{r+d}\) | Stated size condition and \(d\) exact Gaussian-digit divisions |
| \(n=p^k-1\), any \(r\) selected indices | \(2r+1\) | \(p^{r+1}\) | \(r+1\le n\), one exact Gaussian-digit division |
| \(n=p^k+d\), any \(r\) selected indices | At most \(2(r+2d)\) | At most \(p^{r+2d}\) | \(1\le d\le p^k\), extra directions are unchanged actual identities |
| \(n=p^k+d\), last \(r\) indices | \(2(r+d)\) | \(p^{r+d}\) | \(d\le r\le p^k\) |

The proofs are in [fixed deficits](fixed-deficit-identity-cores.md),
[arbitrary one-deficit selections](arbitrary-one-deficit-identity-cores.md),
and [arbitrary surplus selections](surplus-identity-cores.md).
The last two remove an important restriction to end selections.
For a global crossing, the original selected operator still needs the
stated invertible energy-filtered endpoint or the strict identity endpoint.

For one deficit, the contraction has only one non-scalar local
eigenchannel. Its weighted regular eigenbasis need not be a valid
global endpoint change: it mixes actual and formal directions, and
its endpoint determinant vanishes when \(p\mid r+1\).
In particular, the two-mode qutrit case has this singularity.

The geometry checker covers 505 one-deficit selections, including
all 325 pairs at \(n=26\), with 24,883 interpolation identities
and 5,111 compact-versus-adjugate identities. The surplus checker
has complete polynomial identities in characteristics two, three,
and five, including regressions to the separately implemented binary
one-surplus construction. Geometry checks are distinguished from
literal full-tensor controls in the individual notes.

## 3. Off-diagonal operators with complete physical costs

The [binary near-power engine](near-power-identity-operator-reduction.md),
[qutrit deficit engine](qutrit-deficit-operator-search.md), and
[binary above-power engine](above-power-identity-operator-reduction.md)
construct exact operator profiles. They retain factorial units and
signs, local denominator clearing, every physical spectator, both
paired physical partitions, and ordinary polynomial degree.

Completed finite searches include:
Every "no gain" entry below means no detected gain in the stated finite
grids and refinements; it is not an all-scale obstruction.

| Family | Completed scope | Finding |
|---|---|---|
| Binary terminal/split operators at \(n=31\) | 42,560 specified kernels per geometry | No gain in the recorded grids and refinements |
| Binary operators at \(n=33\) | 42,560 kernels, 1,111 exact profile classes | No gain |
| Qutrit terminal/split operators at \(n=26\) | 1,129 strict-endpoint kernels per geometry | No gain |
| Qutrit parent-exponent and energy-conjugated operators | Eight further complete runs, 9,032 kernel/geometry cases | No gain; interior selected modes included |
| Resonant Laurent operators | 11,796,480 kernels in each of two overlapping support windows | No gain in the declared finite schedule |
| Joint binary and prime-field polynomial filters | Eight retained catalogues and 6,899 literal witnesses | No gain in the declared finite schedule |

The [prime-field affine constructor](prime-deficit-affine-operator-synthesis.md)
then allows independent coefficient polynomials in all nineteen
occupation-preserving qutrit entries. Its exact identity-endpoint
constraints use 57 variables on the unit window or 95 variables on
five specified sparse powers. Four completed searches made 4,985
exact affine queries along 242 prescribed greedy paths. The 27 initial
or terminal profile records, counted across the runs, give no detected
gain. These paths do not exhaust the coefficient spaces or every
intermediate feasible witness. A fresh replay checks all 16,721
coefficient equations and all 54 cuts of those saved outputs against
the final source.

The [arbitrary-index API](arbitrary-one-deficit-operator-search.md)
has ten original physical cut controls and thirty-six independent
full-cost comparisons. The
[energy-conjugated family](deficit-energy-conjugated-operators.md)
changes the actual coefficients and charges their complete Laurent
support; its diagonal conjugation is not treated as free.

See the [above-power catalogue](above-power-operator-searches.md),
[resonant catalogues](resonant-laurent-catalogues.md), and
[joint-filter records](joint-filter-search-records.md) for exact
supports, multiplicities, grids, and scope. Grouped coefficient
coverage is not a count of individually constructed dense matrices.

## 4. Every nonnegative physical bulk ordering is available

The [reflection compression](reflection-dual-rail-compression.md)
reduces the physical alphabet from \(2^M\) to \(2^{M/2}\) while
retaining the Gaussian rate. This makes bounded defects cheaper to
represent, without changing the score.

The [all-order reflection construction](reflection-hot-all-bulk-gauges.md)
uses two explicit regular cross pivots and physical-leg swaps. It
supports every ordering of nonnegative physical bulk heights with
the same bounded active core. Hot-state compensation can be negative,
but all cold-bit costs and the exact regularity inequalities remain
in the calculation.

The [separate-leg searches](reflection-hot-leg-bulk-searches.md) and
[all-order searches](reflection-hot-all-bulk-searches.md) completed
five finite schedules. They include 25-candidate fixed-cone runs at
intervals 64 and 128, all 16 cones for 13 tensors at both lengths,
and a further interval-scaled box at 128. The best values remain the
Gaussian values, approximately 0.7075072113 and 0.7072351452.
Both exceed \(1/\sqrt2\).

The saved-point replay reconstructs all 543 candidate/cone best entries,
their spectra, and their physical energies. Some entries describe the
same point. It reproduces all pressures exactly on the recorded
platform and verifies 272 chain-rule comparisons. It does not certify
global minima over continuous parameter boxes.

One optimizer proposal had a negative slack of about \(2.12\cdot10^{-22}\).
The strict regularity gate rejected it before scoring. The successful
fresh run projects proposals into the declared box and records the
points actually evaluated. No regularity inequality was weakened.

## 5. Complete binary endpoint fibers and finite Laurent exclusions

The [two-cut affine construction](binary-two-cut-affine-synthesis.md)
assigns an independent polynomial to every one of the sixteen tensor
entries, subject only to its literal permutation endpoint. The unknown
tensor need not preserve occupation or factor as a Gaussian.
Prescribed local row divisibilities become exact simultaneous affine
equations over \(\mathbb F_4\).

For the fixed frames and physical heights, the completed search uses
all 8,925 constant flags independently at each cut and 51 explicit
targets. Only the reference target is feasible. Its solution is an
affine line over \(\mathbb F_4\), and all four points have the same
two Smith profiles. All fifty other targets are infeasible.
Necessary common-row conditions reject whole groups safely; only
748 flag pairs require direct intersection across the complete run.

A separate support argument proves that every Laurent tensor meeting
those targets lies in the window \([-4,8]\). Clearing by \(w^4\)
therefore reduces the unrestricted Laurent question for those targets
to degree twelve. This is stronger than merely testing a convenient
degree cutoff. It remains a conclusion about the specified targets,
regular frames, and physical height ray.

A larger simultaneous-target schedule then completed 7,079 selected
row-target vectors, chosen from 23,099 raw sorted vectors, with **129
selected F4 flags per cut**. Its five feasible systems are five monomial
shifts of the same four-point affine line. Every point on every line
was reconstructed, with 50 all-minor cut profiles and 400 endpoint
coefficients checked. Necessary common constraints independently exclude
1,430 targets regardless of the constant flags. The other exclusions
retain the selected-flag scope; this larger target experiment does not
inherit the all-8,925-flag coverage of the separate 51-target run.

## 6. Global Gaussian pressure and complete one-mode exclusion

The [global charged Gaussian theorem](global-charged-gaussian-floor.md)
proves

\[
G_n(J,t)>0\qquad(n\ge1,\ J\in\mathbb R,\ t>0)
\]

at the critical parameter \(\rho=1/\sqrt2\), including both cut
costs and the complete paired physical partition. The analytic statement
holds at every integer length. It concerns the aggregate-charge Gaussian
family; independent physical state gauges remain outside its scope.
The digit-factorization identity gives the corresponding analytic floor
for power-of-two digit counts. It does not assert a separate prime-state
construction or an all-charge theorem for odd digit counts.

The proof combines the four-periodic continuation, a convex-tail
comparison, exact charge folding, certified small-length bases, and a
positive scalar kernel for the final edge range. Ordinary analytic
arguments handle the unbounded parameters and scale tails. The compact
certificates use 290 rectangles for the periodic continuation, 1,627
rectangles for the small bases, and 250 intervals for the final scalar
kernel, all with exact outward-integer bounds. These are rigorous
statements about the displayed pressure, not a new factoring score or
new Lean formalization.

The [charged one-mode Laurent theorem](charged-one-mode-binary-laurent-filters.md),
combined with the global Gaussian floor, excludes every admissible
one-mode binary Laurent occupation filter **as a strict witness for
this pressure criterion**, at every aggregate charge and scale. It
permits arbitrary finite Laurent support, arbitrary multilevel
cancellation, every selected direction, and any fixed finite coefficient
field of characteristic two, under the stated endpoint and degree
accounting. The filter construction retains its power-of-two interval
lengths.

For \(|J|\le2n\), the filter costs at least as much as the Gaussian,
with equality only for a common scalar Laurent monomial. The proof uses
an additional coefficient-at-\(x^{nP}\) barrier, beyond the elementary
degree/order bounds. Exterior charges can give genuine relative gains.
A sharp strict loss cap and a stronger exterior Gaussian estimate,
with its own 70-interval certificate, keep the full pressure positive.
Nontrivial exterior filters can also have ratio exactly one; the
common-monomial equality classification is an interior statement.

The correlated [homogeneous Boolean family](homogeneous-boolean-filter-bounds.md)
has a different scope. Its full pressure is positive for every real
charge at or above the explicit threshold \(t_n^*\), uniformly over
all admissible patterns and degrees. At every scale, positivity holds
for each fixed Boolean degree bound and fixed charge bound once the
binary prime-power length is sufficiently large; at zero charge it
holds for every admissible length. The general correlated small-scale
argument still requires its quantitative charge and degree bounds.
The global Gaussian floor alone does not absorb an arbitrary negative
relative correction. Arbitrary joint multilevel filters remain open.

Zero charge is not a global pressure minimizer. The exact certificate
at \(n=32,t=3,J=323/5\) has a smaller charged pressure than at zero
charge, while both are positive. Thus global positivity does not
justify fixing the optimizer at \(J=0\).

## 7. Verification, interruptions, and remaining routes

The repository separates written uniform arguments, exact finite
identities, rational sign certificates, and numerical optimization.
The named Makefile gates replay the focused constructions. Literal
physical tensors check the reduced engines at tractable sizes, while
the written proofs supply their uniform statements.

The final [validation record](../evidence/frontier-session-20261031-validation.json)
retains the commands and their scope. `make verify` passed all 55 unit
tests and the default frontier/rate checks. The complete
`make verify-global-one-mode-filters` dependency chain also passed,
including 1,917 positive rectangles and 320 positive scalar intervals.
The [integration log](../evidence/frontier-session-20261031-verify.log)
and [global-proof log](../evidence/frontier-session-20261031-global-verify.log)
are retained. All 113 tracked Lean source files and the Lean toolchain
match the starting snapshot; no local Lean build was run.

All four complete qutrit affine search drivers were freshly replayed
against their final source. Every saved JSON field matched, apart from
elapsed time. This establishes current-source reproducibility without
inventing historical source hashes for the earlier runs.

An execution outage interrupted the session. A saved nine-candidate
prefix of the interval-128 separate-leg search was validated and
resumed; its record identifies that boundary. A separate raw-core
transition screen retains only its documented 167-of-615 prefix.
Unretained or interrupted explorations are not described as completed
searches. In particular, the earlier partial polynomial-only complete
flag run is superseded by the completed Laurent run above.

Coupling several modes or using off-diagonal operators escapes the
complete one-mode exclusion. Independent physical gauges lie outside
the aggregate-charge floor. Enlarging coefficient support can escape
the finite searches, but not the unrestricted one-mode Laurent theorem.
The new cores make these coupled searches practical on both sides of
prime powers. For complete binary tensor synthesis, different physical
heights, different regular frames, and nonconstant regular flags remain
outside the finite exclusions.

Any candidate improvement must still provide a strict pressure witness
below \(1/\sqrt2\), with complete coefficients, valid endpoints,
both physical costs, and the uniform transfer data. None of the
completed packages in this checkpoint supplies that witness.
