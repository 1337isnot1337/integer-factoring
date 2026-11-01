# Six-state balanced extensions and positive added-symbol types

## Result and historical comparison

The complete one-state extension screen contains 196 exact tensors:
the retained five-state reference, all 65 six-state symmetric
extensions, and both rectangular orientations for each extension.
The strongest natural-gauge result is
`rho = 0.7226932545985394`, with added four-particle mask 85.

This value is already present in the earlier
[236-alphabet reflection screen](reflection-code-alphabets.md).
It is the six-letter reflection alphabet `{0,2,4,8,10,12}`.
The new certificate below identifies the historical profiles exactly;
this is a reproduced benchmark, not a new best crossing ratio.
The earlier full sixteen-letter reflection ratio remains about
`0.7166536597`. The proved factoring exponent remains `1/2`.

The new finite coverage consists of every nonreflection one-state
extension, both rectangular orientations, unrestricted physical-height
optimization of selected extensions, and positive lower bounds on the
added symbol's type mass. These searches found no improvement over
the six-state benchmark or the existing full-alphabet result.

## 1. Exact separate-alphabet construction

Start with the four-particle masks

\[
C=\{15,23,43,51,77\}
\]

inside the 70 four-particle states on eight modes. For each omitted
state `a`, test the three paired-alphabet choices

\[
(C\cup\{a\},C\cup\{a\}),\qquad
(C,C\cup\{a\}),\qquad
(C\cup\{a\},C).
\]

The first alphabet is used on X/V and the second on Y/U. The tensor
dimensions are therefore `(A,B,B,A)`, with `(A,B)` respectively
`(6,6)`, `(5,6)`, or `(6,5)`. The
[constructor](../tools/balanced_particle_extensions.py) evaluates every
selected strict Gaussian minor, compresses the even powers, and divides
the exact common polynomial gcd. The gcd has value one at the endpoint.
Consequently all normalized tensors retain

\[
T_{xyuv}(1)=\delta_{xv}\delta_{yu}.
\]

The polynomial construction and even-power identity are the same as
in [the five-state balanced construction](balanced-four-particle-codes.md).
Both local matrices now have size `A*B`. Their ranks and Smith profiles
are recomputed from all entries, including the added physical state.

Among the 65 omitted states, eleven are reflection-compatible and
54 are not. Thus the schedule includes 162 nonreflection extension
tensors over the three placements. It does not enumerate all
six-state subsets of the 70 states, coherent physical subspaces,
or choices of different added states on the two paired alphabets.

## 2. Complete natural-gauge screen

The physical energy of a mask is its weighted occupation grade minus
the minimum grade in its alphabet. X/V and Y/U may have different
energy vectors. Infinity uses row height `E_Y(u)/2` and column height
`E_X(v)/2`; zero uses row height `(E_X(x)+E_Y(y))/2` and zero column
height. The scalar infinity shift and the normalization gcd are
retained consistently.

The [complete record](../evidence/balanced-particle-extension-natural-search.json)
contains every one of the 196 cases, their literal alphabet recipes,
complete integer profiles, coefficient digests, natural energies,
and numerical scale/Gibbs costs. The search completed in about
43 seconds.

| Placement | Best added mask | Best natural ratio |
|---|---:|---:|
| Symmetric six by six | 85 | 0.7226932545985394 |
| Rectangular five by six | 85 | 0.7247105375670271 |
| Rectangular six by five | 85 | 0.7247105375670271 |

The best nonreflection natural extensions use masks 216 or 228, which
give the same normalized tensor in this screen. Their symmetric ratio
is about `0.7487772058`; the rectangular ratio is about `0.7418532737`.

The [independent replay](../evidence/balanced-particle-extension-natural-replay.json)
reconstructs all 201,865 tensor entries as complete polynomial minors
of the closed strict Gaussian kernel. It also certifies all 392
complete original-matrix profiles with exact regular polynomial frames,
without calling the native Smith oracle. Every check passed. Numerical
scale minimization remains exploratory; exact profile verification does
not certify a global minimum over scale.

## 3. Exact benchmark bracket and positive sixth mass

For mask 85 the canonical physical alphabet is

\[
(15,23,43,51,77,85),\qquad E=(0,2,6,8,10,12).
\]

Its normalized degree is 24 and its common entry gcd in the compressed
variable is `z^4`. The
[standalone checker](../tools/check_six_state_reflection_benchmark.py)
reconstructs all 1,296 Gaussian minors and both complete Smith profiles.
It matches the historical reflection record by the exact identities

\[
\sigma_{0,\mathrm{old}}=2\sigma_{0,\mathrm{new}}+8,
\qquad
\sigma_{\infty,\mathrm{old}}=2\sigma_{\infty,\mathrm{new}}-8,
\]

with the old physical energies equal to twice the new energies after
permutation. The shifts cancel in the product of profile partition
functions. The factor of two is the change of polynomial variable.

At the exact rational scale

\[
t=551453777/10^9,
\]

the [integer interval certificate](../evidence/six-state-reflection-benchmark-certificate.json)
proves

\[
0.7226932545<\rho(t)<0.7226932547.
\]

It also proves that the Gibbs mass of the added symbol on each paired
alphabet lies strictly between `0.00623` and `0.00625`. This symbol
has the largest energy, so all six symbols have positive mass bounded
below by the same lower bound. This is a genuine six-symbol witness,
but it is the previously retained reflection benchmark.

The arithmetic uses outward dyadic intervals with denominator `2^100`.
Exponential and logarithm bounds have explicit Taylor remainders.
Logarithms of larger partition sums are reduced by powers of two before
the short series is applied. The two signs bracket the unique Gibbs
root at the displayed scale. They do not bracket a global minimum
over all scales or heights.

## 4. Physical gauges with an added-symbol mass floor

A paired alphabet of sizes `A,B` has `4(A+B-2)+1` independent search
variables: four physical legs at each of two places, one reference
height removed from every vector, and a positive common scale. This
gives 41 variables for six by six and 37 for the rectangular cases.

Unconstrained types can hide an ineffective new state by sending its
mass toward zero. The
[bounded Gibbs routine](../tools/bounded_paired_gibbs.py) therefore also
imposes a specified lower mass on the added state. All other type
coordinates retain their usual nonnegative constraint. In a rectangular
case, the floor applies only to the side that contains the new symbol.

For fixed physical energies `g`, positive rate `R`, and lower masses
`ell`, the type problem is

\[
\min_{p_i\geq\ell_i,\,\sum p_i=1}
\left(\sum_i g_i p_i+R\sum_i p_i\log p_i\right).
\]

Strict convexity gives the unique solution

\[
p_i=\max\{\ell_i,c e^{-g_i/R}\},
\]

with `c` chosen to normalize the free coordinates. The implementation
fixes coordinates that violate their floor and resolves the remaining
normalization. At most the number of type coordinates can become fixed.
The free coordinates have equal values of `g_i+R log(p_i)`; active
coordinates have a value at least as large. These conditions are also
sufficient for the constrained type optimum.

Let `L` be the logarithm of the complete two-profile partition product.
The paired equation is

\[
L+\phi_X(R)+\phi_Y(R)=0,\qquad \rho=R/2,
\]

where each `phi` is the constrained minimum above. Its derivative in
`R` is minus the paired entropy. The usual envelope gradient therefore
divides the fixed-type gauge derivative by twice that entropy. This
removes the type variables from the numerical height search while
retaining their mass constraints.

The roots, active-set arithmetic, and gradients in these search runs
are floating point. Their mathematical formulas do not turn the
height optimizer into a global certificate.

## 5. Completed full-gauge schedules

The first [gauge schedule](../evidence/balanced-particle-extension-gauge-search.json)
selects the six best distinct natural tensors in each placement,
plus the five-state reference. All six selected added masks in each
placement are reflection-compatible. Each extension is tested at
mass floors `0`, `0.001`, and `0.01`, with natural, interior warm,
and uniform starts. The reference uses three unconstrained starts.

Every SLSQP run has a 250-iteration limit, height bounds `[-30,30]`,
and scale bounds `[0.0001,2]`. The interior warm start embeds the
reference five-state heights and assigns an initial two-percent
Gibbs mass to the new symbol before profile recomputation.

The schedule completed 165 runs on 19 tensors, with 93,997 profile
queries and 527 rejected trial evaluations. It found no improvement
over the retained natural benchmark.

| Placement | Best with no floor | Best with floor 0.001 | Best with floor 0.01 |
|---|---:|---:|---:|
| Six by six | 0.7226932546 | 0.7226932546 | 0.7230736424 |
| Five by six | 0.7247105376 | 0.7247105376 | 0.7252328844 |
| Six by five | 0.7247105376 | 0.7247105376 | 0.7252327300 |

All table entries use mask 85. The one-percent runs attain the imposed
mass floor; the smaller floor is inactive at the natural benchmark.

A separate [nonreflection schedule](../evidence/nonreflection-particle-extension-gauge-search.json)
selects the three best distinct nonreflection tensors in each
placement: masks 216, 184, and 180. It uses the same starts, floors,
bounds, and iteration limits. All 81 runs completed, with 45,969
profile queries and 1,196 rejected trials.

| Placement | Best with no floor | Best with floor 0.001 | Best with floor 0.01 |
|---|---:|---:|---:|
| Six by six | 0.7270575830 | 0.7375497571 | 0.7930359785 |
| Five by six | 0.7266850782 | 0.7320998797 | 0.7652756849 |
| Six by five | 0.7265708831 | 0.7322342607 | 0.7652792742 |

The best unconstrained nonreflection value is `0.7265708831` for
the six-by-five mask-184 extension. Its added-symbol mass is about
`9.90e-13`, so it approaches the five-state boundary. The positive-floor
values are substantially worse in these finite runs. Neither finding
excludes better gauges or other physical subspaces.

## 6. Replay scope and reproduction

The [general checker](../tools/check_balanced_particle_extensions.py)
certifies the complete natural schedule and independently assembles
every retained full-gauge matrix. For a saved full-gauge state, it
checks the recorded weighted orders, total physical energies, type
floors, sufficient Gibbs KKT conditions, and entropy quotient.

The [reflection-extension gauge replay](../evidence/balanced-particle-extension-gauge-replay.json)
passed 19,201 literal Gaussian-minor identities, 165 retained states,
and 330 independently assembled weighted profiles. The
[nonreflection gauge replay](../evidence/nonreflection-particle-extension-gauge-replay.json)
passed 9,288 literal minor identities, 81 retained states, and 162
weighted profiles. These real-height replays remain numerical; the
natural integer-height profiles have the separate exact frame proofs.

Source hashes, literal seeds, parameters, rejected-trial counts, and
optimizer termination messages are retained. A run can retain a valid
evaluated state even if its optimizer reaches the iteration limit.
No retained penalty from a failed oracle call is treated as a witness.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/search_balanced_particle_extensions.py --output /tmp/balanced-particle-extension-natural-search.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_balanced_particle_extensions.py --natural evidence/balanced-particle-extension-natural-search.json --output /tmp/balanced-particle-extension-natural-replay.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_six_state_reflection_benchmark.py --output /tmp/six-state-reflection-benchmark-certificate.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_balanced_particle_extensions.py --gauges evidence/balanced-particle-extension-gauge-search.json --output /tmp/balanced-particle-extension-gauge-replay.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_balanced_particle_extensions.py --gauges evidence/nonreflection-particle-extension-gauge-search.json --output /tmp/nonreflection-particle-extension-gauge-replay.json
```

The two optimization drivers are
[the top-extension driver](../tools/search_balanced_extension_gauges.py)
and [the nonreflection driver](../tools/search_nonreflection_extension_gauges.py).
Their complete invocation arguments are embedded in their saved records.
