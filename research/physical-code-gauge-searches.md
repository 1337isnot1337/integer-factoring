# Coherent physical codes with unrestricted letter gauges

## Result and scope

This checkpoint extends the physical split construction to coherent sums
over disjoint odd fibers. It also optimizes every physical-letter gauge
for thirty five-letter coordinate codes and twenty selected coherent codes.
The original seed is the reflection crossing with interval length eight
and sixteen physical letters.

The completed searches found no improvement over the full sixteen-letter
baseline `rho = 0.7166536597`. The best five-letter candidate remains the
coordinate code `{0,2,4,8,12}`, with `rho` about `0.7263479477` at its
displayed natural heights. These results do not improve the factoring
exponent or exclude other codes, fields, heights, or optimizers.

The construction is a special case of the
[physical split-projection identity](physical-split-crossings.md).
Its distinct feature is coherent mixing of three or more original
letters, followed by complete profile recomputation and optional
optimization of all endpoint heights.

## 1. Odd fibers give exact split maps

Let `G_1,...,G_q` be nonempty, pairwise disjoint subsets of the original
alphabet. Assume every `|G_i|` is odd. Over `F_2`, let `r_i` be the
indicator column of `G_i`, and let `R` have these columns. Then

```math
(R^T R)_{ij}=|G_i\cap G_j|\pmod2=\delta_{ij}.
```

Thus `L=R^T` satisfies `LR=I_q`. Apply the four physical maps
`(L,L,R^T,R^T)` to the original crossing. The result is a strict crossing
on `q` letters, and its polynomial degree does not increase. This works
over any coefficient field of characteristic two.

For a singleton fiber, this is coordinate restriction. A three-element
fiber replaces one physical letter by a coherent sum of three old
letters. Coefficients can cancel in both cuts. The code does not estimate
the new profiles by selecting old Smith orders; it constructs the full
new polynomial tensor and recomputes both profiles.

The original letter energies are sums of the mode energies
`(28,20,12,4)`. A fiber can combine unequal energies. Therefore a chosen
new height is charged as a height of the new seed. A minimum, mean, or
maximum of the old energies is only a search choice; none is treated as
an inherited profile identity.

## 2. All physical heights and the Gibbs root

A five-letter scalar seed has `8(5-1)+1=33` search variables. Each of
the two places has four independent physical-leg height functions.
The reference letter has height zero. The last variable is `log s`.
Heights in the implementation already contain the factor `log s`.

The existing [rectangular Smith oracle](rectangular-convex-pressure.md)
computes complete profiles with exact `F_4` polynomial arithmetic.
Its pivot comparisons, heights, and final pressures are floating point.
For every parameter vector, the paired physical types are optimized
through their Gibbs equation. The objective is the resulting `rho`,
not pressure at a fixed trial rate.

The gradient follows directly from the envelope equation. Write
`C_pi(x)` for the logarithmic profile cost plus the expected physical
heights, and write `H=H(pi_X)+H(pi_Y)`. At the optimizing types,

```math
C_\pi(x)=2\rho H,\qquad
\nabla_x\rho=\frac{\nabla_x C_\pi(x)}{2H}.
```

The driver supplies this gradient to SLSQP. Weighted heights are bounded
by `[-20,20]`, and `log s` lies in `[0.0001,2]`. The cost is nonsmooth
when an attaining minor changes. Each record retains the starting
candidate if the optimizer returns a worse value, and records its
termination status. A failed convergence flag is not an optimality
certificate. Small numerical variations between reruns are expected.

## 3. Completed finite searches

| Family | Exact seeds | Height or optimization choices | Best recorded behavior |
| --- | ---: | --- | --- |
| Five-letter coordinate codes | 30 | All 33 variables, three starts per seed, up to 250 SLSQP iterations | About `0.726348`; the low-energy code remains best |
| Five-letter coherent codes | 158 | Minimum, mean, and zero heights; all 33 variables for the best twenty, with two starts and up to 350 iterations | `0.7263479477`, attained by the coordinate control |
| Fourteen-letter coherent codes | 57 | One three-letter fiber and thirteen singletons; minimum, mean, and maximum heights | `0.7167908019` |

The coordinate list contains five explicit structured codes and twenty-five
distinct pseudorandom codes containing letter zero. The random seed is
`20261031`.

The coherent five-letter list contains the coordinate control, all
fifty-five ways to enlarge its vacuum fiber to `{0,a,b}` using two
omitted letters, two explicit partitions into triples, and one hundred
pseudorandom triple partitions. The driver sorts and removes duplicate
fiber lists. The resulting list has 158 seeds.

The fourteen-letter list uses consecutive triples in energy order,
every triple among the six highest-energy letters, and thirty random
triples. Removing duplicates leaves 57 seeds. Its best case merges
`{7,11,15}` and assigns the merged letter energy 64, the maximum of
those three energies. The best scale is approximately
`log s = 0.2192131259`.

Every one of the `18,750`, `98,750`, and `2,189,712` endpoint entries
in these three searches is checked directly. Profile calculations use
the complete polynomial seed at the displayed heights. The numerical
minimizations remain finite exploratory searches. In particular, the
necessary-profile targets used as some starting points do not assert
that any seed realizes those targets.

## 4. Independent exact controls

The exact checker compares the staged four-leg transformation with a
literal contraction over all original indices. Four small-code cases
also compare both complete weighted profiles with binary
restriction-of-scalars Smith calculations and independent block-Toeplitz
ranks. The cases include a full five-letter coherent code, not only
two-letter examples.

A further fourteen-letter control checks every tensor coefficient by
literal contraction. Across all controls, the checker verifies 39,219
endpoint entries, eight independent profiles, and 74 Toeplitz ranks.
It also rejects empty, even, overlapping, and out-of-range fibers.
These controls establish the exact construction used by the searches;
they do not certify numerical optima.

## 5. Reproduction

The exact controls require the same optional C++17 Smith kernel as the
physical split checks:

```sh
python3 tools/check_odd_fiber_crossings.py \
  --output evidence/odd-fiber-crossings.json
```

The searches additionally use NumPy and SciPy:

```sh
python3 tools/search_physical_code_gauges.py --suite coordinate \
  --output evidence/physical-code-gauges-coordinate.json
python3 tools/search_physical_code_gauges.py --suite coherent \
  --output evidence/physical-code-gauges-coherent.json
python3 tools/search_physical_code_gauges.py --suite compressed \
  --output evidence/physical-code-gauges-compressed.json
```

The output retains the exact fiber lists, every tested height, the
profiles, numerical parameters, optimization status, and source hashes.
The full-gauge records include the optimized physical types and complete
weighted orders. These data make each proposed parameter vector
reviewable without treating a numerical search as a proof.
