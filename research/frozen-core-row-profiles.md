# A fixed-size profile compiler for polynomial frozen rows

## Result and scope

The eleven-parameter family in
[polynomial-frozen-core-crossings.md](polynomial-frozen-core-crossings.md)
now has a reusable profile compiler. For a power-of-two interval of
length `M`, each of its two local matrices reduces to a `64 x 64` active
matrix and the known Gaussian spectators. At `M=4` the active matrix
is `16 x 16`. The matrix size stays fixed as the interval grows.

The compiler constructs twelve source matrices once. Every affine row
`g` with `g[0]+g[10]=1` is then an exact linear combination over `F4`.
This permits searches that impose coefficient cancellations before
computing any Smith profiles. It also keeps the actual non-Gaussian
coefficients, including cancellations between source rows.

The finite search recorded below found no pressure improvement. It is
a search result, not an obstruction theorem for this family. No new
factoring exponent follows from it.

## 1. The active matrices

The compiler uses the regular full-core basis from
[coordinate-core-crossings.md](coordinate-core-crossings.md). It fixes
the input and output physical cores to the low-low state, then records
each of the twelve height-dressed source output covectors separately.
All operations act inside the existing bounded active subspace; the
orthogonal Gaussian spectator factors are unchanged.

In the full-core basis the two low-low coordinates are each labelled
by core state `5`. Restricting to this coordinate leaves active
multiplicity `16` at `M=4` and `64` at `M>=8`. The twelve resulting
matrices are denoted `A_r(z)`, with all Laurent entry shifts retained.
For a row `g`, the active matrix is exactly

\[
A_g(z)=\sum_{r=0}^{11}g_r A_r(z).
\]

The source components are aligned to a common Laurent shift **for
each entry** before addition. A leading term can therefore cancel
between rows without being replaced by a generic valuation estimate.
The finite field representation uses two binary polynomial planes.

The compiler also removes common binary polynomial factors with
constant term one. It may divide a factor shared by every component,
by every component in one row, or by every component in one column.
These divisions commute with every choice of `g`. They are regular
invertible diagonal operations and preserve every weighted Smith
profile, including profiles with nonintegral real row and column
weights.

These unit divisions require complete polynomial coefficient arrays.
The implementation rejects truncated arrays before normalization. It
records every removed factor, so the original twelve arrays can be
reconstructed exactly. The retained-precision convention is the one
in [retained-profile-precision.md](retained-profile-precision.md).

For example, at `M=64, beta=1`, the largest aligned coefficient lengths
drop from `2596` to `725` bits at infinity and from `1582` to `496`
bits at zero. These are measured sizes for this compilation, not
claimed formulas for all intervals.

## 2. Independent exact controls

`tools/check_frozen_core_row_profiles.py` checks three different parts
of the implementation:

- It restores `110592` component entries exactly after the common,
  row, and column unit divisions.
- It compares `36` weighted profiles at `M=4` with direct profiles of
  the original physical coefficient tensor. The comparisons use both
  `beta=0` and `beta=1`, the Gaussian row, the inverse-repair row, a
  dense affine row, both cuts, and three pairs of charge weights.
- It verifies the complete `M=8` Gaussian profiles, including every
  spectator, at three balanced charge gauges against the independent
  interval formula of
  [rectangular-low-background.md](rectangular-low-background.md).

The saved exact controls are in
`evidence/frozen-core-row-profiles.json`. The separate polynomial-family
checker verifies the endpoint condition and the equivalence of all
nonzero `beta` parameters; hence `beta=0,1` cover the two free-row
families.

## 3. A finite cancellation search

The search uses the critical exponent `rho=1/sqrt(2)` and the scale

\[
\tau=\sqrt2\log M,\qquad \ell=\tau/M.
\]

For each pair `(a,b)` in
`{-2,-1,-1/2,0,1/2,1,2}^2`, its balanced charge weights are

\[
\theta=(a\tau,-a\tau,b\tau,-b\tau,\ell).
\]

At each cut, let `m_G` be the least weighted entry order of the
Gaussian row, before multiplication by `ell`. The search imposes the
exact `F4` linear equations that cancel all component coefficients
whose weighted order is below `m_G-M/4`. It adds the affine endpoint
equation `g[0]+g[10]=1`, then solves the resulting system over `F4`.
Every generated row is evaluated using its complete active profiles
and all Gaussian spectators.

The surviving row space is enumerated when its size is at most the
candidate cap. Larger spaces use deterministic basis-direction
candidates and pseudorandom candidates with seed `310827`. The search
also records the dimension, total possible rows, number evaluated,
and best row separately at every gauge pair.

The aggregate record is
`evidence/frozen-core-critical-pressure.json`. Every recorded run uses
complete polynomial matrices and has zero profile-precision failures.
The Gaussian row attains the minimum found in all of them. Its
critical pressures at the displayed scale are approximately:

| `M` | Gaussian critical pressure |
|---:|---:|
| 8 | 0.224610593648 |
| 16 | 0.094771591219 |
| 32 | 0.043705824907 |
| 64 | 0.021942636016 |

The search covers both parameter orbits at `M=32,64`, and the
nonzero orbit at `M=8,16`. In the nonzero orbit the surviving affine
spaces have dimensions `2`, `3`, or `4`. At `M=16,32,64`, the cap `256`
enumerates every surviving row at every selected gauge pair.

The coefficient threshold is a search heuristic. It can exclude a
row whose loss at one cut is offset by an improvement at the other.
The finite gauge grid does not cover arbitrary balanced gauges, and
the chosen scale does not cover every positive scale. Finite sampled
spaces are labelled explicitly in the evidence. These limitations
prevent interpreting the result as a lower bound for the full family.

## 4. Reproduction

Run the exact coefficient and profile controls with:

```sh
python3 tools/check_frozen_core_row_profiles.py
```

For one complete search run, use:

```sh
python3 tools/search_frozen_core_critical_pressure.py \
  --M 32 --beta 1 --cap 256 \
  --output evidence/frozen-core-M32-replay.json
```

The portable command rebuilds the source matrices from the construction;
it does not depend on scratch caches. `--tau`, `--grid`, `--margin`, and
`--cap` expose the finite search choices. The aggregate evidence lists
the exact choices of each saved run. Field operations and coefficient
cancellations are exact; the displayed pressure values use floating
point logarithms and are exploratory numerical evidence.
