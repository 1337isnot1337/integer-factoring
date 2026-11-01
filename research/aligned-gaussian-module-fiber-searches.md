# Exact fibers in the complete aligned Gaussian module

The [aligned Gaussian complement module](aligned-gaussian-complement-module.md)
has 2,120 independent endpoint-zero directions of degree at most eighteen.
This checkpoint searches that complete displayed module around the
retained five-letter crossing. Coefficients may lie in `F_4`.

The best recorded rates remain `0.7263479477044683` at the natural
physical heights and `0.8335422390094076` at the two-zero heights.
No factoring exponent is improved. The exact new information is a
zero-dimensional reference-preserving fiber in each pair of fixed
regular frames, together with finite searches that explicitly trade
one row order against others. This does not exclude improvements
using different frames or other parameter values in the module.

## 1. The prescribed affine space

Let `B` denote the complete degree-eighteen reference tensor. Its
digest is

```
91e0ee1014b79e5e042053201717ccd95ad33a27ec7f1ee9a6c2144b10aa75e9
```

The source record contains 37 binary primitive slabs and all their
allowed monomial shifts and placements on the first physical label.
The directions are independent over `F_2`, hence also over `F_4`.
Every direction has value zero at one. A parameter vector therefore
gives a literal polynomial crossing

\[
 T=B+\sum_j c_j\Delta_j,\qquad c_j\in\mathbb F_4,
\]

with the strict endpoint and the full degree bound preserved.

For each local place, construct regular polynomial frames `L,R` of
the reference. A target vector asks every entry in row `i` of
`L T R` to vanish below a specified order. These conditions are
exact affine linear equations in the 2,120 coefficients. They are
sufficient conditions for a profile constraint, and are tied to
the displayed reference frames.

At both choices of physical heights, the constraints that preserve
all fifty reference row orders have the unique solution `T=B`.
Consequently none of the fifty individual row orders can be raised
while all other reference row orders are retained in those frames.
This is an exact statement about these linear systems. It does not
say that the reference is isolated among all tensors with the same
Smith profile, since another tensor may require different regular
frames.

## 2. Reusing primitive frame products

The first implementation projects the complete entry-coordinate
system onto each prescribed direction. The new
[source module fiber](../tools/source_module_crossing_fiber.py)
instead transforms each primitive slab and placement once. There
are `37*5*2=370` such frame products for the two local places.
Every shifted direction is then recovered by shifting the already
transformed polynomial coefficients.

At zero, multiplication by `z^k` raises the local order by `k`.
At infinity, it lowers the order by `k`, so the primitive product
must be computed through the requested row cap plus the largest
shift. The implementation retains precisely that extra window.

The reference frames themselves only need their original precision.
If a primitive admits a shift by `k`, its degree is at most `D-k`.
After the common infinity clearing shift and the minimum physical
heights are removed, its local matrix is divisible by `t^k`.
Thus omitted frame coefficients cannot enter the requested window
after the direction is shifted back down. This degree observation
is essential to the reuse argument.

The literal parameter reconstruction and module-membership check
are separate from the frame equations. The optimized constructor
requires binary reference frames and binary primitives; the affine
parameter coefficients remain arbitrary elements of `F_4`.

In the recorded runs, constructing the natural fiber and checking
its baseline rank took about 23 seconds. The earlier entry-projection
construction at the two-zero gauge took about 441 seconds. These
are observed timings on this workspace, not a controlled benchmark.

## 3. Completed finite searches

The fixed physical energy arrays are `(0,10,6,2,8)` for the natural
gauge and `(0,12,6,0,8)` for the two-zero gauge. Both use the same
array on the two boundary pairs.

The forced-row search first relaxes all fifty row floors to zero,
then asks one row to exceed its reference order. Its first fifty
nonreference trials cover every row. Later trials force small
sets of rows. Each feasible system is completed greedily in a
recorded randomized order, sampled over `F_4`, and evaluated with
both complete Smith profiles.

| Search | Completed schedules | Distinct evaluated tensors or target vectors | Best recorded natural rate | Best recorded two-zero rate |
|---|---:|---:|---:|---:|
| Full 2,120-direction row fiber, natural gauge | 300 | 266 target vectors | 0.7263479477044683 | — |
| Full 2,120-direction row fiber, two-zero gauge | 200 | 181 target vectors | — | 0.8335422390094076 |
| Degree-at-most-four submodule, singleton perturbations | 120 coefficient-orbit representatives | All 60 independent directions | 0.7263479477044683 | 0.8335422390094076 |
| Coefficient-priority normal forms in the full module | 300 priorities | 111 nonreference tensors | 0.7263479477044683 | 0.8335422390094076 |

The singleton record is complete for each displayed basis direction
and its nonzero `F_4` coefficients, with conjugate coefficients
identified by the exact Frobenius symmetry. It does not cover all
linear combinations of the sixty directions.

The coefficient-priority search uses the 424 independent shifted
slab directions separately at each of the five first labels. It
eliminates coefficients in different chosen priority orders, then
computes the complete profiles of the resulting literal tensors.
Its 300 priorities produced 111 distinct nonreference tensors.
The best retained nonreference rates in this experiment exceed one.
This normal-form screen is finite and does not optimize the entire
2,120-dimensional affine space.

All reported rates use numerical pressure optimization. The exact
part is the displayed tensor, module membership, endpoint, degree,
linear constraints, and complete valuation profiles.

## 4. Independent checks

The [primitive-shift checker](../tools/check_source_module_crossing_fiber.py)
compares every projected affine equation with the original constructor
on smaller complete examples. It also checks 103,908 literal `F_4`
frame coefficients, including every coefficient omitted as zero
from the sparse equation lists. Sixteen target systems agree exactly.
The actual 2,120-dimensional module has the recorded zero-dimensional
baseline fiber, and all fifty baseline-preserving raises are rejected.

The [finite-search replay](../tools/check_source_module_lattice_searches.py)
checks both full-module records, including their exact parameter
spaces, all fifty fully relaxed single-row feasibility tests, all
retained module memberships, and every displayed row target. It
certifies the complete profiles by regular polynomial frames on
the literal matrices, independently of the native profile oracle.

The singleton replay separately certifies 100 complete profiles and
31,250 endpoint entries. The coefficient-normal-form replay reconstructs
all 24 retained tensors from their parameter masks and certifies
52 complete profiles and 15,000 endpoint entries. These totals
include the relevant baseline controls.

## 5. Reproduction and records

```sh
python3 tools/check_source_module_crossing_fiber.py
python3 tools/check_source_module_lattice_searches.py
python3 tools/check_module_coefficient_normalforms.py evidence/aligned-gaussian-module-normalform-search.json
```

The principal search driver is
[search_source_module_lattices.py](../tools/search_source_module_lattices.py).
The completed full-module records are the
[natural-gauge search](../evidence/aligned-gaussian-module-lattice-natural-d18.json)
and [two-zero search](../evidence/aligned-gaussian-module-lattice-two-zero-d18.json),
with [literal replay controls](../evidence/aligned-gaussian-module-lattice-checks.json).
The implementation controls are in
[source-module-crossing-fiber-checks.json](../evidence/source-module-crossing-fiber-checks.json).

The smaller experiment records are
[degree-four directions](../evidence/aligned-gaussian-module-low-degree4-directions.json),
[singleton search](../evidence/aligned-gaussian-module-low-degree4-singletons.json),
[singleton replay](../evidence/aligned-gaussian-module-low-degree4-checks.json),
[coefficient normal forms](../evidence/aligned-gaussian-module-normalform-search.json),
and [normal-form replay](../evidence/aligned-gaussian-module-normalform-checks.json).

This research and its implementation used substantial OpenAI assistant support.
