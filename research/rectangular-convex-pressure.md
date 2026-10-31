# Rectangular paired-type pressure search and its verification boundary

The optional rectangular F4 search tools now expose the nested minors
selected by weighted Smith elimination, a fixed-type convex formulation,
and the derivative after exact Gibbs optimization of the physical types.
These changes make the search reproducible and independently testable.
They do not certify a new rate or a global optimum of the joint search.

## 1. Fixed-type convex formulation

Fix a rational local tensor with crossing-pair alphabet sizes `A,B`.
At each valuation place, let `g` be the physical letter heights multiplied
by `log(s)`, with `s>1`. For any nonzero `k` by `k` minor `J`, its
weighted valuation is an affine function

```math
d_J\log s+\sum_{i\in\mathrm{rows}(J)}g_i
             +\sum_{j\in\mathrm{columns}(J)}g_j.
```

The [weighted Smith oracle](weighted-smith-separation-oracle.md) proves
that the sum of the first `k` ordered Smith values is the minimum of
these affine functions. Introduce sorted auxiliary values `v_1,...,v_n`,
where `n=AB`, and impose every inequality

```math
\sum_{i=1}^k v_i\le d_J\log s+g(J).
```

At `k=n` impose the determinant equality. At both places, minimize
`log(sum_i exp(-v_i))` plus the restored physical height average for the
fixed paired types. This is a convex objective with affine constraints.
For fixed heights, the constraints say that the auxiliary vector is at
least as spread as the actual ordered profile, with the same total sum.
Convexity of the exponential therefore makes the actual profile an
optimizer. Thus the complete constraint formulation is exact.

The implementation starts with all nonzero-entry inequalities and the
determinant equalities. It then asks the native oracle for attaining
prefix minors and adds their inequalities. A minor obtained at one
height vector remains a valid inequality at every other vector. The
implementation retains row and column indices even when a profile is
assembled from disconnected support components.

Only a finite subset of inequalities is present during separation.
A relaxed optimizer can therefore be better than the actual profile;
it is never accepted without another full oracle evaluation. The code
also applies explicit height and series-window bounds. These bounds
describe the searched region, not a mathematical exclusion outside it.

## 2. Gibbs types and the derivative

Let `G_X,G_Y` be the sums of the four relevant physical heights for the
two crossing-pair alphabets, and let `B` be the product of the two
profile sums. The [paired-type theorem](paired-type-crossing-pressure.md)
gives the positive root `R=2*rho` of

```math
\log B=R\left(\log\sum_x e^{-G_X(x)/R}
                    +\log\sum_y e^{-G_Y(y)/R}\right).
```

The minimizing types are the corresponding Gibbs distributions. The
direct search computes those distributions at each oracle call. Its
gradient follows by the envelope theorem: differentiate the physical
height average and the two log profile sums at the current optimal
types, then divide by twice their total entropy.

The selected minimum-valuation minors give a valid branch derivative
away from ties. This does not make the joint type-and-height problem
convex. The direct method uses multiple starts; the convex method
alternates fixed-type optimization and Gibbs updates. Both retain the
best fully evaluated state and record failed profile queries.

## 3. Generic rank and series precision

The [generic matrix interface](../tools/native_matrix_smith_f4.py) also
accepts arbitrary entrywise Laurent shifts and a supplied generic rank.
It is useful for reduced defect matrices with zero rows and columns.
That rank must come from a separate proof. A finite series window cannot
prove that later pivots vanish. If a promised pivot is missing, the
routine raises an error and does not return an incomplete profile.

The arithmetic in F4 is exact. Floating height comparisons and numerical
optimization are exploratory. A claimed strict crossing improvement
still needs rational heights, sufficient precision, independent exact
profile replay, and an outward-rounded pressure inequality.

## 4. Reproduction and independent controls

The search entry points are
[direct paired-type search](../tools/search_rectangular_f4_pressure.py)
and [fixed-type convex separation](../tools/search_rectangular_f4_convex.py).
Both accept an explicit seed and write their parameters and search scope.

```sh
python tools/check_native_rectangular_smith_f4.py
python tools/check_native_matrix_smith_f4.py
python tools/check_rectangular_f4_optimizer.py
python tools/check_matrix_f4_entry_shifts.py
```

The rectangular controls independently recompute 148 selected minors
with fraction-free F4 determinants. The generic controls enumerate all
minors of small matrices with proved deficient ranks. A separate audit
uses genuinely entry-dependent Laurent shifts: 26 cases, 52 profiles,
832 exact minors, and 154 selected-prefix comparisons at both places.

The optimizer audit checks 168 derivatives against finite differences,
5,125 minor inequalities at held-out gauges, the independent Gibbs root,
and disconnected support assembly. Its largest derivative discrepancy
in the recorded run is below `3e-9`. These are implementation controls,
not asymptotic or global-optimality certificates.

Recorded evidence:
[rectangular profiles](../evidence/native-rectangular-smith-f4.json),
[generic profiles](../evidence/native-matrix-smith-f4.json),
[entrywise Laurent shifts](../evidence/matrix-f4-entry-shifts.json), and
[optimizer audit](../evidence/rectangular-f4-optimizer.json).
