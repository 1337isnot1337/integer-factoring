# Complete physical-letter gauges at every fused singular place

## Result and scope

The [complete singular-place catalogue](fused-all-place-profiles.md) now has
an exact weighted oracle. It accepts an independent height for every
physical symbol at all four legs, at any listed place and in either cut.
It uses every nonzero minor of the actual support blocks, so it does not
fix a row frame or restrict heights to a linear occupation charge.

The oracle also supplies an active-minor subgradient for exploratory
optimization. For fixed physical types, the local corrected log profile
sum is convex in the scale-weighted heights and the local valuation
scale. This convexity concerns one local cost; it does not establish
global optimality when physical types or place assignments also vary.

No improving factoring parameter is asserted by this oracle or by its
finite controls. Every retained numerical candidate needs its endpoint,
complete costs, and an exact strict certificate.

## 1. Reweight every actual determinantal prefix

For a cut matrix `M(t)`, give each row and column an additive physical
weight. A `k`-by-`k` minor `Delta_(I,J)` then acquires the same monomial
factor in every determinant term. If `a>0` multiplies the intrinsic
valuation and `h` is the full vector of physical-letter heights, its
weighted valuation is exactly

```math
a\,val(\Delta_{I,J})+\langle c_{I,J},h\rangle,              \tag{1}
```

where `c_(I,J)` counts each physical letter on the selected rows and
columns. Its entries sum to `k` on each of the four physical legs.

The `k`th determinantal prefix is the minimum of (1) over all nonzero
minors of that size. Successive differences give the complete Smith
orders. Diagonal physical weights preserve every support block, so
the block profiles can be merged without loss.

At a nonlinear closed point, `val` is the irreducible-factor valuation.
All its conjugate roots have the same value for every original minor,
including its physical row and column labels. At infinity, the base
valuation of a minor is `kD-deg Delta`, using the complete homogeneous
degree `D`. No infinity denominator correction is omitted.

The exact API accepts integer or rational weights. Float queries are
provided separately for numerical searches. A float query is not an
exact minor or pressure certificate.

## 2. Rational heights can be analyzed independently at each place

Rational local heights do not require a new global parameterization
of the polynomial seed. Let `L` clear their denominators at one place.
In that local calculation, set `t=tau^L` and multiply entries by the
integer powers `tau^(Lh)`. An original Hasse coefficient of index `j`
becomes the coefficient of index `Lj` before the physical shift.

After restriction to one pair of rational physical types, the total
height on every retained entry is the same scalar `mA`. The original
jet is therefore the coefficient of index `Lj+LmA` in the locally
gauged tensor. The type length can be chosen to clear all fixed
denominators. Apply the ordinary coefficient estimate at scale
`s^(1/L)`, and divide the reported local orders by `L`. The resulting
bound contains exactly `s^(j+mA)`.

Each selected place can use its own fixed denominator. This operation
only analyzes the same original reconstruction functionals. It does
not multiply the original global degree budget. The coefficient ranges
remain `O(m)`, with constants depending on the fixed seed and gauges.

## 3. Physical types are part of every local cost

Let `pi_X` be the common type on `x,v`, and `pi_Y` the common type
on `y,u`. For a height vector `h`, the scalar shift is

```math
A(h)=\langle pi_X,h_x+h_v\rangle
    +\langle pi_Y,h_y+h_u\rangle.                          \tag{2}
```

If the heights are already multiplied by the local logarithmic scale,
the corrected log profile sum is

```math
F(a,h)=\log\sum_i e^{-sigma_i(a,h)}+A(h).                  \tag{3}
```

The [mixed-place criterion](mixed-place-hermite-crossings.md) takes
the maximum corrected cost on each side and adds the complete global
degree cost. One pair of types is used at every place. Optimizing a
different type separately at each place would not give one crossing
on a common physical alphabet.

Two exact invariance controls check (2). Adding a constant `c` to all
heights on one leg adds `c` to every Smith order and to `A(h)`, leaving
(3) unchanged. Also, the literal fused support relation `x+y=u+v`
makes the linear charge `(h_x,h_y,h_u,h_v)=(x,y,-u,-v)` identically
zero on every nonzero entry and on every paired type.

## 4. Local convexity and an explicit subgradient

Write `delta_k(a,h)` for the `k`th determinantal prefix in one block,
with `delta_0=0`. It is a minimum of affine functions by (1), hence
concave. Its successive differences `sigma_k` are nondecreasing.

The entropy identity for log-sum-exp gives

```math
\log\sum_k e^{-sigma_k}
=\max_p\left[-\sum_k p_k sigma_k+H(p)\right].             \tag{4}
```

The maximizing probabilities are nonincreasing within each block.
Thus the maximization may be restricted to that convex set of
probabilities, also when several blocks are present. Summation by
parts gives

```math
-\sum_k p_k sigma_k
=\sum_k(p_{k+1}-p_k)delta_k,\qquad p_{n+1}=0.              \tag{5}
```

Every coefficient on the right is nonpositive. Each expression in
(4) is therefore convex in `(a,h)`, and taking its maximum preserves
convexity. Adding the fixed-type linear term (2) proves the assertion
for (3).

At the current point, choose any minimizing actual minor for each
prefix. Let its physical count vector and intrinsic valuation be
`c_k,v_k`. A valid active-minor subgradient of the log profile sum is

```math
-\sum_k p_k(c_k-c_{k-1}),\qquad
-\sum_k p_k(v_k-v_{k-1}),                                \tag{6}
```

where `p_k` is proportional to `exp(-sigma_k)` and the index restarts
at each block. At tied minima, these are supporting subgradients;
differentiability is not assumed. The gradient of (2) adds the physical
type probabilities on the appropriate four legs.

Jointly changing the type and the heights introduces the bilinear
term (2). The local convexity theorem therefore does not certify the
global searches over both variables.

## 5. Exact controls and search boundary

```sh
python3 tools/check_fused_weighted_places.py \
  --output evidence/fused-weighted-places.json
```

The controls cover all four seed alphabets and every closed place.
They compare all zero-height profiles with the complete catalogue,
check exact nonzero integer and rational gauges, verify constant-shift
and ramification identities, and test the conserved charge above.
For selected positive and negative gauges, independent polynomial
multiplication modulo powers of the irreducible factor checks the
weighted orders. Infinity uses independent block-Toeplitz ranks.
Supporting-plane checks and exact-profile comparisons validate the
separate float-gradient API at the recorded points.

The budgeted search driver is
[`search_fused_mixed_weighted.py`](../tools/search_fused_mixed_weighted.py).
It saves completed trials atomically, records optimizer stopping messages,
and marks a budget-interrupted trial separately. It tracks the actual
maximum local cost even when a smooth upper approximation is used for
optimization. It retains full orders, common paired types, allocations,
and every physical mean shift for its best point in each trial.

This note supplies a search and verification mechanism. A numerical
relative gain by itself is not a stronger factoring bound.
