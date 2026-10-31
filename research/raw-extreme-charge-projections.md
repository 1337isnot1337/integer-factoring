# Charge projections on the raw Gaussian's extreme modes

## Result and scope

The original binary Pascal Gaussian admits arbitrary nonlinear charge
projections on its two extreme physical modes, without first applying
the strict Pascal repair. The construction preserves its existing
typed endpoint. Both profiles reduce to ten active fermion modes and
explicit Gaussian spectators, using the proved full-core reduction.

This is a different family from the
[strict low-core projections](strict-core-charge-projection.md): the
underlying tensor and the selected physical modes both differ. No
improved pressure or factoring exponent is asserted here.

## 1. A local mask preserves a pointwise typed endpoint

Let \(G\) be the original binary Gaussian on a full power-of-two
interval of length \(M\). At the local endpoint, its \(Y,U\)
block is an identity, and its \(V,X\) block is the triangular
Pascal matrix. Matching physical mode totals in a tensor power force
all triangular increases to vanish. Thus the retained global endpoint
is the ordinary crossing.

Suppose a local scalar mask \(m(x,y,u,v)\) satisfies

~~~math
m(x,y,y,x)=1\qquad\text{for every }x,y.
\tag{1}
~~~

Put \(T=mG\), entry by entry. In an ordinary tensor power,

~~~math
T^{\otimes N}(x,y,u,v)
=\left(\prod_{\ell=1}^{N}m(x_\ell,y_\ell,u_\ell,v_\ell)\right)
 G^{\otimes N}(x,y,u,v).
\tag{2}
~~~

An off-diagonal typed endpoint coefficient was already zero, so (2)
cannot revive it. On a diagonal endpoint coefficient, every factor
in the product is one by (1). The typed crossing endpoint is therefore
unchanged. This argument uses the pointwise tensor-power endpoint;
it does not assume that the untyped local Gaussian is a strict swap.

For a charge \(f:\{0,1\}^2\to\mathbb Z/3\), evaluate \(f\)
on modes \(0,M-1\) and take

~~~math
m(x,y,u,v)=
\mathbf1_{f(x)+f(u)=f(y)+f(v)}.
\tag{3}
~~~

Condition (1) holds. The mask preserves polynomiality, does not
increase the maximum entry degree, and keeps the full physical
alphabet. The same three-term \(\mathbb F_4\) phase sum as in
the strict construction evaluates (3), and its coefficients descend
exactly to \(\mathbb F_2\).

## 2. Ten active modes and all Gaussian spectators

Every physical extreme-state charge is a polynomial in the two
corresponding number operators. The
[full-core support theorem](rectangular-full-core-factorization.md)
already transports all creation and annihilation operators on the
extreme modes of all four physical legs into ten active Gaussian
directions at each place. Every transformation and its Gaussian
companion is regular. The theorem applies to any polynomial in those
operators, including the phase projection (3).

Thus, for \(M\ge8\), both cut matrices factor into one active
matrix of order \(1024\), tensored with every remaining Gaussian
mode. For \(M=2,4\), keep all \(2M\) one-particle modes instead.
All full physical alphabets are retained here; the rectangular code
restrictions of the earlier application are not imposed.

The generic active rank is its full matrix order. At the endpoint,
the untyped cut is triangular after a fixed permutation of coordinates,
with diagonal entries one. Masking entries preserves that triangular
support and every diagonal entry. This supplies the rank promise
independently of any truncated-series computation.

## 3. One Gaussian matrix supplies every charge projection

Use the regular physical-core coordinates from the same theorem.
Their basis maps intertwine every physical core creation and
annihilation operator with its coordinate counterpart. In particular,
every diagonal core function is diagonal in those coordinates.

Consequently, first construct the unprojected Gaussian matrix in
those coordinates. Every charge projection is then its exact entry
mask, using the core-state labels on its rows and columns. There is
no need to rebuild a separate fermionic program for each charge
function. The unit vacuum minors certify the full coordinate bases,
and all Gaussian spectators remain explicit.

Arbitrary diagonal physical core gauges are ordinary row and column
weights in these bases. The two bulk parameters from the strict
core-gauge implementation remain valid, because the raw full-core
changes preserve the same total infinity occupation and zero charge
difference. The physical core energies now use actual mode indices
\(0,M-1\), rather than the low pair \(0,1\).

The implementation selects this family with
`CoreChargePressure(..., family='raw_extremes')` in
[`strict_core_charge_gauges.py`](../tools/strict_core_charge_gauges.py).
The data adapter reuses the existing full-core operator images and
records the original physical indices explicitly. The
[fourteen-class screen](../tools/search_raw_extreme_charge_projections.py)
computes exact active profiles and optimizes the scale and bulk
gauges, including all spectators. Its floating roots guide searches;
they are not strict pressure certificates.

The [raw-family checker](../tools/check_raw_extreme_charge_projections.py)
verifies matching-mode typed endpoints and compares complete weighted
cut profiles against the original physical tensor, with independent
quarter-integer core and bulk gauges. Its
[evidence](../evidence/raw-extreme-charge-projections.json) records the
endpoint and full-matrix comparisons separately from the exploratory
long-interval pressure screen.
