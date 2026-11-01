# A complete coupled-quadratic F4 comparison below the protected core

## Family and result

On six physical modes, let D_i delete occupation bit i. For nonzero
`a,b,c` in F4, impose the two differential equations

```math
T_8=D_0D_3+aD_1D_2,\qquad
T_{12}=D_0D_5+bD_1D_4+cD_2D_3.                      \tag{1}
```

Their subscripts are their physical heights. These operators share modes
and have different heights. They are not the disjoint, common-height
paired system of the matroid Hilbert formula.

All 27 coefficient triples in (F4 nonzero) cubed were evaluated. Their
coordinate-retraction codes have 28 letters and the same physical height
multiset. Four also have a nondegenerate kernel pairing and therefore an
orthogonal retraction. Complete natural profiles were obtained for all
31 maps.

The binary coordinate control `(a,b,c)=(1,1,1)` has the smallest product
of the two channel counting polynomials at every `0<r<1` among these
31 maps. The comparison is strict for every other map. Its numerical
balanced pressure is approximately 0.7182108251, so it does not improve
the existing crossing bound or factoring exponent.

The height-eight equation changes states in the previously protected
low core. Thus the comparison is separate from the
[protected-core theorem](coherent-protected-height25.md). Its scope is
the complete finite coefficient family (1), over F4, on these six modes.
It does not exclude larger fields, other equations or other physical maps.

## 1. Actual physical maps

Every operator in (1) is homogeneous in height and occupation degree.
All bit deletions commute with them. Their common kernel on the full
six-mode occupation space is therefore deletion closed.

The [existing differential constructor](f4-differential-pair-codes.md)
row-reduces the kernel separately in each homogeneous cell. Its leading
words give a deletion-closed coordinate retraction. The profile driver
checks the actual equations, split identity and both deletion closures.

Every case has occupation dimensions `(1,6,13,8)`, summing to 28. Four
coefficient triples admit an orthogonal inverse-Gram retraction:

```text
(2,1,2), (2,1,3), (3,1,2), (3,1,3).
```

Here 2 and 3 encode omega and omega plus one. These are exact finite
kernel and Gram computations; no general formula for overlapping
differential equations is inferred from them.

The [complete run](../evidence/tensor-coupled-quadratic-profiles.json)
computes both cuts separately for all 27 coordinate maps. For the four
orthogonal maps it computes the zero cut and uses the proved
equal-profile identity for the infinity cut. The record marks those
four infinity profiles as derived. The run contains 45,472 directly
computed Smith orders and 3,136 derived infinity orders.

## 2. An exact all-scale certificate

Let C_0,C_infinity be the complete counting polynomials of the binary
coordinate control. For another map let Z_0,Z_infinity be its two
complete counting polynomials. The
[checker](../tools/check_tensor_coupled_quadratic_dominance.py) factors

```math
Z_0(r)Z_\infty(r)-C_0(r)C_\infty(r)
=r^u(1-r)^v Q(r).                                   \tag{2}
```

The difference is zero only for the control itself. For 28 other maps,
all coefficients of Q in the power basis are nonnegative and at least
one is positive. The remaining two cases, coordinate triples `(1,1,2)`
and `(1,1,3)`, have `u=48`, `v=1` and quotient degree 193. Their exact
expansions

```math
Q(r)=\sum_{j=0}^{193}b_j r^j(1-r)^{193-j}
```

have nonnegative integer b_j, with a positive coefficient. Both the
factorization and each basis conversion are checked coefficient by
coefficient. Thus (2) is strictly positive throughout `0<r<1` for
all thirty alternatives.

The physical height multisets are identical, so matched physical types
have the same entropy and mean-height costs. The strict channel-product
comparison proves domination at each finite natural pressure scale.
It does not rely on numerical minimization or one cut alone.

## 3. Replay

Run

```sh
python3 tools/check_tensor_coupled_quadratic_dominance.py
```

This verifies the integer identities in the
[comparison receipt](../evidence/tensor-coupled-quadratic-dominance.json)
against the full saved profiles. It binds the source artifact by SHA-256
and explicitly records that no Smith orders were reconstructed by this
quick replay.

The [search driver](../tools/tensor_coupled_quadratic_profiles.py)
reconstructs the entire finite list from the explicit operators. It
records source hashes, in-progress cases and completion. The profile
engine and the differential-space construction remain credited to the
preceding notes. Substantial OpenAI assistant assistance was used for
the experiment, exact comparison and text. No new Lean claim is made.
