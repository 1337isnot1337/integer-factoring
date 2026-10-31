# Polynomial changes of basis on complete physical alphabets

## Construction and scope

Polynomial changes of basis give an effective family of strict crossing
paths. The matrices act on complete physical letters; they need not be
exterior powers of one-particle matrices or preserve a parity grading.
This allows full-letter maps beyond the
[Gaussian physical frames](growing-gaussian-physical-frames.md).

The construction below is exact. No member is claimed to improve the
established factoring exponent `1/2`. Its complete polynomial entries,
actual degree, physical types, and both local rank costs must still be
included in any pressure certificate.

## 1. Paired inverse-transpose frames preserve the crossing

Let `K(z;x,y,u,v)` be a polynomial tensor over a fixed field whose value
at `z=1` is `[x=v][y=u]`. Suppose `A(z)` and `B(z)` are polynomial
matrices with polynomial inverses, acting on the `x,v` and `y,u`
alphabets respectively. Apply the four leg maps

```math
A_x(z),\qquad B_y(z),\qquad B_u(z)^{-T},\qquad A_v(z)^{-T}.
                                                               \tag{1}
```

All resulting coefficients are polynomial. At `z=1`, the paired maps
contract to `A(1)A(1)^{-1}=I` and `B(1)B(1)^{-1}=I`. Thus the
transformed tensor is again a strict crossing at one.

A particularly simple generator in characteristic two is

```math
A(z)=I+f(z)E_{ij},\qquad i\ne j,\qquad
A(z)^{-T}=I+f(z)E_{ji}.                              \tag{2}
```

Since `E_ij^2=0`, applying the same paired transvection twice restores
every polynomial coefficient. Any finite word of these generators
therefore supplies an explicit admissible frame. The polynomial `f`
may be a monomial, an endpoint-trivial polynomial such as `1+z^k`,
or any fixed finite coefficient list.

There is no assertion that this transformation reduces rank. It changes
the polynomial path before tensor powers and a common physical type
are chosen. A different frame may be tested as a different seed. One
cannot freely use unrelated physical frames for different reconstruction
jets while retaining the same typed tensor.

## 2. Exact entry computation and degree accounting

For one paired transvection, let `G_x` and `G_v^T` denote the two
single-leg nilpotent maps. The complete transformed tensor is

```math
K'=K+f(G_xK+G_v^T K)+f^2G_xG_v^T K.                \tag{3}
```

The implementation computes the equivalent two successive leg actions.
It uses exact binary polynomial multiplication and exact arithmetic in
`F4`; no leading-term approximation is made. It checks every entry of
the crossing endpoint after each paired transformation and measures the
largest surviving polynomial degree after all cancellations.

The degree can grow by up to `2 deg(f)` for one pair, and cancellations
can make the actual increase smaller. Search records retain the actual
degree. Raw infinity profiles use the complete reversed polynomial,
including every added coefficient. The code never assigns the old
degree to a changed tensor without recomputing it.

## 3. What a polynomial frame can and cannot change locally

In either cost cut, (1) multiplies the flattened matrix on the left and
right by polynomial matrices with polynomial inverses. At every finite
place these factors are invertible over the local ring. Therefore the
complete **unweighted** finite-place Smith profile is unchanged.

Physical-letter gauges need not commute with these matrices. The
weighted local profile can consequently change. At infinity a
nonconstant polynomial frame also need not be locally invertible with
zero valuation, so its full degree and local orders are recomputed.
These are the mechanisms explored by the bounded search. Invariance
of the unweighted finite-place profile is neither a proof that weighted
pressure is unchanged nor permission to omit the infinity cost.

Laurent frames can be more general: their separate legs can have poles
while the fully transformed tensor is polynomial after cancellation.
They require an exact final pole test and are handled by the separate
coupled-frame constructor, rather than by declaring each intermediate
leg polynomial. The polynomial implementation here also provides an
independent reference for positive-exponent special cases of that
constructor.

## 4. Exact controls and bounded searches

Run the exact checks with

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_physical_polynomial_frames.py \
  --output /tmp/physical-polynomial-frames-replay.json
```

The [checker](../tools/check_physical_polynomial_frames.py) verifies
1,024 complete endpoint entries, 1,024 restored coefficients after
inverting frames, and 1,024 coefficients of the constant crossing,
which must stay fixed for every spectral value. It also checks
commutation of the two disjoint physical pairs.

For four nontrivial frames, independent fraction-free elimination
verifies 128 unweighted finite-place channels and 256 complete weighted
channels in both cuts at zero and infinity. These controls include
negative rational heights and nonbinary field coefficients. The
[evidence](../evidence/physical-polynomial-frames.json) records the actual
degrees and every checked profile.

The [bounded search](../tools/search_physical_polynomial_frames.py)
varies words of full-letter polynomial transvections, all physical
heights, common paired types, and asymmetric cut weights. It retains
the best resolved complete-profile call in each timed trial, including
when an optimizer stops unsuccessfully. Missing native pivots are
rejected rather than interpreted as small ranks.

The separate [sixteen-letter basis search](../tools/search_reflection16_local_bases.py)
applies exact four-letter gates to selected pairs of physical modes,
then also tests full-letter neighbours. Its natural starting profile
is independently checked against all 512 channels of the existing
reflection formula. Both searches use floating optimization and need
independent rational profile and pressure certificates before any new
rate is asserted. A completed search record is separate from the
constructor identities proved here.
