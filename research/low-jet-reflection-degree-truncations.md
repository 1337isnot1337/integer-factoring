# Reflection degree truncations that preserve the low jets

## Uniform construction and scope

For every power-of-two reflection interval `M>=4`, there is a strict
crossing of degree `M^2` that agrees with the original reflection tensor
in every coefficient below `M^2`. Its paired traces are unchanged.
The original tensor has degree `3M^2/2`.

More generally, every positive even degree cap `d` gives an explicit
strict crossing that retains all original coefficients below `d`.
For `d<M^2` its paired traces change, so these paths leave the original
trace-preserving affine family. All the paths have scalar content one.

These are ordinary finite coefficient identities. They do not imply
an improved weighted profile, crossing pressure, or factoring exponent.
In particular, preserving low jets at zero gives no free bound on the
new infinity profile. Both complete cuts still require evaluation.

## 1. One endpoint correction

Let `T_M(w)` be the [reflection-pair compressed crossing](reflection-dual-rail-compression.md),
with physical alphabet `q=2^(M/2)` and endpoint

```math
C_{xyuv}=[x=v][y=u].
```

The [literal Laurent Gaussian formula](../tools/reflection_dual_rail_gaussian.py)
has only even powers of `w`: all its exponents, the common monomial,
and the restored physical shifts are even. Its actual entries are
polynomials, so write

```math
T_M(w)=\sum_{r\ge0} T_r w^r,\qquad T_r=0\text{ for odd }r.
```

The original degree is uniformly `3M^2/2`. Indeed, its nonnegative
natural infinity heights have total `M^2`; regularity of the gauged
reversed kernel bounds every matching exponent by that total. Add
the scalar exponent `M^2/2`. Equality holds with empty `x,y` and full
`u,v`: at the largest degree only the within-leg reflection matchings
contribute, each with coefficient one. The cross `u,v` edges have
strictly positive valuation in that infinity gauge, so they cannot
cancel this leading term.

For any positive even `d`, set

```math
P_d(w)=\sum_{r<d}T_r w^r,
\qquad
\boxed{S_d(w)=P_d(w)+w^d\bigl(C+P_d(1)\bigr).}             \tag{1}
```

The coefficient field has characteristic two. Therefore

```math
S_d(1)=P_d(1)+C+P_d(1)=C.
```

Equation (1) has degree at most `d` and agrees with `T_M` in all
coefficients `0,...,d-1`. It remains even, so its derivative is zero
at every point. The construction uses entrywise truncation and one
coefficient evaluation; it does not require finding a physical frame
or solving a span-membership problem.

## 2. Exact traces and scalar content

The [uniform paired-trace identity](reflection-paired-trace-invariant.md)
is

```math
\sum_x T_{xyux}(w)=(1+w^{M^2})\delta_{yu},\qquad
\sum_y T_{xyyv}(w)=(1+w^{M^2})\delta_{xv}.                  \tag{2}
```

Each corresponding contraction of `C` is zero, since `q` is even.
Contracting (1) therefore gives

```math
\sum_x (S_d)_{xyux}(w)=(1+w^{\min(d,M^2)})\delta_{yu},
\qquad
\sum_y (S_d)_{xyyv}(w)=(1+w^{\min(d,M^2)})\delta_{xv}.      \tag{3}
```

For `d<=M^2`, the truncated trace of `P_d` is the identity, and the
endpoint correction contributes `w^d` times the identity. For
`d>M^2`, the full trace in (2) is already retained and its value at
one is zero, so the correction has zero contraction.

At `d=M^2`, equation (3) forces an actual coefficient of degree
`M^2`. Thus the constructed path has exactly that degree and attains
the trace lower bound. At every even `d<M^2`, the same argument proves
that the actual degree is exactly `d`, with a different trace invariant.
For larger caps the actual degree is read from the constructed tensor.

There is also a direct scalar-content proof for every cap, including
caps for which `1+w^d` has roots other than one. The empty Gaussian
coefficient is

```math
(T_M)_{0000}=w^{M^2/2}.
```

Equation (1) consequently gives

```math
(S_d)_{0000}=w^{\min(d,M^2/2)}.                            \tag{4}
```

Any common scalar factor must therefore be a power of `w`. The
constant tensor is unchanged, and (3) has nonzero constant trace, so
some entry has nonzero constant coefficient. No positive power of `w`
divides every entry. The scalar content is exactly one.

## 3. General double-endpoint truncation

The even-source shortcut has a version for any polynomial source with
endpoint `C` and derivative zero at one. Let `d>=2` be even, and put

```math
P(w)=\sum_{r<d-1}T_r w^r,
\quad A=P'(1),
\quad B=C+P(1)+A.
```

Then

```math
S(w)=P(w)+w^{d-1}A+w^dB                                   \tag{5}
```

has degree at most `d`, endpoint `C`, derivative zero at one, and all
source coefficients through degree `d-2`. In characteristic two,
`d-1` is odd and `d` is even, so the derivative identity follows
directly. For an even source, `A=0` and (5) is exactly (1).

An independent formulation explains the affine connection. Choose

```math
U=E_{00}\otimes I+I\otimes E_{00},\qquad
T_{\min}=C+(1+w^r)U,
```

in the paired coordinates `(x,v)|(y,u)`, with even `r<=d`. Both
`T-T_min` and its derivative vanish at one, so it is divisible by
`(w+1)^2=1+w^2`. Retain quotient coefficients of degree at most `d-2`,
multiply by `1+w^2`, and add `T_min`. This recovers (5). For the
trace-preserving reflection construction, take `r=M^2` and `d>=M^2`.

The direct proof of (1)--(5) works at every reflection interval and
does not assume an affine-span saturation theorem. A separate
[explicit affine witness](degree16-affine-reflection-crossing.md)
already attains the degree-sixteen floor at alphabet four.

## 4. Exact construction and checks

The [constructor](../tools/truncate_reflection_crossings.py) exposes
`truncate_even_crossing(seed, bound=None)` and
`truncate_double_endpoint(seed, bound=None)`. The default cap is
`M^2`. It verifies the endpoint and retained coefficients and records
the actual degree and source tensor hash. Stale Gaussian scalar
metadata is removed, since the result need not remain Gaussian.

Run

```sh
python3 tools/check_truncated_reflection_crossings.py \
  --output evidence/truncated-reflection-crossings.json
```

The [checker](../tools/check_truncated_reflection_crossings.py) reconstructs
all twelve four-letter paths with caps `2,4,...,24`, together with
nine sixteen-letter paths with caps `16,32,48,64,68,72,80,88,96`.
It independently checks the full entries using exponent sets and the
quotient construction, verifies every paired trace, and checks (4)
and a nonzero constant entry. The [evidence](../evidence/truncated-reflection-crossings.json)
stores complete tensor hashes and the exact coefficient ranges retained.

The [sixteen-letter numerical screen](../tools/search_truncated_reflection16.py)
uses complete weighted cut profiles at fixed rational multiples of
the natural physical heights. It optimizes the shared type, scale,
and supporting angle, then runs bounded full-height refinements.
It charges each new path's actual coefficients. Its floating results
are numerical candidates, separate from the exact degree theorem.

The [completed four-letter screen](../evidence/truncated-reflection4-pressure-search.json)
checks all twelve caps in both cut orientations. Every case completes
five scalar-height starts, two full-gauge/common-type refinements, and
three asymmetric refinements. The best balanced value remains the original
degree-24 value, approximately `0.731080936`. The degree-16 and degree-22
paths reach approximately `0.760242` and `0.742803`, respectively, after
full-gauge refinement. No numerical threshold candidate is found.
The [replay checker](../tools/check_truncated_reflection4_search.py)
reconstructs all 24 tensors, binds the source hashes, and freshly evaluates
all 48 saved balanced/asymmetric minima. Its [record](../evidence/truncated-reflection4-pressure-replay.json)
distinguishes these numerical replays from the exact coefficient theorem.
