# An affine reflection path with eight additional singular places

## Exact construction and scope

Two polynomial paths in the coupled physical-frame orbit of the
four-letter reflection crossing have an affine combination with eight
nonzero finite singular points. They form four quadratic closed places
over `F4`. Their complete profiles, both cut determinants, and the zero
and infinity profiles are checked exactly below.

This supplies a concrete input for mixed-place reconstruction. It does
not establish an improved crossing pressure or a factoring exponent.
The additional roots have to be used with all physical type, gauge,
degree, and reconstruction costs included.

Let `T(w)` be the interval-four reflection crossing, with alphabet four
and degree 24. Work over

```math
\mathbb F_4=\{0,1,\omega,\omega^2\},\qquad
\omega^2+\omega+1=0.
```

Use the [literal Bell normalization](coupled-rational-physical-frames.md)
before applying a frame, and restore the original physical coordinates
afterwards. For a square-zero physical matrix `G`, let `T[G,k]` denote
the complete action of `I+w^-k G` on `x,u` and its inverse transpose
on `y,v`.

The two specific paths

```math
T_A=T[\omega E_{0,3},4],\qquad
T_B=T[\omega^2 E_{3,0},8]
```

are polynomial. This is an exact finite coefficient check; poles four
and eight at interval four are outside the range of the uniform
bounded-pole theorem. Both paths have the original leading coefficient
and satisfy the strict crossing endpoint.

Define

```math
S(w)=\omega^2 T_A(w)+\omega T_B(w).                \tag{1}
```

Since `omega^2+omega=1`, (1) has the same endpoint and leading
coefficient. Its actual degree is 24. In the frozen coordinate-frame
catalogue this is the affine word

```text
[[19,3],[89,2]].
```

Here a word means `T+sum c_i(T_i-T)`. The two copies of the base
cancel in this particular word. The full tensor coefficient hash is

```text
dd9490bc487c838a1f29363522a9f56ab06f83becdbc4124a85265198da2ef9c.
```

## 1. Exact cut determinants

Write `R` for the cut `(x,y)|(u,v)` and `C` for `(x,u)|(y,v)`.
Put

```math
\begin{aligned}
f_4(w)&=w^2+w+\omega,\quad &f_6(w)&=w^2+w+\omega^2,\\
f_9(w)&=w^2+\omega w+\omega,\quad
&f_{12}(w)&=w^2+\omega^2w+\omega^2.
\end{aligned}
```

These are four distinct irreducible quadratics over `F4`. The complete
polynomial determinants are

```math
\boxed{\det R(S)=w^{160}f_6(w)^{24}f_9(w)^{24}},   \tag{2}
```

and

```math
\boxed{\det C(S)=f_4(w)^{24}f_{12}(w)^{24}}.        \tag{3}
```

The checker derives both determinants by fraction-free elimination on
the literal 16-by-16 polynomial matrices and compares every coefficient
with (2) and (3). It does not infer these identities only from sampled
finite-field evaluations.

Their degrees are 256 and 96. Equations (2) and (3) are coprime, so the
tensor has no nonconstant common scalar factor. They also prove that
the list below contains every nonzero finite singular place of either
cut, over any algebraic extension.

## 2. Complete profiles at the quadratic places

Realize `F16` by adjoining `u` with `u^2+u+omega=0`. The implementation
encodes `a+b*u` as `a+4*b`, with `0,1,2,3` denoting the four elements
of `F4`. The `F4` Frobenius conjugates are `z` and `z^4`.

At a root of the indicated quadratic, the complete local Smith orders
are as follows. An exponent in parentheses denotes multiplicity.

| Quadratic | Encoded conjugate roots | `R` profile | `C` profile |
|---|---|---|---|
| `f4` | 4, 5 | `0 (16)` | `0 (11), 4 (4), 8 (1)` |
| `f6` | 6, 7 | `0 (11), 4 (4), 8 (1)` | `0 (16)` |
| `f9` | 9, 11 | `0 (11), 4 (4), 8 (1)` | `0 (16)` |
| `f12` | 12, 15 | `0 (16)` | `0 (11), 4 (4), 8 (1)` |

Thus each singular cut has rank 11 at that point and determinant order
24. The undeformed source has full rank 16 on both cuts at all eight
points.

There are two independent controls for every displayed nontrivial
profile. Exact local valuation elimination returns all 16 orders.
Local block-Toeplitz ranks at lengths `1,4,5,8,9` are respectively

```text
11, 44, 59, 104, 120.
```

These equal `sum_i max(length-nu_i,0)` for the displayed profile.
The complete generic cut rank is 16 because the endpoint is a
permutation matrix. In particular, the last two controls show that
no omitted channel has a larger order.

## 3. The remaining two places

At zero, the complete profiles are

```math
R:\quad 0^{[2]},8^{[8]},16^{[6]},\qquad
C:\quad 0^{[16]}.                               \tag{4}
```

At infinity use the homogeneous local tensor `tau^24 S(1/tau)`.
Its complete profiles are

```text
R: 0,2,2,4,6,6,8,8,8,8,10,10,12,14,14,16
C: 0,8,10,10,14,14,18,18,20,22,22,24,24,28,28,28.
```

Their sums are 128 and 288. Together with (2) and (3), each total
homogeneous determinant order is `16*24=384`, as required.

The exact polynomial determinants and profiles are stored in
[the replay evidence](../evidence/affine-quadratic-place-witness.json).
Run the [standalone checker](../tools/check_affine_quadratic_witness.py)
to reconstruct the entire path from the two explicit physical frames.

## 4. Full-cost use of the new places

The [mixed-place theorem](mixed-place-hermite-crossings.md) permits
disjoint row and column place sets. For this example a natural set is

```math
\mathcal A=\{0,6,7,9,11\},\qquad
\mathcal B=\{\infty,4,5,12,15\}.                 \tag{5}
```

Every quadratic contributes two distinct labelled places in the one
fixed field `F16`. Conjugate profiles agree under the same physical
gauge, so a group can share its gauge parameters while its capacity
is counted twice. The two roots must not be counted as a single place.

The [explicit polynomial profile adapter](../tools/affine_place_profile_oracle.py)
accepts the complete tensor and recomputes every weighted cut. It uses
direct `F4` matrices at rational places and the checked `F16/F4`
restriction at the quadratic places. Restriction duplicates every
Smith order twice; the logarithmic channel sum removes that factor of
two. The exact rational-height path remains available independently.

The [adapter checker](../tools/check_affine_place_profiles.py) verifies
30 complete rational-height profiles, 480 channels, 768 endpoint
entries, 51 gradient coordinates, and three scalar-content
normalization controls. Its evidence is
[here](../evidence/affine-place-profile-controls.json).

The search implementations
[for affine words](../tools/search_affine_mixed_places.py) and
[for arbitrary explicit seeds](../tools/search_explicit_mixed_places.py)
retain the full degree charge and one common pair of physical types
at every place. They permit independent physical height gauges and
reconstruction allocations. Their floating outputs are proposals for
exact replay, not certificates.

The auxiliary [rational-place screen](../tools/screen_affine_finite_places.py)
and [quadratic-place screen](../tools/screen_affine_extension_places.py)
only filter candidates by exact finite matrix ranks. Such a rank drop
does not by itself bound a typed jet, a reconstruction piece, or the
resulting crossing pressure.
