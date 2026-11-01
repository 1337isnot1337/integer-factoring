# The complete four-letter paired-trace affine space

## Result

For the four-letter reflection tensor, independent constant paired
physical changes generate the entire polynomial affine space permitted
by its two paired traces, its strict crossing endpoint, and its vanishing
first derivative at that endpoint.

Write `F=F4`, let `T` be the original degree-24 reflection tensor, and put

```math
C_{xyuv}=\delta_{xv}\delta_{yu},\qquad
\chi(w)=1+w,\qquad
\mathfrak s=\{A\in\operatorname{Mat}_4(F):\operatorname{tr}A=0\}.
```

Use the physical pairs `(x,v)` and `(y,u)` to identify a four-leg tensor
with `Mat_4(F) tensor Mat_4(F)`. Let `S` range over the source tensor and
the [complete constant-companion catalogue](constant-companion-reflection-search.md).
Let `R` range over independent constant paired actions

```math
(A_x,B_y,B^{-T}_u,A^{-T}_v),\qquad A,B\in\operatorname{GL}_4(F).
```

Define `W` to be the linear span of the differences `R(S)-T`.
Thus `T` serves as the affine basepoint.
Then the exact equality is

```math
\boxed{
W=\chi(w)^2(\mathfrak s\otimes_F\mathfrak s)[w]_{\le22},
\qquad \dim_F W=15^2\cdot23=5175.
}                                                            \tag{1}
```

Consequently `T+W` consists of precisely the degree-at-most-24 tensors
`S` satisfying

```math
\begin{aligned}
S(1)&=C,& S'(1)&=0,\\
\sum_x S_{xyux}(w)&=(1+w^{16})\delta_{yu},&
\sum_y S_{xyyv}(w)&=(1+w^{16})\delta_{xv}.
\end{aligned}                                                 \tag{2}
```

This is a complete affine-space description for this finite four-letter
family. It does not bound arbitrary crossing tensors, and it is not a
crossing-pressure or factoring-exponent certificate.

## 1. The invariant upper bound

The [paired-trace theorem](reflection-paired-trace-invariant.md) proves
both contraction identities in (2) for the source. Independent paired
physical actions preserve those identities. Constant actions preserve
the degree bound.

Every rational coupled frame used in the companion catalogue is regular
at `w=1`, and its physical action `R(w)` fixes the crossing tensor for
all `w`. The source tensor is even, so `T'(1)=0`. Differentiating gives

```math
\left.\frac{d}{dw}(R(w)T(w))\right|_{w=1}
=R'(1)C+R(1)T'(1)=0,
```

because differentiating `R(w)C=C` gives `R'(1)C=0`. Constant companions
and subsequent independent constant paired actions retain this
property. Thus every affine difference has zero value and zero first
derivative at `w=1`. Entrywise, it is divisible by

```math
\chi(w)^2=1+w^2.
```

It remains to identify the simultaneous kernel of the two contractions.
The ordinary trace functional on `Mat_4(F)` is nonzero: it takes value
one on `E_00`. This remains true although `tr(I_4)=0` in characteristic
two. Therefore

```math
\operatorname{Mat}_4(F)=F E_{00}\oplus\mathfrak s,
\qquad \dim_F\mathfrak s=15.
```

Decompose both tensor factors using this direct sum. The conditions
`(tr tensor id)(V)=0` and `(id tensor tr)(V)=0` kill the three summands
containing an `E_00` factor. Their common kernel is exactly
`s tensor s`, of dimension 225.

For a polynomial difference, this statement applies to every
coefficient. Division by the nonzero scalar polynomial `chi^2`
preserves the two zero-trace identities. Degree at most 24 leaves a
quotient of degree at most 22. This proves

```math
W\subseteq\chi^2(\mathfrak s\otimes\mathfrak s)[w]_{\le22},
\qquad \dim_F W\le5175.                           \tag{3}
```

No numerical rank estimate is used in this upper bound.

## 2. Exact attainment of the upper bound

The [degree-sixteen checkpoint](degree16-affine-reflection-crossing.md)
preserves a 3,035-dimensional constructive prefix of `W`. Its complete
polynomial basis is the input to the
[saturation program](../tools/saturate_paired_trace_affine_space.py).
The original prefix and its provenance records remain unchanged.

For an off-diagonal matrix unit `E=E_ij`, the constant transvection
`I+cE` has inverse `I+cE`. On either physical pair its tensor action is

```math
R_{I+cE}=I+cL+c^2Q,\qquad
L=E_{\rm left}+E^T_{\rm right},\qquad
Q=E_{\rm left}E^T_{\rm right}.                    \tag{4}
```

For `c=1,omega`, the two coefficient rows `(c,c^2)` form an invertible
matrix over `F4`. Hence `L` and `Q` are linear combinations of the
two differences `R_(I+cE)-I`.

Apply these operators to the source tensor and to known differences.
Their images are again valid affine directions. There are 48 operators:
two physical pairs, twelve off-diagonal positions, and the two
operators in (4). The exact program stores both field-coefficient
bitplanes of every complete polynomial tensor and uses normalized
finite-field elimination.

The constructed span reaches dimension 5,175. Every retained direction
is checked entrywise for zero endpoint, zero derivative, and both zero
paired traces. The invariant upper bound (3) then proves equality (1).
Once that dimension is reached, processing further closure images is
unnecessary: all of them must lie in the already attained ambient
space.

The frozen replay reaches the ceiling after 230 queue vectors and
11,040 operator images. It checks 1,324,800 tensor entries for both
endpoint conditions and 165,600 polynomial entries of the paired
traces. Its canonical basis occupies 36,710 compressed bytes.

The [saturation record](../evidence/paired-trace-affine-saturation.json)
binds the input prefix, exact source files, constructed dimension, and
canonical output basis. A full replay is

```sh
python3 tools/saturate_paired_trace_affine_space.py \
  evidence/independent-paired-affine-degree-basis.json.gz \
  --output /tmp/paired-trace-affine-saturation.json \
  --basis-output /tmp/paired-trace-affine-canonical-basis.json.gz
```

The packed physical actions in (4) also have independent
exponent-dictionary controls in the degree-sixteen checkpoint. The
dimension statement concerns full polynomials, not evaluations at a
finite collection of places.

## 3. A sparse canonical basis

Take the following fifteen matrices as a basis of `s`:

```math
E_{ij}\quad(i\ne j),\qquad
E_{ii}+E_{00}\quad(i=1,2,3).                       \tag{5}
```

For every ordered pair `A,B` from (5) and every `0<=r<=22`, use

```math
w^r(1+w^2)\,A_{xv}B_{yu}.                        \tag{6}
```

These are 5,175 independent directions. Their span is the right-hand
side of (1). Every vector in (6) has only two polynomial degrees and
at most four nonzero tensor entries. The
[compressed canonical basis](../evidence/paired-trace-affine-canonical-basis.json.gz)
is a convenience copy; formulae (5)--(6) define it completely.

The array order is `x+4y+16u+64v`. A coefficient of degree `r` occupies
bit position `25*(x+4y+16u+64v)+r` in each of two field bitplanes.
The first plane is the coefficient of one and the second the coefficient
of `omega`. Every canonical basis vector lies in the first plane.

Because (1) is the full allowed space, multiplying a direction by `w`
and retaining it only when the degree stays at most 24 adds no new
directions. More general linear constraints on values or Hasse jets
can be imposed directly on this basis.

## 4. The complete degree-sixteen family

Let `P=E_00` and use the paired coordinates to define

```math
U=P\otimes I_4+I_4\otimes P,\qquad
T_{\min}(w)=C+(1+w^{16})U.                        \tag{7}
```

Since `tr(P)=1` and `tr(I_4)=0`, both partial traces of `U` equal
`I_4`. The tensor (7) therefore satisfies (2), has degree sixteen,
and has zero first derivative. It belongs to `T+W` by (1).

All degree-at-most-sixteen solutions of (2) are exactly

```math
T_{\min}+(1+w^2)(\mathfrak s\otimes\mathfrak s)[w]_{\le14}.
                                                               \tag{8}
```

This is an affine space of dimension `225*15=3375`. Each solution has
actual degree sixteen because a paired trace has degree sixteen.
Each has scalar polynomial content one: a common factor would divide
`1+w^16=(1+w)^16`, whereas the strict endpoint has nonzero entries.

There is also a direct solution preserving the low coefficients of the
source. Put

```math
Q(w)=\frac{T(w)-T_{\min}(w)}{1+w^2},
\qquad
S(w)=T_{\min}(w)+(1+w^2)\operatorname{trunc}_{<15} Q(w).
                                                               \tag{9}
```

The discarded quotient begins in degree fifteen. Thus `S` agrees with
`T` in every degree below fifteen, and it satisfies (2) with degree
sixteen. Since the source and `T_min` are even, the quotient and its
truncation are even. Both tensors have zero degree-fifteen coefficient,
so they in fact agree in every degree below sixteen.

Equivalently, if `T_<16` denotes the source truncation below degree
sixteen, then

```math
S(w)=T_{<16}(w)+w^{16}\bigl(C+T_{<16}(1)\bigr).   \tag{10}
```

The [uniform low-jet truncation theorem](low-jet-reflection-degree-truncations.md)
proves (10) directly for all reflection sizes. Equation (1) supplies
the additional finite-family conclusion here: this four-letter
truncation belongs to the complete independent-paired affine family.

## 5. What this description does and does not resolve

The earlier 3,035-dimensional prefix did not contain a degree-sixteen
representative preserving the source constant coefficient. Equations
(9)--(10) show that the complete space does. The prefix infeasibility
was therefore a limitation of that recorded prefix, as its original
statement specified.

This description removes the need to search physical-frame words merely
to obtain a spanning set. It also replaces dense polynomial bases by
the sparse basis (6). It does not determine local Smith profiles or
their weighted pressure. A tensor chosen from (8) still needs complete
cut profiles, compatible common physical types, and its actual degree
charged in the crossing argument.
