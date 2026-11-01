# Direct interpolation with prescribed crossing jets

## Constructive result and scope

Polynomial crossing tensors can be constructed by solving each entry's
Hermite equations first, then coupling the remaining choices through
the two paired traces. At alphabet four and degree cap 24, the coupled
stage has at most 800 displayed field equations, with rank at most 713.
It avoids eliminating a dense 5,175-direction polynomial basis for
every proposed target.

The [exact implementation](../tools/entrywise_reflection_constraints.py)
accepts coefficient protections, complementary cut rectangles, and
complete prescribed finite Hasse jets. It returns the full polynomial
tensor and replays every required identity. In particular, it constructs
a degree-sixteen path whose row cut agrees with a rank-one Bell tensor
through order five at `omega`, while its column cut does so at
`omega^2`. Thus fifteen of the sixteen local Smith orders are at least
six at each respective place.

This enlarges the available exact jet targets. It does not establish
an improved crossing pressure or factoring exponent. The complete
weighted profiles, common physical types, degree, and reconstruction
costs remain necessary.

## 1. Separate the entry constraints

Work over `F4` with the physical order `x,y,u,v`, and degree cap `d`.
For each of the 256 tensor entries, introduce its `d+1` polynomial
coefficients. Impose the following requested linear equations directly
on that one entry:

- its endpoint value `[x=v][y=u]` at one;
- its first derivative zero at one;
- specified coefficients at zero or at the homogeneous infinity end;
- specified Hasse coefficients at finite places.

The coefficient of degree `r` in the Hasse evaluation of order `j`
at `a` is

```math
\binom rj a^{r-j}.
```

In characteristic two, the binomial coefficient is nonzero exactly
when every binary digit of `j` is a digit of `r`. The construction uses
that identity; the output replay uses independent integer binomial
coefficients reduced modulo two.

For entry `e`, exact elimination either proves inconsistency or returns
one particular polynomial `p_e` and a complete basis `b_(e,k)` of its
homogeneous solution space. Every allowed entry is then

```math
p_e+\sum_k c_{e,k}b_{e,k}.                                \tag{1}
```

The degree bound is built into these variables. A protected source
coefficient above that bound is treated as an equation with zero left
side. If its required value is nonzero, the target is rejected.

## 2. Couple only the paired traces

The remaining global equations are

```math
\sum_x S_{xyux}(w)=(1+w^{16})\delta_{yu},\qquad
\sum_y S_{xyyv}(w)=(1+w^{16})\delta_{xv}.                  \tag{2}
```

Each entry appears in a trace only if `x=v` or `y=u`. Exactly
`2*4^3-4^2=112` entries meet that condition. Free variables from the
other 144 entries have zero columns in this second linear system.
The constructor chooses zero for those free variables when producing
one representative, while counting them in the full affine dimension.

For each of the 112 relevant entries, insert the trace images of its
local basis polynomials into a second exact elimination. The right side
is the desired traces minus those of the particular entries. A solution
reconstructs every polynomial in (1). An inconsistent residual proves
that the complete requested system has no solution.

There are `32(d+1)` displayed polynomial coefficient equations. At each
degree, the two matrix traces have the same final scalar contraction,
giving one redundancy. Equivalently, the simultaneous kernel is
`sl_4 tensor sl_4`, of dimension 225. Every local homogeneous variation
also has value and first derivative zero at one. Its trace image is
therefore divisible by `(w+1)^2`, leaving only `d-1` scalar coefficient
positions. The rank is at most `31(d-1)`, or 713 at degree cap 24 and
465 at degree cap 16. Additional entry constraints may lower that rank.

This two-stage elimination is complete: stage one parameterizes every
entrywise solution, and stage two imposes exactly the remaining global
conditions. When it is feasible, the affine solution dimension is
the sum of all local nullities minus the rank of the trace system.

For caps at most 24, every output satisfying (2) belongs to the
[complete independent-paired affine family](complete-paired-trace-affine-space.md).
The direct output identities already prove that it is a strict crossing;
that interpretation does not require recovering physical-frame words.

## 3. Rank-one jets at opposite places

For an invertible four-by-four matrix `A`, define two constant tensors

```math
\begin{aligned}
B_R(A)_{xyuv}&=A_{xy}(A^{-1})_{uv},\\
B_C(A)_{xyuv}&=A_{xu}(A^{-1})_{yv}.
\end{aligned}                                             \tag{3}
```

Their respective row and column cuts have rank one. Both paired
traces of either tensor are the identity matrix. Consequently the
finite targets

```math
S(\omega+t)\equiv\omega^2 B_R(A)\pmod{t^L},\qquad
S(\omega^2+t)\equiv\omega B_C(B)\pmod{t^L}                \tag{4}
```

are compatible with (2) for `L<=16`. Indeed, the trace polynomial
`1+w^16` has constant Hasse jet `omega^2` at `omega`, constant jet
`omega` at `omega^2`, and zero coefficients of orders one through
fifteen at either point.

If (4) holds with `L>=1`, the corresponding flattening has one unit
Smith channel and fifteen channels of order at least `L`. To see this,
use constant invertible row and column operations to send the nonzero
rank-one leading matrix to a single nonzero entry. Every other entry
is divisible by `t^L`. Eliminating the unit pivot leaves a fifteen-by-
fifteen Schur complement divisible by `t^L`. The generic rank is
sixteen because the endpoint flattening is a permutation matrix.

The checker supplies one exact cap-sixteen solution of (4) with
`L=6`, taking `A` to be the source's zero Bell factor and `B=I`.
Its full tensor hash is

```text
a137ce9c8e40b6e8d0fd1807e3d3904cfdc33fea88946d5337e21374fcd37568
```

It also supplies coordinate-rectangle examples. In one degree-sixteen
example, both cuts have value rank four at both `omega` and `omega^2`.
Those examples establish their requested coefficient identities;
stronger complete-profile or pressure claims require separate replay.

### A uniform seventh-order construction

In fact, every pair of invertible matrices `A,B` admits (4) with `L=7`
and degree cap sixteen. This statement follows directly from scalar
Hermite interpolation and does not depend on a numerical search.
Set

```math
K(w)=(w+1)^2(w^2+w+1)^7.
```

Let `e_1,e_omega,e_(omega^2)` be the scalar CRT idempotents of degree
less than sixteen for the respective moduli `(w+1)^2`,
`(w+omega)^7`, and `(w+omega^2)^7`. Put

```math
T_0=e_1 C+\omega^2 e_\omega B_R(A)
                 +\omega e_{\omega^2} B_C(B).
```

Both paired traces of `C` vanish because the alphabet has size four
in characteristic two. The two trace scalars of `T_0` are equal.
They agree with `tau=1+w^16` modulo each of the three coprime
moduli. Their difference from `tau` is monic of degree sixteen,
so that difference is exactly `K`. Consequently

```math
S=T_0+K V                                                     \tag{5}
```

has the required traces, endpoint, first derivative, and two seventh-
order Bell jets whenever the constant tensor `V` has both paired
traces equal to the identity. Such a carrier exists: either `B_R(A)`
or `B_C(B)` works. The carriers form an affine space of dimension 225.
Because their paired traces are nonzero, every carrier is nonzero;
the leading coefficient of (5) is `V`, and its degree is exactly
sixteen.

The [explicit constructor](../tools/bell_jet_crt_crossings.py) implements
(5) and checks its full coefficient identities. Useful identities for
this representative are

```math
e_1=w(w^2+w+1)^7,\qquad
e_\omega(0)=\omega,\qquad e_{\omega^2}(0)=\omega^2.
```

Thus choosing `V=B_R(A)` gives `S(0)=B_C(B)` and leading coefficient
`B_R(A)`; choosing `V=B_C(B)` exchanges those two constant ends.
These value-rank statements do not specify all local Smith orders.
The extra carrier freedom and complete profiles must still be charged
by the mixed-place pressure criterion.

The [separate CRT checker](../tools/check_bell_jet_crt_crossings.py)
uses integer binomial coefficients and an independent field multiplier
to replay 32 complete tensors, including nonsymmetric relative Bell
forms and carriers outside the span of the three pairings. Its
[record](../evidence/bell-jet-crt-controls.json) separates the scalar
Hermite equations, full tensor identities, and rank-31 constant trace
map. Run it with

```sh
python3 tools/check_bell_jet_crt_crossings.py \
  --output evidence/bell-jet-crt-controls.json
```

## 4. Coefficient protections and search interpretation

The numerical searches use the natural physical heights
`h=(0,6,2,8)`. At zero, an entry coefficient of degree `r` has charged
index `r+h_x+h_y`; at the source's homogeneous infinity end it has
index `24-r+h_u+h_v`. Exact protection retains the source coefficient
below a selected cutoff at either end.

A separate, explicitly named zero-pattern mode retains only the source
zeros in that region. Existing nonzero coefficients may change or
disappear. This weaker condition may admit lower degrees. It is not a
Smith-order guarantee: changing coefficients can destroy cancellations
inside minors. The complete profiles are recomputed for every candidate.

The [bounded search](../tools/search_entrywise_reflection_jets.py)
prioritizes full Bell jets, then complementary rectangles. It records
the actual degree, exact feasibility, all finite ranks, source hashes,
and completed targets. Natural two-place pressure is only a diagnostic
for retaining candidates. Promising full seeds are passed to the
[mixed-place](mixed-place-hermite-crossings.md) and
[asymmetric mixed-place](asymmetric-mixed-place-crossings.md) objectives.

## 5. Independent controls

```sh
python3 tools/check_entrywise_reflection_constraints.py \
  --output evidence/entrywise-reflection-constraint-controls.json
```

The [checker](../tools/check_entrywise_reflection_constraints.py)
compares 36 complete small affine solution spaces with all 12,240
possible field assignments in their declared dimensions. It checks
six full tensor targets, including both feasible and infeasible cases.
Every feasible output is replayed entrywise for its endpoint, first
derivative, paired polynomial traces, protected coefficients, and
finite Hasse equations. The [evidence](../evidence/entrywise-reflection-constraint-controls.json)
records the exact source hashes and resulting tensor hashes.
