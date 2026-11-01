# Stable small physical codes from larger reflection intervals

## Result and scope

A small physical alphabet can be constructed directly inside a large
reflection tensor. The ambient alphabet does not need to be expanded.
For bounded-energy labels, the primitive tensor eventually becomes
independent of the reflection interval. This gives a finite, explicit
family of non-Cartesian physical restrictions beyond the earlier
[interval-eight code screen](reflection-code-alphabets.md).

The first sixteen-letter example has degree 88 and numerical natural-gauge
cost `rho = 0.71473663132372`. Its complete natural profiles are computed
from all 256 channels of both cuts. This is below the earlier
sixteen-letter interval-eight natural value `0.71665365972819`, but is
still above `1/sqrt(2)`. No improved factoring exponent is claimed.

The constructor, interval-stability statement, polynomial coefficients,
and endpoint identity below are exact. Pressure minimization and later
gauge searches are numerical diagnostics unless separately certified.

## 1. Restrict the uniform implicit tensor

Let `M>=4` be a power of two and put `m=M/2`. The
[uniform reflection compression](reflection-dual-rail-compression.md)
has physical alphabet `F_2^m` and coefficient formula

```math
T_M(w;x,y,u,v)=w^{M^2/2}\operatorname{Pf}(K_M(w)|_{x,y,u,v}).
```

Here a physical bit indicates whether the corresponding logical
Gaussian variable is present in the principal submatrix. The source
proof shows that every restored coefficient is a polynomial and that

```math
T_M(1;x,y,u,v)=[x=v][y=u].                              \tag{1}
```

Choose any nonempty physical label sets `C_X,C_Y` and retain only
`x,v in C_X` and `y,u in C_Y`. Equation (1) immediately remains the
literal paired crossing on these smaller alphabets. The two sets need
not have the same size and need not be Cartesian sets of logical bits.

[`restrict_implicit_reflection.py`](../tools/restrict_implicit_reflection.py)
evaluates only the requested Pfaffians, in the regular zero gauge.
It restores each physical height separately, rejects any pole, and
checks every endpoint entry. It makes no rank claim by deleting
one-particle Smith channels: all later cut profiles are computed from
the complete restricted tensor.

If a label set is closed under deleting occupied bits, every recursive
Pfaffian submask is another requested physical tuple. Consequently at
most `|C_X|^2 |C_Y|^2` Pfaffians are cached. Each uses at most one sum
over the remaining occupied bits. Natural-energy prefixes have this
closure even when only some labels at their last energy are retained,
because every bit has strictly positive energy.

## 2. Reverse the logical coordinates

Write the hot coordinate as `a=m-1-i`, where `i` is the source logical
index. Its natural physical energy is

```math
E_a=8a+4.                                               \tag{2}
```

Thus a fixed hot state `S`, a finite subset of nonnegative integers,
maps to the physical mask

```math
\phi_M(S)=\sum_{a\in S}2^{m-1-a},
\qquad E(S)=\sum_{a\in S}(8a+4),                       \tag{3}
```

provided every `a<m`. In these coordinates all cross-leg entries of
the source Laurent kernel are independent of `M`:

```math
\begin{aligned}
K_{v_a,y_b}=K_{u_a,x_b}
  &=(1+w^2)w^{4(a\mathbin{\&}\neg b)},\\
K_{v_a,u_b}&=(1+w^2)w^{4(a\mathbin{\&}b)},\\
K_{y_a,x_b}&=(1+w^2)w^{-4(1+(a\mathbin{|}b))},\\
K_{v_a,x_b}=K_{y_a,u_b}&=[a=b].
\end{aligned}                                         \tag{4}
```

The kernel is symmetric with zero diagonal. The bitwise operations
in (4) are on nonnegative integers; the first expression means the
bits in `a` that are absent from `b`.

For example, with `i=m-1-a,j=m-1-b`, the source exponent
`-2M+4(i&j)` equals `-4(1+(a|b))`. The other identities follow from
the same complementation inside the fixed power-of-two bit interval.

The remaining source entries lie within one physical leg. Such an
entry joins the two hot coordinates satisfying

```math
a+b=m-1.
```

Their combined energy is exactly

```math
(8a+4)+(8b+4)=8m=4M.                                  \tag{5}
```

## 3. Exact stabilization theorem

Fix two finite hot-state codes. Suppose they fit in the interval and
every retained state has energy strictly less than `4M`. No retained
state can contain both endpoints of an edge in (5). Therefore no
within-leg kernel edge occurs in any requested principal submatrix.
Every remaining entry is one of the interval-independent entries (4).

There is consequently one Laurent tensor `G_C(w)`, depending only on
the two hot codes, such that

```math
T_{M,C}(w)=w^{M^2/2}G_C(w)                              \tag{6}
```

for every eligible interval. If `M' >= M` is another eligible interval,
the complete coefficient identity is

```math
T_{M',C}(w)=w^{(M'^2-M^2)/2}T_{M,C}(w).                 \tag{7}
```

Both sides are polynomial tensors by the source theorem. After exact
division by their full common scalar polynomial and normalization at
`w=1`, their primitive strict tensors are identical. This division is
legal because a common factor cannot vanish at the strict endpoint.
Over the binary source field its value at that endpoint is one.

The energy inequality is strict. At energy `4M`, the within-leg
reflection pairs can reappear, so (7) is not asserted without checking
their absence. Nor is stabilization a statement that arbitrary
physical restrictions have the same pressure.

## 4. The first sixteen-letter stabilized code

For `M=16`, order physical masks by `(natural energy, integer mask)`
and retain the first sixteen. They are

```text
0, 128, 64, 192, 32, 160, 16, 96,
144, 8, 224, 80, 136, 4, 208, 48.
```

Their energies are

```text
0, 4, 12, 16, 20, 24, 28, 32,
32, 36, 36, 40, 40, 44, 44, 48.
```

The maximum is `48<64=4M`, so the theorem applies. The raw tensor has
degree 176 and common scalar `w^88`. The primitive degree is 88. Its
complete coefficient hash is

```text
3429c10e5cec92f1eb60e5a1c9aa23b1e2b74e9b205f3458bc953dbc82e5f7f1
```

At `M'=32`, the corresponding masks are the listed integers shifted
left by eight bits. The checker verifies every one of the 65,536
coefficients in (7), whose common monomial exponent is 384. Both
constructions cache exactly 65,536 Pfaffians.

Natural physical gauges use the complete label energies in the same
way as [the earlier code pressure](reflection-code-alphabets.md).
The two complete integer profiles are recorded in the evidence.
At a positive scale `t`, write their channel sums as `Z_R(t),Z_C(t)`.
The common paired Gibbs type is determined by the positive root `R` of

```math
\log Z_R(t)+\log Z_C(t)
  =2R\log\sum_{x\in C}e^{-tE(x)/R},\qquad \rho=R/2.    \tag{8}
```

An independent evaluation of (8), using only those complete integer
profiles, gives the value in the opening paragraph. A scalar search
then varies `t`; it does not replace the entropy denominator by a
constant. The finite-profile computation and the transcendental
pressure evaluation are recorded separately.

## 5. Exact controls and ongoing finite search

Run

```sh
OPENBLAS_NUM_THREADS=1 python tools/check_implicit_reflection_codes.py \
  --output evidence/implicit-reflection-code-controls.json
```

The [checker](../tools/check_implicit_reflection_codes.py) passed:

| Exact control | Count |
|---|---:|
| Complete coefficients against existing interval-four/eight tensors | 65,936 |
| Independent original selected-minor coefficients at interval sixteen | 324 |
| Complete interval-sixteen/thirty-two stabilization coefficients | 65,536 |
| Strict endpoint entries | 197,332 |
| Complete integer natural-profile channels through two interfaces | 512 |

The selected-minor controls use the original strict one-particle matrix
and exact polynomial determinants, independently of the Laurent
Pfaffian constructor. The profile comparison uses the rectangular
tensor interface and a separately assembled matrix interface; both
perform exact finite-field arithmetic, with integral pivot weights.
The [evidence](../evidence/implicit-reflection-code-controls.json)
binds their sources and retains the complete profiles.

[`search_implicit_reflection_codes.py`](../tools/search_implicit_reflection_codes.py)
constructs every boundary choice for a declared list of natural-energy
prefix sizes, removes full scalar content, and records complete
natural-gauge pressures. Its finite search record and any full-gauge
refinement are separate numerical evidence. No numerical minimum in
this family is a global exclusion or a factoring exponent certificate.
