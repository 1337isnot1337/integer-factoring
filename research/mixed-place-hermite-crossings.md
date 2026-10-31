# Crossing pieces from opposite flags at several places

## Result and scope

A crossing piece can have one rank estimate from a combination of jets
at several places, and its other rank estimate from jets at a disjoint
set of places. Complete jet flags on the two sets are opposite. A basis
adapted to both flags supplies an effective crossing decomposition.

This changes the reconstruction step of the
[same-place Hermite construction](multipoint-hermite-crossings.md).
For its ternary polynomial seed, an exact mixed-place certificate gives
`rho < 9/10`. Every unweighted same-place certificate for that seed has
`rho >= log(54)/(4 log(3)) > 9/10`. Thus the new reconstruction gives
a strict improvement for this seed, after allowing every same-place
degree allocation and every pair of local scales.

The improved seed parameter is still much larger than `1/sqrt(2)`.
It does not improve the established factoring score `1/2`. The result
is a new sufficient construction and a precise mechanism to search.

## 1. Polynomial seed and disjoint place sets

Let `F` be a fixed finite field and `q>=2`. Let the fixed tensor

```math
K(z;x,y,u,v)\in F[z],\qquad 0<=x,y,u,v<q
```

have degree at most `D` and satisfy

```math
K(1;x,y,u,v)=[x=v][y=u].                                  \tag{1}
```

Choose two nonempty disjoint sets of rational places `A,B` on
the projective line over `F`, neither containing one. At a finite
place `a`, use the section `K_a(t)=K(a+t)`. At infinity use
`K_infinity(t)=t^D K(t^{-1})`. The complete degree `D` is charged;
common scalar factors and surplus degree are not discarded.

For each `a in A`, let `nu_(a,i)` be the local Smith orders in the
first cut `(x,y)|(u,v)`. For each `b in B`, let `mu_(b,i)` be the
orders in the second cut `(x,u)|(y,v)`. Both cuts have generic rank
`q^2`. Define

```math
U_a(s)=\sum_i s^{-\nu_{a,i}},\qquad
V_b(s)=\sum_i s^{-\mu_{b,i}}.                              \tag{2}
```

Choose positive rational numbers `alpha_a,beta_b` with

```math
\sum_{a\in A}alpha_a=\sum_{b\in B}beta_b=1.                 \tag{3}
```

For a real `s>1`, put

```math
R(s)=\max_{a\in A}U_a(s^{1/alpha_a}),\qquad
C(s)=\max_{b\in B}V_b(s^{1/beta_b}).                        \tag{4}
```

**Mixed-place theorem.** If a rational `1/2<=rho<1` satisfies

```math
\boxed{s^D R(s)C(s)<q^{4rho},}                             \tag{5}
```

then the seed supplies effective scale-one crossing pieces with
parameter `rho`. There are `O(lambda)` pieces at binary label width
`lambda`, and their complete tables and rank factorizations are
constructible in `2^O(lambda)` bit operations.

When the scales and powers in (5) are rational, raising the inequality
to a fixed integer clears rational exponents. A final certificate then
uses exact integer comparisons. The theorem does not rely on numerical
optimization of (3)--(4).

## 2. Opposite Hasse-jet flags

Let `P_N` be the `(N+1)`-dimensional space of polynomials of degree at
most `N`. A finite-place jet is the functional

```math
J_{a,k}(P)=[t^k]P(a+t).
```

The infinity jet is `J_(infinity,k)(P)=[t^k]t^N P(t^{-1})`.
These are Hasse coefficients, so every statement below is valid in
positive characteristic without division by factorials.

At the places in `A`, order jets so that a prefix of length `r`
contains an initial segment of jets at every place and

```math
r_a<=alpha_a r+O(1).                                      \tag{6}
```

One explicit order repeats a fixed word containing `L alpha_a` copies
of each place, where `L` clears all denominators. On each visit, append
the next unused jet at that place. The error in (6) is at most `L`.
The first `r` jets span a subspace `R_r` of `P_N^*`. Build `C_c`
from `B` and `beta` in the same way.

For all `0<=r,c<=N+1`,

```math
\dim R_r=r,\qquad \dim C_c=c,\qquad
\dim(R_r\cap C_c)=\max(0,r+c-N-1).                        \tag{7}
```

To prove (7), vanishing of the first `r` jets is precisely divisibility
by the corresponding finite-place factors, together with the infinity
degree cap if applicable. It gives `r` independent conditions while
`r<=N+1`. The two place sets are disjoint, so imposing both jet sets
gives `r+c` conditions, capped by the whole dimension. Equivalently,

```math
\dim(Ann(R_r)\cap Ann(C_c))=\max(0,N+1-r-c).
```

Taking dimensions of the sum of the two annihilators proves the last
identity in (7). Infinity can occur in at most one set, and the same
degree-cap argument applies.

For each `0<=k<=N`, choose a nonzero functional

```math
L_k\in R_{k+1}\cap C_{N-k+1}.                             \tag{8}
```

This intersection has dimension one. Moreover, `L_k` does not belong
to `R_k`, since `R_k` and `C_(N-k+1)` have zero intersection. Hence
`L_0,...,L_N` are a basis of the dual space. In particular, computable
field elements `c_k` satisfy

```math
P(1)=\sum_{k=0}^N c_k L_k(P)\qquad(P\in P_N).              \tag{9}
```

This basis is effective: write the Hasse jets in the monomial basis,
solve the `N+1` homogeneous systems in (8), and solve (9). Straight
Gaussian elimination takes a polynomial number of field operations
in `N`, even without sharing work between systems. All matrices have
dimensions `O(N)`; no changing-parameter oracle is being assumed.

## 3. Two rank estimates for the same piece

Set `N=Dm`. Apply (9) entrywise to `K(z)^(tensor m)`. Its pieces are

```math
A_{m,k}=c_k L_k(K^{\otimes m}),\qquad 0<=k<=Dm.             \tag{10}
```

Their sum is the crossing on `q^m` letters. Each piece has two
representations. Its first-cut representation is a linear combination
of the first `k+1` row-flag jets. Its second-cut representation uses
the first `Dm-k+1` column-flag jets. The two linear combinations need
not have the same coefficients or involve the same places; they define
the same functional by (8).

For a fixed local Smith profile, the standard coefficient estimate is

```math
rank[t^j]M(t)^{\otimes m}
 <=poly(m) u^j\left(\sum_i u^{-\nu_i}\right)^m,
\qquad u>=1,\quad j=O(m).                                \tag{11}
```

At row place `a`, (6) gives `j<=alpha_a k+O(1)`. Insert
`u=s^(1/alpha_a)` in (11). Summing the `O(m)` jet matrices in the
representation of (10) yields

```math
r(A_{m,k})<=poly(m)s^k R(s)^m.                            \tag{12}
```

The same argument for the column flag gives

```math
c(A_{m,k})<=poly(m)s^{Dm-k}C(s)^m.                        \tag{13}
```

The coefficient indices cancel on multiplication:

```math
r(A_{m,k})c(A_{m,k})
 <=poly(m)\left[s^D R(s)C(s)\right]^m.                    \tag{14}
```

The degree factor is `s^D`, not `s^(2D)`. Each flag has a complete
basis, but a single piece only uses complementary prefix lengths.

Under (5), (14) gives the required sum of factor widths with a
logarithmic overhead. For arbitrary `lambda`, take the least `m`
with `q^m>=2^lambda`, restrict to `2^lambda` labels, and use ordinary
rank factorization to cap both widths at `2lambda`. The strict margin
absorbs polynomial factors. Finitely many small widths can be supplied
by a fixed coordinate decomposition.

All `q^(4m)` polynomial entries have degree `O(m)`. Their construction,
the polynomial-cost flag computation, every complete coefficient table,
and exact rank factorization take `exp(O(m))` field operations. The
field and seed are fixed, so this is the required bit cost. The
[scale-one transfer](scale-one-modular-transfer.md) then applies.

## 4. A strict separation on the ternary seed

Use the existing seed over `F_3`:

```math
K(z;x,y,u,v)=[x+y=u+v]\bigl(1-(x-u-z(v-u))^2\bigr).        \tag{15}
```

Its degree is two, and its endpoint at one is the crossing. Its local
profiles, independently reproduced by complete block-Toeplitz ranks,
are:

| Place | First cut | Second cut |
| --- | --- | --- |
| `0` | `0` nine times | `0,1,1,1,1,2,4,4,4` |
| `2` | `0,0,0,1,1,1,2,2,2` | `0` nine times |
| infinity | `0,0,0,1,1,1,2,2,2` | `0` nine times |

Take `A={2,infinity}`, `B={0}`, `alpha=(1/2,1/2)`, `beta=(1)`,
and `s=5/3`. Then

```math
\begin{aligned}
R(s)&=3(1+s^{-2}+s^{-4})=2793/625,\\
C(s)&=1+4s^{-1}+s^{-2}+3s^{-4}=2593/625,\\
s^2R(s)C(s)&=2414083/46875.
\end{aligned}                                                     \tag{16}
```

Exact integer comparison gives

```math
2414083^5<3^{18}46875^5.                                  \tag{17}
```

Thus (5) is strict at `rho=9/10`; by its strict margin, a rational
parameter below `9/10` is also available.

### The whole unweighted same-place criterion is weaker here

Allow any nonnegative allowances `d_0,d_2,d_infinity` summing to two,
and any independent scales in the old same-place criterion. If
`d_0>=1`, the six second-cut channels of orders `0,1,1,1,1,2` have
mean one. For every `r>=1`, their arithmetic-geometric mean gives

```math
\sum_i r^{-\mu_i}r^{d_0}>=6.
```

The first cut at zero has nine zero orders. Its contribution is at
least nine, so `F_0>=log54`.

If `d_0<1`, one of `d_2,d_infinity` is greater than `1/2`. At that
place, the six first-cut orders `0,0,0,1,1,1` have mean `1/2` and
give the same lower bound six at every scale. The other cut again
contributes at least nine. Consequently

```math
\max_j F_j>=\log54,\qquad
rho>=\frac{\log54}{4\log3}>\frac9{10}.                    \tag{18}
```

The last strict inequality is the integer comparison `54^5>3^18`.
This excludes a same-place certificate at `9/10` without assuming that
the exploratory numerical optimum was found.

There are no helpful additional singular places over a field extension.
The first-cut determinant spends its entire homogeneous degree eighteen
at place two and infinity, nine at each. The second-cut determinant
spends its entire degree eighteen at zero. At every other place both
local profiles are zero, giving cost at least `log81` for a same-place
piece. Additional zero-allowance places can simply be omitted.

Numerically optimizing the displayed unweighted formulas gives about
`0.948005` for the old criterion and `0.896878` for the mixed-place
criterion. Those numerical values are exploratory. The strict
separation is proved only with (17)--(18). The comparison does not
claim an optimum for weighted versions of either construction.

## 5. Independent paired types and local physical gauges

The construction also permits rectangular physical alphabets and fixed
rational paired types `pi_X,pi_Y`, with positive entropy sum `H`.
Take `m` through multiples of their fixed type denominators. All local
physical heights in this statement are integers, so every type shift
`m A_a` and `m B_b` is an integer on these admissible lengths.
At a row place `a`, gauge the first cut by independent leg heights and
let its restored profile sum be `U_a`. Let `A_a` be the total mean
height on the chosen paired types. Define `V_b,B_b` in the second cut
at column place `b` in the same way. These are the height conventions
of [the weighted same-place note](weighted-multipoint-hermite-crossings.md).

On every retained tensor entry, a row gauge is the scalar power
`t^(m A_a)`. Thus the original row jet of index `j` has gauged index
`j+m A_a`; the reconstruction functional is unchanged. The corresponding
column shift is `m B_b`. Restoring these shifts in (12)--(13) gives

```math
\begin{aligned}
R_{pi}(s)&=\max_a s^{A_a/alpha_a}U_a(s^{1/alpha_a}),\\
C_{pi}(s)&=\max_b s^{B_b/beta_b}V_b(s^{1/beta_b}).
\end{aligned}                                                     \tag{19}
```

The sufficient fixed-type criterion is

```math
\log(s^D R_{pi}(s)C_{pi}(s))<2rho H.                       \tag{20}
```

Tensor with the pair-swapped copy, keeping the two functional indices
separate, to obtain equal physical alphabets. Its rank-cost calculation
is exactly the [paired-type transfer](paired-type-crossing-pressure.md).
Negative local orders are handled by a fixed scalar clearing power and
its restored coefficient shift.

One common pair of types must satisfy every row and column place cost
in (19). Choosing a separate Gibbs type at every place would not define
one crossing on a fixed alphabet. No such substitution is made here.

The `3/4` floor for same-place Hermite estimates does not apply to
(20): a piece in (10) is generally a combination of several local jets,
with its two cut estimates coming from disjoint flags. This observation
is not a claim that (20) can pass `1/sqrt(2)`.

## 6. Exact reconstruction checks

Run

```sh
python3 tools/check_mixed_place_hermite.py \
  --output evidence/mixed-place-hermite.json
```

The checker builds both Hasse-jet flags and solves the intersections
over characteristics two, three, and five. It checks every defining
linear equation, both representations of every piece, complete
crossing reconstruction, and independent exact cut ranks for small
tensor powers. The ternary profiles and (17)--(18) are reproduced.
One negative control lets the two place sets overlap; their flags then
fail the claimed one-dimensional intersection property.

The seed (15), original Smith estimates, and factoring transfer are
prior work recorded in this repository. The opposite-flag reconstruction,
separation argument, and checker used OpenAI assistant support.
