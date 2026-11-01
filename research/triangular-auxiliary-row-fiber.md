# A triangular auxiliary fiber preserves every word Smith profile

## Result and scope

For the four-letter reflection Gaussian, fix its natural physical gauges
and its two actual local **row lattices**. The endpoint-preserving Laurent
perturbations that remain in both row lattices form a vector space of
dimension `14`. Every such perturbation is a strictly height-increasing
operation on the common physical `x` leg.

This classification also controls auxiliary matrices of arbitrary finite
dimension. Give each of the fourteen directions an arbitrary constant
auxiliary matrix, including matrices that do not commute. Every resulting
ordered word, after auxiliary boundary contraction, is a regular
unitriangular physical change of the scalar tensor power at both endpoints.
Consequently its complete word Smith profiles are unchanged.

Thus this entire fiber cannot improve a certificate obtained by counting
regular rank-one channels with their endpoint orders. The conclusion
includes state-dependent choices, complete DVR basis changes, larger
full-tail state families, and blocking. It concerns the specified fiber
and certificate method. It is not a lower bound on actual coefficient
ranks, a classification of all crossings, or a factoring lower bound.

The scalar seed is the `M=4` instance of
[reflection-pair compression](reflection-dual-rail-compression.md).
The finite interpolation argument is related to
[polynomial crossing interpolation](polynomial-crossing-interpolation.md),
but the constraints here prescribe the actual **row module**. This is
different from prescribing divisibility of transformed matrix rows.

## 1. Seed, gauges, and the exact linear fiber

Remove the common even ramification from the reflection seed. Write the
resulting polynomial tensor as

```math
T(z;x,y,u,v),\qquad x,y,u,v\in\{0,1,2,3\}.
```

Its degree is `12`, and

```math
T(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.
```

The four physical heights are

```math
h=(0,3,1,4).
```

At zero use the `(x,y)|(u,v)` cut and multiply an entry by
`z^(h(x)+h(y))`. At infinity put `tau=1/z`, use the `(x,u)|(y,v)`
cut, and multiply an entry by
`tau^(h(u)+h(v)+4)`. The common extra power `4` makes the infinity
matrix regular. Denote the resulting matrices by `T_0` and `T_inf`.

Their ordered Smith profiles are

```math
\begin{aligned}
\nu_0&=(4,5,7,8,9,10,11,12,12,13,14,15,16,17,19,20),\\
\nu_\infty&=(0,1,3,4,5,6,7,8,8,9,10,11,12,13,15,16).
\end{aligned}
```

In particular, each profile is the multiset of subset sums of
`1,3,5,7`, with an additional `4` at zero.

Let `Delta` be a Laurent tensor with `Delta(1)=0`. Require

```math
\mathrm{row}_{\mathcal O_e}(\Delta_e)
\subseteq
\mathrm{row}_{\mathcal O_e}(T_e),
\qquad e\in\{0,\infty\},                         \tag{1}
```

where `O_e=F_2[[tau]]` is the respective local power-series ring.
The flattenings and gauges of `Delta_e` are the same as those of `T_e`.

Both row modules in (1) are regular. Therefore every monomial `z^k`
that occurs in the entry at `(x,y,u,v)` lies in the forced window

```math
-h(x)-h(y)\le k\le h(u)+h(v)+4.                  \tag{2}
```

There is no separate degree cutoff. The condition `Delta(1)=0` gives
the complete entrywise basis

```math
(1+z)z^k E_{x,y,u,v},\qquad
-h(x)-h(y)\le k<h(u)+h(v)+4.                     \tag{3}
```

This basis contains `3,072` variables.

To impose (1), take regular unit matrices satisfying

```math
L_e T_e R_e=\mathrm{diag}(\tau^{\nu_{e,j}}).
```

The row-module condition is exactly

```math
\tau^{\nu_{e,j}}\mid(\Delta_e R_e)_{i,j}
\quad\text{for every }i,j.                       \tag{4}
```

No left change is applied to `Delta_e` in (4). Right multiplication by
`R_e` identifies the row module with the diagonal row module; the
regular invertible matrix `L_e` merely changes its generators.
Only the finite jets preceding the specified column orders occur in
(4). The checker constructs deterministic unit changes through
precision `32`, which exceeds every required order.
This suffices for the exact row module: a residual term divisible by
`tau^32` remains regular after division of any column by its target
power, and the resulting error to the identity has zero constant term.
It therefore changes the diagonal generators by a regular unit matrix.

The exact binary constraint matrix has rank `3,058`. Its kernel has
dimension `14`.

### Complete description

Every solution of (1) is uniquely of the form

```math
\Delta(z)=N_x(z)T(z),
\qquad
N_x(z)=
\sum_{h(a)>h(b)}\ \sum_{k=h(b)-h(a)}^{-1}
c_{a,b,k}(1+z)z^k E_{a,b},                       \tag{5}
```

with binary coefficients `c_(a,b,k)`. Here `E_(a,b)` sends physical
letter `b` to letter `a` at the `x` leg. There are

```math
\sum_{h(a)>h(b)}(h(a)-h(b))=14
```

displayed directions. The checker proves that all fourteen satisfy
the constraints and are independent. The exact kernel dimension then
proves completeness, including all Laurent degrees allowed by (1).

Each direction strictly increases `h`. At zero its coefficient in
the gauged leg map is

```math
z^{h(a)-h(b)+k}(1+z),
```

which is regular because `k>=h(b)-h(a)`. At infinity its coefficient
is `tau^(-k-1)(1+tau)`, which is regular because `k<=-1`.
Consequently `I+N_x` and its finite geometric-series inverse are
regular at both ends.

## 2. Arbitrary noncommuting auxiliary coefficients

Let `D` be any positive integer. Consider every matrix-valued crossing

```math
K(z;x,y,u,v)=T(z;x,y,u,v)I_D+\Delta(z;x,y,u,v),    \tag{6}
```

where each auxiliary entry of `Delta` satisfies (1) and vanishes at
`z=1`. By (5), this is the full family

```math
K(z;x,y,u,v)=\sum_b A(z;x,b)T(z;b,y,u,v),         \tag{7}
```

with

```math
A(z;a,b)=\mathbf1_{a=b}I_D+
\sum_{k=h(b)-h(a)}^{-1}(1+z)z^k C_{a,b,k}
\quad\text{when }h(a)>h(b),                       \tag{8}
```

and zero off-diagonal entries otherwise. The fourteen matrices
`C_(a,b,k)` are arbitrary elements of `Mat_D(F_2)`. They need not
commute. Thus the family has `14D^2` independent binary parameters.
Every member has the strict matrix endpoint `K(1)=E I_D`.

Take constant auxiliary boundaries `f,e` with `f^T e=1`, as required
to obtain the scalar crossing endpoint. For a word of length `t`, put

```math
\Gamma_t(\mathbf x,\mathbf y,\mathbf u,\mathbf v)
=f^T K(x_1,y_1,u_1,v_1)\cdots
K(x_t,y_t,u_t,v_t)e.                              \tag{9}
```

Expanding (7), with the auxiliary factors kept in their original order,
gives

```math
\Gamma_t=U_t\,T^{\otimes t},
\qquad
(U_t)_{\mathbf x,\mathbf b}
=f^T A(x_1,b_1)\cdots A(x_t,b_t)e.                \tag{10}
```

The operator `U_t` acts only on the physical `x` word. At either
endpoint, every factor in (10) is regular after the corresponding
leg gauge. A nonzero factor has `h(x_i)>=h(b_i)`, and a physical
off-diagonal factor has a strict inequality. It follows that

```math
(U_t)_{\mathbf x,\mathbf b}\ne0,
\quad\mathbf x\ne\mathbf b
\quad\Longrightarrow\quad
\sum_i h(x_i)>\sum_i h(b_i).                     \tag{11}
```

For equal physical words, all local factors are `I_D`, so the diagonal
entry is `f^T e=1`. Order words by total height. Then `U_t` is regular
unitriangular. More explicitly, `(U_t-I)^(4t+1)=0`, and

```math
U_t^{-1}=\sum_{j=0}^{4t}(U_t-I)^j
```

is regular. This argument uses no commutativity of the auxiliary
matrices.

The `x` word is on the row side of both chosen cuts. Therefore (10)
is a regular invertible row operation on each full word flattening.
It preserves the entire Smith profile, at every word length:

```math
\boxed{\quad
\mathrm{Smith}_e(\Gamma_t)
=\mathrm{Smith}_e(T_e^{\otimes t}),
\qquad e=0,\infty.\quad}                        \tag{12}
```

The same argument applies to any seed and any common physical leg
with a strict height order, whenever the matrix-valued leg map is
regular and unitriangular at the two relevant endpoints. The exact
linear classification above establishes that this general lemma
covers the **whole specified fiber** for the reflection seed.

## 3. Consequence for regular-channel pressure

For a regular matrix `M` with finite Smith orders `mu_j`, define

```math
\Phi_s(M)=\sum_j s^{-\mu_j},\qquad s>1.
```

Suppose a regular rank-one decomposition has the form

```math
M=\sum_{a}\tau^{w_a}\ell_a r_a,
```

with `w_a>=0` and regular factors. Then

```math
\Phi_s(M)\le\sum_a s^{-w_a}.                     \tag{13}
```

To see this directly, reduce modulo `tau^P`. The matrix then defines
a binary linear map of rank `sum_j max(P-mu_j,0)`. A rank-one term
of charged order `w_a` factors through a space of dimension at most
`max(P-w_a,0)`. Rank subadditivity gives the corresponding inequality
for each `P`. Multiply by
`(1-s^-1)^2 s^(-(P-1))` and sum over `P>=1`; the identity

```math
(1-r)^2\sum_{P\ge1}(P-w)_+r^{P-1}=r^w
\quad(0<r<1)
```

proves (13). Uniform scalar shifts handle a Laurent presentation.

A regular-channel word certificate, after the auxiliary terminal
factor is included, is exactly such a rank-one decomposition.
By (12), its counted endpoint pressure must therefore be at least

```math
\Phi_s(T_e)^t.
```

Any valid exponential growth rate for those counts is at least the
scalar rate `Phi_s(T_e)`. This remains true when the channel choices
depend on states or when several consecutive sites are refactored
together. Additional states and a more complete local basis oracle
cannot remove the full-word lower bound (13).

For this seed, write

```math
g(s)=(1+s^{-1})(1+s^{-3})(1+s^{-5})(1+s^{-7}).
```

The two regularized scalar rates are `s^-4 g(s)` and `g(s)`.
Undoing the common infinity shift contributes `s^4`, so the raw
two-endpoint product is bounded below by

```math
\boxed{\quad B(s)\ge g(s)^2.\quad}               \tag{14}
```

This is exactly the original scalar profile product at the same
physical gauges. The auxiliary fiber cannot reduce it. The statement
does not exclude other gauges, another row lattice, a deformation
outside this fiber, or a method that uses information beyond regular
channel pressure.

The common-type physical charge must still be restored when translating
this comparison into a crossing certificate. In the present normalization
its four letter energies are `G=(0,12,4,16)`. In particular, `g(s)^2`
alone is not a factoring rate: although it tends to one for large `s`,
the type entropy and physical charge also change. The same-type optimum
is described by

```math
\log B(s)=4\rho\log\!\left(\sum_{a=0}^3
s^{-G(a)/(4\rho)}\right).
```

Inequality (14) says that this auxiliary fiber cannot improve the
original scalar comparison under the same gauges and type transfer.

## 4. Exact checks

Run:

```sh
python3 tools/check_triangular_auxiliary_row_fiber.py
```

The standard-library checker constructs the eight-variable Gaussian
seed directly, checks all `256` endpoint entries, derives both Smith
changes, and verifies the `3,072`-variable kernel classification.
It does not depend on a coefficient table or another branch's helper.

An independent noncommuting control uses `E_01` on the direction
`(a,b,k)=(2,0,-1)` and `E_10` on `(1,2,-2)`. Their products do not
commute and yield a nontrivial auxiliary return path at length two.
At both endpoints the checker verifies the complete word operators
through length three, their polynomial inverses, and all `131,072`
entries of the two-site identity (10). The proof above covers every
dimension, word length, and choice of the fourteen coefficient matrices.

Evidence is saved in
[`evidence/triangular-auxiliary-row-fiber.json`](../evidence/triangular-auxiliary-row-fiber.json).

No improved factoring exponent is claimed.
