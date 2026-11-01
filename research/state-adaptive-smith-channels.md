# Refactor the physical channels after every auxiliary lattice state

## Result and scope

A local channel factorization need not be fixed throughout an ordered
crossing. At each auxiliary row lattice, one can factor the complete
physical action again. The new channels may be linear combinations of
channels used at another state. A finite closed state graph, its exact
local identities, and positive rational potentials give an all-length
coefficient-rank bound.

This strengthens a specific existing certificate. For the binary companion
crossing from [exact rational module orbits](exact-rational-module-orbits.md),
at the same scale `s=2`, the new certificate has rate product

```math
B_{\mathrm{adaptive}}\le\frac{41}{3}
<\frac{96}{7}\le B_{\mathrm{fixed}}.
```

Here `B_fixed` is the exact grouped-word rate for the previously recorded
fixed local factorization. Its infinite-precision lattice orbit is already
closed. The improvement therefore cannot be obtained by merely increasing
the precision of that fixed factorization's module states.

The numerical bound remains worse than the elementary scalar certificate
`256/27`. It does not improve the factoring exponent. It demonstrates a
way to retain cancellation between channels before charging their branches,
which the preceding fixed-channel theorem explicitly left open. There is
no claim of optimality over other factorizations, state bases, evaluation
points, or the actual coefficient ranks.

## 1. A local identity that permits changing channels

Work over a fixed finite field `F`, with `O=F[[tau]]`. Let `K(z;a,b)` be
a polynomial auxiliary matrix of size `D` and degree at most `L`. The
indices `a,b` are the two physical halves of one cost cut, each containing
`q^2` possibilities. The original four labels still have alphabet size
`q`. Constant boundary vectors `f,e` specify the scalar ordered product.

At zero use `K_0(tau)=K(tau)`. At infinity use
`K_infinity(tau)=tau^L K(tau^-1)`. Both local matrices are regular. The
choice of cut and the physical half carrying both auxiliary indices is
fixed separately at each endpoint.

Let a state be a nonzero finitely generated row lattice
`M subset O^(1 by D)`, with a chosen basis matrix `B_M`. Its rank may be
smaller than `D`. Suppose that for every state there are finitely many
identities

```math
B_M K_e(\tau;a,b)
=\sum_j\tau^{\nu_{M,j}}
  H_{M,j}(\tau;a)B_{M_j}R_{M,j}(\tau;b),                 \tag{1}
```

where the orders `nu` are nonnegative integers, every entry of `H` and
`R` is regular, and `M_j` is another state. The `R` factors are scalar.
The number of rows of `H` is the rank of `M`, and the number of its columns
is the rank of `M_j`. No shared channel index between different states is
required.

Assume `f=c B_(M_start)` for a regular coordinate vector `c`. The final
boundary order at state `M` is

```math
\eta_M=\min_i\mathrm{ord}_\tau(B_Me)_i,             \tag{2}
```

with `eta=infinity` if every entry is zero.

Iterating (1) gives an exact decomposition over directed paths in the
state graph. For one path, the product of its `H` factors and the terminal
vector is a scalar function of the left physical word. The product of its
`R` factors is a scalar function of the right physical word. Its total
order is at least the sum of its transition orders plus (2).

The choice of the next factorization depends only on the path's state.
It never depends on the opposite physical word. Thus each path still gives
a separated left/right term. This is why the state-dependent choice does
not invalidate a matrix-rank bound.

## 2. Construct the identities from the full lattice action

Flatten `B_M K_e` into a matrix with rows indexed by
`(input lattice coordinate, left physical half, outgoing auxiliary index)`
and columns indexed by the right physical half. Perform minimum-valuation
Schur elimination over `F(tau)`. Every step has the form

```math
A=\sum_j\tau^{\nu_j}\ell_jr_j,                         \tag{3}
```

where `ell_j` and `r_j` have regular entries and each contains a unit.
A zero residual stops exactly, so its generic rank can be smaller than
`q^2`. The implementation verifies every entry of (3) after the complete
elimination.

Reshape `ell_j` into rows of length `D`, one for each input lattice
coordinate and left physical half. Let `M_j` be their exact `O`-span.
Compute a row Hermite basis `B_(M_j)` and express each reshaped row in
that basis. The resulting coordinate matrices are `H_(M,j)`, and the
scalar vector `r_j` gives `R_(M,j)`. The coordinates are regular because
the rows belong to the generated lattice. Exact span witnesses verify
both the Hermite basis and these coordinate identities.

Some component of `ell_j` is a unit, so every target lattice has order
zero. No common scalar valuation is lost when forming the next state.
This proves (1) and gives a deterministic construction of all its data.

One may instead replace each target by the larger lattice

```math
M_j+\tau^P O^D.                                         \tag{4}
```

The same rows still have regular coordinates in this larger lattice, so
(1) remains an exact identity. Starting with `Of+tau^P O^D` is also valid.
There are finitely many such full-tail states over a finite field. This
provides a finite sufficient certificate at every fixed precision `P`.
It is distinct from silently adding a tail to an exact state.

Changing precision can also change the chosen local factorization.
We do not assert monotonicity or completeness of this new precision
procedure. Exact mode has no tail and can have an infinite state orbit;
a resource-budget failure is reported as inconclusive.

## 3. The all-length rank bound

Fix a rational `s>1`. At zero define the nonnegative transition matrix

```math
(T_0)_{M,N}=\sum_{j:M_j=N}s^{-\nu_{M,j}},
```

and at infinity define

```math
(T_\infty)_{M,N}=\sum_{j:M_j=N}s^{L-\nu_{M,j}}.           \tag{5}
```

The terminal vector in either graph is `r_M=s^(-eta_M)`, with
`s^(-infinity)=0`. Let

```math
Q_e(t)=e_{M_{\mathrm{start}}}^{T}T_e^t r_e.
```

Suppose the scalar ordered crossing is

```math
fK(z;x_1,y_1,u_1,v_1)\cdots K(z;x_t,y_t,u_t,v_t)e
=\sum_{k=0}^{Lt}z^k\Gamma_{t,k}.
```

For opposite cost cuts at the two endpoints, (1) gives

```math
\begin{aligned}
\mathrm{rank}_R\Gamma_{t,k}
 &\le (Lt-k+1)s^{-k}Q_\infty(t),\\
\mathrm{rank}_C\Gamma_{t,k}
 &\le (k+1)s^k Q_0(t).                                  \tag{6}
\end{aligned}
```

For example, at zero a path of total order `v` contributes nothing when
`v>k`. Otherwise its coefficient of order `k` is a sum of at most
`k-v+1` outer products, obtained by splitting the Taylor degree between
its scalar left and right functions. The number of eligible paths is
at most `s^k` times their total `s^(-v)` weight. Summing proves the
second inequality. The first applies the same argument to the coefficient
of order `Lt-k` at infinity; its factor `s^(Lt)` is included in (5).
Rational factors may have poles away from the chosen endpoint. Their
Taylor series are regular at that endpoint, which is all (6) needs.

If positive rational vectors satisfy

```math
T_eh_e\le\beta_e h_e,
```

then `Q_e(t)<=C_e beta_e^t`, where one explicit choice is

```math
C_e=(h_e)_{M_{\mathrm{start}}}\max_M\frac{(r_e)_M}{(h_e)_M}.
```

Consequently the sum of square-root products of coefficient ranks is at
most a polynomial in `t` times
`(beta_infinity beta_0)^(t/2)`. This is the existing crossing transfer
with a new way to certify its complete coefficient ranks. A rate product
strictly below `q^(2 sqrt(2))` would pass the half-rate threshold.

For fixed finite certificate data, the construction remains effective.
There are `O(t)` coefficient pieces. Direct polynomial products construct
their entries in polynomial time per fixed physical word. Their full
tables and rank factorizations can therefore be constructed in
`q^(O(t))=2^(O(lambda))` field operations. The field, auxiliary dimension,
state graph, degree, and rational scale are constants. No input-dependent
advice or unexplained channel oracle is introduced.

## 4. A strict gain on the companion crossing

Use `F=F_2` and binary physical labels. Write

```math
C=[x=u][y=v],\qquad E=[x=v][y=u],\qquad N=C+E.
```

The companion local matrix before a constant change of auxiliary basis is

```math
K(z)=\begin{pmatrix}zN&C\\zE&0\end{pmatrix},\qquad
f=(1,0),\qquad e=(1,1)^T.
```

It satisfies `K(1)e=Ee` and `fe=1`. Conjugation by
`S=((1,0),(1,1))`, for which `S^2=I`, changes both boundaries to the
first coordinate. The checker uses the same conjugated matrix as the
preceding exact-orbit certificate, without changing the physical family.

The exact infinity states, in their generated order, have row bases

```math
\begin{gathered}
(1,0),\quad
\begin{pmatrix}1&0\\0&\tau\end{pmatrix},\quad
(1,1),\quad
\begin{pmatrix}\tau&0\\0&1\end{pmatrix},\quad
I_2,\quad
\begin{pmatrix}1&\tau\\0&\tau^2\end{pmatrix}.
\end{gathered}
```

The exact zero states have row bases

```math
\begin{gathered}
(1,0),\quad
\begin{pmatrix}1&1\\0&\tau\end{pmatrix},\quad
\begin{pmatrix}\tau&0\\0&1\end{pmatrix},\quad
I_2,\quad
\begin{pmatrix}1&1+\tau\\0&\tau^2\end{pmatrix}.
\end{gathered}
```

All transitions are reconstructed from (3). At `s=2` their matrices are

```math
T_\infty=
\begin{pmatrix}
0&2&2&1&0&0\\
0&2&0&0&3&0\\
8&0&0&0&0&0\\
0&0&0&0&0&8\\
0&8&0&0&0&0\\
0&2&0&3&0&0
\end{pmatrix},\qquad
T_0=
\begin{pmatrix}
1&1&1/2&0&0\\
1&0&0&1&0\\
3/4&0&0&0&1\\
1&1&0&1/2&0\\
1&0&0&1&0
\end{pmatrix}.                                          \tag{7}
```

The positive rational potentials

```math
h_\infty=(3,3,4,4,4,3)^T,\qquad
h_0=(1000,911,730,1075,911)^T
```

satisfy exactly

```math
T_\infty h_\infty=6h_\infty,\qquad
T_0h_0\le\frac{41}{18}h_0.                               \tag{8}
```

The terminal orders are respectively `(0,0,0,1,0,0)` and
`(0,0,1,0,0)`. Thus one can take
`C_infinity=1`, `C_0=1000/911`, and

```math
B_{\mathrm{adaptive}}\le6\cdot\frac{41}{18}=\frac{41}{3}.
```

## 5. The old exact word rate is strictly larger at the same point

Rebuild the fixed local factors and their exact normalized lattice orbits
from the preceding note. Their two-state infinity matrix at `s=2` is
`((0,8),(3,2))`, with positive eigenvector `(4,3)` and eigenvalue six.
The four-state quotient of its seven-state zero graph is

```math
W_0=
\begin{pmatrix}
1/2&1&1&0\\
1/2&0&3/2&0\\
0&1&1&1/2\\
0&3/2&1/2&0
\end{pmatrix}.                                          \tag{9}
```

With `g=(1000,840,947,758)^T`, exact multiplication gives

```math
W_0g\ge\frac{16}{7}g.                                   \tag{10}
```

Every component of `g` is positive. The initial state is reachable and
the free-boundary terminal vector is positive, so iterating (10) proves
an exponential rate at least `16/7`. The preceding
[minimal-boundary theorem](minimal-boundary-module-completeness.md) identifies
this exact free rate with the fixed-boundary grouped-word limsup rate for
these minimized factors. Hence

```math
B_{\mathrm{fixed}}\ge6\cdot\frac{16}{7}
=\frac{96}{7}>\frac{41}{3}.                              \tag{11}
```

Equations (8)--(11) compare the same crossing, point, physical cuts, and
auxiliary side. The fixed factorization is the one already recorded in
the repository. The new method changes its channel decomposition at later
states. The comparison is not a claim that every fixed factorization is
worse, nor that the adaptive certificate dominates it at every scale.

## 6. Reproduction and remaining question

Run:

```sh
python3 tools/check_adaptive_smith_lattices.py \
  --output evidence/adaptive-smith-lattices.json
```

The standard-library checker verifies exact closure, every local rational
factorization, regular coordinates in every target lattice, the endpoint,
both upper potentials, and the independent old lower potential. It rebuilds
the old exact module orbits and their quotient. Literal polynomial products
check all 4,368 endpoint entries through length three, and independently
computed coefficient ranks satisfy 81 exact bounds, including the full-tail
controls. Scalar and auxiliary-scalar cases reproduce the elementary
profiles. A deliberately insufficient state budget is rejected.

Files:

- [Adaptive factorization engine](../tools/adaptive_smith_lattices.py)
- [Exact checker](../tools/check_adaptive_smith_lattices.py)
- [Local identities and potentials](../evidence/adaptive-smith-lattices.json)

The next question is constructive: can a different ordered crossing and a
finite graph of identities (1) achieve `B<q^(2 sqrt(2))`? The present result
does not answer it. It does show that a completed fixed-channel module
calculation need not be the final bound on the same local crossing.

The polynomial crossing, rational arithmetic, exact DVR Hermite algorithm,
and companion example come from the linked prior notes and their credited
sources. The state-dependent local identities, transfer proof, implementation,
and strict comparison here used substantial OpenAI assistant support. No
Lean formalization accompanies this result.
