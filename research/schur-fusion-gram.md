# The homogeneous Schur central piece as a fusion Gram matrix

## Result and scope

The full homogeneous Schur construction has an exact central-piece
description in characteristic two, and its two central flattenings
have the same rank. Write

```math
\mu_h:\Lambda_h\otimes\Lambda_h\longrightarrow\Lambda_{2h},
\qquad \mu_h(a,b)=a[h]b,
\qquad q_h=\operatorname{rank}\mu_h.
```

Here `a[r]` is the source's weight-shifting mode, so the displayed map
is the usual vertex-algebra product with field exponent zero. If
`h` is a power of two, the central row matrix `A_0` of the full
homogeneous Schur family satisfies

```math
\max\{0,2q_h-p(2h)-1\}
\ \le\ \operatorname{rank}A_0
\ \le\ \min\{p(h)^2,q_h+1\}.
\tag{1}
```

More precisely, invertible changes of the two label pairs turn its
non-vacuum part into the Hall Gram matrix of the fusion map. The other
source operator order identifies the column flattening with the
transpose of an invertible transform of the row flattening. This
identification is integral before reduction. It uses a nondegenerate
bilinear form over `F_2`; it does not use a positivity argument there.

Exact finite calculations give:

| `h` | `p(h)` | `p(2h)` | Fusion rank `q_h` | Fusion Gram rank | Central rank |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 2 | 1 | 0 | 1 |
| 2 | 2 | 5 | 3 | 1 | 2 |
| 4 | 5 | 22 | 16 | 12 | 12 |
| 8 | 22 | 231 | 229 | 227 | 228 |

The asymptotic size of `q_h` along powers of two remains open in this
note. The table does not prove a new factoring exponent or a barrier
for all crossing constructions. It identifies a specific algebraic
question governing the central piece of this particular full-label
family.

The mode definitions and crossing identity are from
[Appendix B, pinned revision `16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The integral lattice, divided translations, and prime-power duality
are proved in [the integral Schur construction](integral-schur-crossing.md).
The rational Gram operator is developed in
[the homogeneous-dual construction](homogeneous-duals.md).
This note proves the adjoint and fusion reduction for those exact modes.

## 1. An integral adjoint involution

Use the Hall form, with Schur functions orthonormal, and set

```math
\mathcal D(p_1)=0,\quad \mathcal D(p_m)=mp_{m-1}\ (m\ge2),
\qquad \mathcal U(p_m)=mp_{m+1},
```

```math
T=e^{\mathcal D},\qquad
\omega(p_m)=(-1)^{m-1}p_m,\qquad J=T\omega.
```

Both `T` and `omega` preserve the integral Schur lattice. Since
`omega D omega=-D`, the map `J` is an integral involution. In particular,

```math
Jp_m=(-1)^{m-1}\sum_{j=1}^{m}\binom mjp_j,
\qquad J^2=1.
\tag{2}
```

For every polynomial `a` and every integer mode `r`,

```math
\boxed{\ (a[r])^*=(Ja)[-r]\ }.
\tag{3}
```

The star means the adjoint for the Hall form. To prove (3), multiplication
by `p_k` has adjoint `k partial/partial p_k`. The source's creation and
annihilation coefficients for the parameter `p_m` are respectively

```math
\binom{k-1}{m-1},\qquad
(-1)^{m-1}k\binom{k+m-1}{m-1}.
```

Substitute (2) and use the two binomial identities

```math
\sum_{j=1}^{m}\binom mj\binom{k-1}{j-1}
=\binom{k+m-1}{m-1},
```

```math
\sum_{j=1}^{m}(-1)^{m+j}\binom mj\binom{k+j-1}{j-1}
=\binom{k-1}{m-1}.
```

They give (3) for one operator parameter. A source monomial is a
normally ordered product, with all multiplication operators before
all derivatives. Taking its adjoint reverses order and exchanges these
two commuting groups, leaving the result normally ordered. Consequently
the one-parameter identities multiply without additional contractions.
Linearity proves (3).

All entries are integral in Schur coordinates. Thus (3) is an identity
of integral pairings and remains true after reduction modulo any prime.

## 2. Homogenizing the adjoint at the fusion mode

For a homogeneous parameter `c` of weight `w`, the source coefficients
also give the translation identity

```math
(\mathcal Uc)[r]=(r-w)c[r].
```

For one power sum this is the elementary binomial recurrence. For a
product, the derivation rule sums `r_j-w_j` over the individual modes.
Iteration yields

```math
(\mathcal U^{(k)}c)[r]
=\binom{r-w}{k}c[r],
\qquad \mathcal U^{(k)}=\mathcal U^k/k!.
\tag{4}
```

The divided powers preserve the integral lattice. Formula (4) is first
an identity over the rationals and then an identity of integral maps;
it does not instruct us to invert `k!` in characteristic two.

On weight `h`, define

```math
G_h=\sum_{k\ge0}\mathcal U^{(k)}\mathcal D^{(k)},
\qquad B_h=G_h\omega=\omega G_h.
\tag{5}
```

The sum is finite. The derivations are adjoint, so `G_h` is self-adjoint.
It commutes with `omega`, because each summand has equally many upward
and downward translations. The homogeneous vacuum pairing is

```math
[1]a[-h]b=\langle B_ha,b\rangle
\qquad(a,b\in\Lambda_h).
\tag{6}
```

For `a` of weight `h`, the polynomial `D^(k) omega a` has weight `h-k`.
Applying (4) with `r=h` gives coefficient `binom(k,k)=1`. Therefore

```math
(Ja)[h]
=\sum_k(\mathcal D^{(k)}\omega a)[h]
=\sum_k(\mathcal U^{(k)}\mathcal D^{(k)}\omega a)[h]
=(B_ha)[h].
\tag{7}
```

Combining (3) and (7) gives the needed homogeneous adjoint:

```math
\boxed{\ (a[-h])^*=(B_ha)[h]\ },
\qquad a\in\Lambda_h.
\tag{8}
```

The matrix of `B_h` in the Schur basis is symmetric and integral. At
`h=2^r`, it is invertible over `F_2`, by the prime-power determinant
argument in the integral Schur construction. No diagonalization over
`F_2` is required.

## 3. The exact central matrix

Now work over `F_2` and let `h=2^r`, including `r=0`. Index all Schur
functions of weight `h` by `x=1,...,m`, where `m=p(h)`. Take

```math
b_v=s_v,\qquad a_x=B_h^{-1}s_x.
```

By (6), `a_x[-h]s_y=[x=y]`. In the exchanged order of the source's
crossing identity, the homogeneous angle coefficients are `binom(2h,i)`.
Modulo two only `i=0,2h` remain. The central operator is thus

```math
Q'_0(a_x,b_v)=b_v[h]a_x[-h]+b_v[-h]a_x[h].
\tag{9}
```

The first summand, applied to `s_y` and evaluated at coordinate `s_u`,
is `[x=y][v=u]`. Here `b_v[h]1=b_v`, directly from the source's creation
rule. For the second summand, (8) gives

```math
\begin{aligned}
\langle s_u,b_v[-h]a_x[h]s_y\rangle
&=\langle (B_hs_v)[h]s_u,a_x[h]s_y\rangle\\
&=\langle\mu_h(B_hs_v,s_u),
             \mu_h(B_h^{-1}s_x,s_y)\rangle.
\end{aligned}
\tag{10}
```

Let `F_h` have rows indexed by Schur functions of weight `2h` and columns
indexed by `(x,y)`, with that column equal to `mu_h(s_x,s_y)`. Put

```math
K_h=F_h^{\mathsf T}F_h,\quad
L=B_h^{-1}\otimes I_m,\quad R=B_h\otimes I_m,
\quad v_{(x,y)}=[x=y].
```

With rows `(x,y)` and columns `(v,u)`, equations (9)--(10) are exactly

```math
\boxed{\ A_0=L K_hR^{\mathsf T}+vv^{\mathsf T}\ }.
\tag{11}
```

The full weight-`2h` Hall form is nondegenerate over `F_2`, since its
Schur matrix is the identity. If `W=im F_h`, then the radical of its
restricted form is `W intersect W^perp`, of dimension at most
`p(2h)-dim W`. Hence

```math
2q_h-p(2h)\le\operatorname{rank}K_h\le q_h.
```

The invertible label maps in (11) preserve rank, and the remaining
matrix has rank one. This proves (1). The same conclusion also follows
from Sylvester's matrix rank inequality; no real inner product is being
transported to characteristic two.

The original operator order gives the additional identity

```math
Q_0(a_x,b_v)=a_x[-h]b_v[h]+a_x[h]b_v[-h].
```

Its first summand is
`<mu_h(s_x,s_u),mu_h(s_v,s_y)>`, since `B_h a_x=s_x`.
Its second summand is `(B_h^{-1})_(x,u)(B_h)_(v,y)`.
Define `vec(C)_(x,u)=C_(x,u)` for a symmetric matrix `C`. With
column-flattening rows `(x,u)` and columns `(v,y)`, we obtain

```math
C_0=K_h+\operatorname{vec}(B_h^{-1})
              \operatorname{vec}(B_h)^{\mathsf T}.
\tag{12}
```

On the other hand, multiplying (11) by the inverse label maps gives

```math
L^{-1}A_0(R^{\mathsf T})^{-1}
=K_h+\operatorname{vec}(B_h)
             \operatorname{vec}(B_h^{-1})^{\mathsf T}
=C_0^{\mathsf T}.
\tag{13}
```

Here `B_h`, its inverse, and `K_h` are symmetric. Thus the central
row and column ranks agree exactly, and (1) bounds both of them.

## 4. What an asymptotic image theorem would imply

For example, if an independent argument proves

```math
q_h=(1-o(1))p(2h)\qquad(h=2^r\longrightarrow\infty),
\tag{14}
```

then (1) implies `rank A_0=(1-o(1))p(2h)`. Its logarithm would have
the ambient leading term `sqrt(2) log2 p(h)`. This would show that
the central piece of the full homogeneous Schur family retains the
ambient exponent in both flattenings. Their rank product would have
logarithm `2sqrt(2) log2 p(h)+o(log p(h))`. Thus this specific full-label
construction would require `rho>=1/sqrt(2)`, even when actual matrix
ranks replace declared channel counts. Equation (14) is an explicit
open input, not a consequence of the four finite cases above. Label
subspaces or different piece constructions require separate analysis.

Ordinary multiplication cannot simply replace fusion in (14). Already
at `h=2`, their images are different. In decreasing-part notation,
over `F_2`,

```math
\begin{aligned}
\mu_2(s_{11},s_{11})&=s_{1111}+s_{31},\\
\mu_2(s_{11},s_2)=\mu_2(s_2,s_{11})&=s_{22},\\
\mu_2(s_2,s_2)&=s_{211}+s_4.
\end{aligned}
```

The functional `[s_1111]+[s_31]` annihilates all these vectors, but is
one on the ordinary product

```math
s_{11}^2=s_{1111}+s_{211}+s_{22}.
```

Thus `Lambda_h Lambda_h` is not contained in `im mu_h`, even at a
prime-power weight. A characteristic-zero triangular argument in the
power-sum basis does not repair this: that basis is not an integral
Schur basis, and relations such as `p_1^2=p_2` appear modulo two.
An integral image comparison would require a separate proof.

## 5. Reproducible checks

Run the quick suite with

```sh
python3 tools/check_schur_fusion_gram.py
```

To reproduce the recorded table including `h=8`, run

```sh
python3 tools/check_schur_fusion_gram.py --weight-eight \
  --output evidence/schur-fusion-gram.json
```

The checker uses exact integer arithmetic throughout. It independently
computes symmetric-group characters by beta-set rim-hook removal,
checks their integral orthogonality, and compares the small character
tables with the separate Jacobi--Trudi implementation. To compute
fusion, it retains the common denominator `(h!)^2` inside the modulus
`2(h!)^2`, applies the source's integer operator recursion, converts
to Schur coordinates, verifies divisibility, and only then reduces
modulo two. No even denominator is inverted.

It checks the adjoint, the homogenization identity, the vacuum pairing,
and the central formula against the independent complete-piece engine
through `h=4`. It verifies the full row-to-column reshuffling identity
at every reported weight. It also checks the `h=2` multiplication counterexample
and a sign-sensitive adjoint negative control. The optional `h=8`
fusion calculation recovers central rank `228`, independently matching
the earlier direct central-piece calculation recorded in
[the Schur evidence](../evidence/integral-schur-crossing.json).

These finite calculations verify identities and implementations at
the displayed sizes. The proof above establishes (1), (3), (8),
(11), and (13) for all stated weights. Neither part establishes (14).
