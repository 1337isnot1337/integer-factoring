# An all-field binary two-piece inequality and an F4 boundary pencil

This note proves a small rank inequality over every field, then uses it
to resolve a three-piece boundary-stable construction over
`F4=F2[omega]/(omega^2+omega+1)`. The latter has unequal terminal
projections, so the common-terminal reduction does not apply. Its
transition spectral base nevertheless has exact minimum `4`.

The statements concern this construction and the specified spectral
transfer. They do not exclude other boundary-stable constructions.

## 1. The binary two-piece inequality

For a tensor `T(x,y,u,v)` on four binary labels, let `r(T)` and `c(T)`
be its ranks in the cuts `(x,y)|(u,v)` and `(x,u)|(y,v)`. Write

```math
E_2(x,y,u,v)=[x=v][y=u],
\qquad \sigma(T)=\sqrt{r(T)c(T)}.
```

**Lemma.** Over any field, if `A+B=E_2`, then

```math
\boxed{\sigma(A)+\sigma(B)\ge4.}                 \tag{1}
```

**Proof.** If one tensor is zero, the other is the crossing and has
both ranks `4`. Otherwise all four ranks belong to `{1,2,3,4}`.
Rank subadditivity gives

```math
r(A)+r(B)\ge4,\qquad c(A)+c(B)\ge4.             \tag{2}
```

If `r(A)=1`, write `A(x,y,u,v)=U(x,y)V(u,v)`. Its other flattening
is `U tensor V`, hence `c(A)` is `1`, `2`, or `4`. The same statement
holds with the cuts or the two tensors exchanged. In particular,
none of the four rank pairs can be `(1,3)` or `(3,1)`.

The remaining integer possibilities satisfying (2) and having
`sqrt(r(A)c(A))+sqrt(r(B)c(B))<4` are exactly

```math
((1,2),(3,2)),\quad ((2,1),(2,3))
```

and their exchanges. This small enumeration can be checked without
floating point: for positive integers `p,q`, the condition
`sqrt(p)+sqrt(q)<4` is equivalent to
`p+q<16` and `4pq<(16-p-q)^2`.

It remains to rule out the first profile; symmetry handles the rest.
Its rank-one factorization has `rank(U)rank(V)=2`, so one matrix has
rank one and the other rank two. Let `K` be the permutation matrix
that swaps two binary tensor factors. Multiplying the column
flattenings by `K` turns their sum into the identity and gives

```math
T=K(U\otimes V),\qquad \operatorname{rank}T=2,
\qquad \operatorname{rank}(I-T)=2.               \tag{3}
```

The two images in (3) sum to the four-dimensional ambient space,
since every vector is `Tx+(I-T)x`. Their dimensions sum to four,
so they have zero intersection. The image of `T(I-T)=(I-T)T` lies
in both images. Consequently `T^2=T`.

On the other hand,

```math
T^2=(VU)\otimes(UV)
```

has rank at most one, because one of `U,V` has rank one. This
contradicts `rank(T)=2`. This argument uses no assumption on the
characteristic. `QED`.

### A weighted-crossing version

Let `B(y,u)` be an invertible `2 by 2` matrix. A change of basis on
the `y` leg sends `B` to the identity and preserves both flattening
ranks. Therefore any two-piece decomposition of

```math
H(x,y,u,v)=[x=v]B(y,u)
```

satisfies the same lower bound `4`.

When `rank(B)=1`, a related two-piece lower bound is `2`. If both
pieces are nonzero, each has two positive integer ranks and hence
cost at least one. If one piece is zero, the other is `H`, whose
two ranks are both `2`. Thus, for any binary `B`,

```math
H=A+C\quad\Longrightarrow\quad
\sigma(A)+\sigma(C)\ge2\operatorname{rank}(B).   \tag{4}
```

This includes the zero case. It is useful in characteristic two,
where the two pieces may be written `F` and `F+H`.

## 2. A three-piece unequal-terminal construction over F4

Put `omega^2+omega+1=0` and let

```math
P(y,u)=[y=u=0],\quad Q(y,u)=[y=u=1],
\quad (\beta_1,\beta_2,\beta_3)=(P,Q,P+\omega Q).
```

Use an auxiliary space with basis `eta` together with the four
pair vectors `f(x,v)`. Define

```math
a_1(x,v)=f(x,v)+[x=v]\eta,\qquad
a_2(x,v)=\omega f(x,v)+[x=v]\eta,\qquad
a_3(x,v)=f(x,v).                                  \tag{5}
```

Choose arbitrary scalar tensors `F_i(y,u,x,v)`. They define the
three covectors uniquely by

```math
b_i(y,u)^T\eta=\beta_i(y,u),\qquad
b_i(y,u)^T f(x,v)=F_i(y,u,x,v).                    \tag{6}
```

Set `A_i=a_i b_i^T`. Since
`beta_1+omega beta_2+beta_3=0` and
`beta_1+beta_2=[y=u]`,

```math
\left(\sum_i A_i(x,y,u,v)\right)\eta
  =[x=v][y=u]\eta.                               \tag{7}
```

Thus every choice of the three free tensors gives a valid common
right-eigenvector endpoint. With `phi(eta)=1` and `phi(f)=0`, the
initial boundary ranks are `(2,2,0)`. The terminal ranks are
`(1,1,2)`.

Define `D_i(y,u,x,v)=[x=v] beta_i(y,u)`. The three transition
tensors in row `i` are

```math
(T_{i1},T_{i2},T_{i3})
 =(F_i+D_i,\ \omega F_i+D_i,\ F_i).               \tag{8}
```

Their ranks use the two transition cuts
`(y,x)|(u,v)` and `(u,x)|(y,v)`. These are the ordinary two cuts
of Section 1 after permuting the four leg names. Write
`G_ij=sigma(T_ij)`.

The [diagonal rank-one word transfer](diagonal-rank-one-word-transfer.md)
constructs amplified crossing seeds with rank-sum bound
`s^T G^(k-1)t`. The question here is whether the live spectral
radius of `G` can be smaller than four.

## 3. Exact minimum of the transition base

**Theorem.** For every choice of the three tensors `F_i` over F4,
the live transition matrix has spectral radius at least `4`.
The choice `F_1=F_2=F_3=0` attains `4`.

**Proof.** Take the strictly positive vector `h=(1,1,2)`.
For `i=1,2`, each nonzero scalar multiple of `D_i` has both cut
ranks `2`, and therefore cost `2`. If all three tensors in (8)
are nonzero, their costs are at least one and

```math
G_{i1}+G_{i2}+2G_{i3}\ge4.                       \tag{9}
```

If one transition is zero, the other two are nonzero scalar
multiples of `D_i`: the three possible roots for `F_i` are
`0`, `D_i`, and `omega^(-1)D_i`, which are distinct. Their costs
are both two, so (9) still holds.

For the third row, `beta_3=diag(1,omega)` is invertible. Apply the
weighted-crossing version of (1) first to `F_3` and `F_3+D_3`,
and then to `omega F_3` and `omega F_3+D_3`. Scalar multiplication
does not change a nonzero tensor's ranks. Adding the two bounds
gives

```math
G_{31}+G_{32}+2G_{33}\ge8.                       \tag{10}
```

Equations (9)--(10) are exactly `Gh>=4h`. Iteration, or the
elementary positive-vector spectral bound, implies `spr(G)>=4`.

This also holds after removing states unreachable from the initial
boundary. Every such removed outgoing transition from a reachable
state is zero, so the restricted vector satisfies the same
inequality. Every terminal rank is positive, so each reachable
state can terminate; there is no additional coaccessibility
pruning that could invalidate the bound.

When every `F_i` is zero, state three is unreachable. The live
matrix on states one and two is

```math
\begin{pmatrix}2&2\\2&2\end{pmatrix},
```

whose spectral radius is four. `QED`.

The theorem remains valid over every extension field of F4: the
proof only needs an element satisfying `omega^2+omega+1=0` with
`omega` different from zero and one. It is not an inference from
a finite search.

## 4. Exact controls and the next open family

`check_binary_two_piece_f4_boundary.py` checks the four exceptional
integer profiles, the rank-one Kronecker identity, and the
projection contradiction. It also checks exact endpoint and
transition identities for arbitrary sampled F4 tensors, as well
as both cut ranks and the positive-vector lower certificate.
These finite controls support implementation correctness; the
universal claims are proved above.

A four-piece family with two independent free ket functions has
transitions

```math
(F_i+G_i+D_i,\ \omega F_i+\omega^2G_i+D_i,
  F_i,\ G_i),                                    \tag{11}
```

and terminal projections `(P,Q,P+omega Q,P+omega^2 Q)`.
The separate [four-piece pencil proof](four-piece-f4-boundary-pencil.md)
extends the lower-vector argument and gives exact spectral minimum four
for this family as well.
