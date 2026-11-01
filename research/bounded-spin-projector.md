# Constant-state spin projectors and exact cut reciprocity

## Uniform results

The [coherent spin encoders](bounded-spin-crossing.md) admit an exact
orthogonal projector with **thirty auxiliary states**, independent of
length. Restricting to exactly `2^lambda` labels needs sixty states.
These are spatial matrix-product representations of the projector, not
bounds on the two crossing cut ranks used by the factoring evaluator.

There is also an exact symmetry of the
[six-vertex coefficient tensors](bounded-spin-six-vertex.md). After the
specified column normalization, their two cut operators satisfy

```math
\boxed{R(T_{A,B,C})=C(T_{C,B,A}).}                    \tag{1}
```

For the original pencil, writing `C_j=[z^j]G_n(z)`, this becomes

```math
\boxed{R(C_j)=C(C_{n-j}).}                           \tag{2}
```

Thus one cut operator space and an explicit coefficient involution suffice
for both rank queries. These identities hold uniformly; no favorable
rank-growth estimate follows from them alone.

## 1. A thirty-state projector

Put the path vectors into the columns of `F`, and let `Gamma` be their
diagonal Gram matrix. Orthogonality gives

```math
P_n=F\Gamma^{-1}F^T
   =\sum_P\frac{v_Pv_P^T}{\nu_P},\qquad
P_n^2=P_n=P_n^T.                                    \tag{3}
```

All norms are products of 2, 6 and 12, so (3) is defined over
`Z[1/6]`, and hence over every field of characteristic different from
two and three. Its image is freely spanned by the selected path vectors.
Over a field its ordinary matrix rank is the integer number of paths,
even if that number reduces to zero in the field.

To compute the entry `P_n(x,y)`, sum over the encoder path while reading
the two binary words together. A state is

```math
(h,r,s),\qquad 0\le h\le3,\quad0\le r,s\le h.
```

There are `1^2+2^2+3^2+4^2=30` such states. A local transition chooses
the path step up or down, applies the two magnetic transitions of the
original encoder, and multiplies their scalars. At a downward step
from height `h`, also divide by `h(h+1)`. The product of those divisions
is exactly `1/nu_P`.

Both boundary states are `(0,0,0)`. Summing the transition products is
therefore exactly (3). Equivalently there are four fixed thirty-by-thirty
matrices `M_ab`, one per pair of input bits, with

```math
P_n(x,y)=e_0^T M_{x_0y_0}\cdots M_{x_{n-1}y_{n-1}}e_0.       \tag{4}
```

For a fixed word pair, the forward contraction is linear in `n` with a
constant state space. The implementation uses sparse transitions and
does not enumerate the exponentially many paths.

### Exactly a power-of-two number of labels

The selected paths are the first `2^lambda` in up-first order. Compute
their last path `Q` by the existing polynomial-time unranking routine.
While reading a candidate path, retain whether its prefix equals `Q`
or is earlier. A later prefix is rejected. There are only two surviving
comparison states, so the projector needs at most sixty auxiliary states.
The right boundary accepts either comparison state. This selects exactly
the desired initial segment, including its final path.

The cutoff makes the local matrices depend on the known step `Q_i`, but
their size remains constant. It does not require one auxiliary state per
physical label. [`bounded_spin_projector.py`](../tools/bounded_spin_projector.py)
implements both versions over `F_5`.

## 2. An ambient representation preserving both cut ranks

Let `K_n` be any local four-leg kernel tensor power, and let the two
label frames coincide. In the row arrangement, the encoded operator is

```math
G=(F^T\otimes F^T)K_n
            (F\Gamma^{-1}\otimes F\Gamma^{-1}).
```

Its ambient extension is

```math
\widehat G=(P_n\otimes P_n)K_n(P_n\otimes P_n)
 =(F\Gamma^{-1}\otimes F\Gamma^{-1})
       G(F^T\otimes F^T).                            \tag{5}
```

Every external leg map in the final expression has full column or row
rank, as appropriate. Hence (5) preserves the ordinary rank in **both**
crossing cuts: in either arrangement it multiplies the encoded matrix
by an injection and a surjection. The constant-state projectors therefore
give an exact ambient tensor-network model for the same rank problem.

At the crossing endpoint the ambient tensor is `(P_n tensor P_n)S`,
rather than a full crossing on all binary words. Its action is the
crossing on the selected image and zero on the complementary directions.
Consequently, full-alphabet auxiliary-kernel restrictions cannot be
invoked just from the small auxiliary dimension in (4).

## 3. Proof of cut reciprocity

Use `epsilon=[[0,1],[-1,0]]`. Every return-to-zero path vector is fixed
by `epsilon^tensor n`, as proved in the
[spin-zero invariance argument](bounded-spin-three-jets.md).
On a local tensor, exchange the `y` and `u` legs and apply epsilon on
both exchanged legs. Its six weights transform as

```math
(a,f,b,e,c,d)\longmapsto(c,d,-b,-e,a,f).
```

Thus the balanced products change by `(U,V,W)->(W,V,U)`. The two minus
signs cancel because the numbers of `b` and `e` vertices agree. Since
the epsilon action fixes each encoder, the unnormalized coefficient
tensor obeys

```math
T^{\rm raw}_{A,B,C}(P,Q,R,S)
   =T^{\rm raw}_{C,B,A}(P,R,Q,S).                    \tag{6}
```

For the row cut, arrange columns in `S,R` order and divide by
`nu_S nu_R`. For the column cut, use row order `P,R`, column order
`S,Q`, and divide by `nu_S nu_Q`, as in the existing rank checkers.
After substituting the exchanged labels in (6), the two matrices have
identical entries and identical column denominators. This proves (1),
including an operator identity rather than just equality of ranks.

For the original specialization `U=1,V=(1-z)^2,W=z^2`, exchanging
`U,W` is `z^n G_n(1/z)` at the coefficient level: each balanced
monomial has total degree `n/2`. Equating coefficients proves (2).

The induced tensor involution fixes the crossing and exchanges its
two cut ranks. Any favorable decomposition may therefore be combined
with its image, each scaled by one half. This yields a decomposition
closed under the involution, with the same maximum rank product and
only twice the rank sum and piece count. This constant loss does not
change a strict asymptotic saving. Individual pieces need not be fixed
by the involution; forcing that stronger condition could increase ranks.

## 4. What this does and does not resolve

The projector representation permits exact composition of encoded
operators by inserting another copy of `P_n tensor P_n` between their
local tensor powers. For any fixed number of factors, the resulting
auxiliary dimension is bounded independently of `n`. It supplies a
concrete starting point for transfer or recoupling identities without
assuming a categorical quotient interpretation of the encoders.

The thirty-state bond refers to a cut between **positions** in (4).
The ordinary matrix rank separating the two entire binary words is
exponential, equal to the path count. The two crossing cut ranks are
likewise different partitions. A small bond alone does not bound them,
and operator words of growing length need growing composition state
spaces. Uniform favorable ordinary ranks remain the missing target.

## Checks

Run `make verify-bounded-spin-projector` (NumPy required only for the
checker). It compares all 139,793 word-pair entries in the stated full
and binary-selected cases with independent integral path vectors,
checks idempotence, symmetry and exact ranks, and adds 128 ten-site
binary-cutoff controls. Omitting norm division is a negative control:
the raw two-site sum is not idempotent.

`make verify-bounded-spin-six-vertex` now checks (1) on every coefficient
matrix at lengths four, six and eight, and (2) on all nine original
eight-site coefficient matrices. Receipts:
[projector](../evidence/bounded-spin-projector.json) and
[six-vertex](../evidence/bounded-spin-six-vertex.json).
The proofs and implementations used substantial Codex AI assistance.
No improved factoring exponent is asserted.
