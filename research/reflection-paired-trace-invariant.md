# A uniform paired-trace identity for reflection crossings

## Result

Let `M>=4` be a power of two, `m=M/2`, and `q=2^m`. Let
`T_M(w)` be the complete polynomial reflection crossing from
[`reflection_dual_rail_gaussian.py`](../tools/reflection_dual_rail_gaussian.py).
Over every field of characteristic two,

```math
\sum_{x=0}^{q-1}T_{x y u x}(w)
   =(1+w^{M^2})\delta_{y,u},\qquad
\sum_{y=0}^{q-1}T_{x y y v}(w)
   =(1+w^{M^2})\delta_{x,v}.                         \tag{1}
```

Both identities survive arbitrary paired rational physical frames
`(A_x,B_y,B^{-T}_u,A^{-T}_v)`, whenever the resulting tensor is
polynomial. They also survive every affine combination of such tensors
whose coefficients sum to one. Consequently every polynomial strict
crossing in this orbit-affine family has:

- scalar polynomial content exactly one;
- degree at least `M^2`.

The degree lower bound is below the original degree `3M^2/2`. It does
not exclude degree reductions by more general affine combinations.
It excludes removing a common scalar factor, including factors over
extension fields. This is a structural restriction on this stated
family, not a lower bound for arbitrary crossing constructions or
integer factoring.

## 1. The complete Laurent Gaussian

Use commuting square-zero variables `x_i,y_i,u_i,v_i`, with
`0<=i<m`. In characteristic two the Gaussian associated with an
alternating kernel `K` can be defined without division:

```math
\Gamma(K)=\prod_{a<b}(1+K_{ab}z_a z_b).
```

Its square-free coefficient on a set `S` is the Pfaffian of `K_S`.
The complete reflection tensor is the coefficient tensor of
`w^(M^2/2) Gamma(K)`, where `K` is the Laurent kernel in the linked
implementation. The regular gauges used there are removed exactly;
they do not change this coefficient identity.

Put

```math
t=w^4,\quad s=1+w^2,\quad
a=1+w^{-2M},\quad b=1+w^{2M},\quad J_{ij}=[i+j=m-1].
```

The diagonal blocks of `K`, in `x,y,u,v` order, are
`aJ,aJ,bJ,bJ`. Its off-diagonal blocks are

```math
K_{xy}=w^{-2M}sG,\quad K_{uv}=sH,\quad
K_{xu}=sL^T,\quad K_{yv}=sL^T,\quad
K_{xv}=K_{yu}=I,                                  \tag{2}
```

with transposes on the reverse blocks. The entries are

```math
G_{ij}=t^{i\mathbin{\&}j},\qquad
H_{ij}=t^{m-1-(i\mathbin{|}j)},\qquad
L_{ij}=t^{j\mathbin{\&}\neg i}.                    \tag{3}
```

In (3), the complement is restricted to the `log_2(m)` index bits.
The matrices are tensor products of the one-bit matrices

```math
G(z)=\begin{pmatrix}1&1\\1&z\end{pmatrix},\quad
H(z)=\begin{pmatrix}z&1\\1&1\end{pmatrix},\quad
L(z)=\begin{pmatrix}1&z\\1&1\end{pmatrix},
```

using `z=t^(2^r)` on bit `r`. Direct two-by-two multiplication gives,
with `d=(1+t)^(m-1)`,

```math
GJG=HJH=LJL^T=L^TJL=dJ,
\qquad GJL^T=L^TJH=dI.                            \tag{4}
```

Tensor products prove (4) at every size. Also

```math
s^2d=b,\qquad a=w^{-2M}b.                         \tag{5}
```

## 2. Closing a paired alphabet

Multiply `Gamma(K)` by `prod_i(1+x_i v_i)` and take the coefficient
of all the `x_i` and `v_i`. A contributing monomial originally uses
the same subset of the `x` and `v` variables: the extra product supplies
the complementary matching pairs. Thus this operation is exactly
the sum over equal physical letters `x=v` in (1).

Multiplication adds the identity matching to the `x,v` kernel block.
The existing block in (2) is also identity, so it cancels. The block
on the eliminated variables is now

```math
K_{xv,xv}=\operatorname{diag}(aJ,bJ).
```

Work first over the rational-function field. Since `m` is even, `J`
is alternating and invertible. Gaussian elimination gives the
Pfaffian of this eliminated block times the Gaussian on its Schur
complement. One may verify this coefficientwise: the block triangular
congruence has determinant one and turns each relevant principal
submatrix into the direct sum of the eliminated block and its Schur
complement.

On the remaining variables `y,u`, formulas (2)--(5) give

```math
\begin{aligned}
S_{yy}
 &=\left(a+\frac{w^{-4M}s^2d}{a}+\frac{s^2d}{b}\right)J
   =\left(a+w^{-2M}+1\right)J=0,\\
S_{uu}
 &=\left(b+\frac{s^2d}{a}+\frac{s^2d}{b}\right)J
   =\left(b+w^{2M}+1\right)J=0,\\
S_{yu}
 &=I+\left(\frac{w^{-2M}s^2d}{a}+\frac{s^2d}{b}\right)I
   =I.
\end{aligned}                                                    \tag{6}
```

This is precisely the identity-pairing Gaussian `prod_i(1+y_i u_i)`.
The scalar, including the original common monomial, is

```math
w^{M^2/2}\operatorname{Pf}(aJ)\operatorname{Pf}(bJ)
 =w^{M^2/2}(ab)^{m/2}
 =b^m
 =1+w^{M^2}.                                      \tag{7}
```

Equations (6)--(7) prove the first identity in (1). It is a polynomial
identity, so the rational-function derivation imposes no restriction
on the value substituted for `w`. The Laurent kernel is invariant
under simultaneously exchanging `x,y` and `u,v`. This proves the
second identity as well.

## 3. Preservation under physical frames

Apply arbitrary invertible matrices `A(w),B(w)` on the physical
letters, using the paired action

```math
T'=(A_x,B_y,B^{-T}_u,A^{-T}_v)T.
```

In the `x=v` contraction, the sum over the output `x` contracts
`A^T A^{-T}=I`. On the remaining matrix the action is conjugation
by `B`. The right side of (1) is a scalar identity, so it is fixed.
The other contraction similarly conjugates by `A`.

This argument is over the rational-function field and therefore
permits poles in the frames. If the final tensor is polynomial, the
same polynomial identity holds. Constant basis changes, complete
polynomial frames, the coupled pole-cancelling families, and the
independent Laurent-series families are all included.

Linearity preserves (1) under `sum_i c_i T_i` whenever
`sum_i c_i=1`. In particular it includes the affine strict paths used
by the current searches. No individual intermediate regularity is
needed beyond the existence of the final polynomial tensor.

## 4. Scalar content and degree

If a polynomial `h(w)` divides every entry of such a tensor, (1)
shows that it divides `1+w^(M^2)`. Since `M^2` is a power of two,

```math
1+w^{M^2}=(1+w)^{M^2}.
```

The strict crossing endpoint at `w=1` has entries equal to one. Thus
`1+w` cannot divide every entry, and `h` is a unit. This remains true
after any coefficient-field extension.

The same identity contains a nonzero coefficient of degree `M^2`.
Some entry of the tensor must have at least that degree. The bound
does not assert that degree `M^2` is attainable within the family.

## 5. Exact controls and finite span search

[`check_reflection_paired_trace.py`](../tools/check_reflection_paired_trace.py)
independently clears the denominators in the Gaussian Schur complement,
checks every entry through interval sixteen, and contracts the complete
four- and sixteen-letter polynomial tensors. It also checks arbitrary
paired bases, nonconstant Laurent-series frames, and affine combinations.
The evidence is in
[`reflection-paired-trace.json`](../evidence/reflection-paired-trace.json).

A separate exact search rebuilt all 7,176 tensors from the complete
5,355-generator rank-one catalogue. Their affine evaluation spans have
dimensions 102 at zero, 181 at each nontrivial `F4` point, 316 at the
two nontrivial points together, and 415 at all three points together.
None contains the scalar-zero target. The uniform proof above is
stronger than those finite tests; the span computation remains useful
for constructing other finite-rank targets.

The span program is
[`solve_dense_affine_scalar_content.py`](../tools/solve_dense_affine_scalar_content.py).
Its [complete record](../evidence/dense-affine-scalar-content.json) binds
every reconstruction to the exact source catalogue. No rank-pressure
or global optimality conclusion is inferred from these dimensions.
