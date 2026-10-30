# The central homogeneous piece on primitive inputs has no alternating cancellation

## Result and scope

Let `V=Q[t_1,t_2,...]`, with `wt(t_m)=m`, and use the operators `a[r]`
of [Appendix B, equation (47)](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
Write `V_h` for its weight-`h` component. Define derivations

```math
D t_1=0,\quad D t_m=mt_{m-1}\ (m\ge2),\qquad Ut_m=mt_{m+1},
```

and let `W_h=ker(D|V_h)`. These are the primitive spaces from the
[homogeneous-dual construction](homogeneous-duals.md). For each `j>=0`,
let `Pi_j` be the orthogonal projection from `V_j` to `W_j`; `W_0=Q`.

For `a,y in W_h`, define

```math
A_i(a,y)=a[i-h]y,\qquad C_j(a,y)=\Pi_j A_j(a,y).
\tag{1}
```

**Primitive central-rank theorem.** After an invertible change of the
operator label and a swap of the two row labels, the central homogeneous
piece restricted to `W_h` is a positive semidefinite Gram matrix. Its rank
is exactly

```math
\operatorname{rank}\left(\bigoplus_{j=0}^{2h}C_j:
W_h\otimes W_h\longrightarrow\bigoplus_{j=0}^{2h}W_j\right).
\tag{2}
```

Thus alternating cancellation cannot lower the central rank below the
joint image of these primitive fusion maps. This statement holds for all
positive `h` over the rationals. It does not assert a large enough image
to establish an asymptotic rank lower bound, and it proves no lower bound
for factoring algorithms. It concerns the full partition-label tensor;
an arbitrary binary-label restriction requires a separate argument.

The exact checker verifies the identities directly through `h=6`.
This note proves an algebraic restriction, not a new factoring exponent.

## 1. Inner products, adjoints, and a sign change

Use the positive Fock inner product

```math
\langle t^\mu,t^\nu\rangle=z_\mu[\mu=\nu],\qquad
z_\mu=\prod_m m^{\operatorname{mult}_m(\mu)}
                 \operatorname{mult}_m(\mu)!.
```

Then `D*=U` and `[D,U]=2H`, where `H` is the weight operator.
Put

```math
S(t^\mu)=(-1)^{|\mu|-\ell(\mu)}t^\mu,\qquad Y=e^D.
```

Here `S` is an involution, `DS=-SD`, and `S` preserves every `W_h`.
The exact adjoint formula for every polynomial `a` is

```math
(a[r])^*=(YSa)[-r].
\tag{3}
```

For a one-part state it follows from `nabla_k*=nabla_{-k}` and

```math
(-1)^{m-1}\sum_{j=1}^m\binom mj\binom{-r-1}{j-1}
=\binom{r-1}{m-1}.
```

For a monomial, transpose each normally ordered product in equation (47).
Transposition changes annihilations to creations and reverses their
order, so the result remains normally ordered. The one-part identity
applies independently to every part. This proves (3) by linearity.
For `a in W_h`, `YSa=Sa`, hence `(a[r])*=(Sa)[-r]`.

Consider the central bilinear form on row and column pairs

```math
T_h(a\otimes y,b\otimes u)
=\langle u,Q_0(a,Sb)y\rangle.
```

All states have weight `h`, so Appendix B's angle coefficient is
`binom(2h,i)`. Using the exchanged-order form of `Q_0` and (3) gives

```math
T_h(a\otimes y,b\otimes u)
=\sum_{i=0}^{2h}(-1)^i\binom{2h}{i}
       \langle A_i(a,y),A_i(b,u)\rangle.
\tag{4}
```

Changing `b` to `Sb` is invertible on `W_h`. Reading an output with the
Fock inner product is a nondegenerate change of output coordinates.

## 2. Decompose each intermediate vector into primitive components

The normal-order formula gives, for homogeneous `a in V_h`,

```math
[D,a[r]]=(r+h-1)a[r-1]+(Da)[r-1].
\tag{5}
```

One direct verification uses `[D,nabla_k]=-k nabla_{k+1}`, with
`nabla_0=0`, and the coefficient identity

```math
(1-q)\binom{-q}{m-1}
=(m-q)\binom{-q-1}{m-1}
 +m\binom{-q-1}{m-2}.
```

No creation-annihilation crossing occurs when a mode changes sign:
the intervening zero mode contributes zero. Since both `a` and `y`
are primitive, (5) specializes to

```math
D A_i(a,y)=(i-1)A_{i-1}(a,y).
\tag{6}
```

The orthogonal primitive decomposition and its ladder norm are

```math
V_i=\bigoplus_{j=1}^i U^{i-j}W_j,\qquad
\|U^n v\|^2=n!(2j)_n\|v\|^2\quad(v\in W_j),
\tag{7}
```

for `i>=1`. Here `(x)_n=x(x+1)...(x+n-1)`, with `(x)_0=1`.
These follow from `D*=U` and
`D U^n v=n(2j+n-1)U^{n-1}v`. The vacuum contributes only at weight
zero because `U1=0`.

Solving (6) separately in the summands of (7) gives

```math
A_i(a,y)=\sum_{j=1}^i
 \frac{(j)_{i-j}}{(i-j)!(2j)_{i-j}}
 U^{i-j}C_j(a,y),\qquad i\ge1.
\tag{8}
```

Substitute (8) into (4) and use orthogonality. The result is

```math
T_h=\sum_{j=0}^{2h}(-1)^j c_{h,j}\, C_j^*C_j,
\tag{9}
```

where `c_{h,0}=1`, and, for `j>=1`,

```math
c_{h,j}=\sum_{n=0}^{2h-j}(-1)^n
 \binom{2h}{j+n}\frac{(j)_n^2}{n!(2j)_n}.
\tag{10}
```

The adjoint in (9) uses the positive tensor-product inner product on
`W_h tensor W_h` and the Fock inner product on `W_j`.

## 3. Every coefficient magnitude is strictly positive

Let `B(j,j)=integral_0^1 t^(j-1)(1-t)^(j-1) dt`. The beta-integral
identity and an elementary binomial rearrangement give

```math
\frac{(j)_n}{(2j)_n}
=\frac1{B(j,j)}\int_0^1t^{j+n-1}(1-t)^{j-1}\,dt,
```

```math
\binom{2h}{j+n}\frac{(j)_n}{n!}
=\frac{j\binom{2h}{j}}{j+n}\binom{2h-j}{n}.
```

Expanding `(1-tu)^(2h-j)` therefore turns (10) into

```math
c_{h,j}=
\frac{j\binom{2h}{j}}{B(j,j)}
\int_0^1\!\int_0^1
 t^{j-1}(1-t)^{j-1}u^{j-1}(1-tu)^{2h-j}\,du\,dt>0.
\tag{11}
```

The integrand is positive in the open unit square. Thus the signs in
(9) depend only on the parity of the primitive weight `j`.

## 4. Skew symmetry puts the two signs on separate subspaces

The source fields satisfy the elementary identity

```math
Y(a,z)y=e^{zU}Y(y,-z)a.
\tag{12}
```

For completeness, it can be derived directly from equation (47), without
assuming further vertex-algebra results. Use the formal states
`E(x)=exp(sum_m x_m t_m)`. Their normally ordered fields satisfy

```math
Y(E(x),z)E(y)
=\exp(K(x,y;z))
 \exp\left(\sum_m y_mt_m+\sum_m x_m e^{zU}t_m\right),
```

where

```math
K(x,y;z)=\sum_{m,k\ge1}
(-1)^{m-1}\frac{(m+k-1)!}{(m-1)!(k-1)!}
 x_m y_k z^{-m-k}.
```

Now `K(x,y;z)=K(y,x;-z)` and `e^{zU}e^{-zU}=Id`.
The two exponential expressions in (12) consequently agree. Taking
coefficients in the auxiliary variables proves (12) for all polynomials.
All identities are formal; each fixed auxiliary coefficient has a finite
principal Laurent part.

For two states of weight `h`, the coefficient of `z^(j-2h)` in (12) is

```math
A_j(a,y)=\sum_{k=0}^j
 \frac{(-1)^{j-k}}{k!}U^k A_{j-k}(y,a).
```

Apply `Pi_j`. Every term with `k>0` lies in the image of `U`, which is
orthogonal to `W_j`. Therefore

```math
C_j(y,a)=(-1)^j C_j(a,y).
\tag{13}
```

Let `tau(a tensor y)=y tensor a`. Combining (9), (11), and (13) gives

```math
T_h(\tau v,w)=\sum_{j=0}^{2h}c_{h,j}
                    \langle C_jv,C_jw\rangle.
\tag{14}
```

The right side is positive semidefinite, and its kernel is exactly
`intersection_j ker C_j`. Since `tau` is invertible, (14) proves the
rank equality (2). It also shows why a claim that the unswapped central
form itself is positive semidefinite would be incorrect: its odd
primitive channels carry a negative sign.

## 5. Exact finite checks and the unresolved asymptotic question

[check_primary_central_gram.py](../tools/check_primary_central_gram.py)
computes primitive projectors over the rationals and compares the actual
source operator products with (4), (8), (9), (13), and (14). The primary
projector on `V_j` is evaluated as a polynomial in `X_j=UD`:

```math
\Pi_j=\prod_{k=1}^{j-1}
 \frac{X_j-\chi_k}{-\chi_k},\qquad
\chi_k=j(j-1)-k(k-1),
```

with the identity projector at `j=0,1`. This uses the eigenvalues proved
in the homogeneous-dual note. All checks use exact rational arithmetic;
no numerical positivity threshold or modular rank is used.

| `h` | `dim W_h` | Central restricted rank | Stacked fusion rank |
| ---: | ---: | ---: | ---: |
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 |
| 3 | 1 | 1 | 1 |
| 4 | 2 | 4 | 4 |
| 5 | 2 | 4 | 4 |
| 6 | 4 | 16 | 16 |

These finite ranks do not prove a first-order asymptotic profile.
To derive such a profile from (2), one still needs a lower bound on the
joint image of the `C_j` as `h` grows. The full-input channel-span theorem
does not supply that bound for primitive input pairs.

Run `python3 tools/check_primary_central_gram.py`. The recorded output is
[primary-central-gram.json](../evidence/primary-central-gram.json).
The new argument and checker used substantial OpenAI assistant support.
