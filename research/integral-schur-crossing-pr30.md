# The complete partition crossing over the integral Schur lattice

## Result and scope

The complete Appendix B partition construction has an integral form
whose crossing multiplier is one. Its states and operator parameters
use Schur functions, rather than power-sum monomials. It therefore
reduces to every characteristic without deleting states or inverting
the source's factorial multiplier.

The two explicit channel counts remain the ordinary partition counts.
Their declared leading width parameter is still `rho=1/sqrt(2)`.
This note does not prove a stronger factoring exponent. It provides a
complete-state modular construction in which actual ranks may differ
from their characteristic-zero values.

The operators and crossing identities are the prior Appendix B
construction at
[`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Their definitions are in
[`Boson/Defs.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean),
[`Boson/Defs2.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs2.lean),
and
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The contribution here is the integral lattice and its normalization,
with independent finite checks.

## 1. Integral states and normal symbols

Write `Lambda_Z` for the ring of symmetric functions over the integers.
Its complete symmetric functions `h_m` freely generate the ring, and
its Schur functions `s_lambda` form an integral basis. Work first over
the rationals, where the power sums `p_m` freely generate the ring.
These are the variables denoted `t_m` in the source.

The Hall pairing and adjoints satisfy

```math
\langle p_\lambda,p_\mu\rangle=z_\lambda[\lambda=\mu],
\qquad \langle s_\lambda,s_\mu\rangle=[\lambda=\mu],
\qquad p_k^\perp=k\frac{\partial}{\partial p_k},                 \tag{1}
```

where `z_lambda=product_k k^(m_k) m_k!`. Multiplication by an element
of `Lambda_Z` preserves the integral lattice. Its Hall adjoint does
too: every matrix entry in the Schur basis is the transpose of an
integer multiplication matrix.

For a homogeneous parameter `f` of weight `m`, package the source modes
as

```math
Y(f,z)=\sum_r f[r]z^{r-m}.                                    \tag{2}
```

The normally ordered symbol of `Y(p_m,z)` is the sum of

```math
\begin{aligned}
C_m(z)&=\sum_{\ell\ge m}\binom{\ell-1}{m-1}p_\ell z^{\ell-m},\\
D_m(z)&=\sum_{k\ge1}(-1)^{m-1}\binom{k+m-1}{m-1}
                   p_k z^{-k-m}.                             \tag{3}
\end{aligned}
```

Here the first summand belongs to the first tensor factor of
`Lambda_Q tensor Lambda_Q`, and the second to the second factor.
The symbol `a tensor b` acts as multiplication by `a`, followed by
`b^perp`. Multiplication of symbols is commutative normal-symbol
multiplication. It is not ordinary composition of their represented
operators. The source's recursion is precisely the normal ordering
obtained by substituting (3) into polynomial parameters.

Define

```math
H(u)=\sum_{m\ge0}h_m u^m
    =\exp\left(\sum_{m\ge1}\frac{p_m u^m}{m}\right),
\qquad H(u)^{-1}=E(-u).                                      \tag{4}
```

The normal-symbol generating series for `Y(H(u),z)` is

```math
\boxed{\frac{H(z+u)}{H(z)}\ \otimes\
       \frac{H(z^{-1})}{H((z+u)^{-1})}.}                       \tag{5}
```

To verify it, the creation exponent contributed by `p_ell` is

```math
\sum_{m\ge1}\frac{u^m}{m}\binom{\ell-1}{m-1}z^{\ell-m}
       =\frac{(z+u)^\ell-z^\ell}{\ell}.
```

The second-symbol exponent contributed by `p_k` is

```math
\sum_{m\ge1}\frac{u^m}{m}(-1)^{m-1}
          \binom{k+m-1}{m-1}z^{-k-m}
       =\frac{z^{-k}-(z+u)^{-k}}{k}.
```

These are the logarithms of the two ratios in (5).

Every coefficient of (5) is integral in the relevant completed tensor
ring: `H` and its reciprocal have integral coefficients, and the
expansion of `(z+u)^(-k)` uses integer binomial coefficients. No
convergence or rearrangement of an infinite numerical sum is needed.
For a fixed input weight `a` and output weight `b`, only second-factor
weights at most `a` can act, and the first-factor weight is then fixed
by `b-a`. Each resulting matrix coefficient uses a finite sum.

It follows that every mode of `h_m` preserves `Lambda_Z`. The normal
symbol map is multiplicative in its parameter, and
`Lambda_Z=Z[h_1,h_2,...]`. Therefore every mode of every integral
symmetric-function parameter preserves `Lambda_Z`. In particular,
every mode of a Schur function, and of any integral linear combination
of Schur functions of different weights, has an integer matrix.

## 2. The source dual transform is an integral involution

Let the rational algebra map `T` be specified by the source duals:

```math
T(p_m)=(-1)^{m-1}\sum_{j=1}^{m}\binom mjp_j.                  \tag{6}
```

Its generating series has the particularly simple form

```math
\boxed{T(H(t))=E\!\left(\frac{t}{1+t}\right).}                \tag{7}
```

Indeed, substitute (6) into the logarithm of `H(t)` and use
`binom(m,j)/m=binom(m-1,j-1)/j`. Equating coefficients gives, for
`m>=1`,

```math
T(h_m)=\sum_{j=1}^{m}(-1)^{m-j}\binom{m-1}{j-1}e_j.          \tag{8}
```

Consequently `T` preserves `Lambda_Z` and its weight filtration. Every
homogeneous component of `T(f)` is integral when `f` is integral.
Also `T(E(t))=H(t/(1-t))`, so applying (7) twice gives `T^2=id`.

The source's power-sum vacuum identity, after cancelling the nonzero
integer factorial in characteristic zero, is

```math
[1]T(p_\lambda)[-h]p_\mu
       =z_\mu[\lambda=\mu]\qquad(|\lambda|=|\mu|=h).          \tag{9}
```

By rational bilinearity this is the Hall pairing on all homogeneous
parameters of weight `h`. In particular,

```math
[1]T(s_\lambda)[-h]s_\mu=[\lambda=\mu].                      \tag{10}
```

All entries in (10) are integers by Section 1. Equation (10) is thus an
integral identity and reduces directly modulo any integer. It does not
cancel a vanishing factorial in the target characteristic.

## 3. Complete crossing pieces with multiplier one

Choose any `d` distinct partitions of weight `h`, and set

```math
a^\lambda=T(s_\lambda),\qquad b^\lambda=s_\lambda.
```

States use the full Schur basis at every intermediate weight. For a
parameter `a` supported in weights at most `h`, retain the source angle
operation

```math
a^{\langle i\rangle}
   =\sum_{w=0}^{h}\binom{h+w}{i}a_w,
\quad 0\le i\le2h.                                          \tag{11}
```

It is integral by (8). Define, as in the source,

```math
Q_s(a,b)=\sum_{i=0}^{2h}(-1)^i
             a^{\langle i\rangle}[s-h+i]b[h-s-i],
\quad -2h\le s\le2h,                                        \tag{12}
```

and let `A_s(x,y,u,v)` be the coefficient of `s_u` in
`Q_s(a^x,b^v)s_y`. The source exchange identity reverses the two
operators in every summand. The source summation identity and (10)
give

```math
\sum_{s=-2h}^{2h}A_s(x,y,u,v)=[x=v][y=u].                    \tag{13}
```

The source identities extend to these parameters by rational
bilinearity. Section 1 shows that every matrix coefficient on both
sides, including the two factorizations, is integral. Thus they are
identities over the integers, and scalar extension gives the same
identities over every commutative unital coefficient ring.

In the exchanged order, the intermediate weight is `s+i`. In the
original order it is `2h-s-i`. Inserting the complete Schur basis gives
explicit integral factors of lengths

```math
R_s=\sum_{i=0}^{2h}p(s+i),\qquad
C_s=\sum_{i=0}^{2h}p(2h-s-i),                                \tag{14}
```

with `p(j)=0` for `j<0`. These are factorization lengths over a ring;
their reductions bound ordinary matrix ranks over every field.
The construction retains all `p(h)` states at its main weight.

For binary labels choose `d=2^lambda` and a suitable `h` with
`p(h)>=d`. The source partition estimate applied to (14) gives the
same declared bound

```math
\log_2R_s+\log_2C_s
       \le 2\sqrt2\,\lambda+O(\log(\lambda+1)),               \tag{15}
```

using the least such `h`. The present normalization changes neither
the number of shifts nor that asymptotic estimate. It does permit
direct investigation of the complete reduced matrices over small
fields, where their actual ranks can be smaller.

## 4. An exact basis-change audit

Let `C_(lambda,mu)=chi^lambda(mu)` be the ordinary character table of
the symmetric group at weight `h`, and put
`S_(lambda,mu)=C_(lambda,mu)/z_mu`. Thus
`s_lambda=sum_mu S_(lambda,mu)p_mu` and `S C^T=I`.

Let `A_s^old` be the source integer matrix whose parameter families
are `T(p_x)` and `(h!/z_v)p_v`, with power-sum input and output states.
Rows are `(x,y)` and columns are `(v,u)`. The Schur matrix is exactly

```math
A_s^{\rm Schur}
 =\frac1{h!}(S\otimes S)A_s^{\rm old}(C^T\otimes C^T).        \tag{16}
```

On the output coordinate, extracting the Schur coefficient uses the
Hall pairing, which contributes `chi^u(mu)`. On the `b` parameter,
the coefficient `chi^v(mu)/z_mu` cancels the old factor `h!/z_mu`,
leaving `chi^v(mu)/h!`. These explain both right-hand factors and the
single division in (16). The result is an integer matrix by the proof
above, even though this formula computes it through rational values.

At `h=2`, in partition order `((1,1),(2))`, the central matrix is

```math
\begin{pmatrix}
13&4&5&18\\9&6&3&12\\12&3&6&9\\18&5&4&13
\end{pmatrix}.                                               \tag{17}
```

Modulo two its two crossing flattenings have ranks two and one.
This is a complete-state rank drop, but the nine raw shift pieces
have square-root rank sum `4+7 sqrt(2)`, which does not give an
improving seed. The finite observation is evidence for a new modular
search family, not an asymptotic compression theorem.

## 5. Full homogeneous duals at prime-power weights

There is also a complete-state version of the homogeneous modular
construction. This corollary uses the eigenvalue calculation already
proved in [homogeneous-duals.md](homogeneous-duals.md), Sections 1--2,
and used for the smaller regular quotient in
[modular-crossing.md](modular-crossing.md), Section 6.

Let `omega(p_m)=(-1)^(m-1)p_m`, and let
`Y(p_m)=sum_(j=1)^m binom(m,j)p_j`. Both preserve the integral Schur
lattice: `omega` conjugates Schur labels and `Y=T omega`. Let `G_h`
be the positive Hall Gram operator on weight `h`, specified by

```math
\langle G_h f,g\rangle=\langle Yf,Yg\rangle.
```

Projection onto a homogeneous weight is integral, and Hall adjoints
of integral maps are integral. Thus `G_h` has an integer Schur matrix.
The homogeneous vacuum pairing is

```math
\mathcal P_h(f,g)=[1]f[-h]g
      =\langle G_h\omega f,g\rangle.                          \tag{18}
```

The cited characteristic-zero eigenvalue calculation gives

```math
\det G_h=\prod_{j=1}^h
       \binom{2h-1}{h-j}^{m_j},\qquad
m_1=1,\quad m_j=p(j)-p(j-1)\ (j\ge2).                       \tag{19}
```

If `h=p^r` for a prime `p`, every displayed binomial is a unit
modulo `p`. The lower `r` base-`p` digits of `2h-1` are all `p-1`,
whereas `0<=h-j<h`. Lucas's digit identity therefore gives a product
of nonzero residues. Equations (18)--(19), together with the integral
unimodularity of `omega`, show that the **full** homogeneous vacuum
matrix on all `p(h)` Schur states is nonsingular modulo `p`.

Invert this matrix over `F_p` and use its rows to choose homogeneous
Schur-linear parameters `a^lambda` dual to `b^mu=s_mu`. Their vacuum
pairing is the identity, so the same crossing construction applies.
Every angle factor is now the scalar `binom(2h,i)`. For `p=2`, only
`i=0,2h` survive. For odd `p`, only `i=0,h,2h` survive, with
coefficients `1,2,1`. Suppressing the other terms before reduction is
valid because each fixed-intermediate-weight composition is integral
by Section 1.

This is an algebraic corollary. Restricting `h` to powers of `p` can
overshoot the smallest weight needed for a binary label width by a
factor approaching `p`. Its logarithmic state dimension can then
overshoot by a factor approaching `sqrt(p)`. Sparse angles alone do
not therefore give the all-width bound (15), much less an improved
factoring exponent.

## 6. Independent finite checks and remaining work

`check_integral_schur_modes.py` constructs Schur functions from the
Jacobi--Trudi determinant over the power sums. This is independent of
a character implementation based on the Murnaghan--Nakayama rule.
It checks the source modes, the transform (7), the vacuum duals,
both orders of every small crossing piece, integral entries,
the exact summed crossing, and small-field ranks.

The default mode checks use operator weights through four, input
weights through four, and output weights through six. They examine
3,960 Schur matrix entries. The crossing checks use weights one
through three. A negative control records that reducing the old
factorial-scaled crossing before changing its integral lattice gives
the zero endpoint at weight two, whereas the new integral crossing
reduces to the identity.

The companion [effective construction audit](effective-integral-schur-transfer.md)
charges Schur conversion, exact rational arithmetic, factor-table output,
and coefficient bit sizes. It supplies the uniform `2^O(lambda)`
construction bound. The finite checker by itself establishes only the
reported finite identities; neither note improves the factoring exponent.

### Larger complete-state profiles

The optional NumPy implementations compute source power-sum factors, make
the exact Schur change of basis, and reduce each resulting integral entry.
They retain sufficient prime-power precision before dividing known factors
of the factorial, then check the complete crossing endpoint. Integer matrix
products have an explicit overflow bound. Schur conversion uses the
division-free constructor from the companion audit.

| Family over `F_2` | Weight | External dimension | Central ranks | Raw sum of square-root rank products |
| --- | ---: | ---: | --- | ---: |
| Inhomogeneous | 2 | 2 | `(2,1)` | `13.899495` |
| Homogeneous | 2 | 2 | `(2,2)` | `6` |
| Homogeneous | 4 | 5 | `(12,12)` | `44.324555` |
| Homogeneous | 8 | 22 | `(228,228)` | `1777.546262` |

The homogeneous weight-eight run checks 7,730,448 entrywise divisibilities,
the invertibility of the full 22-dimensional vacuum Gram matrix, and every
entry of the normalized crossing. Odd shift pieces vanish at weights two
and four, but some are nonzero at weight eight. Therefore the first two
weights do not support a general even-shift assertion.

The [complete weight-four span search](full-schur-span-recombination.md)
also tests every binary linear recombination, allowing an arbitrary number
of summands. Its certified minimum cost is still too large. None of the
listed raw profiles gives an improved finite seed.

```sh
python3 tools/full_schur_crossing.py --weights 1 2 3 4 5 --prime 2 \
  --output evidence/full-schur-f2.json
python3 tools/homogeneous_full_schur.py --weights 2 4 8 \
  --output evidence/homogeneous-full-schur-f2.json
```

Add `--pieces PATH` to export the final weight's normalized binary matrices.
For the two recombination inputs use weight four alone and the corresponding
`evidence/full-h4-pieces.json` or `evidence/homogeneous-h4-pieces.json` path.
