# A binary site-symmetry barrier for actual ranks

## Claim and consequence

Let all four physical word alphabets be `{0,1}^m`. Suppose

```math
\sum_i A_i=\kappa E_m,\qquad \kappa\ne0,
```

where every `A_i` is invariant under simultaneously permuting the `m`
positions in its four words. Write `r(A)` and `c(A)` for ranks in cuts
`(x,y)|(u,v)` and `(x,u)|(y,v)`. Over any fixed field,

```math
\boxed{\max_i r(A_i)c(A_i)\ge 2^{3m-o(m)}.}          \tag{1}
```

In characteristic zero the loss in the exponent is `O(log m)`; in fixed
positive characteristic `p` it is `O_p(log^2 m)`. No bound on the number
of summands is needed. Consequently a family with
`r(A_i)c(A_i)<=2^(4 rho m+o(m))` must have **`rho>=3/4`**, strictly above
the post-half target `1/sqrt(2)`.

This ordinary proof covers every scalar binary polynomial seed, at every
degree, with arbitrary coefficient recombinations, several seed curves,
and parameters that vary with `m`, provided the pieces retain this common
site symmetry. It is an actual-rank obstruction, stronger in scope than
the [bare Smith criterion](multipoint-smith-floor.md). It does not rule
out arbitrary factoring algorithms, genuinely coupled larger alphabets,
physical code restrictions, or pieces that break the stated symmetry.
Section 8 extends the same obstruction to any common partition into
`o(m)` position classes; merely using a slowly growing menu of scalar
local seeds is therefore insufficient.
The established factoring saving exponent remains **1/2**. There is no
new Lean formalization or independent mathematical review of this proof.

## 1. Replace the seed curve by the universal tensor

Let `T` be the four-leg binary tensor with 16 independent entries
`a_(xyuv)`, numbered `8x+4y+2u+v`. Every site-invariant tensor has the form

```math
A=L(T^{\otimes m})                                 \tag{2}
```

for a linear functional on degree-`m` polynomials in the 16 variables.
Indeed, the coefficient of each monomial in `T^(tensor m)` is exactly
the sum of the distinct tensor coordinates in one permutation orbit.
These orbit sums form a basis of the invariant space in every
characteristic. No division by an orbit size or by `m!` is involved.

Thus it suffices to construct one common polynomial covector that
detects `E_m` and forces high ranks in both cuts of the same summand.
For the crossing, the nonzero local coordinates are `a_0,a_6,a_9,a_15`,
all equal to one.

## 2. A common degree-four polynomial

Regard `T` as a two-by-two array of matrices
`F_(xv)=(a_(xyuv))_(y,u)` and put `D_(xv)=det F_(xv)`. Define

```math
P=D_{00}D_{11}-D_{01}D_{10}.                        \tag{3}
```

Its eight integer monomials are

```text
 a0 a6 a9 a15 - a0 a6 a11 a13 - a2 a4 a9 a15 + a2 a4 a11 a13
-a1 a7 a8 a14 + a1 a7 a10 a12 + a3 a5 a8 a14 - a3 a5 a10 a12.
```

In particular, `P(E)=1`. For **either** four-by-four cut matrix `M`,
with row and column indices starting at zero, direct integer expansion
gives

```math
P=M_{01}\operatorname{Cof}_{01}(M)
 +M_{03}\operatorname{Cof}_{03}(M)
 -M_{10}\operatorname{Cof}_{10}(M)
 -M_{12}\operatorname{Cof}_{12}(M).                 \tag{4}
```

Let `A=diag(1,1,0,0)` and `B=diag(0,1,0,1)`. Equivalently,

```math
P=\operatorname{tr}(A M B\operatorname{adj}M)-\det M.\tag{5}
```

The same `P` therefore belongs to a useful polynomial representation
coefficient space in both cuts, even for a seed whose ordinary cofactor
spans are disjoint. The degree-eleven control in
[matching-curve-controls.md](matching-curve-controls.md) illustrates
why looking only for common minors misses this certificate.

## 3. The representation, including characteristic two

On the 15-dimensional space of trace-zero matrices consider

```math
R(M)(X)=M X\operatorname{adj}(M).
```

For invertible `M` this is conjugation tensored with the determinant;
the displayed degree-four polynomial formula defines it for all `M`.
It is the Weyl module `Delta(lambda)` with

```math
\lambda=(2,1,1,0).                                 \tag{6}
```

One way to identify the module in every characteristic is to use the
integral highest-root Weyl module of `SL_4`, which is `Lie(SL_4)`, and
then twist the `GL_4` action by `det`. The highest-root identification
over the integers is stated in Garibaldi and Nakano,
[*Bilinear and Quadratic Forms on Rational Modules of Split Reductive Groups*](https://doi.org/10.4153/CJM-2015-042-5),
section 4.5, page 404. In characteristic two this module contains the
scalar matrix `I`; semisimplicity is neither asserted nor needed.

Here is explicit coefficient membership of `P`:

- In odd characteristic, set `A*=2A-I`, `B*=2B-I`. Both are trace zero,
  and `tr(A* M B* adj(M))=4P`. Thus `P` is a matrix coefficient of `R`.
- In characteristic two, `B` and `I` are both trace zero. Formula (5)
  writes `P` as the sum of two matrix coefficients of `R`, because
  `R(M)(I)=det(M) I`. No division by two or four is used.

For any representation `W`, write `C(W)` for the linear span of its
matrix coefficients. We have proved
`P in C(Delta(lambda))` in both cuts over every field.

## 4. Amplification without semisimplicity

Extend the field to its algebraic closure; ranks do not change. Tensor
powers of the natural four-dimensional module have Weyl filtrations.
The multiplicity of `Delta(mu)` in degree `d` is the integer number
`f^mu` of standard tableaux of shape `mu`. See Doty,
[*Schur-Weyl duality in positive characteristic*](https://arxiv.org/pdf/math/0610591v3),
Lemma 13 and the formal-character argument on pages 13--14.
The [determinant note](joint-determinant-barrier.md) explains the same
filtration-to-rank step for rectangular weights.

In particular `Delta(t lambda)` occurs in `V^(tensor 4t)` with
multiplicity

```math
h(t)=f^{(2t,t,t,0)}
 =\binom{4t}{2t}\binom{2t}{t}
   \frac{t+2}{(2t+1)(2t+2)}
 \ge\frac{2^{6t}}{2(4t+1)(2t+1)^2}.                \tag{7}
```

These multiplicities count separate diagonal filtration factors; they
are not field scalars reduced modulo `p`. In a fixed filtration basis
the polynomial matrix `M^(tensor 4t)` is block triangular. Whenever
`L` is nonzero on some coefficient of `Delta(t lambda)`, each of its
`h(t)` identical diagonal factors is a nonzero matrix after applying
`L`. Selecting one nonzero entry from each yields a triangular minor,
so the rank is at least `h(t)`.

Suppose now the characteristic is `p>0` and `t=p^j`. The Frobenius twist
`Delta(lambda)^(j)` is a cyclic highest-weight module of weight
`t lambda`: the Frobenius map on the algebraically closed group is
surjective, so the original highest-weight generator is still cyclic.
The universal property of a Weyl module therefore gives a surjection

```math
\Delta(t\lambda)\twoheadrightarrow\Delta(\lambda)^{(j)}.
```

Coefficient spaces of quotients are contained in coefficient spaces
of the original modules. The freshman's-dream identity now gives
`P^t in C(Delta(t lambda))`. This remains true when `P` was a sum of
coefficients, as in characteristic two.

## 5. Base-p digits give every tensor length

A subsequence of lengths is insufficient when a construction may
choose its own lengths. Write

```math
m=4t+s,\quad 0\le s<4,\qquad
t=\sum_{j\ge0}d_jp^j,\quad 0\le d_j<p.
```

Split the `m` tensor positions into `d_j` chunks of length `4p^j`,
and `s` spare positions. Tensor the fixed Weyl filtrations of the
chunks. The product module `tensor_j Delta(p^j lambda)^(tensor d_j)`
occurs as a diagonal factor at least

```math
B_p(t)=\prod_j h(p^j)^{d_j}                         \tag{8}
```

times. Its quotient is the tensor product of the Frobenius twists
used above. Consequently its coefficient space contains `P^t`.
For the spare positions choose the entry `ell=a_0`, with `ell(E)=1`.
The selected diagonal factors, including `V^(tensor s)`, contain the
common coefficient polynomial

```math
Q=P^t\ell^s,\qquad Q(E)=1.                          \tag{9}
```

By (7), the denominator loss in (8) has logarithm
`O_p(sum_(j<=log_p t) (j+1))=O_p(log^2(t+2))`. Thus

```math
B_p(t)\ge2^{6t-O_p(\log^2(t+2))}.                 \tag{10}
```

Now apply the physical tensor covector corresponding to the polynomial
`Q` to `sum_i A_i=kappa E_m`. In (2) this says
`sum_i L_i(Q)=kappa`, so some **same** `i` has `L_i(Q)!=0` in both cut
descriptions. The selected-factor argument gives both ranks at least
`B_p(t)`. Their product satisfies (1), uniformly at every `m`.

## 6. Characteristic zero

In characteristic zero `Delta(lambda)` is irreducible. Every irreducible
constituent `Delta(mu)` of `Delta(lambda)^(tensor t)` has highest weight
dominated by `t lambda`. Since `P^t` is a coefficient of that tensor
product, `L(P^t ell^s)!=0` detects at least one such constituent,
tensored with `V^(tensor s)`, in each cut. The detected constituents
may differ between the cuts; the following bound holds for all of them.

The hook formula gives, for `|mu|=4t`,

```math
f^\mu=\frac{(4t)!\prod_{i<j}(\mu_i-\mu_j+j-i)}
                  {\prod_{i=1}^4(\mu_i+4-i)!}
 \ge\frac{(4t)!/\prod_i\mu_i!}{(4t+4)^6}
 \ge\frac{2^{6t}}{(4t+1)(2t+1)(4t+4)^6}.           \tag{11}
```

The second inequality uses convexity of `log(k!)`: among shapes
dominated by `(2t,t,t,0)`, the factorial product is at most
`(2t)! t! t!`. Central binomial coefficients give the last bound.
The semisimple decomposition of `V^(tensor 4t)` has `f^mu` copies of
each constituent. Detecting one therefore bounds each cut rank below
by (11). The remainder `s<4` changes only a constant factor in the
exponent, proving (1) with an `O(log m)` loss.

## 7. Exact checks, limits, and the next construction question

Run

```sh
make verify-binary-site-symmetry
```

The standard-library checker verifies the two cut identities over
`Z`, the odd- and even-characteristic coefficient expressions,
`P(E)=1`, the hook formula through `t=64`, all dominated shapes through
`t=12`, and 24 exact base-prime amplification controls. It also checks
literal four-copy matrices for 128 binary, 46 ternary and 54 quinary
nonzero-witness choices in an eight-orbit subspace. Those finite checks
support the algebra; they do not independently establish the cited
representation-theoretic results or replace a proof review.

A single nonsymmetric coordinate tensor has both ranks one and is
detected by the same polynomial covector. The negative control confirms
why common site symmetry is essential. In particular, the theorem
does not justify symmetrizing an arbitrary proposed construction:
symmetrization can increase its ranks or vanish in small characteristic.

Increasing scalar degree, mixing scalar curves, or using one common
word relabeling preserves the obstruction. A useful next candidate
must change that structural assumption: for example, a physical code
whose restrictions alter the available covectors, a growing collection
of position-dependent local operators, or a genuinely coupled larger
alphabet. Any such construction still needs a uniform rank **upper**
bound, a complete crossing identity, and full evaluator cost accounting
before it can improve the factoring exponent.

## 8. Sublinearly many position classes still do not suffice

Suppose there is a common partition of the `m` positions into `k`
nonempty classes of sizes `m_1,...,m_k`, and every piece is invariant
under permutations **within** each class. No invariance between
different classes is assumed. The pieces need not factor across
classes: arbitrary joint coefficient functionals are allowed.

Replace the single universal tensor by independent universal tensors
`T_1,...,T_k`. Each piece has a representation

```math
A_i=L_i\left(\bigotimes_{b=1}^k T_b^{\otimes m_b}\right).
```

For each class set `t_b=floor(m_b/4)`, `s_b=m_b-4t_b`, and use the
common witness

```math
Q=\prod_{b=1}^k P(T_b)^{t_b}\ell(T_b)^{s_b}.
```

It again has value one at `E_m`. Tensoring the separate class
filtrations gives the exact positive-characteristic bound

```math
\max_i r(A_i)c(A_i)\ge
 \left(\prod_{b=1}^k B_p(t_b)\right)^2.             \tag{12}
```

The same summand is detected in both cuts, as before. With constants
depending only on the fixed characteristic, the logarithm of (12) is
at least

```math
3m-O_p\left(\sum_b\log^2(m_b+2)\right)
\ge3m-O_p\left(k\log^2(m/k+2)\right).              \tag{13}
```

For the last step, `log^2(x+2)` is concave for `x>=1`, since its
second derivative has the sign of `1-log(x+2)<0`; apply Jensen's
inequality. In characteristic zero, (11) gives the still smaller
loss `O(k log(m/k+4))`.

If `k=o(m)`, put `L=m/k`. The loss divided by `m` is
`O_p(log^2(L+2)/L)=o(1)`. Thus the leading floor remains **3/4**
even for a growing number of classes. This includes periodic scalar
seeds with a fixed period, bounded menus of local polynomials with
arbitrary position assignments, and sublinearly growing menus. Their
degrees and the joint coefficient recombinations may vary freely.

The conclusion concerns one common position partition for all pieces.
It does not cover unrelated piece-dependent partitions, or a
noncommuting auxiliary product whose value changes when positions
of the same nominal local type are exchanged. Nor does a linear
number of classes by itself provide a useful upper bound: it only
removes this particular asymptotic obstruction.
