# Almost full intermediate spans for the full restricted label spaces

## Result and scope

Restricting the labels to partitions with small parts preserves the original
Appendix B crossing identity. It does not, by itself, confine the intermediate
vectors to the small-part partition space. This note gives a quantitative
obstruction to that proposed span bound.

Let `p_L(h)` count partitions of `h` whose parts are at most `L`, and put
`D=p_L(h)`. For `L>=4` and `h>=2L^2`, consider the source's row factors on
the **full** restricted label set. Already among their fixed `i=2h` blocks,
one block has rational rank at least

```math
\frac{D}{4h+1}
\left(\left\lfloor\frac{h}{8L^2}\right\rfloor+1\right)^{L-2}
\ \ge\
\frac{D^2}{2^{8L}h(4h+1)}.
\tag{1}
```

The original operators are unchanged. The proof observes a subset of output
coordinates to witness a lower bound on the original intermediate rank.
It does not insert projections into a piece or alter its crossing identity.

In particular, if `D>=2^lambda`, `L=o(lambda)`, and `log h=o(lambda)`, then
some such block has width at least

```math
2\log_2D-o(\lambda).
```

In the proposed regime `h=lambda^c`, `c>2`, and
`L=Theta(lambda/log lambda)`, this is almost the entire two-label dimension.
It rules out a uniform bound of the form `p_L(2h)*2^o(lambda)` on these
full-domain intermediate spans.

There are two essential limits to this conclusion. It is not a lower bound
on the final piece matrices `A_s`: composing the two factors and summing
over `i` can reduce rank. Also, a chosen binary label subset can omit the
short partitions used in the proof. Section 6 states the exact retained-subset
bound; no obstruction for every binary encoding is claimed.

The operator definitions and row factorization are from Claude (Anthropic),
*A 2^o(n) Time Algorithm for 3-SAT*, at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The source definitions used below are `term`, `opMon`, `angle`, `aFam`,
and `PieceFacts.row` in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
This note is an ordinary mathematical argument, with no new Lean proof.

## 1. The precise block being bounded

Write `ell(mu)` for the number of parts of a partition `mu`. Let
`P_L(h)` be the set counted by `p_L(h)`. For `-2h<=s<=2h`, define the matrix

```math
F_s(\mu,\nu;\tau)
=(-1)^{h-\ell(\mu)}
 [t^\tau]\,(t^\mu)[s+h]t^\nu,
\qquad
\mu,\nu\in\mathcal P_L(h),\quad |\tau|=s+2h.
\tag{2}
```

Its rows are the ordered label pairs `(mu,nu)`. Its columns include every
partition of `s+2h`, without a restriction on its largest part.

This is exactly the `i=2h` block of the source's row factor, up to the
displayed row signs. Indeed,

```math
a^\mu=\prod_{n\in\mu}(-1)^{n-1}
       \left(\sum_{j=1}^{n}\binom nj t_j\right).
```

Every monomial in this polynomial has weight at most `h`. The only
weight-`h` term is `(-1)^(h-ell(mu))*t^mu`. At `i=2h`, the factor
`binom(h+j,2h)` kills every homogeneous component of weight `j<h` and is
one for `j=h`. The outer sign `(-1)^(2h)` is one. The row-factor operator
shift is `s-h+i=s+h`, and its intermediate energy is `s+i=s+2h`.

Thus a lower bound on (2) is a lower bound on an actual intermediate
factor block. It does not identify that block with the full piece matrix.

## 2. An injective map from the normally ordered operators

Let `P_L=Q[t_1,...,t_L]`. Define a final coordinate observation `Pi_(2L)`
by setting `t_j=0` for `j>2L`. For `a,b` in `P_L`, put

```math
\mathcal H_L(a\otimes b)
=\sum_{r\in\mathbb Z}\Pi_{2L}(a[r]b).
\tag{3}
```

This sum is finite. If a monomial of `a` has `d` factors and a monomial
of `b` has weight `E`, a surviving term has at most `d` created factors,
each of weight at most `2L`. Differentiations are applied before any
multiplications, as specified by the source's `term`. The output weight
therefore lies between zero and `E+2Ld`, and the mode lies between
`-E` and `2Ld`.

We claim that (3), extended linearly to the tensor product, is injective.
Use independent formal variables

```math
X_1,\ldots,X_L,\qquad Y_1,\ldots,Y_L.
```

The creator and annihilator coefficients from the source are

```math
C_n=\sum_{l=1}^{2L}\binom{l-1}{n-1}t_l,
\qquad
\Gamma_{nk}=(-1)^{n-1}\binom{k+n-1}{n-1}k.
\tag{4}
```

Here `Gamma_(nk)=k*binom(-k-1,n-1)`. Define the constant-coefficient
differential operator

```math
\Delta=\sum_{n=1}^{L}\sum_{k=1}^{L}
        \Gamma_{nk}\,\partial_{X_n}\partial_{Y_k}.
\tag{5}
```

On each polynomial, only finitely many powers of `Delta` are nonzero.
Consequently `exp(Delta)` is a well-defined linear automorphism of
`Q[X,Y]`, with inverse `exp(-Delta)`.

Let `T` be the substitution

```math
X_n\longmapsto C_n,\qquad Y_k\longmapsto t_k.
```

The original normal ordering gives exactly

```math
\mathcal H_L(a\otimes b)
=T\left(\exp(\Delta)\bigl(a(X)b(Y)\bigr)\right).
\tag{6}
```

To verify the coefficients, put
`D_n=sum_k Gamma_(nk)*partial_(Y_k)`. The differential operators commute,
so the expression before substitution is

```math
\sum_{\gamma\in\mathbb N^L}
\frac{\partial_X^\gamma a(X)}{\gamma!}
\,D_1^{\gamma_1}\cdots D_L^{\gamma_L}b(Y).
```

For a monomial `a`, the derivative divided by `gamma!` chooses which
operator factors are annihilators. The remaining factors are creators.
All annihilators act on `b` before the creators are substituted, exactly
as in the source. Only derivatives in `Y_1,...,Y_L` can act nontrivially.
Observing the final modes through `2L` replaces each creator by (4).
This proves (6). Equivalently, one can retain a mode variable `z` and
then set `z=1` after the final observation has made the Laurent sum finite.

The order in (6) matters: apply `exp(Delta)` first and substitute second.
Differentiating after substitution would also differentiate the creator
forms and would not be the source's normal ordering.

Finally, `T` is invertible. Modulo the low variables `t_1,...,t_L`, its
matrix on the creator forms is

```math
B_{nk}=\binom{L+k-1}{n-1},\qquad 1\le n,k\le L.
```

The polynomials `binom(x,n-1)` have successive degrees `0,...,L-1` and
leading coefficients `1/(n-1)!`. Their evaluation determinant at
`x=L,...,2L-1` is

```math
\det B
=\frac{\prod_{1\le j<k\le L}(k-j)}
        {\prod_{n=1}^{L}(n-1)!}=1.
\tag{7}
```

Together with the unchanged low variables, the creator forms therefore
give an invertible linear change of the `2L` polynomial variables.
Equation (6) is a composition of two linear isomorphisms. This proves
the injectivity of (3).

The coordinate observation here is only a rank witness. If observed
vectors are independent, their original unobserved vectors are independent.
There is no claim that inserting this observation between two operators
preserves a crossing identity.

## 3. Keep the injective image within the source's energy window

Let `S_L(h)` be the partitions in `P_L(h)` with at most `3h/(2L)` parts,
and let `S` be their number. Restrict the operator input in (3) to their
span, and allow every state monomial in `P_L(h)`.

The domain has dimension `SD`. Injectivity gives an image of that same
dimension. Every observed output has weight at most

```math
h+2L\frac{3h}{2L}=4h.
```

Decompose the output into its `4h+1` possible weights `E=0,...,4h`.
The rank of the whole map is at most the sum of the component ranks.
Some component therefore has rank at least `SD/(4h+1)`.

That component is the observed coefficient map for the single mode
`r=E-h`. Set `s=E-2h`. Then `-2h<=s<=2h`, `r=s+h`, and the component
is a submatrix of (2), after harmless row signs. Its columns are the
weight-`E` partitions with parts at most `2L`; (2) includes those
columns and possibly more. Hence

```math
\max_{|s|\le2h}\operatorname{rank}_{\mathbb Q}F_s
\ge\frac{SD}{4h+1}.
\tag{8}
```

The maximizing `s` can depend on `h,L`. This argument does not claim
that the central block `s=0` attains the bound.

## 4. There are many short operator partitions

Put `m=floor(h/(8L^2))`. Choose the multiplicities of parts
`1,...,L-2` independently in `{0,...,m}`. Their total number of parts
is at most `h/(8L)`, and their total weight `w` is at most `h/16`.
There remain `R=h-w>=15h/16` units of weight.

Choose

```math
u=(-R)\bmod L,\qquad
v=\frac{R-(L-1)u}{L},
```

and use `u` copies of `L-1` and `v` copies of `L`. The numerator defining
`v` is divisible by `L`. It is nonnegative because `u<=L-1` and
`R>=15h/16>=(15/8)L^2>(L-1)^2`.

These added parts have total number

```math
u+v=R/L+u/L\le h/L+1.
```

The complete partition therefore has at most

```math
\frac{9h}{8L}+1\le\frac{3h}{2L}
```

parts, using `h>=2L^2` and `L>=4`. Different choices of the first
`L-2` multiplicities give different partitions. Thus

```math
S\ge(m+1)^{L-2}.
\tag{9}
```

Equations (8)--(9) give the first bound in (1).

## 5. Exact comparison with the small-part partition counts

A partition counted by `D=p_L(h)` is determined by its multiplicities
`x_2,...,x_L`, with `sum_(j=2)^L j*x_j<=h`; the multiplicity of part one
fills the remaining weight. Comparing the disjoint unit cubes at these
lattice points with weighted simplices gives

```math
\frac{h^{L-1}}{(L-1)!L!}
\le D\le
\frac{(h+\sum_{j=2}^{L}j)^{L-1}}{(L-1)!L!}.
\tag{10}
```

For the lower bound, round each point of the simplex down coordinatewise.
For the upper bound, every cube at a feasible lattice point lies in the
simplex with the displayed larger threshold. This proves both bounds
without a fixed-`L` partition asymptotic or a period-dependent error term.

Since `sum_(j=2)^L j<=L^2<=h/2`, the upper bound in (10) is at most
`(3h/2)^(L-1)/((L-1)!L!)`. The elementary inequality `n!>=(n/3)^n`
and `L-1>=L/2` give

```math
D\le
\frac{3^{3L-2}h^{L-1}}{L^{2L-1}}.
```

Also, (9) gives `S>=(h/(8L^2))^(L-2)`. Consequently

```math
\frac DS
\le 3^{3L-2}8^{L-2}\frac h{L^3}
\le 216^L h
\le2^{8L}h.
\tag{11}
```

Substituting (11) in (8) proves the second bound in (1).

The same volume comparison at `2h` and `h` yields

```math
\frac{p_L(2h)}D
\le\left(2+\frac{\sum_{j=2}^L j}{h}\right)^{L-1}
\le(5/2)^{L-1}\le2^{2L}.
\tag{12}
```

Combining (1) and (12),

```math
\frac{\max_{|s|\le2h}\operatorname{rank}_{\mathbb Q}F_s}{p_L(2h)}
\ge\frac{D}{2^{10L}h(4h+1)}.
\tag{13}
```

If `D>=2^lambda`, `L=o(lambda)`, and `log h=o(lambda)`, the right side
is at least `2^(lambda-o(lambda))`. A factor `2^o(lambda)` cannot cover
this gap. The restricted labels therefore do not justify replacing the
original channel spaces by `p_L(2h)` times a subexponential factor,
uniformly over the pieces, on their full label domain.

For the suggested `h=lambda^c`, `c>2`,
`L=Theta(lambda/log lambda)` regime, the volume bounds also give

```math
\log D=L\log(h/L^2)+O(L+\log h)=\Theta(\lambda).
```

Thus (1) is `D^(2-o(1))`, while (12) is `p_L(2h)=D^(1+o(1))`.
The tensor product of the two label spaces, rather than the ordinary
product space alone, is reflected in the intermediate span.

## 6. What follows for a retained binary subset

The injectivity statement (6) applies to any operator and state subspaces.
Let a selected label set `B` have size `d`, and let

```math
s_B=\#\{\mu\in B:\ell(\mu)\le3h/(2L)\}.
```

The energy-window argument gives the exact conditional bound

```math
\max_{|s|\le2h}\operatorname{rank}_{\mathbb Q}(F_s|_{B\times B})
\ge\frac{d\,s_B}{4h+1}.
\tag{14}
```

If every selected label is short, this becomes `d^2/(4h+1)`.

If `d=2^lambda`, `log h=o(lambda)`, and the encoding retains
`2^(lambda-o(lambda))` short operator labels, (14) again gives almost
full intermediate rank.
If it excludes those short labels, the count (9) for the full space
does not transfer to the retained set. This proof does not settle that
case. An encoding can also be studied through different intermediate
factorizations; (1) does not rule out cancellations or compression after
the two factors have been composed into `A_s`.

This is a restricted obstruction and a guide to the remaining problem.
It establishes no new factoring rate and no lower bound on the final
half-piece tensors or on factoring algorithms.
