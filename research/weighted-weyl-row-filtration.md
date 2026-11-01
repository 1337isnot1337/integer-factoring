# An all-energy row filtration for the binary weighted Weyl curve

## Result and scope

Let `h=2^R`, `R>=1`, and let `d=b(h)` be the number of partitions
of `h` into powers of two. Work over any characteristic-two field.
For the weighted Weyl curve

```math
P(z)(\alpha,y,u,\beta)
=[\alpha+u=\beta+y]\,[t^\alpha]
\prod_{j=0}^R(1+t_j)^{\beta_j}
 (1+(1+z)^{2^j}t_j)^{y_j},
\qquad P(1)=E_d,
\tag{1}
```

every tensor `T=L(P(z))` in its coefficient span has row
characteristic polynomial

```math
\det(tI-T_{\rm row})=(t-L(1))^{d^2}.
\tag{2}
```

Thus `L(1)!=0` forces full row rank `d^2`. Under this same
condition, the column rank is at least

```math
b(2h)-R-2.
\tag{3}
```

Every decomposition of `E_d` into tensors from this coefficient
span has a summand with `L(1)!=0`. For `R>=2`, the bound in (3)
is at least `d`, so that summand alone contributes at least
`d^(3/2)` to the rank sum. Since `3/2>sqrt(2)`, this excludes the
strict finite-seed criterion at every power-of-two energy `h>=4`.
Asymptotically, this family cannot have parameter below `rho=3/4`.

For arbitrary energies tending to infinity, Section 5 gives a
parity-sector restriction to a power-of-two shell whose dimension
is `d^(1-o(1))`. Consequently every original coefficient-span
crossing decomposition has rank sum at least `d^(3/2-o(1))`,
without a power-of-two assumption on the energy sequence.

These statements concern the original coefficient span.
They do not automatically include sums of arbitrary external
symmetry images. The separate
[protected-minor theorem](weighted-weyl-protected-minor.md) closes
that larger symmetry span at energy 32. Neither result changes the
best factoring exponent.

The later [exact multiplication law](weighted-weyl-additive-row-law.md)
gives an elementary alternative proof of (2): at power-of-two energy,
every coefficient-span tensor satisfies `T^2=L(1)^2 I`. The module
argument below remains useful for its structural interpretation, and
the column minor and arbitrary-energy restriction remain unchanged.

## 1. Realizing each row block as a group-module quotient

The row matrix has blocks indexed by the integer vector

```math
\delta=\alpha-y=\beta-u,\qquad \sum_j2^j\delta_j=0.
```

Fix one such block. Put

```math
w_j=2^j,\qquad K_j=2h/w_j-1,\qquad N_j=K_j+\delta_j.
\tag{4}
```

All `N_j` are nonnegative. Indeed, `|delta_j|<=h/w_j`, so
`N_j>=h/w_j-1`.

For one coordinate, an elementary coefficient change gives

```math
[t^\alpha](1+t)^\beta(1+qt)^y
=[s^\alpha](1+(1-q)s)^\beta(1-qs)^{-u-1},
\quad \alpha-y=\beta-u.
\tag{5}
```

To prove (5), write `y=alpha-delta` and substitute
`t=s/(1-qs)` in the formal coefficient residue. The derivative
contributes `(1-qs)^(-2)`, leaving exponent
`delta-beta-1=-u-1` on the second factor. This identity is valid
before reduction to characteristic two.

Take `q=(1+z)^w=1+z^w`. Since `K+1=2h/w` is a power of two and
`alpha<=h/w<K+1`, the negative exponent in (5) can be replaced,
for this coefficient, by `K-u`. Hence the local entry is

```math
[s^\alpha](1+z^w s)^\beta
               (1+(1+z)^w s)^{K-u}.
\tag{6}
```

Let `V` be the natural two-dimensional module for `SL_2`, with
coordinates `x,y`. Formula (6) is a matrix entry of the
Frobenius-twisted action of

```math
g(z)=\begin{pmatrix}z&1+z\\1&1\end{pmatrix}\in SL_2(F[z])
\tag{7}
```

on `Sym^(N_j)(V)^(j)`. The input monomial is
`x^(beta_j)y^(K_j-u_j)`, and the output monomial is
`x^(alpha_j)y^(K_j-y_j)`.

There is an actual invariant quotient giving the needed occupation
range. If `delta_j<=0`, use the whole symmetric power. If
`delta_j>0`, quotient by the span of monomials whose first exponent
is outside the interval `[delta_j,K_j]`.

Here is the invariance check. The Weyl element interchanges the two
outside intervals. The diagonal torus preserves each monomial.
For an upper unipotent, a coefficient from an outside exponent
`a<delta` to an inside exponent `b` is

```math
\binom{K+\delta-a}{b-a}.
```

Write `K=2^B-1` and `r=delta-a>0`. The upper integer is
`2^B+(r-1)`, while `r<=b-a<2^B`. Lucas' identity makes the
coefficient zero: a number at least `r` cannot be a bit-subset of
`r-1`. An exponent above `K` cannot enter the interval under this
upper unipotent. The Weyl element gives the corresponding assertion
for the lower unipotent. These subgroups generate `SL_2`, proving
invariance. This is a group-module quotient, not just a matching
of dimensions or characters.

Let `U_delta` be the tensor product of these local quotients with
the indicated Frobenius twists. It is a quotient of

```math
W_\delta=\bigotimes_{j=0}^R
             \operatorname{Sym}^{N_j}(V)^{(j)}.
\tag{8}
```

Its subspace of total weighted first exponent `h` has exactly the
original row-block basis. The quotient intervals ensure that both
occupation vectors are nonnegative; their equal weighted energies
then impose all remaining bounds. Conversely, every original label
pair lies in those intervals because each occupation is at most
`h/w_j<=K_j`.

Thus the row block of (1) is precisely
`pi_h g(z)|_(U_delta,h)`, where `pi_h` projects to this torus
weight space.

## 2. The composition factors that can reach this weight

We may extend the field to an algebraic closure; this does not
change matrix rank. The only standard representation-theory input
used here is the highest-weight classification of simple rational
`SL_2` modules. No complete reducibility is assumed.

For the standard classification and Frobenius tensor description, see
J. C. Jantzen, *Representations of Algebraic Groups*, second edition,
Mathematical Surveys and Monographs 107, AMS (2003),
[publisher record](https://doi.org/10.1090/surv/107).
The actual invariant quotient, weight exclusion, and common filtration
needed in this application are proved below.

The particular simple modules needed have an elementary explicit
description. If `m` has binary digits `m_j`, then

```math
L(m)=\bigotimes_{j:m_j=1}V^{(j)}.
\tag{9}
```

Its weights are `m-2e`, one-dimensional and indexed by bit-subsets
`e` of `m`. To see simplicity directly, its torus weights are all
distinct. A nonzero submodule contains a weight line. Coefficients
of the upper and lower unipotent actions first reach an extremal
line, and then every basis line; distinct subsets give distinct
powers of the unipotent parameter. Hence (9) is simple of highest
weight `m`, and highest-weight classification identifies it with
the corresponding simple module. In particular,
`V tensor L(a)^(1)=L(2a+1)`.

The total polynomial degree of (8) is independent of `delta`:

```math
N=\sum_jw_jN_j
 =\sum_{j=0}^R(2h-w_j)
 =2Rh+1.
\tag{10}
```

Because `sum w_j delta_j=0`, the integer `delta_0` is even.
Therefore `N_0=2h-1+delta_0` is odd. The equivariant monomial map

```math
V\otimes\operatorname{Sym}^{a}(V)^{(1)}
   \longrightarrow\operatorname{Sym}^{2a+1}(V),
\qquad v\otimes f^{(1)}\longmapsto vf^2
```

is an isomorphism: each odd-degree monomial has exactly one odd
exponent. Consequently,

```math
W_\delta=V\otimes (W'_\delta)^{(1)},
\qquad \deg W'_\delta=(N-1)/2=Rh.
```

The latter degree is even. All simple factors of `W'_delta` have
even highest weight, since all its torus weights are even. By (9),
all simple factors of `W_delta`, and hence of `U_delta`, have
highest weight

```math
m\equiv1\pmod4.
\tag{11}
```

Use the torus convention that the weight of a monomial is its
second exponent minus its first. The desired space has weight

```math
\mu=N-2h=2(R-1)h+1\equiv1\pmod{2h}.
\tag{12}
```

If a factor `L(m)` contains this weight, then
`e=(m-mu)/2` is a bit-subset of `m`. Its bounds are `0<=e<=h`,
because `mu<=m<=N`, and (11)--(12) make `e` even.

There is no positive such `e`. Its least set bit would be `2^k`
with `1<=k<=R`. But `m=mu+2e`, while `mu` has no set bits in
positions `1,...,R`. Since `2e` has no set bit below position
`k+1`, the bit `2^k` is absent from `m`, contradicting the
bit-subset condition. This also covers `e=h`.

Thus only `m=mu` contributes, and the relevant weight space is
its extremal line with first exponent zero. On that line the
projected action of (7) is the scalar `g_22(z)^m=1`.

## 3. The common strict filtration

Choose a composition series of `U_delta`. Torus weight spaces
are exact, so its weight-`mu` parts give a filtration of the row
block. Every group-invariant term is preserved by `g(z)` and by
the torus projection. Each nonzero induced quotient is the
one-dimensional line from Section 2, on which the projected
polynomial acts as the constant one.

In a basis adapted to this filtration, `P(z)` is upper triangular
with every diagonal entry equal to one. The basis is independent
of `z`. All positive-degree coefficient matrices are strictly
upper triangular in this same basis, while the constant
coefficient has diagonal one. Applying a coefficient functional
`L` gives diagonal `L(1)` throughout. Summing the block dimensions
gives `d^2`, proving (2).

This establishes the needed nilpotent filtration. Small exact
calculations also find the stronger divided-power multiplication
law for the normalized row coefficients, but that additional law
is not needed or asserted by this proof.

## 4. The column rank and the forced summand

Column blocks are indexed by occupation sums
`sigma=alpha+u=beta+y` of weighted energy `2h`.
The diagonal identity from the protected-minor note is

```math
P(z)(\alpha,u,u,\alpha)
=z^{\sum_j w_j(\alpha_j\mathbin{\&}u_j)}.
\tag{13}
```

For every sum `sigma` whose binary occupation expansion has no
single bit of weighted size `2h`, its bits have sizes at most
`h` and total size `2h`. Such dyadic bits split into two groups
of total size `h`: repeatedly merge equal smaller bits until a
group of size `h` is obtained. This gives `alpha+u=sigma` with
`alpha_j & u_j=0` for every coordinate. Equation (13) then gives
the constant polynomial one at a diagonal entry of that column
block. If `L(1)!=0`, the block is nonzero.

There are `b(2h)-1` sums using only modes of weights at most `h`.
Exactly `R+1` are the excluded single-bit sums, one for each mode.
All other sums admit the displayed split. This proves (3).

For `R>=2`, the elementary recurrence
`b(2h)=b(2h-1)+b(h)` shows `b(2h)-R-2>=b(h)`. Indeed,
`b(2h-1)>=R+2`: use the all-ones partition, one partition with
a single part `2^j` for each `1<=j<=R` and ones otherwise, and
the distinct partition with three parts of size two and ones
otherwise.

Finally, the common diagonal scalar from Section 3 is an intrinsic
linear functional on the tensor span. It takes value one on
`P(1)=E`. In any decomposition of `E`, some summand has nonzero
scalar, hence row rank `d^2` and column rank at least `d` for
`h>=4`. This proves the announced rank-sum bound without any
restriction on the number of summands or on extension-field
coefficients.

## 5. A parity-sector restriction for every energy

Now let `h>=8` be arbitrary. Choose `h'` to be the largest power
of two at most `h/2`, and write

```math
h=2h'+e,\qquad 0\le e<2h',\qquad
e=\sum_j2^j p_j,\quad p_j\in\{0,1\}.
\tag{14}
```

Embed every degree-`h'` occupation label `a` as
`iota(a)=p+2a`. Pad lower labels with zero coordinates when
necessary. This is an injection into the degree-`h` labels.
Restrict all four external labels using this same injection.
The resulting polynomial tensor is exactly

```math
P_h(z)\big|_{\iota^4}=z^e P_{h'}(z^2).
\tag{15}
```

Here is a coordinatewise proof. Put `q=1+z^{2^j}`. If `p_j=0`,
the local polynomial is the square of the smaller polynomial,
and extracting an even exponent replaces `q` by `q^2`. If
`p_j=1`, factor out

```math
(1+t)(1+qt)=1+(1+q)t+qt^2.
```

The remaining factor is a square. Its odd coefficient therefore
selects only `(1+q)t=z^{2^j}t`; the even coefficient left behind
is again the smaller entry with `q` replaced by `q^2`. Multiplying
over coordinates gives the scalar `z^e` and replaces `z` by
`z^2`. The occupation-conservation equation becomes the smaller
one after cancelling the common `p` and dividing by two. This
also verifies zero entries outside the conservation support.

Restriction and coefficient extraction commute. Every restricted
tensor from the coefficient span of `P_h` therefore belongs to
the coefficient span of `P_{h'}`. Restriction sends `E_{b(h)}`
to `E_{b(h')}` and cannot increase either flattening rank.
Sections 1--4 force one restricted summand to have row rank at
least `b(h')^2` and column rank at least `b(h')`. Consequently
every coefficient-span decomposition at energy `h` satisfies

```math
S\ge b(h')^{3/2}.
\tag{16}
```

This is an exact finite lower bound. To express its asymptotic
size, use the binary-partition identity

```math
b(2m)=\sum_{j=0}^m b(j)\le(m+1)b(m).
```

Monotonicity and `h<4h'` give

```math
b(h)\le b(4h')
\le(2h'+1)(h'+1)b(h').
\tag{17}
```

The elementary bound `log b(h)=Omega((log h)^2)` implies
`poly(h)=b(h)^{o(1)}`. For completeness, choose
`L=floor(log_2 h)` and independently set the occupations in
modes `1<=j<=floor(L/2)` to integers below
`2^{floor(L/2)}/L`; their total energy is at most `h`, and the
remaining energy is filled with ones. This gives
`2^{Omega(L^2)}` distinct degree-`h` labels. Equations (16)--(17)
prove

```math
S\ge b(h)^{3/2-o(1)},\qquad \rho\ge\frac34
```

along every sequence of energies tending to infinity, over every
characteristic-two coefficient field.

The finite conclusion (16) need not exceed `b(h)^sqrt(2)` at
small arbitrary energies, so it is not being used as a uniform
finite-seed exclusion. Also, label restriction generally does
not commute with row-matrix multiplication. This corollary does
not extend the separate constant-word-algebra theorem to
arbitrary energies, and it does not cover the full parameterized
algebra generated by both orientations.

## Verification

```sh
python3 tools/check_weighted_weyl_filtration.py --output evidence/weighted-weyl-row-filtration.json
```

The standard-library checker verifies the local invariant-quotient
coefficients, the source-to-`SL_2` coefficient change, the binary
simple-weight exclusion, and explicit no-carry column witnesses.
At energies 2, 4, and 8 it independently constructs the row blocks
and computes a common nilpotent image filtration for every positive
coefficient together with `A_0-I`. Energy 6 is a negative control:
the analogous full-shell filtration fails there. The checker also
verifies the local odd/even coefficient identity, the full
four-label parity restriction at energies 2 through 31, and the
dimension comparison through energy 1024. The all-energy theorem
uses the exact restriction (15), not a full-shell nilpotent
filtration at non-power-of-two energies or an inference from
finite rank tables.
