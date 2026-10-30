# Crossing families from coded central quotients

## Result and scope

This note constructs a family of characteristic-two crossing decompositions
whose noncommutativity is controlled by a binary linear code. It proves an
exact constant-cut reduction, two upper bounds for that cut, a quotient
obstruction, and an upper bound for every row-cut coefficient when central
generators have weight two. The complete column-cut profiles still require
separate estimates. These results do not establish a stronger factoring
bound.

The construction interpolates between extraspecial groups with one central
bit and products of eight-element Heisenberg groups with independent central
bits. It has a useful parameter range: both the constant-cut upper bound and
the known quotient obstruction permit normalized width below
`1/sqrt(2)`. This is a necessary screening condition, not a proof that the
middle coefficients satisfy that width.

## 1. The coded group and its exact crossing family

Let `L:F_2^m -> F_2^t` be surjective, with nonzero columns
`v_1,...,v_m`. Define a group on triples

```math
G_L=\{(a,b,c):a,b\in\mathbb F_2^m,\ c\in\mathbb F_2^t\}
```

by

```math
(a,b,c)(A,B,C)=
(a+A,b+B,c+C+L(a\mathbin{\&}B)).                 \tag{1}
```

Here `&` is coordinatewise multiplication. The bilinearity of the final
term proves associativity. Write

```math
n=2m,\qquad \lambda=n+t,\qquad d=2^\lambda.
```

The horizontal generators `X_i=(e_i,0,0)` and `Y_i=(0,e_i,0)` and
central generators `Z_j=(0,0,e_j)` are involutions. Distinct pairs commute,
and

```math
[X_i,Y_i]=Z^{v_i}.
```

Every element has a unique normal word, with the `X` generators first,
then the `Y` generators, and then the `Z` generators. The normal-word
bits of `(a,b,c)` are `(a,b,c+L(a&b))`.

Choose positive integer generator weights `w_j` and put

```math
F_S(z)=\prod_{j=1}^{\lambda}
       \bigl(g_j+z^{w_j}(1+g_j)\bigr)^{s_j},
\qquad K_z(S,a)=[a]F_S(z),                       \tag{2}
```

where `s` is the normal word of `S`, and the product follows the same
fixed order. Define the scalar four-leg tensor

```math
T_z(x,y,u,v)=[xy^{-1}=vu^{-1}]K_z(xy^{-1},y^{-1}u).
                                                        \tag{3}
```

At `z=1`, every active factor in (2) is the identity element. Therefore

```math
T_1(x,y,u,v)=[x=v][y=u]=E_d(x,y,u,v).              \tag{4}
```

Writing `T_z=sum_i z^i A_i` gives an exact crossing decomposition
`E_d=sum_i A_i`. The number of pieces is at most `1+sum_j w_j`.
Unit weights and weights one on horizontal generators and two on central
generators both give `O(lambda)` pieces.

For a uniformly computable sequence of matrices `L`, all entries of this
family and exact rank factorizations of its scalar cuts can be constructed
in `2^O(lambda)` time. This follows by constructing matrices of order
`d^2` and performing exact elimination. It does not require an uncharged
search for a good code. A concrete uniform code candidate can, for example,
use columns `1,alpha,...,alpha^(m-1)` in `F_(2^t)`, where `alpha` is a
deterministically selected primitive element and `m>=t`. Selecting the
field and primitive element by exhaustive finite-field arithmetic costs
`2^O(t)`. The first `t` columns already span. No complete profile bound is
claimed for this candidate.

## 2. Exact constant-cut compression

At `z=0`, (2) is the normal word for `S`, so `K_0(S,a)=[a=S]`.
The row cut is a direct sum of `d` permutation matrices of order `d`:

```math
r_0=d^2.                                        \tag{5}
```

For a horizontal vector `h=(a,b)`, define the square and polar maps

```math
Q(h)=L(a\mathbin{\&}b),\qquad
B(h,H)=L(a\mathbin{\&}B+A\mathbin{\&}b),
\quad H=(A,B).                                  \tag{6}
```

The square of `(a,b,c)` is `(0,0,Q(h))`; conjugation changes only its
central coordinate, by `B(h,H)`.

**Proposition.** The constant column-cut rank is exactly

```math
c_0=\sum_{z\in\mathbb F_2^t}\operatorname{rank}M_z,
\qquad
M_z((h,\delta),H)=[\delta=B(h,H)],\quad Q(h)=z.    \tag{7}
```

To prove this, the constant support in (3) is equivalently
`u=yxy^(-1)` and `v=x^2y^(-1)`. In the column cut `(x,u)|(y,v)`, the
central value `z=x^2=yv` indexes disjoint blocks. Within a block, only
the horizontal coordinates of `x,y` and the central difference between
`x,u` matter. The unused central coordinates duplicate rows and columns.
Deleting duplicates leaves precisely (7).

Each row of `M_z`, as a function of `H`, is the Boolean polynomial

```math
\prod_{j=1}^t(1+\delta_j+B_j(h,H)).                \tag{8}
```

It has degree at most `t` in `n` variables. Consequently

```math
c_0\le 2^t\sum_{j=0}^t\binom{n}{j}.              \tag{9}
```

The polynomial representation is over the quotient by `H_i^2=H_i`;
these multilinear monomials are a basis for all Boolean functions.

## 3. A sharper bound from the column matroid

Let

```math
I(L)=\sum_{J\subseteq\{1,\ldots,m\}\atop
             \{v_j:j\in J\}\ \mathrm{independent}}2^{|J|}.
                                                        \tag{10}
```

This is the number of independent subsets of the column matroid after
replacing each column by two parallel copies, one for each horizontal
coordinate.

**Proposition.**

```math
c_0\le 2^t I(L)
   \le 2^t\sum_{j=0}^{t}\binom{m}{j}2^j.          \tag{11}
```

For fixed `h`, every row in (8) is a Boolean function of
`sum_i v_i(a_i B_i+b_i A_i)`. Consider the coefficient of a horizontal
monomial supported on coordinates `P`. Set all other coordinates to
zero. If the corresponding column vectors have rank less than `|P|`,
the restricted function factors through fewer than `|P|` independent
linear forms. Its algebraic degree is therefore less than `|P|`, so
that coefficient is zero. Thus only independent sets of the parallel
extension can occur. Their monomials span a space of dimension `I(L)`
in each fiber `z`, proving (11).

This proof also gives a smaller exact rank computation. As `delta`
varies, the rows in (8) span all products
`prod_(j in J) B_j(h,H)`, `J subseteq {1,...,t}`. Expand these products
in the `I(L)` permitted monomials and eliminate their coefficient rows.
The resulting rank equals the original evaluation-matrix rank. It is
not merely an upper estimate.

For a single central bit, all columns are one. Then `I(L)=1+2m`.
For `m>=2`, every square fiber spans the full horizontal space, and
the two ranks in (7) both equal `2m+1`. Hence

```math
c_0=4m+2\quad(m\ge2),\qquad c_0=5\quad(m=1).       \tag{12}
```

For the zero square fiber, the standard basis vectors suffice. For
the one square fiber and any coordinate in pair `i`, choose another
pair `j`. Both `e_(a_j)+e_(b_j)` and that vector plus the chosen
coordinate have square one, and their sum is the chosen coordinate.
This proves the spanning assertion without a finite extrapolation.

## 4. What the abelian quotient forces

The quotient of `G_L` by its center is the binary vector space
`V=F_2^n`. There is a four-leg linear projection from every coefficient
in (3) to the corresponding abelian word coefficient on `V`.

Fix the central coordinates of `x,y` to zero, and sum independently
over the central coordinates of `u,v`. When horizontal conservation
holds, choosing the center of `u` determines the center of `v`.
The resulting kernel is the group-algebra projection

```math
\overline K_z(\overline S,\overline a)
   =\sum_{a:\ a\mapsto\overline a}K_z(S,a).       \tag{13}
```

Every central factor in (2) projects to one. The horizontal control
bits depend only on `overline S`, so (13) is independent of the lift
`S`. It is the ordinary abelian binary-word family with the chosen
horizontal weights. The projected endpoint is `E_(2^n)`.

Neither flattening rank increases under a product of four leg maps.
When all horizontal weights are equal, the quotient pieces are
invariant under all `n` binary-position permutations. Arbitrary
linear recombinations of the original coefficients retain this
property after projection. The
[binary-position invariant theorem, pinned in PR #26](https://github.com/1337isnot1337/integer-factoring/blob/023e4fc4482d72f25738f6b1d0f9cd31583baed1/research/modular-binary-word-barrier.md)
therefore forces one common piece to have both ranks at least

```math
2^{3n/2-O((\log(n+1))^2)}.                         \tag{14}
```

Thus any such complete family has normalized bottleneck

```math
\rho\ge\frac{3n}{4(n+t)}-o(1).                   \tag{15}
```

In particular, a single central bit cannot produce an asymptotic
improvement through this interface, despite the polynomial constant
rank in (12). A linear number of central bits is needed to leave this
specific obstruction behind.

For `t/n -> epsilon` with `0<epsilon<1/2`, (9) gives the constant
piece upper bound

```math
\rho_0\le\frac12+
 \frac{\epsilon+H_2(\epsilon)}{4(1+\epsilon)}+o(1).
                                                        \tag{16}
```

The interval approximately `0.0606602<epsilon<0.2357917` makes this
upper bound and the lower bound (15) both less than `1/sqrt(2)`.
For example, `t=m/4`, so `epsilon=1/8`, gives the respective values
`0.6485698763` and `2/3`. The matroid bound can improve the first
number. None of these comparisons bounds the complete profile.

## 5. An upper bound for every row coefficient

Give horizontal generators weight one and central generators weight
two. Fix a normal control word `S`. Let `k` be its number of active
horizontal generators, and let `r` be the rank of the span of:

- The standard central vectors explicitly active in that word.
- The columns `v_i` for pairs in which both `X_i` and `Y_i` are active.

The subgroup `H_S` generated by the active generators has order
`2^(k+r)`. The group algebra of `G_L` is a direct sum of
`d/2^(k+r)` copies of its regular module for `H_S`.

Choose a basis of the intersection of `H_S` with the ambient central
coordinate subgroup. This has dimension `r`; unpaired horizontal
generators can also be central in `H_S`, but they are already among
the `k` horizontal generators. The ordered monomials in
`1+g` for its `k` horizontal generators and `r` central basis
generators form a vector-space basis. Give their exponents weights
one and two, respectively. Monomials of weight at least `a` form
a multiplicative descending filtration. Indeed, commuting a
horizontal pair introduces a factor `1+Z^v`, whose central expansion
has weight at least two; all other correction terms have still
higher weight. Squares of the individual basis factors vanish.

The Hilbert coefficients of this filtered space are

```math
h_l=[x^l](1+x)^k(1+x^2)^r.                        \tag{17}
```

The coefficient `[z^i]F_S(z)` belongs to filtration degree at least
`i`. Multiplication by it maps degree at least `a` into degree at
least `a+i`. Separating the inputs of degree below `a` proves the
rank upper bound

```math
R_i(k,r)=\min_{a\ge0}
 \left(\sum_{l<a}h_l+\sum_{l\ge a+i}h_l\right).    \tag{18}
```

If the weighted degree `W(S)` of the active control word is less
than `i`, that coefficient is zero. Using the exact row-block
formula for orbit kernels therefore gives

```math
r_i\le
 \sum_{S:\ W(S)\ge i}
  \frac{d}{2^{k(S)+r(S)}}\,R_i(k(S),r(S)).          \tag{19}
```

This bound is computable from active-set statistics of `L`, without
constructing the full crossing matrices. It is exact at the two
ends: it gives `r_0=d^2` and `r_(n+2t)<=1`. It is an upper bound in
the middle, not an assertion that the relevant multiplication maps
have maximal rank.

The same weighted monomial argument applies directly to the graded
restricted enveloping algebra with square-zero primitive generators
and relation `[y_i,x_i]=sum_j L_(j,i) z_j`. Its multiplication has
no higher-degree group-algebra correction terms. Turning that
simplification into complete crossing-rank estimates is a separate
step.

## 6. Reproduction

The checker compares the compressed polynomial matrices with
independently constructed function-evaluation matrices in small
cases. It also compares (19) with actual convolution-matrix ranks.

```sh
python3 tools/check_noncommutative_crossing_coded.py --larger \
  --output evidence/noncommutative-crossing-coded.json
```

The evidence contains the code matrices, exact fiber ranks, matroid
bounds, and the finite row checks. Larger constant-cut cases can
stop once every fiber spans the proved ambient monomial space;
this is a complete rank certificate, rather than a partial-rank
estimate.
