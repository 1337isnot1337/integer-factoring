# Polynomial crossings from an arbitrary frozen-core output row

## Result and scope

There is an eleven-parameter affine family of polynomial crossing tensors
obtained by freezing all four extreme cores low and varying a single
height-dressed output covector. The only endpoint constraint is

\[
g_0+g_{10}=1.
\tag{1}
\]

The construction works over every field of characteristic two and for
every finite mode set in the existing rectangular construction. It
contains the Gaussian all-low background and the all-low restriction of
the full inverse repair. A small exact coefficient identity shows that
the family also contains tensors outside the standard affine Gaussian
coefficient model defined in Section 3, even over the rational function
field.

This gives a concrete way to leave the fixed-background hypothesis of
the [all-scale obstruction](all-scale-balanced-obstruction.md). It does
not establish a pressure below `1/sqrt(2)` or a new factoring exponent.
The family has only the bulk alphabet, so retaining a useful entropy
gain remains a separate problem.

The construction uses the one-particle matrix, residual endpoint and
height convention of
[coordinate-core-crossings.md](coordinate-core-crossings.md) and
[rectangular-residual-repair.md](rectangular-residual-repair.md).
The affine-row observation and its exact checks used OpenAI assistant
support.

## 1. The polynomial family

Let `S` be a finite mode set containing `0` and `T=max S>0`, and choose a
power of two `M>T`. Put `c=1+w`. Let `G` be the Gaussian tensor formed
from all minors of the direct one-particle matrix, and set

\[
Q=x_0y_Tu_0v_T,\qquad K=G(1+Q/c^2).
\tag{2}
\]

On the source output cores, use the rectangular states

\[
v:\quad (\mathrm{low},\varnothing,\mathrm{top},\mathrm{both}),
\qquad
y:\quad (\mathrm{low},\varnothing,
         \mathrm{top}+\beta\,\mathrm{both}).
\tag{3}
\]

Both input cores `x,u` are fixed low. The mixed input parameter `alpha`
therefore has no effect. All bulk modes in `S\setminus\{0,T\}` remain
free on all four legs. Each leg has alphabet size

\[
b=2^{|S|-2}.
\]

Index a source output core pair by `r=v+4y`. Its height in units of `T`
is

\[
h(r)=h_X(v)+h_U(y),\qquad
h_X=(0,0,1,1),\quad h_U=(0,0,1).
\tag{4}
\]

Write `K_{r,0}` for the bulk tensor with this source output pair and
with the low-low input pair. For arbitrary field coefficients
`g=(g_0,...,g_11)`, define

\[
\boxed{\quad L_g=\sum_{r=0}^{11}g_r c^{h(r)}K_{r,0}.\quad}
\tag{5}
\]

Here the output logical core is also declared low-low. The summands in
(5) are components of one output covector; no source core remains in
the physical alphabet.

**Polynomiality.** Every coefficient of `L_g` is a polynomial in `w`,
without any type restriction and without condition (1).

The Gaussian part of each summand has the nonnegative multiplier
`c^{h(r)}`. A coefficient of `QG` can contribute only when both source
outputs contain their top mode. Such a source has `h(r)=2`, including
the `both` term of the mixed `y` covector. Its factor `c^2` cancels the
only denominator in (2). All Gaussian minors are polynomial, so every
summand is polynomial separately. This argument does not require a
full inverse matrix or an invertible completion of the chosen row.

## 2. The complete endpoint condition

The extreme residual matrix of the rectangular construction is
`F(alpha,beta)`. Its low-low input column satisfies the exact identity

\[
F(\alpha,\beta)e_0=e_0+e_{10}.
\tag{6}
\]

This follows either from the seven literal exterior factors in the
residual endpoint or from its coefficient matrices. None of the
`alpha`, `beta`, or `alpha beta` coefficient matrices has a nonzero
entry in column zero.

Make the ramification `w=1+t^T`. On each bulk physical letter `a`, put

\[
H(a)=\sum_{j\in S\setminus\{0,T\}}j\,n_j(a).
\]

Multiply a coefficient of `L_g` by

\[
t^{H(y)+H(v)-H(x)-H(u)}.
\tag{7}
\]

The factor `c^{h(r)}=t^{T h(r)}` in (5) supplies exactly the source
core part of the original regularizing index gauge. Thus the gauged
summands are regular in `t`, and their constant tensors are the
corresponding entries of `F e_0` times the bulk crossing. Equation (6)
gives

\[
\widehat L_g(0)=(g_0+g_{10})\,[x=v][u=y].
\tag{8}
\]

Consequently (1) is necessary and sufficient for this gauged local
endpoint to equal the crossing. It is one nontrivial affine equation
on twelve coefficients, leaving eleven free parameters. No constraints
on the unused columns of a full residual repair are needed.

Now take any positive tensor power and restrict `x,v` to one common
full-letter word type and `u,y` to an independent common full-letter
word type. The total exponent in (7) cancels on each crossing pair.
The resulting tensor is already polynomial by Section 1, and evaluating
it at `w=1` gives the exact rectangular crossing on the two type classes.

The ungauged local value at `w=1` need not itself be the identity
crossing; bulk Pascal transitions can remain before the type
restriction. Equation (8) is the local gauged statement used in the
proof. It also proves full generic rank of both local flattenings when
(1) holds, because their regular gauged specializations are invertible.

The polynomial degree is finite for every fixed seed. The standard
coefficient count and the paired-type rank-cost transfer therefore
apply. Profiles must be computed for the actual tensor `L_g`.

## 3. Included backgrounds and a non-Gaussian member

The row `g=e_0` is the existing all-low Gaussian background. A quartic
source cannot occur in that row. For a full interval, its two profiles
are those of a Gaussian interval of length `M-2`, as proved in
[rectangular-low-background.md](rectangular-low-background.md).

Taking `g` to be row zero of `F(alpha,beta)^{-1}` gives the all-low
restriction of the full inverse repair. Its target core has height
zero, so the inverse repair's height dressing is precisely (5).
Equation (6) also immediately verifies (1) for this row.

For example, take `S={0,1,2,3}`, `M=4`, `beta=1`, and

\[
g=e_3+e_9+e_{10}.
\tag{9}
\]

This is row zero of the inverse at `alpha=beta=1`, where the partial
repair is unavailable. It already has coefficients in the binary
field. Let `a_I(w)` denote its coefficient of the bulk exterior
monomial indexed by `I`. Its vacuum coefficient is

\[
a_\varnothing=w^2.
\tag{10}
\]

More generally this vacuum calculation gives
`1+c^2 w^{2(M-1-T)}` for the same row at arbitrary ambient `M>T`;
equation (10) uses `T=M-1`.

Consider the affine Gaussian coefficient model

\[
a\prod_i(1+\ell_i z_i)\prod_{i<j}(1+q_{ij}z_i z_j),
\]

where the `z_i` are square-zero exterior variables and all coefficients
may lie in the rational function field. Every tensor in this model must
satisfy, for any three distinct variables `i,j,k`,

\[
a_\varnothing a_{ijk}
 +a_{ij}a_k+a_{ik}a_j+a_{jk}a_i=0.
\tag{11}
\]

To verify necessity, expand the factors on these three variables.
Terms using variables outside the triple cannot contribute. The three
copies of the product of the linear coefficients equal one copy in
characteristic two, giving (11).

For the row (9), direct Gaussian minors give the exact polynomial
identity

\[
\boxed{\quad
a_\varnothing a_{x_1x_2y_1}
 +a_{x_1x_2}a_{y_1}
 +a_{x_1y_1}a_{x_2}
 +a_{x_2y_1}a_{x_1}
 =w^6(1+w)^{14}\ne0.
\quad}
\tag{12}
\]

Thus this member is not of that affine Gaussian form, even after
allowing rational linear and quadratic coefficients. The Gaussian
row `e_0` gives zero in the same control. The nonzero identity is an
algebraic distinction; it does not imply a favorable pressure.

## 4. Exact controls and remaining search

The constructor is
[rectangular_frozen_core_rows.py](../tools/rectangular_frozen_core_rows.py).
It constructs twelve polynomial source components by direct minors,
then combines them over `F4`. The implementation enforces (1) when a
crossing tensor is requested.

Run the independent finite controls with

```sh
python3 tools/check_rectangular_frozen_core_rows.py \
  --output evidence/polynomial-frozen-core-crossings.json
```

The [evidence](../evidence/polynomial-frozen-core-crossings.json) checks
the residual column for all sixteen parameter pairs, every source
component's regular ramified endpoint, independent two-copy types for
dense rows, and coefficient equality with separately constructed
partial and inverse seeds. It includes a gapped mode set and ambient
parameters larger than the maximum mode. It also checks (12) and
rejects a row with zero endpoint.

The finite controls support the ordinary proofs in Sections 1 and 2;
they do not replace those proofs. Initial floating searches have not
found an improvement over the included Gaussian member. The open task
is to identify a row, mode set, and gauge whose complete pressure is
strictly below the critical value, or to prove an obstruction for this
larger family.
