# Original-field multiplication blocks and a bundled middle-layer algorithm

## Result

The [monoid block theorem](multiplicative-monoid-blocks.md) admits a
stronger form over the **original** coefficient field. Odd characters
need not be split into scalar characters. Their polynomial coordinate
can be retained inside the `M`-dimensional message space, with a binary
change of basis and no extra coefficient-field degree.

For the adaptive modulus `M=2^m+1`, let `omega(M)` count its distinct
prime factors. There is an explicit decomposition of the residue
multiplication tensor using at most

```math
R'\le M8^{\omega(M)}
 \le M^{1+3/\log_2(2D+1)}                            \tag{1}
```

bilinear products over the existing root field `F=F_(2^(2m))`.
A power-of-two grid version, including its linear evaluation and
interpolation work, computes a complete multiplicative convolution in

```math
M64^{\omega(M)}\operatorname{polylog}M
 \le M^{1+6/\log_2(2D+1)}\operatorname{polylog}M       \tag{2}
```

coefficient operations after the explicit field and modulus
preparation. The usual `poly(m,ell)` bit cost per operation is charged
when using the recovery algebra `F tensor F_(2^ell)`.

There is also a precise consequence for the **whole** source `MID`
vertex. Treat its `A` first-operand inputs, `A` second-operand inputs,
and `A` outputs as three bundles. Its three-way tensor rank satisfies

```math
M^A\le\mathrm{rank}_{\mathrm{CP}}(\mathrm{MID})
       \le(R')^A=M^{A(1+O(1/\log D))}.              \tag{3}
```

Its contraction against three dense bundle vectors has the analogous
bound from (2), with `M` replaced by the bundle dimension `M^A`
and `omega(M)` replaced by `A omega(M)`.

These are ordinary constructive results. A bundle contains all of its
`A` residue coordinates; the factors in (3) need not separate among
the individual neighboring gates. Section 7 explicitly accounts for
that distinction. The complete factoring score remains `1/2`.

## 1. Descent without a scalar extension

Fix `q=p^e`, with `p` odd, and write `p-1=dh`, where `d` is a
power of two and `h` is odd. Define binary polynomials

```math
\begin{aligned}
\Psi_0(Z)&=Z^h-1,\\
\Psi_r(Z)&=
 \frac{Z^{hp^r}-1}{Z^{hp^{r-1}}-1}
 =\sum_{j=0}^{p-1}Z^{jhp^{r-1}},\qquad r\ge1,\\
c_r&=\deg\Psi_r=h\varphi(p^r).
\end{aligned}                                                    \tag{4}
```

The polynomials for `0<=r<=e-1` are pairwise coprime and their
product is `Z^(hp^(e-1))-1`. This follows directly from the
squarefreeness of odd-order cyclotomic polynomials in characteristic
two. A root of `Psi_r` has `p`-primary order exactly `p^r`, while
its `h`-part is unrestricted.

The isomorphism below already exists over `F_2`, and therefore over
every coefficient algebra of characteristic two:

```math
F[\mathbb Z/p^e]\simeq F\ \times
\prod_{r=0}^{e-1}
 F[Z,X,Y]/(\Psi_r(Z),X^d,Y^{e-r}).                  \tag{5}
```

The algebra on the left uses multiplication of residues. Choose
the same unit generators `u=g^(hp^(e-1))`, `v=g^d` as in the
previous note. For `a=p^j a_0!=0`, write `a_0=u^s v^t`. Its
component in factor `r` is

```math
[a]\longmapsto
\begin{cases}
 Z^t(1+X)^sY^j\pmod{\Psi_r,X^d,Y^{e-r}},&j<e-r,\\
0,&j\ge e-r.
\end{cases}                                                    \tag{6}
```

The first component is the augmentation, and `[0]` has all other
components zero. This is the same homomorphism as the split
character map, with all characters of `p`-primary order `p^r`
retained in the single algebra `F[Z]/Psi_r`.

One can prove (5) directly by CRT in the odd unit-group algebra.
Equivalently, extend faithfully to the explicitly split algebra of
the previous note and use its proven isomorphism there. Faithful
extension shows that the already binary map (6) is an isomorphism
before extension as well. Its dimension identity is

```math
1+\sum_{r=0}^{e-1}c_r d(e-r)=p^e.                  \tag{7}
```

Crucially, `c_r` is now a coordinate-axis length inside (7).
It is not a new factor in the bit length of every scalar. The
coefficient field for (1)--(3) is still `F`.

The checker includes `p=7,23,29` over fields containing the
required additive roots but not all odd multiplicative characters.
The binary change of basis remains invertible in each case.

## 2. The basis maps have small linear circuits

Here is a complete algorithm for the basis maps, not only their
existence. On valuation layer `j`, let
`H_j=hp^(e-1-j)`. Arrange the unit coefficients by their even
exponent `s<d` and odd exponent `t<H_j`.

Apply the Boolean zeta transform in `s`. It has `O(d log d)`
additions per odd exponent, and its inverse is the same transform
in characteristic two. The resulting coefficient of `X^k Z^t`
is the required binomial moment. Reduce its polynomial in `Z`
modulo every `Psi_r` with `r<=e-1-j`.

These reductions have a simple linear-time implementation. For
`r>=1`, put `P=hp^r`, `Q=hp^(r-1)`, so `c_r=P-Q`.
First fold exponents modulo `P`, since `Psi_r` divides `Z^P-1`.
The remaining `Q` coefficients in degrees `c_r,...,P-1` are
removed with the `p` terms of `Psi_r`. This costs `O(P)` additions,
and `P<=3c_r/2`. For `r=0`, merely fold modulo `h`.

The inverse CRT is equally explicit. In the cyclic algebra of
length `H_j`, define

```math
E_{\le r}(Z)=\sum_{a=0}^{H_j/(hp^r)-1}Z^{a hp^r},
\qquad
E_0=E_{\le0},\quad E_r=E_{\le r}+E_{\le r-1}.       \tag{8}
```

Each averaging length is odd, so no nonunit scalar is inverted.
The element `E_r` is one on the roots of `Psi_r` and zero on
the other factors. A representative of degree below `c_r` is
embedded by multiplying by `E_r`. This multiplication consists
only of a periodic repetition and, for `r>0`, a second periodic
repetition after folding modulo `hp^(r-1)`. It costs `O(H_j)`
additions. Sum these embedded components and invert the Boolean
zeta transform.

This gives both local basis maps in
`p^e polylog(p^e)` binary-coefficient additions. Recover the
coefficient of `[0]` from the augmentation. For a general odd
`M`, use CRT to permute the original residue coordinates and
apply the local maps along each tensor axis. The total cost is
`M polylog M`, up to a fixed power of the logarithm. Coordinate
permutations and table lookups can be implemented by sorting;
they do not require unit-cost random access to exponentially
long words on a finite-alphabet machine.

The exact checker implements these periodic CRT maps and compares
them with an independent full binary basis matrix and its inverse.

## 3. Exact bilinear counts over the original field

A local block in (5) has three coordinate lengths
`c_r,d,e-r`. Evaluate the two input polynomials at grids of
sizes `2c_r-1`, `2d-1`, and `2(e-r)-1`. Multiply pointwise,
interpolate, reduce in `Z` using `Psi_r`, and truncate in `X,Y`.
The exact bilinear count of this construction is

```math
R'_{p^e}
=1+\sum_{r=0}^{e-1}
 (2c_r-1)(2d-1)(2(e-r)-1).                         \tag{9}
```

CRT gives `R'=prod_(p^e||M) R'_(p^e)`. In particular,

```math
R'_{p^e}\le8p^e,\qquad R'\le M8^{\omega(M)}.        \tag{10}
```

The first bound follows by comparing each summand with eight
times its block dimension. One can sharpen it to
`R'_(p^e)<=4(1+1/p)p^e` using the sum in the earlier split count.

All interpolation nodes lie in the original root field. An axis
has length at most `M-1`. In the source family, `M=2^m+1` and
`m>=2`, so the `2c-1` distinct nodes needed for any one axis
fit inside `F_(2^(2m))`. Their differences are nonzero field
elements and hence units after tensoring with the recovery field.

As proved in the previous note, every prime divisor of the source
modulus is at least `2D+1`. Thus
`omega(M)<=log2 M/log2(2D+1)`, giving (1). The flattening with
one multiplicative input fixed to the identity has rank `M`.
Hence this product count is optimal to leading order in
`log_M R'` for the scalar multiplication tensor.

## 4. Charge evaluation, interpolation, and reduction

The arbitrary-node count (9) is useful as a tensor decomposition.
To give a fast algorithm for all its linear operations as well,
round each grid size up to a power of two. This section gives an
explicit additive-subspace algorithm with no assumption about
power-of-two multiplicative roots in characteristic two.

Take linearly independent `alpha_0,...,alpha_(b-1)` in `F`
over `F_2`, and let `V_j` be their first `j`-dimensional span.
The subspace polynomial

```math
s_j(X)=\prod_{v\in V_j}(X-v)
```

is linearized, with at most `j+1` monomials. Starting with
`s_0(X)=X`, it obeys

```math
s_{j+1}(X)=s_j(X)^2+s_j(\alpha_j)s_j(X),
\qquad s_j(\alpha_j)\ne0.                           \tag{11}
```

To evaluate a polynomial `f` of degree below `K=2^j` on
the affine subspace `a+V_j`, divide it by the sparse polynomial
`s_(j-1)(X-a)`:

```math
f=q\,s_{j-1}(X-a)+r,\qquad\deg q,\deg r<K/2.
```

On the first half of the subspace, evaluate `r`. On the second
half, evaluate `r+cq`, where `c=s_(j-1)(alpha_(j-1))` is a
known nonzero field element. Sparse long division costs
`O(K log K)` operations at this level, so the whole recursion
costs `O(K log^2 K)`.

Interpolation reverses these steps. Interpolate the two halves
to `r` and `r+cq`, subtract them, multiply by `c^(-1)`, and
reconstruct `f`. Its cost is the same. This proof applies to
coefficients in any commutative `F`-algebra because every scalar
division is by a nonzero element of `F`.

For an axis of length `c`, choose

```math
K_c=2^{\lceil\log_2(2c-1)\rceil}<4c.
```

All such subspaces fit in the existing field: their dimensions
are at most `m+1<=2m`. Apply the evaluation/interpolation
recursion along every tensor-grid axis. A block with `a` axes
has fewer than `4^a` grid points per input coordinate. Its
complete operation cost is that grid size times a fixed
polynomial in `log M`.

Each global block has at most `3 omega(M)` axes. Consequently,
the sum of the rounded grid sizes is at most

```math
\widehat R\le M64^{\omega(M)}.                      \tag{12}
```

The `Z` reductions are the linear-time folds from Section 2;
the `X,Y` truncations just discard coefficients. Thus (12)
charges evaluation, pointwise multiplication, interpolation,
and reduction. Together with the binary basis maps, it proves
(2), not merely its nonlinear multiplication count.

The coefficient field and a primitive additive `M`th root are
the already explicit objects from the
[carry Fourier construction](carry-fourier-elimination.md).
If constructing them afresh, retain its `2^O(m)` preparation
term. Trial factoring of `M`, unit generators, and all residue
permutations are also explicit and fit `2^O(m)` preparation.
The bound (2) is the subsequent cost for arbitrary dense input
vectors; it does not erase that preparation term.

The subspace-polynomial method is standard additive FFT
technology. For primary context, see Shuhong Gao and Todd
Mateer, [*Additive Fast Fourier Transforms over Finite Fields*](https://www.math.clemson.edu/~sgao/papers/GM10.pdf),
IEEE Transactions on Information Theory 56 (2010), 6265--6272,
DOI `10.1109/TIT.2010.2079016`. The elementary recursion and
the cost needed here have been given above in full.

## 5. The Fourier closure descends too

Write `T_(r,v)(Z,X)` for the coefficient of `Y^v` in (5),
and let `z=f(0)`. Define

```math
\iota_X(X)=\frac X{1+X},\qquad
\iota_Z(Z)=Z^{-1}\pmod{\Psi_r},\qquad
J(Z)=1+Z+\cdots+Z^{h-1}\pmod{\Psi_0}.               \tag{13}
```

Inversion preserves each `Psi_r`. The polynomial `J` is the
trivial-character projector inside `F[Z]/Psi_0`: it is one at
`Z=1` and zero at every other `h`th root. No scalar character
extension is needed to construct or use it.

Let `G_r(Z,X)` be the unsplit Gauss polynomial

```math
G_r(Z,X)=
\sum_{a\in(\mathbb Z/p^{r+1})^\times}
 \zeta_{p^{r+1}}^a Z^{t(a)}(1+X)^{s(a)}
             \pmod{\Psi_r,X^d},                   \tag{14}
```

where `a=u^s v^t` uses the same unit coordinates as (6).
For `r>=1`, Fourier has the exact formula

```math
T_{r,w}(\mathcal F f)
=G_r(Z,X)\,
 \iota_Z\iota_X\bigl(T_{r,e-r-1-w}(f)\bigr).        \tag{15}
```

For `r=0`, the formula includes all tame characters, and the
trivial one has the following explicit correction:

```math
\begin{aligned}
T_{0,w}(\mathcal F f)
={}&G_0\,\iota_Z\iota_X\bigl(T_{0,e-1-w}(f)\bigr)\\
&+J(Z)X^{d-1}
 \left(z+\sum_{v=e-w}^{e-1}T_{0,v}(f)(1,0)\right),\\
z(\mathcal F f)
={}&z+\sum_{v=0}^{e-1}T_{0,v}(f)(1,0).
\end{aligned}                                                    \tag{16}
```

These identities follow either directly from the unit-kernel
argument in the previous note or by applying its character
formulas after a faithful split extension. All coefficients and
operations in (15)--(16) already lie over `F`, so the equality
descends. The precise unit identity is

```math
G_r\,\iota_Z\iota_X(G_r)
=(1+X)^{d/2}+[r=0]J(Z)X^{d-1}.                    \tag{17}
```

The right side squares to one. This gives the inverse of every
Gauss block over the original coefficient field, including when
`F[Z]/Psi_r` is not a field.

These linear operations also have explicit small circuits.
For `Z` inversion, reverse exponents modulo `hp^r` and use
the linear reduction in Section 2. For `X` inversion, change
from `X` to `U=1+X` using the Boolean zeta transform, permute
the exponents of `U` by `j -> -j mod d`, and change back.
The valuation reversal is a permutation. Multiplication by a
Gauss polynomial uses the two-axis case of Section 4.

The Gauss polynomials themselves can be prepared by arranging
their root powers in unit-coordinate arrays, applying the same
Boolean transform, and reducing in `Z`. Their total size is
polynomially bounded by the original `M`-coordinate space.
The local Fourier circuit therefore takes
`p^e polylog(p^e)` coefficient operations; CRT tensoring gives
an `M polylog M` circuit up to a fixed logarithmic power.
This is a closed original-field version of the earlier
Fourier--Mellin interface.

## 6. The entire MID tensor, with its factors and costs

The source's middle vertex at fixed `k_1` is

```math
\mathrm{MID}(a,b)
=L_{\mathrm{out}}\bigl(L_P(a)\mathbin{\odot}L_Q(b)\bigr),
                                                               \tag{18}
```

where `a,b` are `A`-vectors of residues, `L_P=L_Q` is the
length-`A` DFT, `odot` is coordinatewise residue multiplication,
and `L_out` is the inverse-stage DFT with its diagonal twiddle.
The principal-root argument gives invertibility of all three
linear maps over the actual composite residue ring. A linear
bijection of residue vectors permutes the corresponding
indicator basis over the coefficient field.

Consequently the indicator tensor of (18) is the `A`th tensor
power of the scalar residue multiplication tensor, followed by
three explicit bundle permutations. Tensoring the decomposition
from (9) proves the upper bound in (3).

For the lower bound, fix the second transformed operand to the
all-ones residue vector. Then the first input maps bijectively
to the output. The resulting flattening contains an identity
matrix of order `M^A`. Thus its rank is at least `M^A`.

The decomposition is effective. Prepare the scalar input and
output factor arrays using the binary basis maps and Lagrange
interpolation in Section 3. Their complete written size is
`2^O(m)`. A bundle factor is a product of `A` scalar factors
after applying one of the residue-vector permutations in (18).
Hence its queries require the charged scalar-table access and
polynomial work in `A,m`; there is no implicit tensor oracle.

Given three arbitrary dense bundle vectors, permute their
entries according to (18), perform one multiplicative
convolution on `(Z/M)^A`, and take the remaining dot product.
The tensor-product version of Sections 2--4 gives operation
count

```math
M^A64^{A\omega(M)}\operatorname{poly}(A,\log M)
=M^{A(1+O(1/\log D))}\operatorname{poly}(A,\log M).  \tag{19}
```

This is a concrete improvement over enumerating the `M^(2A)`
input settings for that **three-bundle contraction problem**.
The checker constructs every coefficient of the actual
`M=5,A=2` middle tensor at both `k_1` values and compares its
explicit bundled decomposition against the source formula.

## 7. Consequence and limit for the complete constructor

The input size in (19) includes three dense vectors of length
`M^A`. Neighboring vertices in the factoring network do not
generally supply three independent vectors with that grouping:
one neighboring gate is shared among many middle vertices.
The resulting joint arrays must still be represented and
contracted.

A literal insertion of the CP factors makes the cost visible.
Replace one middle vertex by three bundle factors and use `A`
copy vertices for the scalar rank indices. Each bundle factor
has `A` original residue coordinates and `A` rank coordinates.
The original `3A log M` wire bits remain, and the new three
rank incidences add `3A log R'` wire bits. The setting bits of
these factors and copy vertices total

```math
3A\log_2M+4A\log_2R',                              \tag{20}
```

whereas the original deterministic middle setting used
`2A log2 M` bits. Rounding each individual alphabet upward
adds only its explicit ceiling, not a hidden saving. Thus
inserting this decomposition into the old equality-message
network does not by itself improve its written-table bound.

Equation (19) gives a faster dense **bundle contraction**;
equation (20) explains why it is not automatically a smaller
ordinary network constructor. A recursive algorithm would
need to preserve or exploit these coupled bundle functions
after each contraction. No such all-level invariant is
assumed here, and no stronger factoring exponent is claimed.

## 8. Exact checks and attribution

Run:

```sh
python3 tools/check_unsplit_monoid_contraction.py \
  --output evidence/unsplit-monoid-contraction.json
```

The checker includes genuinely unsplit coefficient fields,
complete binary basis matrices, fast periodic CRT maps,
literal polynomial block products, affine-subspace evaluation
and interpolation, all reductions, complete residue
convolutions, descended Fourier transforms, scalar CP tensors,
and both complete two-frequency middle tensors. The published
JSON records the actual finite counts and fields. These tests
supplement the ordinary uniform proof.

The source transform formulas are `hatOut`, `faOut`, and
`midOut` in
[`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).
The previous note gives the semigroup-algebra context and
primary reference. The unsplit refinement arose from the
independent audit's observation that local multiplicative
convolution need not split odd characters. The direct block
descent, complete original-field algorithms, and bundled
middle-layer bound are proved here.
