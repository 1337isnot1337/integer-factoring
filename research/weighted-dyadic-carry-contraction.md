# Weighted dyadic carries and one joint Frobenius quotient

## Result and scope

The whole carry layer of the exact multiplication network can be
replaced, on its relevant inputs, by a product of dyadic Fourier
weights. Nonzero carry weights make every such factor invertible.
The sparse recovery theorem permits these weights even when their
values and the digit radix change with the recovery field.

After summing one operand, the large digit boxes and the high parts
of all carries combine into **one common Frobenius quotient**. The
remaining low carry factor has at most
`exp(O(D log D))` monomials. Under the adaptive dimensions this is
`exp(O(n^(9/20) log n))`, below the cost of the sparse recovery family.

This is an exact reformulation of the coupled factoring contraction.
It removes the clipping conditions from the relevant Fourier carry
tensor, permits invertible factors at every frequency, and specifies
a more structured remaining coefficient problem. It does not supply
a faster method for that coefficient problem. The established score
remains `1/2`.

The result uses the ordinary proofs in
[adaptive-transform-constructor.md](adaptive-transform-constructor.md),
[operand-length-transform-constructor.md](operand-length-transform-constructor.md),
and [frobenius-compatible-sparse-recovery.md](frobenius-compatible-sparse-recovery.md).
It is a different extension of the carry function from the one in
[carry-fourier-elimination.md](carry-fourier-elimination.md). The two
extensions agree on every convolution produced by the input boxes;
their values elsewhere need not agree.

## 1. A dyadic carry domain with no new factor assignments

Let `D` and `b` be powers of two, with `D>=2` and `b>=2`. Write

```math
b=2^w,\qquad D=2^t,\qquad
K=D(b-1),\qquad R=Db=2^{w+t},
\qquad M=2^m+1>Db^2.                            \tag{1}
```

Let `T` be the number of coefficient positions, and pad the target
integer with radix-`b` digits `n_0,...,n_(T-1)`. Operand digit lengths
`d_p,d_q` satisfy

```math
d_p+d_q\le T,\qquad \min(d_p,d_q)\le D.
```

Consequently, their convolution coefficients obey

```math
a_j=\sum_{i+k=j}p_iq_k,\qquad
0\le a_j\le D(b-1)^2,\qquad a_{T-1}=0.          \tag{2}
```

Replace each internal carry domain `[0,K]` by the full dyadic interval
`[0,R-1]`. Set the two endpoint carries to zero. Impose

```math
a_j+r_j\equiv n_j+b r_{j+1}\pmod M,
\qquad r_0=r_T=0.                              \tag{3}
```

On the domain (2), these congruences are exact integer equations.
Indeed, every permitted carry satisfies

```math
0\le a_j+r_j
\le D(b-1)^2+Db-1
=Db^2-Db+D-1<Db^2<M,
```

and

```math
0\le n_j+b r_{j+1}\le b-1+b(Db-1)=Db^2-1<M.
\tag{4}
```

Multiplying the integer equations by `b^j` and summing proves `pq=N`.
Conversely, the genuine multiplication carries satisfy `r_j<=K`:
if `r_j<=K`, then
`a_j+r_j<=D(b-1)^2+D(b-1)=bK`, so the next carry is at most `K`.
The equations determine each successive carry uniquely from `r_0=0`.

Thus padding to `[0,R-1]` introduces no new accepted assignment and
does not change the unique factor-assignment correspondence.
The carry bit width stays `w+t`, exactly the width already used by
the adaptive constructors. The existing table-size bounds used that
full bit width, so they remain valid.

In a directed local table, one may retain the ordinary quotient
`floor((a+r)/b)`, clip an out-of-domain quotient to a fixed encoded
value, and reject it by a separate local weight. The proof above is
about all globally accepted assignments, where the convolution bound
(2) holds. It does not treat arbitrary modular residues as bounded
integer coefficients.

## 2. Nonzero carry weights preserve deterministic recovery

Choose a permitted odd-prime recovery degree `ell`, round `w` to a
multiple of `ell` as in the sparse recovery note, and construct

```math
E=\mathbb F_{2^\ell},\qquad
F=\mathbb F_{2^{2m}},\qquad
\mathcal A=F\otimes_{\mathbb F_2}E.
```

Let `zeta in F` have order `M`. Fix **one** element
`gamma in E\{0,1}` before scanning the evaluation points
`alpha,beta in E\{0,1}`. Give internal carry `r_j` the weight
`gamma^(r_j)` exactly once. For a divisor `d`, put

```math
\rho_\ell(d)=\sum_{j=1}^{T-1}r_j(d,N/d),
\qquad
c_{\ell,d}=\gamma^{\rho_\ell(d)}\ne0.           \tag{5}
```

The weighted prefix value is the evaluation of

```math
\mathcal D^{\gamma}_{N,k,x}(X,Y)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
 c_{\ell,d}X^dY^{N/d}.                         \tag{6}
```

Its exponent support is exactly the divisor support. The coefficients
are units in `E` and hence in `A`. They are constant while the whole
punctured grid for this degree is scanned. They may vary with `ell`,
the chosen radix, and the fixed choice of `gamma`.

The support-universal corollary in Section 1 of the sparse recovery
note applies. A completely vanishing grid still forces a collision
of an anchored exponent pair modulo `2^ell-1`, independent of all
these coefficient choices. The same prime-degree product therefore
detects every nonempty divisor prefix.

The divisor-one contribution remains `X Y^N`: multiplying one by
`N` produces no carries, so `rho_ell(1)=0`. Removing that monomial
and using the existing prefix search gives the smallest prime factor
and then complete factorization. The cost bounds and admissible
degree family from the sparse recovery note are unchanged.

It would be invalid to vary `gamma` with `alpha` or `beta` during this
argument without a further proof. The coefficients in (6) must stay
fixed within each grid.

### A smaller query family with aligned carry weights

There is a stronger specialization that avoids that issue entirely.
Fix `beta in E\{0,1}` once per degree and scan only
`alpha in E\{0,1}`. Choose

```math
a_\ell=2^{(-t\bmod\ell)},\qquad
\gamma=\beta^{a_\ell},\qquad \gamma^D=\beta.
\tag{6a}
```

The inverse-Frobenius exponent in (6a) is a power of two, so `gamma`
is nontrivial. The univariate query is the evaluation of

```math
F_{N,k,x,\ell}(X)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
 \left[\beta^{N/d}\gamma^{\rho_\ell(d)}\right]X^d.
\tag{6b}
```

Every coefficient in brackets is nonzero, and the distinct first
exponents `d` are fixed across all degrees. The univariate form of
the support-universal recovery theorem therefore applies directly.
Its common anchored product is
`product_(d!=d_0)|d-d_0|<=2^(n(K-1))`.

The query family now has only `2^ell-2` field points per degree,
instead of its square. Both operand boxes and all carry boxes still
have unit Fourier weights. The divisor-one term to remove is
`beta^N X`. No reconstruction of the coefficients in brackets is
required: only their proved nonzeroness matters.

Alignment (6a) also makes every large-box binomial in the joint
quotient below have exactly the same scalar `beta`.

## 3. The whole carry Fourier transform is a product

Extend the weighted carry function to all `a in (Z/MZ)^T` by

```math
W^\gamma_N(a)=
\sum_{0\le r_1,\ldots,r_{T-1}<R}
 \gamma^{\sum_jr_j}
 \prod_{j=0}^{T-1}
 [a_j+r_j=n_j+b r_{j+1}\text{ in }\mathbb Z/M\mathbb Z],
\quad r_0=r_T=0.                               \tag{7}
```

Section 1 proves that this extension gives the required weight on
every valid convolution. Define `S_L(Y)=sum_(0<=v<L)Y^v` and use the
Fourier convention

```math
\widehat W^\gamma_N(s)=\sum_a W^\gamma_N(a)
 \zeta^{\sum_j s_ja_j}.
```

Each choice of carries fixes one modular coefficient vector. Therefore

```math
\widehat W^\gamma_N(s)
=\zeta^{\sum_js_jn_j}
 \prod_{j=1}^{T-1}
 S_R\!\left(\gamma\zeta^{b s_{j-1}-s_j}\right).
\tag{8}
```

There are no nested inequalities or endpoint clipping terms in (8).
They were removed by the legitimate extension (7), not discarded
from an equality on arbitrary integer coefficient vectors.

Every factor in (8) is a unit. For any frequency `c`, set
`Y=gamma zeta^c`. Since `ell` is odd,
`gcd(M,2^ell-1)=1`, so

```math
(1+Y)^{-1}=
 \frac{S_M(Y)}{1+\gamma^M}.
```

Because `R` is a power of two,

```math
S_R(Y)=\frac{1+Y^R}{1+Y}.                       \tag{9}
```

The numerator is also a unit: `gamma^R` is a nontrivial element of
`E`. Thus the **entire Fourier transform** (8) is a unit at every
frequency. In particular it has an explicitly evaluable pointwise
inverse.

Constructing `F`, `E`, and `zeta` costs `2^O(m+ell)` bit operations,
as already charged in the preceding notes. Evaluating (8), or its
inverse, at a specified frequency vector costs `O(Tm)` algebra
operations using binary powering and geometric sums. Each operation
has polynomial bit cost in `m ell`. This statement does not enumerate
all `M^T` frequency vectors for free.

## 4. An exact coupled Fourier formula with no zero box factors

First suppose the second operand has full digit boxes `0<=q_j<b`.
The first operand may have the dyadic digit restrictions used for
prefix recovery. Write `P` for its permitted digit vectors. Since
`ell|w`, both `alpha^b=alpha` and `beta^b=beta`.
Write `Z_N^gamma(alpha,beta)` for this operand-box instance's
contribution to (6). In the operand-length constructor, summing its
at most `n` disjoint bit-length instances gives (6).

For a fixed first operand and frequency vector put

```math
L_j(p,s)=\sum_i p_i s_{i+j}\pmod M.
```

The entire weighted contraction is

```math
Z_N^\gamma(\alpha,\beta)
=\frac1{M^T}\sum_{s\in(\mathbb Z/M\mathbb Z)^T}
 \widehat W^\gamma_N(s)
 \sum_{p\in P}\alpha^{\sum_i p_i}
 \prod_{j=0}^{d_q-1}
 S_b\!\left(\beta\zeta^{-L_j(p,s)}\right).
\tag{10}
```

The usual orthogonality identity is valid because `M` is odd. Its
scalar image is one in characteristic two, so the normalization is
a unit. The same formula follows by expanding all boxes and then
using (7).

Every individual factor in every summand of (10) is a unit: the
carry factors by Section 3, the second-operand boxes by the sparse
recovery note, and the remaining monomial and character phases by
construction. Zero outputs arise from cancellation between summands.
This is fully compatible with empty divisor prefixes.

A partial second-operand digit interval `[a,a+2^r-1]` replaces its
factor by `Y^a S_(2^r)(Y)`, which is still a unit. For the stronger
common-Frobenius identity below, the second operand can in fact use
full boxes in both adaptive constructors.

For the operand-length constructor, condition the first operand on
bit length `s`. Its second-operand digit count is
`ceil((n-s+1)/w)`. Filling the partial top second-operand digit keeps
that digit count unchanged. It introduces no accepted product:
`p>=2^(s-1)` and `pq=N<2^n` force `q<2^(n-s+1)` already.
The degree and coefficient bounds in (2) still hold. The first
adaptive constructor already permits the full digit boxes.

## 5. Combine the large boxes into one Frobenius quotient

Work in the finite cyclic group algebra

```math
\mathcal R=
 \mathcal A[Z_0^{\pm1},\ldots,Z_{T-1}^{\pm1}]
 /(Z_0^M-1,\ldots,Z_{T-1}^M-1).
```

Let `e_j` be its coordinate exponent vectors and set

```math
v_j=b e_{j-1}-e_j\quad(1\le j<T),\qquad
u_j(p)=\sum_i p_i e_{i+j}.
```

For a fixed first operand, define

```math
G_p(Z)=\prod_{j=0}^{d_q-1}(1+\beta Z^{-u_j(p)}),
\qquad
H_D(Z)=\prod_{j=1}^{T-1}(1+\gamma^D Z^{D v_j}),
\tag{11}
```

and the low carry polynomial

```math
B_\gamma(Z)=\prod_{j=1}^{T-1}S_D(\gamma Z^{v_j}).
\tag{12}
```

Every binomial in (11) is a unit in the group algebra. For example,
if `Y=beta Z^u`, then `Y^M=beta^M` and
`(1+Y)S_M(Y)=1+beta^M`, a nonzero scalar in `E`.
The same proof applies to `gamma^D`.

Let `sigma` be the Frobenius `f->f^b` of `R`. It fixes `E`, and on
polynomials with coefficients in `E` it substitutes `Z_j^b` for
`Z_j`. Split every carry uniquely as

```math
r_j=e_j'+D z_j,\qquad 0\le e_j'<D,\quad0\le z_j<b.
```

The identity `S_(Db)(Y)=S_D(Y)S_b(Y^D)` and (9) now give the **joint**
identity

```math
\left[\prod_{j=1}^{T-1}S_{Db}(\gamma Z^{v_j})\right]
\left[\prod_{j=0}^{d_q-1}S_b(\beta Z^{-u_j(p)})\right]
=B_\gamma(Z)\,
 \frac{\sigma\!\left(H_D(Z)G_p(Z)\right)}
      {H_D(Z)G_p(Z)}.
\tag{13}
```

This combines all the large carry and second-operand boxes, retaining
their joint frequency dependence. It does not rely on independent
copies telescoping or on averaging separate Frobenius orbits.

For a group-algebra element `f`, write `[Z^a]f` for the coefficient
of the residue vector `a`. Fourier orthogonality in (10) yields

```math
Z_N^\gamma(\alpha,\beta)
=[Z^{-n}]\left[
 B_\gamma(Z)
 \sum_{p\in P}\alpha^{\sum_i p_i}
 \frac{\sigma(H_DG_p)}{H_DG_p}
 \right].                                      \tag{14}
```

The sign matches the convention in (8): the input box exponents are
`-u_j(p)`, and the external character is `zeta^(s dot n)`.

Polynomial (12) has a direct expansion with at most `D^(T-1)`
monomials. In the adaptive network `T=2D` and `D=O(P^18)` with
`P^40<=n`, so

```math
\log_2(D^{T-1})=O(D\log D)
=O(n^{9/20}\log n).                            \tag{15}
```

This low part can be expanded and charged explicitly without changing
any target cost `exp(O(n/(log n)^alpha))` for fixed `alpha<1`.
The remaining task is to obtain the specified cyclic coefficients of
the first-operand sum of the common quotients in (14).

### One common scalar denominator

Use the fixed-`beta` alignment (6a). Let `nu=T-1+d_q`, and list the
`nu` exponent columns consisting of the `D v_j` and the `-u_j(p)`.
Then

```math
H_DG_p=\prod_{h=1}^{\nu}(1+\beta Z^{c_h(p)}),
```

and the binomial inverse formula gives

```math
\frac{\sigma(H_DG_p)}{H_DG_p}
=\frac{
 \sigma(H_DG_p)\prod_{h=1}^{\nu}S_M(\beta Z^{c_h(p)})
 }{(1+\beta^M)^{\nu}}.                         \tag{16}
```

The denominator is a nonzero scalar. It is independent of the first
operand, frequency, and queried prefix within this operand-box
instance. Thus it can be cleared explicitly.

For the bit-length decomposition, one can even multiply each instance
by its own scalar `(1+beta^M)^nu` before adding the values. Each
divisor belongs to exactly one bit-length instance, so (6b) merely
acquires a known nonzero multiplier on each coefficient. The same
univariate recovery proof still applies. This permits a polynomial
numerator formulation of the remaining cyclic coefficient problem,
without a formal division gate or a shared denominator assumption
across different instances.

After this rescaling, the known divisor-one term to remove is
`(1+beta^M)^(nu for the bit-length-one instance) beta^N X`.
The earlier unscaled subtraction is used only for the unscaled query.

The concise rational expression does not give those coefficients at
unit cost. A single binomial inverse may have `M` monomials, and the
ambient group algebra has dimension `M^T`. The original first-operand
sum is also still present. Any exponent improvement must give a
uniform algorithm for these coupled operations and charge its bit
cost; none is asserted here.

## 6. Reproducible exact checks

Run:

```sh
python3 tools/check_weighted_dyadic_carry.py \
  --output evidence/weighted-dyadic-carry.json
```

The checker compares padded modular assignments with exact products,
enumerates literal weighted carry transforms, verifies full coupled
Fourier sums, and compares the common-Frobenius formula as an exact
cyclic group-algebra identity. It checks coefficients of nonempty and
empty factor sets, genuine nonzero carries, endpoint conventions,
and the everywhere-unit box assertions. Small exhaustive models and
parameters with `ell|w` are identified separately in the evidence.
The aligned tests also recover factors from univariate queries with
fixed `beta`, check the inverse-Frobenius choice of `gamma`, and repeat
recovery after the separate bit-length-instance scalars are cleared.
