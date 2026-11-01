# An exact spectator factor at every complete root-of-unity block

## Result and scope

At the root orders three, five, and fifteen used by the fusion constructor,
every alphabet whose size is a multiple of the root order has an exact
tensor-product description. The larger crossing consists of the base
crossing and a full crossing on the quotient labels. Its constant degree
therefore does not provide free additional physical letters.

Use the fixed `F16` encoding in
[explicit-fused-row-formula.md](explicit-fused-row-formula.md). The pairs
`(theta,ell)` below are `(2,3)`, `(8,5)`, and `(4,15)`, where `ell` is
the exact multiplicative order of `theta`. This note proves the following
statement for **every integer `t>=1`**, using two exact finite polynomial
identities for each fixed pair.

**Theorem.** Let `m=t ell-1`. After complete scalar gcd removal and
the explicit joint endpoint normalization, the fused matrix satisfies

```math
\widehat F_m[(c,d),(a,b)]
=\widehat F_{\ell-1}[(c_0,d_0),(a_0,b_0)]
  [C=A][D=B],                                         \tag{1}
```

where `a=ell A+a_0`, `b=ell B+b_0`, `c=ell C+c_0`, and
`d=ell D+d_0`, with every remainder in `0,...,ell-1`. Thus the strict
crossing tensor is, after the corresponding physical relabeling,

```math
T_m=T_{\ell-1}\otimes\mathrm{Crossing}_t.                \tag{2}
```

The base endpoints are invertible and their primitive degrees are
`ell-1`, as checked directly. Hence (1) supplies strict crossings of
degree `ell-1` for all these multiples, rather than only the previously
checked finite cases.

At **every local place, in both cost cuts**, the complete Smith profile
of the large tensor is the base profile with each order repeated `t^2`
times. This is an exact consequence of (2), including after the joint
endpoint normalization. It is not a claim about arbitrary additional
physical gauges, changes of basis, or different decompositions of the
spectator crossing.

## 1. The divided row has identical invariant blocks

Work first with the unscaled degree-one row `L_m` from the explicit
formula. For `m=t ell-1`, every coefficient depends only on `m mod ell`
and the middle weight `b mod ell`. The coefficient decreasing the
middle weight vanishes at `b=j ell`; the coefficient increasing it
vanishes at `b=(j+1)ell-1`. These are the identities

```math
\theta+\theta^{2b-2m-1}=0\quad(b=j\ell),
\qquad
\theta^m+\theta^{m-2b-2}=0\quad(b=(j+1)\ell-1).
\tag{3}
```

Consequently each block of middle weights
`j ell,...,(j+1)ell-1` is invariant. On every block, the row is
literally `L_(ell-1)` on the remainder, with identical scalar and
spectral-variable conventions. The coarse middle label never changes.

This row statement alone does not determine the full fusion matrix:
the `m` input carrier bits still belong to one symmetric block. The
following generating-variable calculation accounts for that block.

## 2. Generate the complete input symmetric vector

Number carrier positions from `0` to `m-1`, with the lower position
written first. For a word of weight `a`, its inversion count is

```math
\mathrm{inv}(w)
=\sum_{i:w_i=1}(m-1-i)-\frac{a(a-1)}2.
\tag{4}
```

Introduce a formal variable `X`. The inversion-weighted symmetric
vector therefore has the exact expression

```math
S_{m,a}=\theta^{a(a-1)/2}[X^a]
 \bigotimes_{i=0}^{m-1}
 \left(|0\rangle+X\theta^{i-(m-1)}|1\rangle\right).
\tag{5}
```

If `m=t ell-1`, the local factor is `theta^(i+2)`, independent of `t`.
The prefactor in (5) also depends only on `a mod ell`. Indeed, for
`a=j ell+a_0`, the difference of its exponents is

```math
\frac{a(a-1)-a_0(a_0-1)}2
=\ell\,\frac{j(j\ell+2a_0-1)}2\in\ell\mathbb Z,
\tag{6}
```

because `ell` is odd. No division in the coefficient field occurs in
(4)--(6); these are integer exponents.

## 3. Two finite transfer-period identities

On the fine middle weights `b=0,...,ell-1`, define an `ell`-by-`ell`
matrix over `F16[X,z]` for each output carrier bit `e`:

```math
T_i^e[b+a-e,b]
=X^a\theta^{(i+2)a}
 L_{\ell-1}[b+a-e,e;a,b](z\theta^{2i}),
\qquad a\in\{0,1\}.
\tag{7}
```

Entries outside the fine middle block vanish by (3). The matrices are
periodic in `i` with period `ell`. In explicit polynomial form, their
only possible entries are

```math
\begin{aligned}
T_i^0[b,b]&=\theta^{b-1}+z\theta^{1-b+2i},\\
T_i^0[b+1,b]&=X\theta^{i+2}
                   (\theta^{-1}+\theta^{-2b-3}),\\
T_i^1[b,b]&=X(\theta^{i-b}+z\theta^{3i+b+4}),\\
T_i^1[b-1,b]&=z\theta^{2i}(\theta+\theta^{2b+1}).
\end{aligned}                                                    \tag{8}
```

Put `h(z)=1+z^ell`. The finite identities are

```math
\boxed{\quad
T_0^0T_1^0\cdots T_{\ell-1}^0=h(z)I,
\qquad
T_0^1T_1^1\cdots T_{\ell-1}^1=X^\ell h(z)I.
\quad}                                                          \tag{9}
```

The exact checker multiplies the displayed bivariate polynomial
matrices for each of the three fixed field parameters. Every
off-diagonal coefficient vanishes, and every diagonal polynomial is
exactly the one in (9). It also matches (8) against the independently
certified explicit row constructor. The record retains the full
resulting coefficient lists; no floating arithmetic or matrix-rank
estimate is involved.

Equation (9) holds starting at **any** cyclic shift of the index. To
see this, each `T_i^e` is triangular with nonzero polynomial diagonal,
and is therefore invertible over the fraction field `F16(X,z)`.
Cyclically rotating a product equal to a scalar identity conjugates
that identity and changes nothing. Periodicity then covers all integer
starting positions. The resulting equality is again an identity of
polynomial matrices.

Only the six finite identities in (9) need checking to support every
`t` in the theorem. No unproved formula for arbitrary odd root orders
is being asserted.

## 4. Strip complete periods from every carrier output

The canonical carrier output word of weight `d` is
`0^(m-d)1^d`. Its coefficient in `S_(m,d)` is one. The all-degree
fusion invariance theorem permits extracting that coefficient to
recover the full matrix entry.

Rows process the carrier positions from high to low. Thus the matrix
acting on a fine middle block, for this canonical output, is

```math
P_{m,d}(X,z)=
 \left(\prod_{i=0}^{m-d-1}T_i^0\right)
 \left(\prod_{i=m-d}^{m-1}T_i^1\right),                \tag{10}
```

where products are written in increasing index order from left to
right. This convention makes the highest-index row act first.

Write `d=j ell+r`, with `0<=r<ell`. The zero segment in (10) has
`t-j-1` complete periods followed by `ell-1-r` positions. The one
segment has `j` complete periods and `r` positions. Removing those
periods by (9), including its cyclic-shift version, leaves

```math
P_{t\ell-1,j\ell+r}(X,z)
=h(z)^{t-1}X^{j\ell}P_{\ell-1,r}(X,z).
\tag{11}
```

The remaining product is exactly the base word, with the same
spectral positions. Its degree in `X` is at most `ell-1`.
Taking the coefficient `X^a` in (5) therefore forces
`floor(a/ell)=j=floor(d/ell)`. Equations (5)--(6) give precisely
the base fine matrix entry. The coarse middle label is unchanged
by Section 1.

Let `H_m` denote the full fused matrix formed from the unscaled
degree-one rows, before its final common gcd is removed. We have proved
the exact polynomial identity

```math
H_{t\ell-1}[(c,d),(a,b)]
=h(z)^{t-1}H_{\ell-1}[(c_0,d_0),(a_0,b_0)]
  [C=B][D=A].                                         \tag{12}
```

Thus the coarse part before endpoint normalization is a swap, not an
identity. This is the source of the nonscalar raw endpoints at these
resonant parameters.

## 5. Scalar gcd removal and joint endpoint normalization

Every fine entry of the base matrix occurs in (12), and all other
entries are zero. Its complete monic gcd is consequently
`h(z)^(t-1)` times the base matrix's complete monic gcd. The original
raw fusion word additionally contains the product of the explicitly
removed row factors. Passing between monic row factors and the
unscaled rows contributes only a nonzero field scalar.

Therefore the primitive raw matrices obey

```math
F_{t\ell-1}=s_t\,F_{\ell-1}\otimes P_t,
\qquad s_t\in F16^\times,                            \tag{13}
```

under the quotient/remainder ordering; `P_t` swaps the two coarse
labels. This statement accounts for every removed polynomial factor
and every row-normalization scalar.

Let `E=F_(ell-1)(1)`. The base endpoint is invertible in each of the
three certified cases. Hence the larger raw endpoint is
`s_t E tensor P_t` and is invertible for every `t`, even when the
characteristic divides `t`. Since `P_t^2=I`, explicit right
normalization gives

```math
(s_tF_{\ell-1}\otimes P_t)(s_tE\otimes P_t)^{-1}
=(F_{\ell-1}E^{-1})\otimes I.
\tag{14}
```

This proves (1). The base primitive degree is `ell-1`; multiplication
by an invertible constant endpoint matrix preserves it. Thus every
normalized matrix in (1) has that same exact degree.

## 6. Complete local profiles and the unweighted pressure cost

In either cost cut, the spectator crossing in (2) is an invertible
constant matrix of size `t^2`. Tensoring with it repeats every local
Smith channel of the base by exactly `t^2`. This holds at finite
places, closed places over the base field, and infinity with the
unchanged homogeneous degree.

For tensor length `r`, every coefficient functional of the larger
tensor is the corresponding base coefficient tensor times a spectator
crossing on `t^r` letters. Both cut ranks multiply by `t^(2r)`.
In the exact adaptive implementation, channel counts and all jet
costs likewise multiply by `t^(2r)`; the same is true of the full
dimension cap. The greedy jet order is unchanged, and every rank-cost
product is multiplied by `t^(4r)`.

Consequently, for the same unweighted coefficient-functional method,
a base exponential rank parameter `rho_base` becomes

```math
\rho_t=
\frac{\log t+\rho_{\rm base}\log\ell}{\log t+\log\ell}.
\tag{15}
```

If `rho_base<1`, this is at least `rho_base` and tends to one with
`t`. The constant degree in this family therefore cannot improve
that unweighted method by simply adding quotient labels.

Equation (15) is deliberately scoped. It does not rule out physical
gauges coupling fine and coarse letters, correlated type restrictions,
additional crossing decompositions of the spectator, or nonsplit
alphabet sizes. Those constructions must charge their own full ranks
and physical alphabets.

## 7. Reproduction

```sh
python3 tools/check_fused_frobenius_spectators.py \
  --output evidence/fused-frobenius-spectators.json
```

The default checks certify the six period identities, their literal
row interpretation, the generating-vector coefficients, the three base
endpoints and degrees, and several complete larger matrices. For the
larger previously observed controls, add

```sh
python3 tools/check_fused_frobenius_spectators.py \
  --literal-multiple 2:4 --literal-multiple 2:5 \
  --literal-multiple 4:2 \
  --output evidence/fused-frobenius-spectators.json
```

An argument `theta:t` requests a full independent compressed construction
at `m=t ell-1`, followed by complete coefficient and raw-endpoint
comparisons. These controls test the transcription of the theorem;
the proof for arbitrary `t` is (3)--(14).

The elementary fusion and explicit-row formulas are prior repository
work cited above. The transfer-period certificate, all-multiple
factorization, and profile accounting were developed with OpenAI
assistant support. No Lean formalization is added.
