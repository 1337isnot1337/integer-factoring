# Streamed finite seeds and exact polynomial certificates

## Result and scope

A finite crossing certificate can be amplified with **polynomial-time
factor queries and polynomial working space**. The full factor tables
need not be constructed or stored. This strengthens the effective
interface in [finite-seed amplification](finite-seed-amplification.md).

There are two sufficient certificates. Throughout, the seed and its
field are fixed, and the external label length tends to infinity.

1. An exact decomposition on `d` labels with rank sum
   `S=sum_i sqrt(r_i c_i)` gives width parameter
   `rho=log(S)/(2 log(d))`.
2. An exact polynomial crossing of degree `D`, with the two polynomial
   factor certificates below, gives

   ```math
   \rho=\frac{\log B}{4\log d},\qquad
   B=\frac{F_a(x)F_b(x)}{x^D},\quad 0<x\le1.
   ```

The second certificate can use a **rational** `x`. Its strict inequality
can therefore be verified with integers. No numerical optimization or
asymptotic extrapolation is needed for a positive result.

Neither theorem supplies a seed that improves the current time
exponent. The worked binary example has `rho=H_2(1/4)>1/sqrt(2)` and is a
negative control. These theorems provide constructive interfaces for a
future seed, including polynomial-space factor recovery through the
[fixed-field](scale-one-modular-transfer.md) and
[bit-tagged recovery](bit-tagged-witness-recovery.md) results.

## 1. Streamed factors for a rank-sum certificate

Suppose `E_d=sum_{i=1}^k A_i` over a fixed finite field, with exact seed
factorizations

```math
A_i^{\rm row}=L_iR_i,\qquad
A_i^{\rm col}=\widetilde L_i\widetilde R_i.
```

Their inner dimensions are `r_i` and `c_i`. Remove zero pieces. The seed
matrices, their factors, and `k` are constants.

For a count vector `n=(n_1,...,n_k)`, `sum n_i=t`, let `B_n` be the sum
of all tensor products with exactly `n_i` copies of `A_i`. A row channel
is a pair consisting of a piece word `w` of that type and a local rank
index `s_j in [r_{w_j}]` at each of its `t` positions. Its two factors are

```math
\prod_{j=1}^t L_{w_j}[(x_j,y_j),s_j],
\qquad
\prod_{j=1}^t R_{w_j}[s_j,(u_j,v_j)].
```

Summing over the local indices gives the chosen tensor product.
Summing over piece words gives `B_n`. Thus this is an exact factorization
with

```math
R_n={t\choose n_1,\ldots,n_k}\prod_i r_i^{n_i}
```

channels. The other cut has the analogous factorization with
`C_n=multinom(t;n) prod_i c_i^{n_i}` channels. The multinomial theorem
gives

```math
R_nC_n\le S^{2t},\qquad S=\sum_i\sqrt{r_ic_i}.           \tag{1}
```

All words of a fixed type have the same number `prod_i r_i^{n_i}` of
local index choices. Unrank the word using multinomial counts, then
unrank its local indices in mixed radix. Each counter has `O(t)` bits.
The query uses polynomial time and space in `t`. It does not enumerate
earlier channels.

An entry of `B_n` is a coefficient of a product of `t` linear forms in
the fixed number `k` of variables. A coefficient dynamic program uses
at most `(t+1)^k` field elements. This also gives polynomial-time entry
queries.

Take the least `t` with `d^t>=2^lambda`, and restrict each external leg
to the first `2^lambda` words. Restricting the displayed factors gives
the restricted tensor without changing the channel counts. If a channel
count exceeds `2^(2lambda)`, replace that cut's factorization by an
identity matrix times the restricted piece matrix. Its entries are
available from the preceding dynamic program. This independently caps
each width at `2lambda`; it cannot increase the product in (1).

Pad each channel index set with zero factors to the next power of two.
Since `t log_2 d=lambda+O(1)`, the sum of the two binary widths is

```math
\le2t\log_2 S+2=4\rho\lambda+O(1),\qquad
\rho=\frac{\log S}{2\log d}.                          \tag{2}
```

There are `binom(t+k-1,k-1)` pieces. All piece, entry, and factor queries
use polynomial time and space in `lambda`. This proves the stronger
effective interface without global Gaussian elimination.

## 2. A polynomial certificate with a rational tilt

Let `P(z)` be a four-tensor on `d` labels, of degree at most `D`, with
`P(1)=E_d`. The degree bound may exceed the actual degree. Suppose the
following identities of polynomial matrices are supplied exactly:

```math
P_{\rm row}(z)
 =L(z)\operatorname{diag}(z^{a_1},\ldots,z^{a_s})R(z),
                                                               \tag{3}
```

```math
z^D P_{\rm col}(1/z)
 =\widetilde L(z)\operatorname{diag}(z^{b_1},\ldots,z^{b_{s'}})
    \widetilde R(z).                                           \tag{4}
```

All exponents are nonnegative integers, and all outer matrices have
polynomial entries over the fixed field. The matrices may be rectangular.
No minimality of the supplied factorizations is required. Polynomial
Smith form supplies one way to obtain them: factor each invariant factor
as a power of `z` times a polynomial with nonzero constant term, then
absorb that latter polynomial into an outer matrix. Exact verification
of (3)--(4) is sufficient even without verifying Smith minimality.

Define the polynomials with **integer** coefficients

```math
F_a(X)=\sum_{i=1}^s X^{a_i},\qquad
F_b(X)=\sum_{i=1}^{s'} X^{b_i}.
```

For any fixed rational `0<x<=1`, put

```math
B=F_a(x)F_b(x)x^{-D}.                                 \tag{5}
```

For every `t`, the coefficient tensors

```math
Q_j=[z^j]P(z)^{\otimes t},\qquad 0\le j\le Dt,
```

sum to `E_(d^t)`. They have explicit two-cut factorizations whose channel
counts satisfy

```math
R_jC_j\le(Dt+1)^2 B^t.                                \tag{6}
```

### Proof of the channel bound

For the row cut, choose a state word `i=(i_1,...,i_t)`, let
`a(i)=sum_l a_{i_l}`, and choose `0<=k<=j-a(i)`. The corresponding factors
are

```math
[z^k]\prod_{l=1}^t L(z)[(x_l,y_l),i_l],
\qquad
[z^{j-a(i)-k}]\prod_{l=1}^t R(z)[i_l,(u_l,v_l)].          \tag{7}
```

The coefficient convolution proves the factorization. It has

```math
R_j=\sum_i\max(0,j-a(i)+1)
\le(j+1)\#\{i:a(i)\le j\}.
```

Since `0<x<=1`, every eligible word has `x^{a(i)}>=x^j`. Therefore

```math
R_j\le(j+1)x^{-j}F_a(x)^t.                              \tag{8}
```

In the reversed column polynomial, the same tensor `Q_j` is coefficient
`Dt-j`. Applying the identical argument to (4) gives

```math
C_j\le(Dt-j+1)x^{-(Dt-j)}F_b(x)^t.                       \tag{9}
```

Multiplying (8) and (9) proves (6). Using the same tilt in both cuts is
what makes the bound uniform in `j`.

## 3. Uniform factor queries

For a nonnegative budget `v`, define

```math
H_0(v)=v+1,\quad H_0(v)=0\ (v<0),\qquad
H_m(v)=\sum_{i=1}^s H_{m-1}(v-a_i).                     \tag{10}
```

Then `H_t(j)=R_j`. The dynamic program has `O(t^2)` states because
`j<=Dt` and the seed is fixed. Its counts have `O(t+log(t+1))` bits.
At each position, the intervals of lengths `H_{m-1}(v-a_i)` identify
the next state from a channel index. After the last state, the remaining
index is `k`. This unranking takes polynomial time and space.

Each coefficient query in (7) multiplies `t` polynomials of bounded
degree. Truncation at the requested degree yields a polynomial-time
dynamic program over the fixed field. Entries of `Q_j` are computed in
the same way directly from the seed tensor. Repeat this construction for
the column certificate.

As in Section 1, restrict to a binary prefix, cap each channel count by
an identity factorization when necessary, and pad with zeros. The number
of pieces is `Dt+1=O(lambda)`. Equation (6) yields

```math
\text{width}_{\rm row}(Q_j)+\text{width}_{\rm col}(Q_j)
\le t\log_2B+O(\log(t+1))
=4\rho\lambda+O(\log(\lambda+1)),
\quad\rho=\frac{\log B}{4\log d}.                       \tag{11}
```

All factor and entry queries use polynomial time and space. The
sublinear width slack is allowed by the existing transfer. If `B<d^4`,
the family has contraction `rho<1`. The stronger condition

```math
B<d^{2\sqrt2}
```

would give a factoring exponent strictly above one half. A sufficient
fully rational certificate is `B^q<d^p`, with positive integers `p,q`
such that `p/(4q)<1/sqrt(2)`. For example,

```math
B^5<d^{14}
```

certifies `rho<7/10<1/sqrt(2)`. The certificate identities, factor
queries, and rational inequality are all finite exact checks.

## 4. Binary example and a negative control

Over `F_2`, set

```math
E(x,y,u,v)=[x=v][y=u],\quad
R(x,y,u,v)=[x=y][u=v],\quad
P(z)=zE+(1+z)R.
```

Use lexicographic pairs in each matrix flattening. In the row cut,
`R=vv^T` with `v=(1,0,0,1)^T`, `R^2=0`, and `RE=R`. In a constant
basis with columns `(v,e_0,e_1,e_2)`, the normalized row matrix `P(z)E`
is the direct sum of

```math
\begin{pmatrix}z&1+z\\0&z\end{pmatrix},\quad z,\quad z.
```

The two-by-two block has the exact factorization

```math
\begin{pmatrix}1+z&1\\z&1\end{pmatrix}
\begin{pmatrix}1&0\\0&z^2\end{pmatrix}
\begin{pmatrix}z&1\\1&0\end{pmatrix}
=\begin{pmatrix}z&1+z\\0&z\end{pmatrix}.                \tag{12}
```

The reversed column matrix is `zI+(E+I)`. Here `E+I` is rank one and
square zero. A constant basis with columns `(e_1+e_2,e_1,e_0,e_3)`
reduces it to a two-by-two block with upper-right entry one. Replace the
left factor in (12) by `[[1,0],[z,1]]` to obtain that block. Both
valuation lists are therefore `(0,2,1,1)`, so

```math
F_a(x)=F_b(x)=(1+x)^2,
\qquad B=(1+x)^4/x.
```

At the rational minimizer `x=1/3`, this is `B=256/27`. The resulting
parameter is

```math
\rho=\frac{\log_2(256/27)}4=H_2(1/4)=0.811278\ldots.
```

It does **not** improve the half-rate construction. The checker verifies
the polynomial identities, complete tensor coefficient factorizations
through three copies, the budget unranking, and the rational count
bounds. It separately checks type-grouped seed factors, including a
binary prefix of a nonbinary label space.

## 5. Four labels still cannot pass the usual Smith estimate

The determinant-degree argument in
[finite-seed amplification, Section 4](finite-seed-amplification.md)
gives the following lower bound on the *usual cumulative Smith
envelope*, with `N=d^2`:

```math
b(d)=\frac{2\log(2N)-\frac{N+1}{N}\log(N+1)}{4\log d}.
```

It applies when a common power of `z` is first removed and the actual
degree is used. Those normalizations cannot worsen the envelope. The
previous note excludes dimensions two and three. Dimension four is
also excluded, despite the small numerical gap:

```math
b(4)>\frac{577}{816}>\frac1{\sqrt2}.
```

Both strict inequalities have short integer certificates:

```math
2^{3544}>17^{867},\qquad 2\cdot577^2>816^2.
```

For the first, multiply the defining logarithmic inequality for `b(4)`
by `64 log 4`, then clear the denominator `51`. It becomes
`32^1632>17^867 4^2308`, which is exactly the displayed power inequality.
The second comparison is immediate by squaring positive numbers.

Thus this Smith estimate needs at least five labels to have any chance
of certifying a half-exponent improvement. This is a limitation of that
upper-bound method, not a lower bound on the coefficient ranks. It does
not limit a more accurate cancellation analysis or a different label
restriction. The rank-sum certificate has a separate
[six-label obstruction](six-label-rank-sum-obstruction.md).

## Reproduction

```sh
python tools/check_streamed_seed_transfer.py \
  --output evidence/streamed-seed-transfer.json
```
