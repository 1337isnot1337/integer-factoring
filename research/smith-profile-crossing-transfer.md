# A finite Smith-profile certificate for polynomial crossing families

## Result

Two local Smith spectra suffice to certify an all-length rank bound for a
polynomial crossing family. The family may use an ordered auxiliary matrix
product. All auxiliary indices are placed in one factor of each separate
cut, so its matrices are never commuted through the opposite factor.

For the binary elementary pencil the certificate gives exactly
`rho=H_2(1/4)=0.811278...`, which is weaker than the established half-rate
factoring construction. No stronger factoring exponent is proved here.
The purpose is an exact finite search criterion for new polynomial
families, including order-dependent ones.

The proof uses elementary Smith reduction over the power-series valuation
ring and coefficient extraction. It does not assume semisimplicity of a
symmetric-group representation and is valid in every characteristic.

## 1. Polynomial matrix-product crossings

Fix an external alphabet size `q>=2`, an auxiliary dimension `D`, a finite
field `F`, and a polynomial auxiliary matrix

```math
K(z;x,y,u,v)\in\operatorname{Mat}_D(F[z])
```

of degree at most `L`. All of these data are fixed independently of the
word length `t`. Suppose vectors `e,f` satisfy `f^T e=1` and

```math
K(1;x,y,u,v)e=[x=v][y=u]e.                            \tag{1}
```

The scalar-identity endpoint `K(1)=E_q I_D` is a special case. Define

```math
P_t(z;\mathbf x,\mathbf y,\mathbf u,\mathbf v)
 =f^T K(z;x_1,y_1,u_1,v_1)\cdots
       K(z;x_t,y_t,u_t,v_t)e
 =\sum_{k=0}^{Lt}z^k A_{t,k}.                         \tag{2}
```

The product keeps its position order. Applying (1) successively gives
`P_t(1)=E_(q^t)`, so the coefficients in (2) are exact crossing pieces.

For the row cut, flatten the local tensor with row coordinates
`(x,y,a,b)` and column coordinates `(u,v)`, where `a,b` are both
auxiliary indices. This is a `q^2 D^2 by q^2` polynomial matrix. For
the column cut use rows `(x,u,a,b)` and columns `(y,v)`. Either cut
may instead put both auxiliary indices on its column side and transpose
the resulting rectangular matrix for the Smith calculation. One may
choose the auxiliary side independently for the two cuts.

The local matrices have generic rank `n=q^2`. In fact contracting their
auxiliary entries with `f_a e_b` at `z=1` produces the invertible
flattening of the crossing. A rank deficiency before this linear row
contraction would be impossible.

## 2. The two finite spectra

Choose one of the two cuts for the calculation at infinity, and use the
other cut at zero. Call them the first and second cuts in this section.

For the first cut, let `M(z)` be its rectangular local matrix. Put

```math
\widehat M(w)=w^L M(w^{-1}).
```

Let `mu_1,...,mu_n` be the Smith orders of this polynomial matrix at
`w=0`, and define its degree spectrum

```math
d_i=L-\mu_i.                                        \tag{3}
```

Negative `d_i` are permitted. Each `d_i<=L`.

For the second local cut, let `nu_1,...,nu_n` be its Smith orders at
`z=0`. They are all nonnegative. Since a nonzero maximal minor has
degree at most `nL`, every finite Smith order in these computations is
at most `nL`. Consequently exact arithmetic modulo `z^(nL+1)` suffices
to recover every order. Full generic rank was proved above.

For a real number `s>=1`, define the positive Laurent sums

```math
U(s)=\sum_{i=1}^n s^{d_i},\qquad
V(s)=\sum_{i=1}^n s^{-\nu_i},\qquad B(s)=U(s)V(s).    \tag{4}
```

**Smith-profile bound.** Every coefficient in (2) satisfies

```math
r(A_{t,k})c(A_{t,k})\le (Lt+1)^2 B(s)^t,             \tag{5}
```

and the complete rank sum satisfies

```math
\mathcal S_t\le (Lt+1)^2 B(s)^{t/2}.                 \tag{6}
```

The roles of the two cuts may be exchanged. A reciprocal change of the
polynomial parameter gives the same option, since it reverses the
coefficient list without changing the crossing endpoint at one.

## 3. Coefficient extraction and ordered auxiliary factors

Smith reduction over `F[[w]]` gives a factorization

```math
M(w^{-1})=U_\infty(w)
             \operatorname{diag}(w^{-d_1},...,w^{-d_n})
             V_\infty(w),                           \tag{7}
```

where both outside factors have power-series entries with nonnegative
orders. The rectangular factor carrying both auxiliary indices is
interpreted as an auxiliary matrix for each physical half-label and
Smith-channel index. The opposite factor has only scalar field entries.

After `t` copies, a channel word `i=(i_1,...,i_t)` contributes the scalar
power `w^(-sum_j d_(i_j))`. Its first factor is the original ordered
auxiliary matrix product, with `f^T` and `e` applied at the ends. Its
second factor is the product of scalar functions of the other physical
half-labels. The scalar factors can be collected without changing the
order of any auxiliary matrices. This gives one rank-one term per
channel word before extracting series coefficients.

Let `D(i)=sum_j d_(i_j)`. To obtain coefficient `w^(-k)`, the two
outside series must contribute total nonnegative degree `D(i)-k`.
There are exactly `D(i)-k+1` possible degree splits when `D(i)>=k`,
and no possible splits otherwise. Each split is a rank-one outer
product across the first cut, even though its auxiliary factor may
depend on every position in its fixed order. Therefore

```math
r(A_{t,k})\le (Lt-k+1)N_\infty(t,k),\qquad
N_\infty(t,k)=\#\{i\in[n]^t:D(i)\ge k\}.             \tag{8}
```

The same argument at zero for the second cut gives

```math
c(A_{t,k})\le(k+1)N_0(t,k),\qquad
N_0(t,k)=\#\{i\in[n]^t:\sum_j\nu_{i_j}\le k\}.       \tag{9}
```

For example, the coefficient split in (8) is made after the entire
ordered auxiliary product is formed. It is not a choice of independent
split degrees at every site, which would introduce an incorrect
exponential overhead. Only the two full-word factor degrees are split.

For `s>=1`, direct nonnegative counting gives

```math
N_\infty(t,k)\le s^{-k}U(s)^t,\qquad
N_0(t,k)\le s^k V(s)^t.                              \tag{10}
```

Multiplying (8)--(10) proves (5). There are at most `Lt+1` coefficient
pieces, so summing their square-root rank products proves (6).
No analytic convergence of the series is used: every needed
coefficient is a finite formal-power-series calculation.

## 4. A numerical factoring criterion with finite exact data

Put

```math
\rho=\frac{\log B(s)}{4\log q}.                      \tag{11}
```

For fixed local data and fixed rational `s>=1`, equations (5)--(6)
give a polynomial-piece crossing family of total channel width
`4rho lambda+O(log(lambda+2))`, with `lambda=t log2 q`.
For an arbitrary binary width, take the least `t` with `q^t>=2^lambda`
and restrict each leg to the first `2^lambda` words. The excess
`t log2 q-lambda` is bounded by the fixed constant `log2 q`.

Each entry of (2) is computed by multiplying `t` fixed-size polynomial
matrices and extracting one coefficient. All coefficient tables and both
exact rank factorizations cost `2^O(lambda)` field operations; over the
fixed field this is also `2^O(lambda)` bit operations. Thus the usual
effective crossing-to-factoring transfer applies whenever `rho<1`.

In particular a verified local polynomial tensor and rational `s` with

```math
B(s)<q^{2\sqrt2}                                    \tag{12}
```

would prove an exponent above `1/2`. A stricter rational exponent can
avoid irrational comparisons in a finite certificate. For instance,
`B(s)^5<q^14` gives `rho<7/10<1/sqrt(2)`.

The theorem does not assert the existence of such a local tensor.
Large Smith valuations guessed from a truncated matrix are not enough:
the truncation must meet the proven `nL+1` precision bound, or a full
factorization certificate must be supplied.

## 5. The elementary binary pencil as an exact control

Over `F_2`, let `E_2` be the one-bit swap, let `I_4` be the row
identity, and take

```math
P(z)=I_4+z(E_2+I_4).
```

Its row degree spectrum is `(1,0,0,-1)`. Its column orders at zero
are `(0,1,1,2)`. Therefore

```math
B(s)=(s+2+s^{-1})(1+2s^{-1}+s^{-2})
     =\frac{(s+1)^4}{s^3}.                           \tag{13}
```

The exact identity

```math
B(s)-\frac{256}{27}
 =\frac{(s-3)^2(27s^2+14s+3)}{27s^3}                 \tag{14}
```

shows that the best bound is attained at `s=3`, and

```math
\rho=\frac14\log_2(256/27)
     =2-\frac34\log_2 3
     =H_2(1/4).
```

This recovers the known two-pair exponent through a different, entirely
local valuation calculation. It does not improve the half-rate theorem.

## 6. A degree-budget limit of this certificate

The basic local certificate cannot improve the factoring score when the
external alphabet has `q<=4`, regardless of the auxiliary dimension or
the ordering of its matrix products. This limits the sufficient bound
(4), not the actual tensor ranks or every MPO construction.

First remove any scalar power of `z` common to every entry of `K`.
It shifts all `d_i` and all `nu_i` by the same amount and leaves
`B(s)` unchanged. Now some constant entry is nonzero. Let `L` be the
actual polynomial degree and assume `L>0`; a constant family already
has `rho=1`. Put `n=q^2` as above.

Some entry has degree `L`, so one infinity Smith order is zero and
`max_i d_i=L`. The sum of the degree spectrum is the largest degree
of a nonzero maximal minor, hence `sum_i d_i>=0`. At zero, the nonzero
constant entry gives `min_i nu_i=0`, while maximal-minor degree gives
`sum_i nu_i<=nL`.

Apply arithmetic--geometric mean to the other `n-1` terms in each
Laurent sum. For `s>=1`,

```math
U(s)\ge s^L+(n-1)s^{-L/(n-1)},\qquad
V(s)\ge 1+(n-1)s^{-nL/(n-1)}.                        \tag{15}
```

Writing `w=s^(Ln/(n-1))` therefore gives

```math
B(s)\ge\frac{(w+n-1)^2}{w^{(n+1)/n}}
 \ge\frac{4n^2}{(n+1)^{1+1/n}}.                      \tag{16}
```

The last positive function decreases until `w=n+1` and increases
afterward, as follows by differentiating its logarithm. Thus every
certificate of the form (4) has

```math
\rho\ge\beta(q):=
\frac{\log\!\left(4q^4/(q^2+1)^{1+1/q^2}\right)}
     {4\log q}.                                    \tag{17}
```

| Local alphabet `q` | Certificate floor `beta(q)` |
| ---: | ---: |
| 2 | 0.7743974703 |
| 3 | 0.7332695228 |
| 4 | 0.7071338414 |
| 5 | 0.6890022872 |
| 8 | 0.6569610575 |

The four-letter floor is just above `1/sqrt(2)`. The comparison can be
certified with small exact integers: `17^13872<2^56704` proves
`beta(4)>577/816`, and `2*577^2>816^2` proves
`577/816>1/sqrt(2)`. The checker includes exact comparisons for `q=2,3`
as well. Consequently a search using only this Smith envelope must use
at least five external letters to have any possibility of beating the
half-rate threshold.

The abstract spectra attaining equality in the degree constraints need
not come from a crossing tensor. For example, ideal profiles with
`L=n-1` would be `(n-1,-1,...,-1)` at infinity and `(0,n,...,n)` at
zero. Equations (15)--(17) do not construct those profiles, or any local
tensor realizing them. Ordered auxiliary refinements that retain more
than the two scalar spectra are outside this limitation.

## 7. Search status

The accompanying exact Smith engine currently screens two concrete
degree-two/three local families: the binary three-pairing auxiliary
matrices and rank-one nilpotent transvection loops on physical-pair
times auxiliary space. Initial screens find no certificate satisfying
(12). Those screens test the sufficient upper bound (4); failure of
that bound is not a lower bound on the actual coefficient ranks.

Independent verification checks the Smith orders against block-Toeplitz
ranks, computes full small polynomial tensors from ordered auxiliary
matrix products, and compares every small coefficient rank with
(8)--(9). The remaining constructive target is a local polynomial tensor
whose two spectra satisfy (12), or a refinement that counts only
channel words with nonzero auxiliary products.
