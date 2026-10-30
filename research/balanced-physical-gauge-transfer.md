# Balanced label types and independent physical degree gauges

## Result and scope

The [Smith-profile construction](smith-profile-crossing-transfer.md) can
use independent Laurent monomial weights on the four physical labels in
its two cuts. Restrict every amplified label to the same balanced type.
Each whole-word weight then becomes a known scalar power, which shifts
the coefficient index without changing the tensor.

The balanced label set has size `q^t/poly(t)` for fixed local alphabet
size `q`. This polynomial loss does not change the leading factoring
parameter. The resulting finite criterion strictly enlarges the allowed
certificate data, although no improving example is supplied here.

This is an ordinary conditional construction theorem. It does not change
the proved exponent `1/2`. In particular, optimizing these weights on an
unrestricted label set and omitting their mean shifts is invalid.

## 1. Polynomial crossings and balanced words

Use the same fixed field, alphabet, polynomial tensor, and boundaries as
in the [auxiliary-gauge theorem](auxiliary-gauge-smith-transfer.md):

```math
K(z;x,y,u,v)\in\operatorname{Mat}_D(F[z]),\quad
\deg K\le L,\quad f^Te=1,\quad K(1;x,y,u,v)e=[x=v][y=u]e.
```

Write the scalar ordered product as

```math
P_t(z)=f^T\prod_{j=1}^t K(z;x_j,y_j,u_j,v_j)e
      =\sum_{k=0}^{Lt}z^k A_{t,k}.                         \tag{1}
```

Take `t` divisible by `q`. Let `T_t` consist of the words of length `t`
in which each letter occurs exactly `t/q` times. Its cardinality is

```math
D_t=\frac{t!}{((t/q)!)^q},\qquad
\frac{q^t}{(t+1)^q}\le D_t\le q^t.                        \tag{2}
```

For the lower bound, there are at most `(t+1)^q` letter-count types,
and their multinomial coefficients sum to `q^t`. The balanced coefficient
is a maximum: moving one unit from a count at least two larger than
another increases the multinomial coefficient. Since `q|t`, a maximizing
type can be balanced.

Restrict each of the four word labels in (1) to `T_t`. Denote the
restricted coefficients by `A^bal_(t,k)`. Their sum is the crossing on
`D_t` labels. The same set must be used on each of its paired legs.

## 2. Physical weights and their exact mean

For each of the two cuts choose four integer functions on the alphabet,
independently of the functions chosen for the other cut. At infinity
write them as `a_infinity,b_infinity,c_infinity,d_infinity`, and define

```math
h_\infty(x,y,u,v)=a_\infty(x)+b_\infty(y)
                       +c_\infty(u)+d_\infty(v),
\qquad
H_\infty=\frac1q\sum_{j=0}^{q-1}
 (a_\infty(j)+b_\infty(j)+c_\infty(j)+d_\infty(j)).     \tag{3}
```

Define `h_0,H_0` in the same way for the other cut. The means may be
negative or fractional; `tH_infinity` and `tH_zero` are integers.
For four balanced words,

```math
\sum_{j=1}^t h_\infty(x_j,y_j,u_j,v_j)=tH_\infty,
\qquad
\sum_{j=1}^t h_0(x_j,y_j,u_j,v_j)=tH_0.                \tag{4}
```

Choose optional rational auxiliary gauges `G_infinity(w)` and `G_zero(z)`
as in the auxiliary-gauge theorem. Form the two local rational tensors

```math
\begin{aligned}
\widetilde K_\infty(w;x,y,u,v)
 &=w^{h_\infty(x,y,u,v)}G_\infty(w)K(w^{-1};x,y,u,v)G_\infty(w)^{-1},\\
\widetilde K_0(z;x,y,u,v)
 &=z^{h_0(x,y,u,v)}G_0(z)K(z;x,y,u,v)G_0(z)^{-1}.
\end{aligned}                                                       \tag{5}
```

In either cut, a sum of four single-leg functions is a sum of a function
on its row half and a function on its column half. Its monomial weights
therefore multiply the flattened matrix by invertible diagonal matrices
over the rational-function field. Auxiliary conjugation is invertible on
the side carrying both auxiliary indices. Both local matrices still
have generic column rank `n=q^2`.

Let their respective Smith orders at the origin be `-d_i` and `nu_i`,
with negative orders permitted. Every normalization of Laurent poles
must be restored in the recorded spectra.

For a fixed real `s>=1`, set

```math
B_{\rm bal}(s)
 =s^{H_\infty+H_0}
       \left(\sum_{i=1}^n s^{d_i}\right)
       \left(\sum_{i=1}^n s^{-\nu_i}\right).            \tag{6}
```

**Balanced physical-gauge theorem.** For some fixed constant `C`, every
admissible `t` and every coefficient satisfy

```math
r(A^{\rm bal}_{t,k})c(A^{\rm bal}_{t,k})
 \le C(t+1)^2 B_{\rm bal}(s)^t,
\qquad
\sum_k\sqrt{r(A^{\rm bal}_{t,k})c(A^{\rm bal}_{t,k})}
 \le C(t+1)^2 B_{\rm bal}(s)^{t/2}.                    \tag{7}
```

The two physical weight choices need not agree. They are tools for
bounding the same restricted coefficients, not two different families
being added together.

## 3. The coefficient shifts in the two cuts

First take identity auxiliary gauges. On balanced inputs, the scalar
product of the first local tensor in (5) is exactly
`w^(tH_infinity) P_t(w^(-1))`. Hence the original coefficient `k`
is the coefficient of `w^(tH_infinity-k)` in that product.

A first-cut channel word contributes `w^(-D(i))` times two regular
whole-word series, where `D(i)=sum_j d_(i_j)`. It can contribute only if

```math
D(i)\ge k-tH_\infty.                                  \tag{8}
```

The number of full-word coefficient splits is
`D(i)+tH_infinity-k+1`, bounded by `C_infinity(t+1)` for a fixed
constant. Scalar physical factors and auxiliary matrices have not been
reordered improperly: both auxiliary indices remain in one half.

At zero the product is `z^(tH_zero)P_t(z)`. For a channel word with
`N(i)=sum_j nu_(i_j)`, the required coefficient is `k+tH_zero`.
It can contribute only if

```math
N(i)\le k+tH_0,                                        \tag{9}
```

and there are `k+tH_zero-N(i)+1<=C_zero(t+1)` possible degree splits.
The original range `0<=k<=Lt` supplies this linear bound even when
individual weighted local entries have poles.

Weighted counting of (8)--(9) gives

```math
\begin{aligned}
r(A^{\rm bal}_{t,k})&\le C_\infty(t+1)
       s^{tH_\infty-k}\left(\sum_i s^{d_i}\right)^t,\\
c(A^{\rm bal}_{t,k})&\le C_0(t+1)
       s^{tH_0+k}\left(\sum_i s^{-\nu_i}\right)^t.
\end{aligned}                                                       \tag{10}
```

Their product proves the first assertion of (7); summing over `Lt+1`
coefficients proves the second after changing the constant.

With auxiliary gauges, regularize the two boundary vectors exactly as
in the auxiliary-gauge proof. Their additional valuation shifts are
fixed integers independent of `t`, hence change only `C` in (7).
Equation (4) and the factors `s^(tH_infinity)` and `s^(tH_zero)` remain
unchanged. Thus (6), including both mean terms, is still the exponential
base.

## 4. Effective binary-width transfer

For a desired binary width `lambda`, choose the least multiple `t` of
`q` with `D_t>=2^lambda`. Equation (2) implies

```math
t\log_2 q=\lambda+O(\log(\lambda+2)).                 \tag{11}
```

Balanced words can be enumerated, ranked, and unranked using exact
multinomial coefficients with `O(t)` bits. These operations are
polynomial in `t` for fixed `q`. Restrict to the first `2^lambda`
balanced words. Complete coefficient tables and their two field-rank
factorizations take `2^O(lambda)` bit operations, as in the original
Smith transfer. All field and certificate parameters are fixed.

Consequently (7) gives total channel width

```math
4\rho\lambda+O(\log(\lambda+2)),\qquad
\rho=\frac{\log B_{\rm bal}(s)}{4\log q}.              \tag{12}
```

The [fixed-field transfer](scale-one-modular-transfer.md) and the
[fixed-characteristic recovery](fixed-characteristic-factor-recovery.md)
then apply. A verified `B_bal(s)<q^(2sqrt(2))` would improve the factoring
exponent. No such certificate is claimed here.

The fractional means in (6) cause no verification problem. Choose a
fixed rational `r>=1` and put `s=r^q`. Then

```math
B_{\rm bal}(r^q)=
 \left(\sum_i r^{q d_i+qH_\infty}\right)
 \left(\sum_i r^{-q\nu_i+qH_0}\right),                 \tag{13}
```

which is a positive rational number with integer exponents. Any strict
real improvement has a nearby rational choice of `r`, by continuity.
Alternatively, ramify the polynomial parameter to clear rational height
denominators. The ramification degree is fixed certificate data.

## 5. Scope controls

A gauge can have zero mean on balanced words while taking very different
values on arbitrary words. The restriction in Section 1 is essential.
There is no claim that all unrestricted coefficients preserve their ranks.

The mean corrections are also essential. Multiplying the zero-cut tensor
by `z` raises every zero-cut Smith order by one. Omitting the resulting
`H_zero=1` would falsely divide the bound by `s`, allowing an arbitrarily
small apparent base without changing the restricted crossing.

A weight depending on a pair of physical labels, such as `h(x,y)`, is
outside this theorem unless its whole-word sum is proved constant on the
chosen common external label set. Balanced counts of individual letters
do not fix joint pair counts. An arbitrary pair-weight transformation
can also fail to separate in the other cut.

## Verification and provenance

`check_balanced_physical_gauge.py` checks exact local Laurent spectra
against block-Toeplitz ranks, restricts explicit ordered crossings to
balanced words, compares transformed products with the correct shifted
coefficients, and verifies the admitted-word rank bounds with rational
arithmetic. Negative controls detect omitted mean shifts and the use of
unbalanced words. These are finite transcription checks; the uniform
statement is proved above.

The construction combines the local Smith and rational auxiliary-gauge
arguments with an elementary multinomial type restriction. It was
developed with substantial OpenAI assistant support. No Lean proof or
new unconditional factoring exponent accompanies it.
