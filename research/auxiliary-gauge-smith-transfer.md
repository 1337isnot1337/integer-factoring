# Auxiliary rational gauges for polynomial crossing certificates

## Result and scope

The [local Smith-profile certificate](smith-profile-crossing-transfer.md)
can use a different rational change of auxiliary basis in each cut.
The changes telescope along an ordered matrix product. They alter its
two local Smith spectra, while their boundary costs remain independent
of the word length.

This gives a larger finite search interface. It is an ordinary conditional
transfer theorem, not a new crossing construction or a factoring exponent
above `1/2`. The checker below verifies the transformations and rank
bounds on exact finite examples.

## 1. The unchanged crossing family

Fix a finite field `F`, external alphabet size `q>=2`, auxiliary dimension
`D`, and a polynomial matrix tensor of degree at most `L`:

```math
K(z;x,y,u,v)\in\operatorname{Mat}_D(F[z]).
```

Suppose fixed column vectors `e,f` satisfy

```math
f^T e=1,\qquad K(1;x,y,u,v)e=[x=v][y=u]e.
```

Its ordered scalar product is

```math
P_t(z)=f^T K(z;x_1,y_1,u_1,v_1)\cdots
             K(z;x_t,y_t,u_t,v_t)e
       =\sum_{k=0}^{Lt}z^k A_{t,k}.                 \tag{1}
```

The coefficients sum to the crossing on `q^t` labels. All coefficient
tables in the eventual algorithm are still computed from this polynomial
formula. The rational matrices introduced next are proof and certificate
data; they do not change (1).

## 2. Two independently chosen rational gauges

Choose invertible auxiliary matrices

```math
G_\infty(w)\in\operatorname{GL}_D(F(w)),\qquad
G_0(z)\in\operatorname{GL}_D(F(z)).
```

Use the first in one cut and the second in the other. Define

```math
K_\infty(w)=G_\infty(w)K(w^{-1})G_\infty(w)^{-1},\qquad
K_0(z)=G_0(z)K(z)G_0(z)^{-1}.                       \tag{2}
```

In each chosen cut, put both auxiliary indices on the same physical
half. The opposite factor then has scalar entries. The two choices of
physical half need not agree. Flattening gives rectangular rational
matrices with `n=q^2` columns and full generic column rank. Indeed the
original local matrix has full column rank at one after contraction by
`f_a e_b`, and auxiliary conjugation is an invertible operation on its
rows over the rational-function field.

At `w=0`, let the Smith orders of the first matrix be `-d_1,...,-d_n`.
At `z=0`, let those of the second be `nu_1,...,nu_n`. These are integers;
negative orders are allowed. A common scalar power clears negative
orders before ordinary power-series Smith reduction. Restore that power
after computing the spectrum.

Choose nonnegative integers `a_infinity,b_infinity,a_zero,b_zero` so that
all four boundary vectors below are regular at the corresponding origin:

```math
\begin{array}{ll}
w^{a_\infty} f^T G_\infty(w)^{-1},&
w^{b_\infty} G_\infty(w)e,\\
z^{a_0} f^T G_0(z)^{-1},&
z^{b_0} G_0(z)e.
\end{array}
```

Such integers exist because these vectors are rational. Set
`H_infinity=a_infinity+b_infinity` and `H_zero=a_zero+b_zero`.
They depend only on the fixed certificate.

For a fixed real `s>=1`, put

```math
U(s)=\sum_{i=1}^n s^{d_i},\qquad
V(s)=\sum_{i=1}^n s^{-\nu_i},\qquad B(s)=U(s)V(s).   \tag{3}
```

**Gauge transfer theorem.** There is a fixed positive constant `C` such
that, at every word length and coefficient,

```math
r(A_{t,k})c(A_{t,k})\le C(t+1)^2 B(s)^t,
\qquad
\sum_k\sqrt{r(A_{t,k})c(A_{t,k})}
 \le C(t+1)^2 B(s)^{t/2}.                          \tag{4}
```

The constants may depend on the local tensor, both gauges, the field,
and `s`. They are independent of `t` and `k`.

## 3. Telescoping and the complete coefficient budget

For either gauge, the product of the conjugated local matrices is
`G` times the original product times `G^{-1}`. Consequently the original
scalar product is recovered by the boundary vectors `f^T G^{-1}` and
`G e`. No auxiliary matrix has been commuted through another one.

For the first cut, a Smith channel word `i=(i_1,...,i_t)` has scalar
power `w^{-D(i)}`, where `D(i)=sum_j d_(i_j)`. After the boundary
regularization, its contribution to (1) is

```math
w^{-D(i)-H_\infty} L_i(w)R_i(w),                   \tag{5}
```

where both factors are regular series. The first factor includes the
whole ordered auxiliary product and its two boundary vectors. The second
is scalar and depends only on the other physical half-labels.

The coefficient of `w^{-k}` in (5) is zero unless
`D(i)>=k-H_infinity`. Otherwise, splitting the required nonnegative
degree between the two whole-word factors uses
`D(i)+H_infinity-k+1` rank-one terms. Since `k>=0`, this is at most

```math
C_\infty(t+1),\qquad
C_\infty=1+H_\infty+\max(0,\max_i d_i).             \tag{6}
```

For the second cut, a word has power `z^{N(i)-H_zero}`, where
`N(i)=sum_j nu_(i_j)`. It contributes to coefficient `k` only when
`N(i)<=k+H_zero`. Its number of full-word coefficient splits is at most

```math
C_0(t+1),\qquad
C_0=1+H_0+L+\max(0,-\min_i\nu_i),                  \tag{7}
```

because the original coefficient index is between zero and `Lt`.
This remains valid when the gauged local matrices have Laurent poles or
larger individual degrees than the original polynomial tensor.

Counting the admitted words with the same `s` in the two cuts gives

```math
\begin{aligned}
r(A_{t,k})&\le C_\infty(t+1)s^{H_\infty-k}U(s)^t,\\
c(A_{t,k})&\le C_0(t+1)s^{H_0+k}V(s)^t.
\end{aligned}                                     \tag{8}
```

The coefficient index cancels in their product. Absorb
`C_infinity C_zero s^(H_infinity+H_zero)` into a constant. There are
`Lt+1` coefficients in (1), proving both assertions in (4).

In particular, no factor exponential in `t` is charged for the gauges.
Counting independent coefficient splits at each site would give an
incorrect extra exponential cost. The splits in (6)--(7) are made after
the entire ordered product is formed.

## 4. Effective factoring criterion

For fixed rational gauges and a fixed rational `s>=1`, all data in
(2)--(3) have finite exact descriptions. Rational matrix inversion,
clearing poles, and polynomial Smith reduction give their integer
spectra by finite algebraic computation. These computations concern the
fixed certificate, not an input-dependent oracle.

As in the basic Smith transfer, complete coefficient tables of (1) and
their two rank factorizations cost `2^O(lambda)` bit operations over the
fixed field. Restrict from the least word length with `q^t>=2^lambda`
to binary labels. Equation (4) gives total channel width

```math
4\rho\lambda+O(\log(\lambda+2)),\qquad
\rho=\frac{\log B(s)}{4\log q}.                    \tag{9}
```

Thus a verified certificate with `B(s)<q^(2sqrt(2))` would yield an
exponent above one half through the fixed-field crossing transfer and
the [fixed-characteristic recovery](fixed-characteristic-factor-recovery.md).
For example the stronger exact inequality `B(s)^5<q^14` suffices.
The constructor and evaluator costs remain part of that transfer.
No such improving certificate is supplied here.

## 5. Cheap diagonal gauges

An especially simple search uses

```math
G(z)=\operatorname{diag}(z^{h_1},\ldots,z^{h_D}),
\qquad h_a\in\mathbb Z.                             \tag{10}
```

An entry with auxiliary indices `(a,b)` gains exponent `h_a-h_b`.
List these Laurent exponents exactly, subtract the common minimum to
obtain a polynomial matrix, and compute its Smith orders. Add the
common minimum back afterward. At infinity negate the restored orders
to get the degrees `d_i`. Forgetting to restore this scalar shift can
create a false improvement.

If both boundaries are the first coordinate vector and `h_1=0`, then
both boundary vectors are unchanged. In that case `H=0` and the
telescoping identity preserves the scalar polynomial entry by entry,
although internal auxiliary transitions have changed degree.

The two cuts may use different height vectors. A restriction proved for
the two spectra of one polynomial local matrix does not automatically
apply to this pair of separately gauged matrices. This observation does
not assert that the more general certificate will improve the exponent.

## Verification and attribution

Run `python3 tools/check_auxiliary_gauge_smith.py`. The checker uses
integer polynomials over `F_2`, explicitly multiplies small Laurent
matrix products, and compares them with the unchanged polynomial
crossings. It independently checks Smith orders through block-Toeplitz
ranks, restores every common Laurent shift, and tests the resulting
coefficient-rank bounds with exact rational arithmetic.

The basic coefficient-splitting argument comes from the linked local
Smith-profile note. The rational telescoping refinement and its checker
were developed with substantial OpenAI assistant support. Neither adds
Lean formalization or a new unconditional factoring exponent.
