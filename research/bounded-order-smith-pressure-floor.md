# Smith orders bounded by the degree exclude alphabets through ten

## Result and scope

Let a scalar polynomial crossing have actual degree `M` after removal
of its common scalar monomial, and alphabet `2<=q<=10`. Suppose every
Smith order of its row-infinity reversal and its column-zero matrix is
at most `M`. Then its bare local pressure satisfies

```math
\boxed{\rho_{\rm bare}>\frac{442}{625}
       >\frac1{\sqrt2}.}                         \tag{1}
```

This is an all-degree, all-field theorem under the displayed order
hypothesis. In particular, an eight-letter bare crossing that improves
the target must have a Smith order larger than its actual polynomial
degree in at least one of the two cuts.

The hypothesis holds when each of the two cut matrices admits a diagonal
form by constant invertible row and column changes. Consequently (1)
closes the polynomial order/degree pressure method in the complete
eight-letter [Bell diagonal space](bell-diagonal-crossings.md), at every
polynomial degree and in every characteristic different from two. It does
not close arbitrary collective coefficient recombinations in that space,
whose actual coefficient ranks can be better than the order/degree bound.
Independent physical gauges and matrices with larger Smith orders are
also outside the result. No factoring exponent is improved.

## 1. Only a prefix budget and upper bounds are needed

Put `n=q^2`. Let `mu` and `nu` be the increasing row-infinity and
column-zero Smith profiles. The
[shared-minor theorem](shared-minor-smith-floor.md) gives

```math
\mu_1=\nu_1=0,\qquad
\sum_{i=1}^q(\mu_i+\nu_i)\le qM.                 \tag{2}
```

The additional hypothesis is

```math
0\le\mu_i,\nu_i\le M\quad(1\le i\le n).        \tag{3}
```

For `M=0` the pressure is `q^4`, so assume `M>0`. Define

```math
z_i=\frac{\mu_i+\nu_i}{2M},\qquad t=M\log s.
```

Then `z_1=0`, every `z_i<=1`, and `sum_(i=1)^q z_i<=q/2`.
The restored product is

```math
B(s)=s^M\left(\sum_i s^{-\mu_i}\right)
            \left(\sum_i s^{-\nu_i}\right).
```

Cauchy--Schwarz gives

```math
\sqrt{B(s)}\ge\sum_i\exp\bigl(t(1/2-z_i)\bigr). \tag{4}
```

There is no assertion that a general degree-`M` polynomial matrix
satisfies (3). A small entry degree can coexist with much larger Smith
orders through cancellation in its minors. The condition is used
explicitly here.

## 2. One rational probability vector works through ten letters

Set `D=20q(q-1)` and choose a probability vector whose integer counts
with denominator `D` are

```math
\left((7q+6)(q-1),\quad
      [6(q-1)]^{\times(q-1)},\quad
      7^{\times(q^2-q)}\right).                  \tag{5}
```

The first coordinate multiplies `z_1=0`. The next `q-1` coordinates
each have probability `3/(10q)`, and the last `q^2-q` coordinates
each have probability `7/(20q(q-1))`. Therefore

```math
\sum_i p_i z_i
\le\frac3{10q}\frac q2+
   \frac7{20q(q-1)}(q^2-q)
=\frac12.                                       \tag{6}
```

Weighted AM--GM applied to (4), or the finite Gibbs inequality, yields

```math
\log\sqrt{B(s)}\ge H(p)+t\left(\frac12-\sum_i p_i z_i\right)
\ge H(p).
```

Hence `rho=log B/(4 log q)>=H(p)/(2 log q)` at every scale.
For each integer `2<=q<=10`, the exact replay proves

```math
H(p)>\frac{884}{625}\log q.                      \tag{7}
```

At `q=10`, the counts are
`(684,54 repeated9,7 repeated90)` with denominator `1800`. The
required comparison is exactly

```math
1800^{45000}
>10^{63648}\,684^{25\cdot684}
  54^{25\cdot54\cdot9}\,7^{25\cdot7\cdot90}.     \tag{8}
```

The left and right integers have respectively 486,621 and 486,618
bits. The checker compares them directly, without logarithm estimates.
The other alphabet sizes use the same count formula and clear the
rational entropy inequality in the same way. Finally
`2*442^2>625^2`, establishing (1).

## 3. Why the Bell polynomial interface satisfies the hypothesis

Suppose a cut has constant invertible matrices `U,V` for which

```math
U K(z)V=\operatorname{diag}(f_1(z),\ldots,f_n(z)).
```

Every `f_i` has degree at most `M`, because the changes are constant.
The invertible crossing endpoint ensures `f_i(1)!=0`. Thus each
`ord_0 f_i` lies between zero and `M`, and each reversal order
`M-deg f_i` does too.

This argument must hold for both physical cuts, with their own
constant changes. The Bell construction supplies those two changes
explicitly: the row cut is diagonal with eigenpolynomials `s_j(z)`,
and the column cut is diagonal with eigenpolynomials `(Ts(z))_j`.
Both bases are constant and invertible in characteristic different from two;
the Bell transform `T` is also constant. At the crossing endpoint,
both eigenvalue vectors are the all-ones vector. Therefore both
profiles obey (3), at arbitrary polynomial degree. Exchanging the two
physical cuts gives the same argument for the opposite orientation.

The [polynomial synthesis bound](bell-polynomial-synthesis.md) is
exactly this order/degree pressure method, so the theorem rules out
using a longer polynomial path to cross the target at physical
dimension eight. The theorem does not replace actual coefficient-rank
calculations by order intervals when those calculations exploit extra
cancellation.

## 4. These inequalities alone stop excluding at eleven letters

A formal profile pair shows that the prefix-and-upper-bound relaxation
cannot simply be extended to `q=11`. Take `M=20` and

```math
\mu=\nu=(0,11^{\times10},20^{\times110}).         \tag{9}
```

The profiles satisfy (2)--(3), as well as the ordinary determinant
budgets. At `s=5/4`, put

```math
P=5^{20}+10\,4^{11}5^9+110\,4^{20}.
```

Their formal pressure is `P^2/(5^20 4^20)`. Exact arithmetic gives

```math
P^{500}<11^{707}5^{5000}4^{5000},
```

so the formal exponent is below `707/1000<1/sqrt(2)`. This is a
feasible profile witness, not a polynomial crossing construction.

## 5. Replay

```sh
python3 tools/check_bounded_order_smith_floor.py \
  --output evidence/bounded-order-smith-floor.json
```

The [checker](../tools/check_bounded_order_smith_floor.py) verifies the
probability coefficients and all nine entropy comparisons exactly,
then checks the eleven-letter formal witness. Its
[evidence](../evidence/bounded-order-smith-floor.json) records integer
sizes and hashes. No floating arithmetic enters any certificate.
