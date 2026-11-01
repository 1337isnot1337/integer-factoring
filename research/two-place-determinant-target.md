# A determinant target for asymmetric improvement, and this pencil's exact ceiling

## Result and boundary

At the retained five-letter physical heights, an opposite-cut determinant
of degree **at most 69** would suffice for a supporting line strictly
below the half-score threshold. The complete two-place pencil from
[PR57](https://github.com/1337isnot1337/integer-factoring/pull/57), source
`8a3ea3e6aeacba7d0a53764de5579e3707bb2ca7`, has determinant degree
**at least 78 over every characteristic-two field**.

The first statement is a sufficient target for other valid kernels at
these heights. The second excludes that particular escape for the whole
eleven-parameter pencil. It does not exclude every asymmetric pressure
improvement, other physical heights, or other generators. No factoring
exponent above **1/2** is established.

## 1. The sufficient mean condition

Let `nu_i` be the 25 raw infinity orders and `mu_i` the 25 zero orders
of a strict five-letter path. Keep the physical pair energies

```math
 e=(0,10,6,2,8)
```

on both matching pairs, and suppose every `mu_i>=8`. Put
`Z_infinity(t)=sum_i exp(-t nu_i)`, `Z_0(t)=sum_i exp(-t mu_i)`
and `F(t)=sum_j exp(-t e_j)`.

For positive `a,b` with `a^2+b^2=1`, the
[supporting-line criterion](asymmetric-support-line-crossings.md) at
target `R=1` is strict negativity of

```math
 M(a,b,t)=a\log Z_\infty(t/a)+b\log Z_0(t/b)-2\log F(t).
```

Write `S=sum_i nu_i`. If **`S>60`**, take
`x=100/(S-60)`, let `b` tend to zero, put `a=sqrt(1-b^2)` and
`t=bx`. Taylor expansion at the finite uniform distributions gives

```math
 \lim_{b\to0+}\frac{M(a,b,bx)}b
 =\log Z_0(x)+x\left(\frac{52}5-\frac S{25}\right)
 \le\log25-x\left(\frac S{25}-\frac{12}5\right)
 =\log25-4<0.
```

The bound uses `Z_0(x)<=25 exp(-8x)`. Thus `M<0` for sufficiently
small positive `b`. Rational Pythagorean points on the unit circle are
dense at this endpoint, so rational `a,b,t` can retain the strict gap.
Continuity then gives a rational `R<1`; rational type approximation
retains the strict paired Gibbs inequality. The
[complete rectangular transfer](asymmetric-modular-transfer.md) supplies
a uniform factoring exponent strictly above `1/2` once such an actual
effective path has been supplied. The mean condition alone does not
construct that path.

The infinity row-plus-column height mass is 130 here. The determinant
identity for Smith orders therefore gives

```math
 S=130-\deg\det C(t).
```

Consequently `S>60` is equivalent to `deg det C<70`. This condition
is sufficient, not necessary, for an asymmetric improvement.

## 2. Exact full determinant retained by the predecessor

The predecessor uses

```math
 A_q=A+q(t)(K+t^4L),\qquad \deg q\le11,\quad q(1)=0.
```

Its exact commutator and sandwich identities are reconstructed in
[the input constructor](../tools/two_place_f4_census.py). They hold over
every characteristic-two field. After reversing the variable, write

```math
 g(z)=\frac{z^{12}q(z^{-1})}{(1+z)^4}.
```

The normalized infinity matrix is `C(z)+g(z) Delta(z)`, with common
pole eight removed. The full determinant coefficient identity already
retained in PR57 is

```math
 D(z,g)=z^{244}P_0(z)+g^4z^{204}P_4(z),
```

where

```math
 P_0=1+z^4+z^6+z^{10}+z^{12}+z^{14}+z^{16},
 \qquad
 P_4=1+z^6+z^8+z^{16}+z^{22}+z^{24}.
```

The [new checker](../tools/check_two_place_determinant.py) recomputes this
whole determinant by fresh exact fraction-free elimination, instead of
assuming the saved coefficients. It substitutes `g=z^512`. Each
coefficient of every relevant minor has local degree less than 512,
so the result decodes without aliasing. This verifies a polynomial
identity over `F2`, and hence over every characteristic-two field.

## 3. Complete all-field valuation classification

The raw infinity sum is `ord_z D-200`. For `q=0` it is 44. For nonzero
`q` of degree `d`, the two determinant terms have orders 244 and
`252-4d`. They can tie only when `d=2`. Therefore degree one gives
sum 44, and every degree `3<=d<=11` gives **`52-4d`**.

For degree two, write `q=a t^2+b t+(a+b)`, with `a!=0`. After
clearing the unit denominator, the low coefficients in the basis
`1,a^4,b^4` are

| Local degree | Coefficient |
|---:|---|
| 0 | `1+a^4` |
| 4 | `1+b^4` |
| 6 | `1+a^4` |
| 8 | `b^4` |

All other coefficients below degree eight vanish. Frobenius is injective
in every field of characteristic two. Thus `a!=1` gives sum 44;
`a=1,b!=1` gives sum 48; and `a=b=1`, namely `q=t^2+t`, gives sum 52.
This proves the exact universal maximum **52** and minimum determinant
degree **78**. No field sampling is used in this classification.

The same constructor proves that the weighted zero perturbation is
`t^8 q(t)` times a constant matrix, whose `(0,0)` entry is zero.
The reference weighted matrix has minimum order eight and unchanged
unit leading `(0,0)` entry. Hence its minimum remains exactly eight
for every member, justifying the sufficient-target comparison above.

Run `make verify-two-place-determinant`. The
[receipt](../evidence/two-place-determinant.json) retains the full exact
coefficients, the coefficient-degree bound, formal quadratic checks,
and source hashes. This result extracts a new consequence of the
predecessor's determinant identity and retains its credit; it does not
claim that the determinant itself was absent from PR57.

Substantial OpenAI assistant support was used for the proof, fresh
verification, code and text. No new Lean formalization is claimed.
