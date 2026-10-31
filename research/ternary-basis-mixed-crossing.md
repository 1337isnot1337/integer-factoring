# A ternary mixed crossing below 821/1000

## Result

A physical basis change and three explicit local gauges improve the
ternary mixed-place construction from `rho < 9/10` to

```math
\rho<\frac{821}{1000}.
```

The endpoint, all three complete local Smith profiles, and the strict
pressure inequality have exact certificates. The construction works at
every label width by the existing fixed-type transfer. The parameter
remains above `1/sqrt(2)`, so it does not raise the established factoring
exponent `1/2`.

## 1. An endpoint-preserving basis change

Start with the ternary seed

```math
K(z;x,y,u,v)=[x+y=u+v]\bigl(1-(x-u-z(v-u))^2\bigr)
```

over `F3`, with degree two and `K(1;x,y,u,v)=[x=v][y=u]`.
The arithmetic inside the formula is in `F3`. Let

```math
A=\begin{pmatrix}1&0&0\\2&1&0\\1&1&1\end{pmatrix}.
```

Apply `A` to the physical `x,y` legs and `(A^-1)^T` to `u,v`.
Call the resulting tensor `K'`. Invertibility and the inverse transpose
pairing prove its literal crossing endpoint. Its degree is still two.
This change preserves the unweighted profiles in both cuts, but it
changes the available physical-coordinate gauges.

Use the uniform type on all four physical legs. At a tensor length
divisible by three, each retained word contains every symbol equally
often. The sum of its symbols, viewed as ordinary integers `0,1,2`,
is exactly the tensor length. The common paired-type entropy sum is
`H=2 log3`.

## 2. Exact weighted profiles

Choose row places `2,infinity`, each with allowance `1/2`, and the
column place `0` with allowance one. All three are distinct and avoid
the endpoint. At a selected place the height of symbol `i` on a leg
is its displayed slope multiplied by the integer `i`.

| Place and cut | Slopes on `x,y,u,v` | Mean shift | Complete local orders after `t=tau^2` |
| --- | --- | --- | --- |
| `2`, first cut | `0,0,1/2,1/2` | `1` | `2,3,3,4,4,4,5,5,6` |
| infinity, first cut | `-1/2,-1/2,0,0` | `-1` | `-2,-1,-1,0,0,0,1,1,2` |
| `0`, second cut | `-1,-1,1,1` | `0` | `0,2,2,4,4,4,6,6,8` |

The profiles contain every one of the nine generic channels. Negative
orders at infinity are retained. The exact checker obtains these lists
from complete block-Toeplitz ranks after a scalar clearing power; it
does not infer them from a selected determinant or minor frame.

The local ramification is only an analysis of the original Hasse jets.
It does not change the polynomial degree budget. Under `t=tau^2`, an
original jet of order `j` is a coefficient of order `2j`, and after the
physical type restriction its gauged index is `2j+2mA`. Use scale
`s^(1/2)` in the ramified coefficient bound. This restores the precise
factor `s^(j+mA)` in the original local coordinate. Each selected
place can be analyzed separately in this way.

After dividing the displayed orders by two and subtracting the mean
shift, either row place has effective orders

```math
0,\tfrac12,\tfrac12,1,1,1,\tfrac32,\tfrac32,2.
```

The column place has effective orders

```math
0,1,1,2,2,2,3,3,4.
```

Thus the [weighted mixed-place criterion](mixed-place-hermite-crossings.md#5-independent-paired-types-and-local-physical-gauges)
has the especially simple exact expressions

```math
R(s)=C(s)=(1+s^{-1}+s^{-2})^2,
\qquad B(s)=s^2(1+s^{-1}+s^{-2})^4.
```

The two row allowances account for the factor of two in their local
scale. The mean-height correction is essential: omitting it would
incorrectly treat the profiles at the two row places differently.

## 3. A strict rational certificate

Take `s=23/10`. Then

```math
B(s)=\frac{859^4}{100\,23^6},
\qquad
(859^4)^{250}<3^{821}(100\,23^6)^{250}.
```

The inequality is checked directly with integers. It is exactly
`B(s)<3^(4*821/1000)`, the sufficient fixed-type crossing inequality.
Its strict margin also permits a rational parameter slightly below
`821/1000`. The opposite-flag reconstruction, type-class padding,
pair-swapping, and complete `2^O(lambda)` table construction are those
proved in the linked mixed-place transfer.

Within this explicit one-parameter formula, differentiation gives the
unique minimizing scale

```math
s_*=(1+\sqrt{13})/2,
\qquad
\frac{\log B(s_*)}{4\log3}\approx0.8203391769.
```

The decimal describes this fixed formula. It is not a claim of global
optimality over bases, types, local heights, or polynomial seeds. The
proved witness uses only the rational scale and integer comparison.

## 4. Reproduction and limitations

Run

```sh
python3 tools/check_ternary_basis_mixed.py \
  --output evidence/ternary-basis-mixed.json
```

Besides the literal endpoint and complete local profiles, the checker
builds the full uniform-type restriction at length three. It checks all
1,296 endpoint entries, the common physical-height shift at every entry,
and 21 exact coefficient ranks against the channel-count upper bound.
The uniform theorem uses the strict pressure certificate, rather than
extrapolating those finite checks.

This result improves the previous mixed-place ternary witness. It does
not claim a separation from every weighted same-place construction,
nor a new factoring exponent. The seed and transfer are prior work
linked above; the basis/gauge witness and independent exact checks were
developed with OpenAI assistant support.
