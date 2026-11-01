# All-scale dominance in the explicit permutation cone

## Result and scope

The [arbitrary-permutation constructor](permuted-pascal-zero-spectra.md)
has an explicit two-cut cone. Inside that cone, every permutation is
weakly dominated by the original ordered Pascal kernel, at every
common scale and every proposed contraction `rho>1/2`.

This is a comparison of complete pressures with both physical pair
charges included. It is not a claim that individual zero spectra
cannot increase. It removes this particular explicit cone from the
search for a better Pascal pressure. Arbitrary physical gauges,
changes of flags outside the cone, and restrictions of the physical
alphabet remain outside the theorem.

The proof is analytic and holds at all sizes. The
[checker](../tools/check_permuted_pascal_dominance.py) also performs
exact rational pressure comparisons; its
[receipt](../evidence/permuted-pascal-cone-dominance.json) uses no
floating logarithms as proof certificates. No factoring score changes.

## 1. Write the entire pressure residual

Put `L=2^r`, `t>0`, `q=2rho>1`, and `ell(z)=log(1+exp(z))`.
For an exponent permutation `pi`, write

```math
 a_i=E_0(i),\quad e_i=E_\infty(i)=2\pi(i)+1,
 \quad b_i=\beta(i).
```

All `a_i,e_i` lie strictly between zero and `2L`. The explicit
infinity cone requires `b_i<=b_j` under Boolean containment, but
we will first relax that condition entirely.

The logarithmic residual whose nonpositivity would certify the
proposed pressure is the sum over `i` of

```math
 \begin{aligned}
 F(a,e,b)={}&2\ell(-ta)
 +\ell(t(2L-\min\{b,e\}))
 +\ell(-t\max\{b,e\})\\
 &-q\ell\!\left(-\frac{t(a+b-2L)}q\right)
 -q\ell\!\left(-\frac{t(a+e-2L)}q\right).
 \end{aligned}
 \tag{1}
```

The two first-order zero energies are both `a`; the infinity orders
are `min(b,e)-2L` and `max(b,e)`. The two physical pair energies
are `a+b-2L` and `a+e-2L`. Thus (1) includes every one-particle
mode and both physical type partition functions of the preceding
construction. No charge is dropped when `a` changes.

## 2. Even without the cone, the best infinity height is `b=e`

Write `sigma(z)=exp(z)/(1+exp(z))`, a strictly increasing function.
For `b<e`, differentiation of (1) gives

```math
 \partial_b F
 =-t\sigma(t(2L-b))
  +t\sigma\!\left(-\frac{t(a+b-2L)}q\right)<0.
 \tag{2}
```

Indeed `b<e<2L` and `a>0` imply

```math
 \frac{a+b-2L}{q}>b-2L.
```

For `b>e`, the derivative is

```math
 \partial_b F
 =-t\sigma(-tb)
  +t\sigma\!\left(-\frac{t(a+b-2L)}q\right)>0,
 \tag{3}
```

because `b>e>0`, `a<2L`, and `q>1` imply
`(a+b-2L)/q<b`. The function is continuous at `b=e`. Therefore

```math
 F(a,e,b)\ge F(a,e,e)
 \quad\text{for every real }b.
 \tag{4}
```

The point `b=e` need not satisfy the original containment constraints.
Using it only makes the lower bound more favorable. Set

```math
 \Phi(a,e)=F(a,e,e)
 =2\ell(-ta)+\ell(t(2L-e))+\ell(-te)
  -2q\ell\!\left(-\frac{t(a+e-2L)}q\right).
 \tag{5}
```

## 3. The sorted zero energies cannot be smaller than the odd list

The polynomials

```math
 F_i(T)=\sum_{j\supseteq i}T^{\pi(j)},\quad 0\le i<L,
```

form a basis of all degree-below-`L` polynomials: the ordinary upset
matrix is invertible, and the permutation just reorders its coefficient
coordinates. The subspace of such polynomials with order at least `k`
at `T=1` has dimension `L-k`.

Hence at most `L-k` of these independent basis vectors can have
valuation at least `k`. If their valuations are ordered from largest
to smallest, the `j`th one is at most `L-j`, for `1<=j<=L`.
Since `a_i=2L-1-2v(F_i)`, their increasing rearrangement satisfies

```math
 \boxed{a_{(j)}\ge2j-1\quad(1\le j\le L).}
 \tag{6}
```

This dimension argument needs no monotonicity or linearity of `pi`.

## 4. Larger zero energies worsen the fully charged residual

For `0<a,e<2L` and `q>1`, differentiating (5) gives

```math
 \partial_a\Phi
 =-2t\sigma(-ta)
  +2t\sigma\!\left(-\frac{t(a+e-2L)}q\right)>0,
 \tag{7}
```

since `(a+e-2L)/q<a`. Thus, after sorting the `a_i`, replacing
them coordinate by coordinate by `2j-1` can only lower the residual.
The exponent energies remain a permutation of the same odd list.

The mixed derivative is

```math
 \partial_a\partial_e\Phi
 =-\frac{2t^2}{q}\sigma(z)(1-\sigma(z))<0,
 \qquad z=-t(a+e-2L)/q.
 \tag{8}
```

Integrating (8) over a rectangle shows that pairing small `a` with
small `e` and large `a` with large `e` minimizes the sum. Explicitly,
for `a<=a'` and `e<=e'`,

```math
 \Phi(a,e)+\Phi(a',e')
 \le\Phi(a,e')+\Phi(a',e).
 \tag{9}
```

Successive exchanges therefore put the two odd lists in the same
order. Combining (4), (6), (7), and (9) proves

```math
 \boxed{\sum_i F(a_i,e_i,b_i)
 \ge\sum_{j=1}^L\Phi(2j-1,2j-1).}
 \tag{10}
```

The right side is exactly the original ordered Pascal residual at
that same `t,q`. Its matching infinity vector is admissible under
Boolean containment. Thus the comparison does not substitute an
unrealizable object for the final reference kernel.

The pressure residual decreases with `q` wherever the physical
alphabet has positive entropy. Equation (10), valid at every scale,
therefore gives the stated rate comparison. It transfers any proved
floor for the ordered Pascal pressure to this whole explicit cone;
it does not assert a floor for more general gauges or restrictions.

## 5. Exact finite checks of the complete comparison

The checker verifies (6) for all 40,320 permutations at length eight
and checks the affine inequalities underlying every derivative sign.
For complete pressure comparisons it uses rational `q=A/B>1` and
`t=A log u`, where `u>1` is rational. Then

```math
 \exp\bigl(B\sum_i F_i\bigr)
 =\frac{(Z_0 Z_\infty)^B}{Z_{\rm physical}^{A}}
```

is an exact rational number, since all exponents in the three
partition functions are integers. The checker compares each stage
of (10) by integer cross multiplication, including nonmonotone
permutations and height vectors of both signs. These finite controls
check the formula and charge bookkeeping; the derivative and
filtration arguments above establish the all-scale theorem.

```sh
python3 tools/check_permuted_pascal_dominance.py \
  --output /tmp/permuted-pascal-cone-dominance.json
cmp /tmp/permuted-pascal-cone-dominance.json \
  evidence/permuted-pascal-cone-dominance.json
```
