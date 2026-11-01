# A common coordinate core rejects the coupled quadratic coefficient search

## Result and search consequence

For the full six-mode kernels of

```math
T_8=D_0D_3+aD_1D_2,\qquad
T_{12}=D_0D_5+bD_1D_4+cD_2D_3,
\qquad a,b,c\ne0,                                    \tag{1}
```

the natural balanced pressure is at least **5/7**, over every
characteristic-two field, for both the canonical coordinate retraction
and the orthogonal retraction when it exists. The same bound applies
to any homogeneous split whose maps act literally as coordinate maps
on the protected space identified below.

The value 5/7 exceeds the post-half threshold `1/sqrt(2)`, since
`2*5^2-7^2=1`. Thus extending the coefficients of this particular family
from F4 to F16, or to a larger field, cannot improve the factoring exponent
through its inherited natural pressure. This exact screening result
prevented a full 851-class F16 profile scan after a 24-class pilot.
It does not cover arbitrary retractions that modify the protected maps,
other physical gauges, other kernels, or a different transfer criterion.

This note uses the protected-submatrix principle already established in
[the height-cutoff obstruction](coherent-protected-height25.md).
The new finite input is a smaller coordinate core tailored to (1).
There is no new factoring exponent or Lean formalization.

## 1. The height budget is independent of the coefficients

Rescale the six deletion variables by

```math
(d_0,d_1,d_2,d_3,d_4,d_5)
=\left(1,\frac ca,\frac1c,1,\frac a{bc},1\right).       \tag{2}
```

Every coefficient in the two operators becomes one. On the occupation
basis, the corresponding diagonal change multiplies the word w by
the product of d_i over its occupied modes. Conjugating a deletion
operator multiplies D_i by d_i. Therefore (2) is an invertible change
between the kernel in (1) and the kernel with `a=b=c=1`.

This diagonal change preserves physical height and occupation number.
Every nonzero coefficient choice consequently has the same dimension
in every graded cell. The full kernel has 28 letters, with occupation
counts `(1,6,13,8)`. Its physical height polynomial A has coefficients,
starting at height zero,

```text
1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,1,2,1,1,1,1,2,0,1,0,1,0,1.
```

The checker reconstructs this budget from the binary coefficient kernel.
Equation (2), rather than a finite sample of field values, proves that
the budget is valid over every field in the theorem.

## 2. Twenty-one coordinate states are retained literally

Let U be the coordinate span of words containing none of the five pairs

```math
\{0,3\},\{1,2\},\{0,5\},\{1,4\},\{2,3\}.            \tag{3}
```

These are the independent sets of the path with vertex order
`5,0,3,2,1,4`. There are 21 such words. Each individual monomial in
(1) kills every word of U, so U is contained in every kernel. It is
also closed under all bit deletions.

The canonical reduced kernel basis contains these coordinate vectors
unchanged. Its other basis vectors can be reduced to have no coordinates
in U. The canonical coordinate retraction therefore acts as the identity
on U, and all its other rows vanish there.

The same holds for the orthogonal retraction. In this basis the kernel
Gram matrix is the identity on U, with no off-diagonal block joining U
to its coordinate complement. Inverting the Gram matrix preserves that
block structure. Thus both P and L restrict literally to the coordinate
maps on U for the two retractions under consideration.

Restricting all logical labels of either regular cut to U gives the
same protected coordinate tensor. Its two complete profiles agree;
let their counting polynomial be Z. It has 441 terms counted with
multiplicity, degree 92, constant coefficient one, and coefficient one
at degree one.

The protected-submatrix principle gives

```math
Z_0(r)Z_\infty(r)\ge Z(r)^2\qquad(0<r<1).             \tag{4}
```

The fixed coordinate profile is binary, so it remains unchanged after
every characteristic-two field extension. This argument protects the
actual polynomial matrices, not guessed high Smith orders.

## 3. An integer certificate gives the pressure floor

The exact scalar certificate is

```math
Z(s^5)^7>A(s^7)^{10}\qquad(0<s\le1).                  \tag{5}
```

Set `rho_0=5/7` and `r=s^5`. Equation (5), together with (4), implies

```math
\log Z_0(r)+\log Z_\infty(r)
>4\rho_0\log A(r^{1/\rho_0}).                         \tag{6}
```

For any physical type with entropy H and expected height E, its Gibbs
inequality is

```math
H+\frac{\log r}{\rho_0}E
\le\log A(r^{1/\rho_0}).
```

Substitution into the natural pressure expression proves that every
finite scale and every positive-entropy type have pressure strictly
greater than 5/7. The infimum is therefore at least 5/7. This does not
promote pointwise strictness to a strict infimum bound. Applying the
two Gibbs inequalities separately also allows independent physical types
on the two paired alphabets, with the same protected-core hypotheses.

### The interval near zero

For `0<s<=1/2`, write `t=s^7`. The height budget gives

```math
A(t)\le1+t+26t^3\le1+2t,
```

because `26t^2<=26/2^14<1`. A binomial bound and a geometric series give

```math
A(s^7)^{10}\le(1+2s^7)^{10}
\le1+\frac{20s^7}{1-20s^7}.
```

On the other hand, `Z(s^5)^7>=1+7s^5`. The strict comparison follows
from

```math
7(1-20s^7)-20s^2
\ge7(1-20/128)-5=116/128>0.
```

### The remaining compact interval

On `[1/2,1]`, both sides of (5) are increasing functions of s. The
checker gives 124 rational intervals, with maximum dyadic depth ten.
On each interval the left side at its lower endpoint strictly exceeds
the right side at its upper endpoint. Exact endpoint adjacency and
coverage are checked. Every comparison uses integers.

## 4. Replay and computational provenance

Run

```sh
python3 tools/check_tensor_coupled_quadratic_core_floor.py --recompute-profiles
```

The [checker](../tools/check_tensor_coupled_quadratic_core_floor.py) and
[receipt](../evidence/tensor-coupled-quadratic-core-floor.json) retain
the protected words, both complete profiles, the height polynomial,
all rational intervals and source hashes. The recomputation mode
checks both cuts using the established F4 implementation and the new
F16 implementation, for 1764 directly compared Smith orders.
Omitting the option replays the integer certificate from its saved input.

The numerical optimum of the optimistic protected-core expression is
approximately `0.7182021973`; this value was used only to choose the
simple rational target in (5). The certified theorem is the exact
lower bound 5/7.
