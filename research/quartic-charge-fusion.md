# Quartic additive-charge fusion in a homogeneous cube

## Result and remaining cost

There is a polynomial cube tensor that breaks separate color parity,
preserves the antipodal endpoint, and has fewer order-four channels than
the one-hot Gaussian. Let physical colors be `F_2^s`, including vacuum
zero, and put `q=2^s>=4`, `d=q-1`. The coefficient field has characteristic
two; it need not contain a field of order `q`. The color operation is XOR.

The three generic axis Smith profiles begin exactly

```math
0^1,\quad 2^d,\quad 4^{2d},\qquad
\text{all other orders are at least five}.             \tag{1}
```

For `q=8`, this changes the Gaussian prefix `0^1,2^7,4^28` to
`0^1,2^7,4^14`. A 136-dimensional charge block verifies the change by
exact Toeplitz ranks, independently implemented over binary matrices.
The general statement follows from the leading-coefficient argument below.

The construction pays **degree four per occupied vertex**, rather than
three. No complete generic profile, favorable normalized pressure, full
coefficient-rank bound or stronger factoring exponent is established.
The complete score remains `1/2`. This is a positive construction and a
prefix theorem, not a new factoring bound.

The starting kernel and one-hot construction are the
[homogeneous Gaussian cube](cube-homogeneous-profile-budget.md) and the
[worker-4 color-class model](color-class-cube-profile.md), retained at
`1772970a7dbaab83eeb0607014a1111e46718d86`. The new interactions change
its parity hypothesis. The [color-parity floor](color-distinct-sector-obstruction.md)
therefore does not apply. This eight-leg three-axis model is separate
from the already excluded bare four-leg criterion through twelve letters.

## 1. Uniform creation-operator construction

For a vertex set `S`, a creation operator replaces vacuum on all vertices
of `S` by specified nonzero colors and is zero on every other input there.
Two such operators on overlapping supports have product zero. On disjoint
supports they commute. Each has square zero, so `I+f C` is invertible in
characteristic two. Sums of creations on the same support have the same
properties.

Put `ell=x+y+z`, and replace every Gaussian pair kernel by

```math
\widetilde K_{uv}=\ell^2 K_{uv}.                         \tag{2}
```

Starting from the vacuum, apply the pair-creation gates for every vertex
pair and every equal nonzero color. Their product is the color-class
Gaussian with entries multiplied by `ell^w`, where `w` is occupied
vertex count. Pair weights have degree eight, so this base has degree `4w`.

For every four-vertex set `S`, let `C_S` sum the creations of four
**distinct** nonzero colors with XOR zero, counting each assignment to
the four fixed vertices exactly once. Add the gates

```math
\prod_{|S|=4}(I+H C_S),\qquad
H=(xyz)^4\bigl(\ell^4+xyz\ell\bigr).              \tag{3}
```

The polynomial `H` has degree 16 and six monomials. Creating four vertices
therefore preserves degree `4w`. All gates commute; there are at most two
disjoint quartics in any term. This defines a uniform finite tensor, with
no quotient or uncharged auxiliary algebra. For fixed `q`, its complete
local polynomial tensor can be constructed by finite pair/quartic
expansion; the jet checker implements only the restricted calculation
needed below, not the full profile.

At `(1,1,1)`, `ell=1` and `H=0`. The pair construction retains its original
antipodal endpoint and every quartic gate becomes the identity. A labeling
with colors `1,2,4,7` on four vertices and vacuum elsewhere has entry
exactly `H`, so separate color parity is genuinely broken. The full tensor
has actual total degree 32: the all-one-color antipodal entry is nonzero
at the endpoint and homogeneous of that degree.

## 2. The leading coefficient matches in every axis

Fix the x cut and work over `F_2(y,z)[[x]]`. Every cross-pair entry has
valuation two and the same leading coefficient

```math
A=y^2z^2(y+z)^2.
```

The quartic has valuation four and

```math
[x^4]H=y^4z^4(y+z)^4=A^2.                              \tag{4}
```

The cyclic identities hold for y and z. With the original degree-three
grading, a quartic of degree twelve with these three leading coefficients
would have to be the monomial `(xyz)^4`, whose endpoint is one. Degree
four per occupied vertex is the smallest integer degree for this
simultaneous-leading-coefficient construction with a vanishing quartic
endpoint. Thus the cancellation is generic
and simultaneous in the three axes, not a consequence of a finite-field
specialization.

All creation gates supported wholly on one side of the cut can be removed
by invertible row or column operations. Since the gates commute, this
removal is exact even in the presence of crossing gates. Afterwards every
cross-pair contributes at least two powers of x, and every crossing quartic
at least four. An entry with more than four occupied vertices has valuation
above four. This bounds the entire low-jet calculation without enumerating
the full `q^4 by q^4` matrix.

## 3. Exact prefix proof

The all-vacuum entry is a disjoint unit. The one-occupied-vertex rows and
columns have one spatial `4 by 4` block per color. Its first two Smith
orders are `2,4`, with no other order at most four. The checker verifies
this generically in all axes: every entry has order two; the least
nonzero two-minor order is six; all three-minors have order at least eleven.
Multiplication by `ell^2` is a local unit and preserves these orders.

Now consider two distinct colors on each side. For an unordered pair
`{a,b}`, the Gaussian order-four coefficient is `A^2` when the column
has the same unordered pair, and zero otherwise. Within each such pair,
all spatial placements have the same leading entry.

The quartic adds `A^2` precisely when the two pairs differ but have the
same XOR. Such pairs are disjoint: sharing one color and having the same
XOR would force the other color to agree too. Consequently the new
order-four matrix is all-one on every nonzero-charge block. There are
`d` such charges. For each charge there are `(q-2)/2` unordered color
pairs; its old rank was `(q-2)/2` and its new rank is one. Thus the
color-pair rank falls from `binom(d,2)` to `d`.

Repeated equal colors have charge zero. Their two cross-matchings have
the same leading coefficient and cancel in characteristic two, and the
quartic requires distinct colors. Hence they contribute no order-four
channel.

A crossing quartic split `3+1` or `1+3` contributes at order four, but is
constant across the position of its single color. This vector belongs to
the existing order-two spatial span. A row or column operation of order
two removes it through order four. The resulting `3+3` Schur correction
has order at least six and does not change this prefix. The two-occupied
sector does not couple to the one-occupied sector at these orders.

It follows that order two has multiplicity `d`, while order four has
`d` spatial channels plus `d` charge-fused color-pair channels. There
are no odd-order channels below five. This proves (1) for every `q>=4`.
The first strict improvement occurs at `q=8`; at `q=4` each charge already
has only one unordered color pair.

## 4. Exact finite controls and continuation test

```sh
make verify-quartic-charge-fusion
```

The [checker](../tools/check_quartic_charge_fusion.py) uses only the
standard library and the existing polynomial kernel. It checks the
endpoint, homogeneous degrees, all 48 cross-pair leading identities,
108 two-minors and 48 three-minors. Its q8 nonzero-charge block has
4 one-vertex, 36 two-vertex and 96 three-vertex words. At precisions
one through five its Gaussian Toeplitz ranks are `0,0,1,2,7`, and its
fusion ranks are `0,0,1,2,5`.

The coefficient field is `F_256` with modulus `0x11B`, `y=2,z=3`.
The checker verifies that all nonzero coefficients are units, represents
multiplication by each coefficient as an eight-bit binary map, and divides
the exact binary Toeplitz rank by eight. A separate ordinary binary
color-pair matrix verifies ranks 21 and 7. The
[receipt](../evidence/quartic-charge-fusion.json) pins this source.
An earlier independent NumPy F256 calculation gave the same low-jet ranks;
it is not required by the committed replay.

The next question is whether the higher channels move far enough to pay
the extra degree. For a common physical type of entropy `H(p)>0` and mean
occupancy `E`, the correctly restored objective is

```math
\Phi_p(t)=
 \frac{\frac13\sum_a\log Z_a(t)-\frac{32}{3}E\log t}{4H(p)}.
                                                               \tag{5}
```

The three-axis target remains `rho<1/sqrt(3)`, with a valid complete
transfer. Substituting the old degree-three restoration into (5) would
be incorrect. A prefix alone cannot evaluate (5). Further work must give
an explicit growing-color argument for higher orders or a complete profile
with the restored cost; a lower first-layer rank is not a score gain.
No larger matrix or q/field/weight sweep is justified solely by this control.
