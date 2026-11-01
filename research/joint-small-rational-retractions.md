# A rational joint row module and a fixed third Smith order

The complete polynomial-output model for rational joint left retractions
of the five-letter source onto its retained four-letter grid has nine
row generators. Its polynomial right inverse has degree seven, so the
model permits arbitrary poles in the internal maps.

Within this full rational class, impose the natural degree bounds, freeze
the first two complete raw matrix jets at both cuts, and require no order
two at either cut. Then the **third zero-cut Smith order is exactly three**
over every characteristic-two field. A particular normalized 3-by-3 minor
always has coefficient one at degree four. This is an ordinary symbolic
proof, without a finite-field solver assumption.

This obstructs the attempted simultaneous lift of the third order while
keeping the earlier raw jets. It does not exclude improvements in later
orders, changing earlier jets, different embeddings or other weights.
No stronger factoring exponent follows.

A simpler invariant already holds when only the first raw matrix jet
(the leading matrix) is fixed at both cuts. At each cut the principal
2-by-2 minor on indices `(0,1)` has coefficient one at degree one.
Consequently the second Smith order is one at both cuts. The complete
479-dimensional correction space is checked direction by direction;
two changes beginning at degree one cannot affect this coefficient,
so this linear calculation proves the identity over every
characteristic-two field. Higher orders and changed leading factors
remain open.

## Complete rational model

Use the same source `A(z)`, embedding `P0` and retraction `L0` as the
[moving joint construction](joint-moving-retractions.md). Let `Q` insert
the nine discarded coordinates and set `B=Q^T A P0`. Every rational left
inverse of `P0` is `L=L0+T Q^T`, with output `C=C0+T B`.

The [implementation](../tools/joint_small_rational_retractions.py) adapts
and credits the existing [rational row-module method](joint-hot-rational-retraction-module.md).
This is a different, smaller source and retained space. The
[receipt](../evidence/joint-small-rational-retractions.json) supplies
polynomial matrices satisfying `S R=I_9`, `B=U S`, and `det(U)!=0`.
The nine row degrees of `S` are `3,3,0,0,0,9,2,12,9`. Literal multiplication
checks its degree-seven right inverse. Invertibility of `U=B R` at
`z=omega` certifies a nonzero binary determinant polynomial; this rational
rank persists under every characteristic-two extension. The leading row
vectors of `S` are independent.

A polynomial row in the rational span has polynomial coefficients in
`S`: multiply it by `R`. If the row vanishes at one, those coefficients
are divisible by `z+1`. Thus the complete strict polynomial-output model is

```
C = C0 + (z+1) G(z) S(z),       G any polynomial 16-by-9 matrix.
```

Conversely `T=(z+1)G U^-1` realizes every such output. No finite pole
location or pole-order cutoff on `T` is used.

## Complete bounded spaces

Let `H_(x,y)=h_x+h_y`, where `h=(0,1,3,4)`. The regular matrices are
`M_zero=z^-8 D_H C` and
`M_infinity=reshuffle(z^8 C(z^-1) D_H)`. Their regularity gives raw bounds
`val C[i,j]>=8-H_i` and `deg C[i,j]<=8+H_j`. Fixing the first `k` raw jets
at both cuts makes the correction satisfy

```
val delta[i,j] >= 8+k-H_i,       deg delta[i,j] <= 8+H_j-k.
```

Independent leading rows give the predictable degree of a row combination.
Enumerating all bounded coefficient shifts and taking an exact binary
kernel therefore computes the full correction space. These linear spaces
extend to every characteristic-two field.

The checker reproduces dimensions for every row height zero through eight
and every `k=0,...,6`. Across the sixteen actual row heights, fixing two
raw jets leaves dimension 284. Excluding order two at both cuts imposes
rank 140, leaving dimension 144. Fixing three raw jets gives dimension 138.
Exactly **six** quotient directions remain in the second jet.

## A bordered minor is fixed

The baseline has a unit pivot at `(0,0)` and a valuation-one pivot at
`(1,1)` at each cut. Fixing the first two raw jets preserves these and
excludes additional orders zero or one. For remaining row `i` and column
`j`, use the minor on rows `(0,1,i)` and columns `(0,1,j)`. Its coefficient
at degree three is linear in the second-jet change, since two such changes
start at degree four. The checker constructs this functional from literal
order-one cofactors and separately checks agreement with series
elimination. Vanishing of all these coefficients is exactly the absence
of order two.

Split the 144-dimensional kernel into the 138-dimensional higher-jet
subspace and six explicit second-jet directions. Literal expansion over
the polynomial ring in six independent variables gives

```
[z^4] det M_zero[{0,1,2},{0,1,2}] = 1.
```

Repeated variables remain powers; no identity such as `x^2=x` is imposed.
Every higher-jet basis direction has zero effect on this coefficient,
checked by a separate determinant calculation. That effect is linear:
two changes starting at degree three, or one at degree three and one at
degree two, cannot contribute at degree four. Hence the identity holds
for arbitrary coefficients over every characteristic-two field.

The first two orders are zero and one; no order two makes the third at
least three. This minor of valuation four bounds the sum of the first
three orders by four. The third order is therefore three.

Run the exact replay with `make verify-joint-small-rational`. It checks
the inverse, source containment, bounded bases, both bordered-minor systems,
the quotient directions and the symbolic degree-four coefficient. No
finite-field search is used to infer the symbolic identity.

The [clean full verification](../evidence/joint-moving-full-verification.json)
passed at exact `f2a16c0`: defaults and three named joint targets, 344 Python
command blocks, 66 unit tests, five crossing tests and 979 Markdown files.
This rational-module result postdates that source and has its own focused
check. There is no new Lean work or independent mathematical review claim;
the factoring score remains `1/2`.
