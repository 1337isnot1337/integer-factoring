# Jordan complementary maps: a triangular all-power barrier

## Result

Let \(J_{ij}=1\) when \(j=i+1\) on labels \(0,\ldots,6\), with no wraparound edge.
Let \(\mathcal V_J\) be the full formal coefficient space of the cubic output maps \(C=\lambda J\).
All scalar output parameters and all recombinations of their coefficient directions are allowed.

For every field of characteristic two and every positive integer \(n\), any decomposition of \(E^{\otimes n}\)
into tensors from \(\mathcal V_J^{\otimes n}\) has rank cost at least
\[
 \boxed{\quad S_n\geq392^{n/2}>7^{\sqrt2 n}.\quad}
\]
The proof uses a 14-dimensional constant row compression and a 28-dimensional column restriction with a common triangular form.
The result strengthens the Jordan case of [the finite mixed-space search](cubic-mixed-F8-spaces.md).

A separate exact support calculation shows that a nonnegative diagonal gauge on all four tensor legs cannot give a nontrivial degeneration of the entire space.
The successful column filtration uses a restricted matrix space.

## Formal space and field scope

Use the cubic construction of [the complementary-output identity](cubic-complementary-output-maps.md), with the same fixed cyclic output selector.
Retain the five raw coefficients
\[
 (r^3,\;pr^2,\;p^2r,\;qr^2,\;q^3).
\]
The combined coefficient spaces at \(C=J\) and \(C=zJ\), where \(z^3+z+1=0\), have dimension eight.
Selecting indices \((0,1,2,3,4,6,7,8)\) from their concatenated raw lists gives a basis \(G_0,\ldots,G_7\) with \(G_0=E\).
Write \(e(\sum_iw_iG_i)=w_0\).

The checker reconstructs the formal parameter coefficients from \(\lambda=0,1,z\).
The linear parameter coefficients vanish; the constant and quadratic coefficient matrices have entries in F2.
Every formal direction belongs to the selected space, in both flattenings with the same coordinates.
Thus the eight-dimensional formal space descends to F2 and contains every \(\lambda\) over any characteristic-two field.
One can extend the field to contain F8 for the coordinate calculation without changing matrix rank.

## Restrictions

Let \(\mathcal P\) be the fourteen unordered nonedges of the seven-cycle:
\[
 \mathcal P=\{(a,b):0\leq a<b\leq6,\ b-a\notin\{1,6\}\}.
\]
Let \(\mathcal O\) contain both orientations of every pair in \(\mathcal P\), so \(|\mathcal O|=28\).
The row flattening uses ordered indices \((x,y)\) and \((v,u)\), while the column flattening uses \((x,u)\) and \((v,y)\).

For every \(T\in\mathcal V_J\), the principal row restriction on \(\mathcal P\) is exactly
\[
 R(T)\big|_{\mathcal P\times\mathcal P}=e(T)I_{14}.
\]
This is an identity of coefficient matrices, with no specialization of the parameters.

On the column restriction \(\mathcal O\times\mathcal O\), every diagonal entry is \(e(T)\).
Every supported off-diagonal entry goes from an ordered pair \((a,b)\) to a pair \((c,d)\) satisfying
\[
 c+d=a+b+2.
\]
There are exactly twenty such supported positions in the complete eight-direction space.
Order \(\mathcal O\) by increasing label sum, with a fixed tie-breaking order. The restriction is then upper triangular for every \(T\), with all 28 diagonal entries equal to \(e(T)\).

This also gives a constant weight degeneration of the restricted column space:
conjugate by the diagonal powers with exponent minus the label sum.
Every off-diagonal entry acquires positive exponent two, while the diagonal has exponent zero.
The limit is \(e(T)I_{28}\).

## Tensor powers

Tensoring the row compression gives
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}R(T)\geq14^n.
\]
Tensor products of the restricted column matrices are upper triangular, and all diagonal entries equal \(e^{\otimes n}(T)\).
These identities persist under arbitrary linear recombinations in \(\mathcal V_J^{\otimes n}\).
Therefore
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}C(T)\geq28^n.
\]

Since \(e(E)=1\), at least one piece of any decomposition of \(E^{\otimes n}\) has nonzero \(e^{\otimes n}\).
That piece alone has cost at least \(\sqrt{14^n28^n}\).
The exact comparisons \(392>7^3\) and \(3^2>8\) give the displayed strict barrier.

## Why the full four-leg weight degeneration fails

Consider arbitrary diagonal powers on the four tensor legs that preserve the crossing entrywise.
After removing irrelevant constants, every tensor entry has weight
\[
 w(x,y,v,u)=\alpha_x+\beta_y-\alpha_v-\beta_u.
\]
A well-defined direct degeneration of the whole coefficient space requires \(w\geq0\) on every supported entry.
The full support has 121 entries.

For each \(j=1,\ldots,6\), both \((0,j,j,0)\) and \((j,0,0,j)\) occur in the support.
Their weights are negatives of one another, so
\[
 \alpha_j-\beta_j=\alpha_0-\beta_0.
\]
Write \(\alpha_i=\beta_i+c\). The remaining entry weights become
\(\beta_x+\beta_y-\beta_v-\beta_u\).

For \(j=1,\ldots,5\), the supported entries \((0,1,j+1,j)\) and \((j,j+1,1,0)\) again have opposite weights.
They force
\[
 \beta_j+\beta_{j+1}=\beta_0+\beta_1.
\]
Hence \(\beta_i\) alternates between one value on even labels and another on odd labels.
An exact check of all 121 supported entries shows that every such remaining weight vanishes on every entry.
The support-weight equations have rank eleven in the fourteen leg variables, leaving precisely these three gauge parameters.

Thus every feasible nonnegative full-leg weighting acts trivially on the support.
This calculation is specific to the actual cyclic output selector; no arbitrary-label symmetry is assumed.
It does not prevent a degeneration after taking a matrix restriction, as the successful column proof shows.

## Reproduction

The [checker](../tools/check_cubic_jordan_filtration.py) verifies 1,568 row-compression coefficients,
6,272 column-triangular coefficients, all twenty strict upper edges, 72,030 formal parameter entries,
and the complete failed full-leg feasibility certificate.
The [evidence](../evidence/cubic-jordan-filtration.json) records the selected coordinates, ordering, support witnesses, and source hashes.

```bash
python tools/check_cubic_jordan_filtration.py --output /tmp/cubic-jordan-filtration.json
```

The theorem concerns the full scalar Jordan-output space and its tensor powers.
General complementary matrices and coefficient directions from outside this space remain outside its scope.
