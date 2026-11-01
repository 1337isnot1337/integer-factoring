# The two-product theorem excludes rank-sum seeds through 54 labels

## Result

For every field and every exact crossing decomposition on `d` labels,

```math
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}\qquad(2\le d\le54).
```

Consequently a seed satisfying the strict sufficient finite rank-sum
criterion must have at least **55 labels**. This does not assert that
such a seed exists. The factoring exponent remains one half.

The [positive-dual certificate](weighted-cut-dual-certificate.md)
already proves the assertion through 52. The new
[two-product complement theorem](two-product-complement.md)
strengthens the elementary residual certificate enough to cover both
remaining dimensions. No new numerical optimizer or stored dual is
needed for this step.

## Exact strengthened residual bound

Use the actual extremal-piece classification from the
[residual-budget proof](weighted-cut-residual-budget.md). After a
global cut exchange, the selected piece has ranks `(a,b)` with

```math
a_0=(k-1)^2b+1\le a\le a_1=\min(k^2b,d^2).
```

The remaining row and column rank totals are at least `R_0,C_0`.
All preceding formulas remain valid. When `b=2`, strengthen them to

```math
R'_0=\max(R_0,d(d-1)/2),\qquad C'_0=C_0,
```

and strengthen the total-rank surplus to

```math
e'=\max(e,a_0+d(d-1)/2-d^2).
```

These follow directly by applying the two-product theorem to the
selected tensor and rank subadditivity to its complement. The
surplus inequality adds the selected piece's rank `a>=a_0` and the
ordinary column total `d^2`. When `b!=2`, retain the earlier bounds,
including the sharper rank-one/Fitting formula.

For the same positive interpolant weights `w_s` and endpoint-chord
coefficients `mu_s,nu_s` used in the residual proof, put

```math
J_2(d,k,b)=\frac{2d^2+e'}{2k-1}
+\sum_{s=1}^{k-1} w_s\Bigl(
\max\{sd,\min(a_0,s^2b)+\mu_sR'_0+\nu_sC'_0\}
+\max\{sd,b+\nu_sR'_0+\mu_sC'_0\}\Bigr).
```

Every replacement is a stronger valid lower bound multiplied by a
nonnegative coefficient. The identical interpolant argument therefore
proves `sum sqrt(r_i c_i)>=J_2(d,k,b)` in each class.

For `d=53,54`, the checker enumerates **every**
`2<=k<=d` and `1<=b<=floor((d^2-1)/(k-1)^2)`. It also checks the
equal-rank class `k=1`, where the lower bound is `d^2`. All
calculations use integers and fractions. For each resulting rational
bound `L`, it verifies

```math
L^{408}>d^{577},\qquad577^2-2\cdot408^2=1.
```

Thus `L>d^(577/408)>d^sqrt(2)`. This gives the asserted strict
exclusion in both new dimensions, together with the earlier result.

At `d=55`, the scalar class `(k,b)=(51,1)` has a bound below the
target. The checker verifies that fact using the lower Pell
approximation `1393/985<sqrt(2)`. This is a limit of the present
certificate; it supplies neither a tensor decomposition nor an upper
bound on its possible quality.

## Replay

```sh
python3 tools/check_two_product_seed_cutoff.py \
  --output evidence/two-product-seed-cutoff.json
```

The [certificate](../evidence/two-product-seed-cutoff.json) retains
every per-`k` minimum, the exact worst rational values, source hashes,
and the next-dimension control. The all-field proof uses the ordinary
theorems cited above; no Lean formalization is added. Substantial
OpenAI assistant support was used.
