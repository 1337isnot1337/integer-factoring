# Balanced-charge obstruction for the full-alphabet quartic

## Result and scope

The full four-by-four physical alphabet does not remove the balanced
asymptotic obstruction for the explicit partial quartic repair.
There are uniform algebraic entry certificates for both conjugate
parameters `gamma=omega,omega^2` in `F4`.

**Theorem.** Fix `c>1/sqrt(2)`. For either conjugate parameter, all
sufficiently large power-of-two intervals have strictly positive
critical pressure for every real 28-component core/charge gauge vector
`theta` and every positive `ell=log s` satisfying

\[
\theta_0+\theta_1=\theta_2+\theta_3=0,
\qquad M\ell\ge c\log M.                                  \tag{1}
\]

This is the full 29-parameter family, including `ell`, with the
[explicit partial correction](full-alphabet-quartic-repair.md) and
original quartic numerator. It allows arbitrary changes of all four
extreme-state weight functions as `M` varies.

The theorem does not cover unbalanced aggregate bulk charges, scales
below (1), other endpoint repairs or perturbations, or isolated small
interval lengths. It does not improve the factoring exponent.

## 1. Uniform entries for the full alphabet

The [ten-mode geometric concentration](rectangular-full-core-factorization.md)
depends on the positions of the extreme CAR operators, so it applies
without change when the rectangular projections are removed. The
resulting active matrix now has full rank `1024` at each place.
There are three nonempty-state weights on each of `v,y,u,x`, giving
12 per place and four bulk charges.

Use the formal variable `Z=z^M` and the closed sequence calculation
from the [uniform coefficient proof](rectangular-uniform-ray-obstruction.md).
All denominators retain their lowest `Z` order when `z=omega`.
The extra denominators in the full partial correction are powers of
`1+gamma`; these are nonzero constants for the two stated parameters.
Thus the same specialization inequality proves uniform entry-slope
upper bounds.

The certificate uses 34 selected regular vacuum row/column maps. At
`gamma=omega`, all have a nonzero specialized coefficient at the stated
slope. At `gamma=omega^2`, specialize `z=omega^2`: every coefficient
is the Frobenius conjugate of the first calculation. The checker
verifies all 68 coefficients directly.

The [input certificate](../evidence/full-alphabet-leading-inputs.json)
records the four physical core masks, ordinary vacuum coordinates,
28-component weight features, slopes, and leading coefficients.
The features are independently reconstructed from the homogeneous
CAR operators. The selected maps need only be regular intertwiners;
no uniform choice of a full vacuum basis is assumed.

## 2. Sixteen physical minima and six gaps

At balanced charges, the physical energy slopes are

\[
\begin{aligned}
E_X&=(0,-2+v_l+x_l,\ 2+v_t+x_t,\ v_b+x_b),\\
E_Y&=(0,-2+y_l+u_l,\ 2+y_t+u_t,\ y_b+u_b),
\end{aligned}                                             \tag{2}
\]

where every weight is summed over both places. If `gamma_e` denotes
the actual minimum weighted entry slope at place `e`, the leading
critical-pressure coefficient is

\[
C=\min E_X+\min E_Y-\gamma_\infty-\gamma_0.                \tag{3}
\]

The 34 certified entries support exact rational duals for all
16 choices of physical minima. Their lower-bound table is

| `X` minimum | `Y` empty | `Y` low | `Y` top | `Y` both |
|---|---:|---:|---:|---:|
| Empty | 2 | 2 | 4 | 4 |
| Low | 2 | 0 | 4 | 2 |
| Top | 2 | 2 | 8 | 4 |
| Both | 2 | 2 | 4 | 4 |

Only the all-low minimum can have `C=0`. Six further exact duals
bound every other physical state on that face:

| Core | Non-low state | Gap lower bound | Off-face `C` multiplier |
|---|---|---:|---:|
| `X` | Empty | 2 | -1 |
| `X` | Top | 4 | -2 |
| `X` | Both | 2 | -1 |
| `Y` | Empty | 2 | -1 |
| `Y` | Top | 4 | -1 |
| `Y` | Both | 2 | -1 |

The last column means that the dual remains valid when `C` is not
zero: its lower bound is the displayed constant plus that multiplier
times `C`, provided the selected minima are low. These are exact dual
identities, with no numerical tolerance.

The checker also verifies primal points for the 34-constraint model.
Their attainment is a statement about that smaller finite LP. The
uniform conclusions use only the dual lower bounds and the actual
nonzero coefficient witnesses.

## 3. The all-low background is preserved

The full partial correction has diagonal entries one except at
rows `10,11,14,15`, and its five off-diagonal entries occur only in
those rows. Its row `0` is exactly the row of the identity. Under the
physical state convention, row `0` is the all-low `v,y` output.

After fixing all four physical extreme states to low, this output
row leaves the Gaussian unchanged. The quartic contains the top
creations `y_T,v_T`, so its coefficient vanishes under the same
restriction. Hence the conditional tensor is exactly the original
Gaussian with occupied low and empty top on all four legs.

The [all-low background calculation](rectangular-low-background.md)
therefore applies: the conditional spectra are those of the interval
of length `L=M-2`, with the scalar and frozen physical energies
cancelling in pressure. All four core weight functions become scalar
factors and cancel as well.

## 4. Uniformity over changing gauge vectors

The proof of the
[balanced-gauge theorem](rectangular-balanced-gauge-obstruction.md)
now applies with two changes of constants.

Choose a common `B>=0` so that the 34 unweighted coefficients satisfy
`v(entry)<=aM+B` for every sufficiently large power-of-two `M`.
Put `tau=M ell`; let `lambda_e` be the smallest weighted active order
in logarithmic units; remove the `O(ell)` endpoint offsets from the
physical energies; and set

\[
T=m_X^0+m_Y^0-\lambda_\infty-\lambda_0+2B\ell.             \tag{4}
\]

The dual table gives `T>=0`, and gives `T>=2tau` whenever a non-low
state minimizes one of the offset-free physical energies. The
Gaussian comparison still has 20 active one-particle factors and
negative-order sum `-4M`. Both physical cores now have four states,
so the coarse pressure bound becomes

\[
F\ge T-(2B+8)\ell-K_0,
\qquad K_0=20\log2+\sqrt2\log16.                          \tag{5}
\]

If `F<=0`, (1) and (5) eventually force the selected minima to be
low. Every off-face multiplier in the gap table has absolute value
at most two. The actual gaps are consequently at least

\[
2\tau-2K_0-(4B+20)\ell.                                  \tag{6}
\]

There are six non-low states in total. The preserved background and
the elementary physical partition estimate give

\[
F\ge F_{G,M-2}
 -6\sqrt2\,e^{2K_0/\sqrt2}
   \exp\!\left(-\frac{(2M-4B-20)\ell}{\sqrt2}\right).
\tag{7}
\]

The same all-scale Gaussian estimate gives
`F_{G,M-2}>=c_0 min(ell,exp(-ell))` uniformly under (1), for all
large `M`. The error in (7) is smaller than this lower bound,
uniformly in every permitted scale, because `sqrt(2)c>1`.
This contradicts `F<=0` and proves the theorem.

## 5. Reproduction

Run

```sh
python tools/check_full_alphabet_leading_obstruction.py \
  --output evidence/full-alphabet-leading-obstruction.json
```

The [checker](../tools/check_full_alphabet_leading_obstruction.py)
recomputes the 34 formal coefficients for each conjugate parameter,
checks every physical weight feature, and verifies all 22 exact
primal/dual certificates. It uses standard-library rational and
finite-field arithmetic. No fit in the interval length, numerical
optimizer, or native Smith calculation is used by the checker.

The [evidence](../evidence/full-alphabet-leading-obstruction.json)
records all coefficient values and precision boundaries, the dual
bounds, the constants needed for changing gauges, and source hashes.
The exploratory search that found the witnesses is separate from
this finite exact verification.
