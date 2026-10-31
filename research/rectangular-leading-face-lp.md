# Exact linear certificates for the fitted leading core model

## Result and scope

The natural entry valuations of the seven/nine-mode rectangular core
were computed at \(M=16,32,64\). Every nonzero entry fits an affine
function of \(M\), using 16 and 32 and checking 64. This gives a finite
candidate model for the leading pressure coefficient
\(C\) in the [bivariate limit](bivariate-core-pressure-limit.md).

The resulting linear programs have exact rational primal and dual
certificates, with no bound on the gauge variables. Their only
zero-cost physical minimum is the all-low core. On that face, every
other physical core state has a strictly positive energy gap.

These are **exact conclusions about the fitted finite model**. The
continuation of its entry slopes for all large \(M\) has not yet been
certified symbolically. The result also concerns the smaller
19-parameter family at zero aggregate bulk charges. Section 5 records
the separate full 25-parameter model, with the same certified LP values
but only two sizes of entry data. Neither result covers other repairs
or all rectangular crossings.

## 1. The entry model

Write \(\theta=A\tau+B\), where \(\tau=M\log s\). There are 18
coordinates in \(A\): four row/column bulk charges, seven infinity
core weights, and seven zero core weights. The profile scale is
separate. The core weights are ordered

```text
v_low, v_top, v_both, y_low, y_top, u_low, u_top
```

at each place. Impose zero aggregate bulk charges:

~~~math
A_0+A_1=0,\qquad A_2+A_3=0.                           \tag{1}
~~~

For a nonzero entry, let \(a_{ij}\) be its fitted \(M\)-slope,
and let \(f_{ij}\) be its integer row-plus-column feature vector.
The smallest slope at place \(e\) is

~~~math
\gamma_e=\min_{i,j}(a_{ij}+f_{ij}\cdot A).             \tag{2}
~~~

The 1,224 infinity entries and 17,739 zero entries have unchanged
support and matching affine valuations at all three sizes. Entries
with the same feature vector give redundant inequalities in (2).
Keeping the smallest slope for each feature leaves 217 infinity and
331 zero constraints.

For the physical cores, let \(v_l,v_t,v_b\) be the sums of the two
places' \(v\)-weights, and let \(u_l,u_t\) be the sums of the
\(y\)- and \(u\)-weights at both places. Their actual energy slopes
are

~~~math
E_X=(0,-2+v_l,2+v_t,v_b),\qquad
E_Y=(0,-2+u_l,2+u_t).                                 \tag{3}
~~~

The Gaussian active slopes and Gaussian physical minimum energies
cancel in this normalization. The leading coefficient is therefore

~~~math
C=\min E_X+\min E_Y-\gamma_0-\gamma_\infty.            \tag{4}
~~~

There is no missing factor \(\sqrt2\) in (4): the factor multiplying
the physical logarithmic partition cancels the division of its
energies by that factor.

## 2. Twelve unbounded linear programs

Select one minimizing state on each physical core. Add inequalities
saying that its energy is no larger than every other state. Introduce
the two variables \(\gamma_e\), with one inequality

~~~math
\gamma_e-f_{ij}\cdot A\le a_{ij}                     \tag{5}
~~~

for each compressed entry constraint. Minimizing (4) is now a linear
program, with all coordinates of \(A\) and both \(\gamma_e\)
unrestricted in sign and size.

The exact optimum table is:

| Selected \(X\) minimum | \(Y\) empty | \(Y\) low | \(Y\) top |
|---|---:|---:|---:|
| Empty | 2 | 2 | 4 |
| Low | 2 | 0 | 4 |
| Top | 4 | 4 | 8 |
| Both | 4 | 2 | 4 |

The entries are attained primal values and certified dual lower
bounds. Ties are allowed in the selected-minimum constraints, so the
table also rules out any non-low state tying for the minimum on a
\(C=0\) point of this model.

## 3. Quantitative gaps on the zero-cost face

Impose the all-low minimum constraints and the equality \(C=0\).
Five more linear programs minimize each non-low state's energy gap.
Their exact optima are:

| Core | State | Minimum energy-slope gap above low |
|---|---|---:|
| \(X\) | Empty | 2 |
| \(X\) | Top | 4 |
| \(X\) | Both | 2 |
| \(Y\) | Empty | 2 |
| \(Y\) | Top | 4 |

If the entry model continues to all large sizes, every such path has
physical non-low probabilities tending to zero at least as a fixed
power of \(M\), for \(\tau=c\log M\). Fixed finite offsets \(B\)
cannot close a positive slope gap.

This does not itself prove a nonnegative constant pressure term
\(D\). The cut profiles can retain contributions that the physical
type has suppressed. The exact preserved all-low Gaussian submatrix
is a possible source of a further lower bound; its limiting profile
must be analyzed separately.

## 4. Certificates and the remaining symbolic task

The [certificate data](../evidence/rectangular-leading-face-lp.json)
contain the compressed constraints, one original entry realizing each
constraint, all primal points, and sparse rational dual multipliers.
The [checker](../tools/check_rectangular_leading_face_lp.py) verifies
primal feasibility, the dual signs, every dual coefficient identity,
and equality of the primal and dual objectives with `Fraction`
arithmetic. It uses no numerical optimizer.

```sh
python tools/check_rectangular_leading_face_lp.py \
  --output evidence/rectangular-leading-face-lp-check.json
```

Across all 17 lower-bound certificates, only 28 original entry
constraints have nonzero dual multipliers. Their coordinates and
fitted affine forms are included separately. To promote the displayed
lower bounds to all sufficiently large power-of-two intervals, it is
enough to prove the corresponding upper bounds on these 28 entry
slopes. Proving their complete exact valuations is stronger than
necessary. The primal attainability claims would additionally require
the other constraints, but the exclusion conclusions use only the
dual lower bounds.

## 5. Independent weights on all four cores

The full ten-mode construction adds three \(x\)-weights at each
place. Its 73,521 nonzero infinity entries and 57,013 nonzero zero
entries have matching support and integer affine fits between
\(M=16\) and \(M=32\). These compress to 2,718 distinct entry
constraints. A third size has not yet been used to validate this fit.

The same unbounded LP construction, now with 24 gauge-slope variables,
gives exactly the table in Section 2 and the five gaps in Section 3.
The [full-model certificates](../evidence/rectangular-full-leading-face-lp.json)
have separate exact primal and dual witnesses. They also use only 28
entry constraints in their lower-bound certificates.

```sh
python tools/check_rectangular_leading_face_lp.py \
  --input evidence/rectangular-full-leading-face-lp.json \
  --output evidence/rectangular-full-leading-face-lp-check.json
```

Thus independent \(x\)-weights do not open a zero-cost non-low
minimum in this fitted model. Changing the repair's preserved all-low
row remains a distinct way to change the leading problem.

## 6. A conditional asymptotic consequence of the preserved background

Assume that the 28 entry-slope upper bounds used in one of these
models hold for all sufficiently large power-of-two \(M\). Fix a
ray \(\theta=A\tau+B\), with zero aggregate bulk charges,
\(\tau=c\log M\), and \(c>1/\sqrt2\).

The dual bounds imply \(C\ge0\). If \(C>0\), the critical
pressure is positive for all sufficiently large \(M\). If \(C=0\),
the five gap certificates force the physical low states to dominate,
and the sum of the two physical core logarithmic corrections is

~~~math
O(e^{-2\tau/\sqrt2})=O(M^{-\sqrt2 c}).                 \tag{6}
~~~

The [all-low background theorem](rectangular-low-background.md)
identifies the restricted tensor's conditional pressure exactly with
the Gaussian interval of length \(L=M-2\). Its fixed extreme-state
weights cancel as a scalar. A matrix's profile partition dominates
that of any row/column submatrix: the submatrix's sorted Smith orders
are componentwise no smaller, and it has no more nonzero orders.
Consequently

~~~math
F_{\mathrm{full}}\ge F_{G,L}-O(M^{-\sqrt2 c}).         \tag{7}
~~~

The interval expansion gives \(F_{G,L}\ge\ell/12\) eventually,
where \(\ell=\tau/M\). Its positive finite tail dominates its
negative tail, and its infinite-product term is \(\ell/6\).
Because \(\sqrt2 c>1\), the error in (7) is \(o(\ell)\).
Thus the full critical pressure is eventually strictly positive on
every such fixed ray.

This conditional consequence excludes a negative limiting constant
and smaller asymptotic gains on the specified rays. It does not rule
out an isolated finite-size witness, rays with changing coefficients,
nonzero aggregate bulk charges, or changes to the seed that remove
the all-low Gaussian background. The symbolic entry bounds remain
the missing hypothesis for an unconditional statement about these
families.
