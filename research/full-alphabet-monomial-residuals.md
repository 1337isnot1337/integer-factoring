# Full-alphabet crossings restored by critical cubic pairs

## Result and scope

Four pairs of extreme cubic monomials restore the full sixteen-state
ramified residual at canonical quartic coefficient one. This removes
the three-state coding loss at that coefficient. The repair can be
chosen block lower triangular in physical height, or as the complete
residual inverse. Both choices give polynomial crossing tensors at
every independently matching type length.

The determinant identity is uniform over every field of characteristic
two. For pair coefficients `a,b` and canonical quartic coefficient
`gamma`, it is

```math
\det F_{\gamma,a,b}=(1+\gamma+ab)^4.
\tag{1}
```

At `gamma=1`, every nonzero product `ab` restores full rank. The
binary choice `a=b=1` is already sufficient. The resulting examples
have nonzero odd-degree coefficients, so they are not Gaussian compound
tensors. This is a new admissible family, not a certified pressure
improvement or a factoring exponent.

## 1. Critical coefficients change the residual

Use extreme modes `0,T`, put `c=1+w`, and let `G` be the original
Gaussian. For a monomial `Q`, write `Delta(Q)` for output index sum
minus input index sum, as in
[the vanishing-monomial construction](rectangular-core-monomials.md).
Under `w=1+t^T` and the physical index gauge, the coefficient
`P(w)Q/c^r` has ramified excess

```math
e(Q,P,r)=\Delta(Q)+
      \left(\mathrm{ord}_{w=1}P-r\right)T.
\tag{2}
```

Positive excess gives the previously studied endpoint-vanishing terms.
Zero excess is also regular, but contributes a monomial to the residual.
If `P=c^k R` with `R(1)` nonzero, that contribution is `R(1)Q`.
Negative excess is excluded by this construction.

Let `Q0=x0*yT*u0*vT`. The raw seed is the affine sum

```math
G(w)\left(1+\frac{\gamma Q_0}{c^2}
              +\frac{a Q_-}{c^{r_-}}
              +\frac{b Q_+}{c^{r_+}}\right).
\tag{3}
```

The four pairs below have particle-count changes `-1` and `+1`.
Their displayed poles make both terms critical.

| Pair | `Q_-` | `r_-` | `Q_+` | `r_+` |
| --- | --- | ---: | --- | ---: |
| A | `x0*u0*vT` | 1 | `yT*u0*vT` | 2 |
| B | `x0*u0*vT` | 1 | `xT*yT*vT` | 1 |
| C | `x0*yT*u0` | 1 | `yT*uT*vT` | 1 |
| D | `x0*yT*u0` | 1 | `x0*yT*vT` | 2 |

The coefficients can include independent Laurent powers of `w`, since
those powers equal one at the crossing point. More general numerators
are allowed by (2), with their first nonzero local coefficient used in
the residual.

## 2. An explicit finite residual and a polynomial identity

Order each extreme alphabet as `(low,empty,top,both)`, represented by
the two-bit masks `(1,0,2,3)`. The core height vector is
`h=(0,0,1,1)`. Rows of `F` are `v+4y`; columns are `x+4u`.

The full Gaussian ramified residual is the coefficient tensor of

```math
\prod_{a\in\{0,T\}}(1+x_av_a)(1+y_au_a)
          (1+u_Tv_0)(1+x_Ty_0).
\tag{4}
```

Multiply (4) by the affine factor
`1+gamma*Q0+a*Q_-+b*Q_+` in the characteristic-two exterior algebra.
Reading its coefficients in the stated order gives `F`. This defines
every entry without a numerical limiting operation. The bulk endpoint
remains the crossing for any finite mode set containing the extremes.

Sort the sixteen indices by the total core height
`h(i)=h(i mod 4)+h(floor(i/4))`. The cumulative height blocks have
orders `4,12,16`. For each of the four pairs, their determinants are

```math
\det F_{\le0}=1,\qquad
\det F_{\le1}=1,\qquad
\det F=1+\gamma^4+a^4b^4=(1+\gamma+ab)^4.
\tag{5}
```

The [symbolic checker](../tools/check_full_alphabet_monomial_residual.py)
expands these determinants over the polynomial ring
`F2[a,b,gamma]`. It represents each polynomial as its set of nonzero
monomials and uses subset determinant expansion. This is an exact
polynomial calculation, not interpolation from the four field elements.
The reported exponent lists are a directly replayable certificate of
(5).

Therefore every leading height block is invertible whenever
`1+gamma+ab` is nonzero. At `gamma=1`, condition (5) is simply
`ab` nonzero. The singular locus is explicitly rejected by the
constructor. Outside this four-pair subfamily, the generic residual
routine still accepts any list of nonnegative-excess monomials whose
required residual blocks are invertible.

## 3. Constructing the partial correction

There is a unique block lower triangular correction `C` such that

```math
CF=U,\qquad U_{ii}=1,\qquad
U_{ij}\ne0\ (i\ne j)\Longrightarrow h(i)<h(j).
\tag{6}
```

It has an explicit block formula. For each height `s`, let `I_s`
contain all core indices of height at most `s`. For a row `i` of
height exactly `s`, set

```math
C_{i,j}=
\begin{cases}
\left(F_{I_s,I_s}^{-1}\right)_{i,j},&j\in I_s,\\
0,&j\notin I_s.
\end{cases}
\tag{7}
```

The selected row of the product with `F` is the corresponding identity
row on all columns in `I_s`. Every remaining column has greater
height. This proves (6) and the lower block support of `C`.
Its diagonal blocks are invertible by (5). The complete inverse
`C=F^-1` is also valid and gives endpoint identity before type
restriction, although it generally has more nonzero entries.

The partial correction from the earlier full-alphabet quartic family
is recovered when the extra monomial list is empty. The implementation
checks that agreement exactly for all allowed nonbinary base parameters.

## 4. All-length polynomial crossing

On output core pairs define `H_ii=c^h(i)` and repair the raw tensor by

```math
\widetilde K=H^{-1}CH K_{\rm raw}.
\tag{8}
```

If the largest raw local denominator is `c^r`, the numerator of
`c^(r+2)*Ktilde` is polynomial after a fixed global power of `w`
clears Laurent terms. This follows because the height differences
in (8) lie between minus two and two. The constructor cancels any
common local factor exactly; the unreduced bound suffices for the proof.

After ramification and the actual occupied-index gauge, (8) is simply
the constant correction `C` applied to a regular tensor. Its endpoint
is `U` on the core and the crossing on the bulk. Independent matching
full-letter types on `(X,V)` and `(Y,U)` cancel the physical energy
gauges in every tensor-power entry. Equal total height excludes all
strict height decreases from (6), leaving exactly

```math
[x=v][y=u].
\tag{9}
```

Regularity after ramification removes the possible pole at `w=1`.
There are no other finite poles after the fixed `w` clearing. Thus
the typed tensor powers are polynomial, and their degrees grow at most
linearly with word length. This proves the ordinary all-length crossing
interface. The familiar type-counting argument then proves full generic
rank of either local square flattening, independently of truncated
series calculations.

The correction in this family can mix different particle counts, so the
earlier single-sign profile obstruction does not apply. For the four
displayed pairs, however, the height-zero residual block is exactly the
identity. The partial correction therefore fixes every height-zero row.
Each raw cubic contains an output top variable, and is absent from the
all-low output row. The all-low Gaussian submatrix consequently persists
for this partial repair, with the same interval of `M-2` remaining modes.
This is a submatrix lower bound; the existing baseline leading-pressure
certificate has not been extended to these new corrected tensors.
The complete inverse is not asserted to preserve that submatrix.

Field-unit particle scalings do conjugate the residual and its unique partial correction,
so `(a,b)` can be normalized to `(1,ab)` without changing weighted
minor valuations. A Laurent particle scaling does not commute with
the new correction. The two Laurent powers are therefore kept
independent in the present constant-correction search family.

## 5. Fixed active spaces and independent controls

Every new term and every correction entry is an operator on the same
extreme physical variables. Pairs A, C, and D use the existing seven
active modes at infinity and nine at zero. Pair B uses `xT` and the
ten-mode extension. The full alphabet retains every core coordinate
state: the long-interval matrix orders are `128,512` in the smaller
support, and `1024` in either place for the full support. The scalar
spectators remain unchanged.

`prepare_residual_repair` derives and checks `F,C,U` before the native
builder accepts a candidate. The coordinate compiler and local cache
accept the additional keyword `residual_monomials`. The ordinary
`extra_monomials` argument continues to require strictly positive
ramified excess. The gamma-one family can be represented by a gamma-zero
base plus the critical canonical quartic and its critical cubic pair.

The [original-tensor constructor](../tools/check_full_alphabet_monomial_residual.py)
uses Gaussian minors and separately checks the complete exterior
expansion. The [recorded constructor evidence](../evidence/full-alphabet-monomial-residual.json)
contains 27 cases: every rank-restoring pair, each nonzero `F4`
coefficient product, both repairs, and larger or nonconsecutive mode
sets with Laurent powers. It checks 79,872 exterior coefficients,
79,872 ramified endpoint entries, and 293,632 two-copy typed entries.
Every small repaired pair has a nonzero odd-degree coefficient witness.
All three tested singular parameter choices are rejected.

The [profile audit](../tools/check_full_alphabet_monomial_profiles.py)
compares the native program with original four-mode tensors, including
the binary pair, nonbinary pairs, the full-support case, the complete
inverse, Laurent powers, and a canceled local denominator.
Its finite entry and quarter-weight comparisons audit implementation;
Sections 1--4 supply the uniform proof.

```sh
python tools/check_full_alphabet_monomial_residual.py \
  --output evidence/full-alphabet-monomial-residual.json
python tools/check_full_alphabet_monomial_profiles.py \
  --output evidence/full-alphabet-monomial-profiles.json
python tools/screen_full_alphabet_monomial_pressure.py --modes 2 \
  --output /tmp/critical-cubic-screen.json
```

The small-mode screen optimizes both physical types and diagonal
weights at each listed candidate. It records every seed hash, coefficient
choice, Laurent powers, repair, start, and numerical outcome. These
search results do not establish a global optimum or an exact factoring
bound.
