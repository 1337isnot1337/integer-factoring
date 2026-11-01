# Adjugate compatibility closes the bare seven-letter interface

**Later extension:** [rectangular subalphabet entropy witnesses](subalphabet-bare-nine-letter-floor.md)
improve the seven-letter floor to `719/1000` and establish bare floors above
the half-rate threshold at eight and nine letters over every field and degree.
The proof and seven-letter certificate below remain valid.

## Result and scope

Every polynomial crossing on seven physical letters, over every field,
has bare local Smith pressure

```math
\boxed{\rho_{\rm bare}>\frac{709}{1000}
       >\frac1{\sqrt2}.}                         \tag{1}
```

There is no degree cutoff. The pressure conclusion also holds for
arbitrary fixed auxiliary dimension, when a constant auxiliary
contraction gives the crossing endpoint and the two auxiliary indices
stay together on one side of each cut.

Together with the existing
[shared-minor floor through six letters](shared-minor-smith-floor.md),
this closes the bare, ungauged local Smith-pressure interface through
seven letters. Independent endpoint gauges, adaptive channel or
cancellation constructions, and direct coefficient-rank methods remain
outside this result. It is not a lower bound on integer factoring.

The new algebraic ingredient is a quantitative extension of
[inverse-degree compatibility](inverse-degree-crossing-compatibility.md).
Using the polynomial adjugate removes the monomial-determinant
hypothesis. Applying the shared-minor bound to that adjugate supplies
an additional upper-tail constraint. A rational entropy witness then
proves (1) by a single exact integer comparison.

## 1. Profiles and their existing constraints

First take a scalar polynomial crossing `K(z)` with alphabet `q>=2`.
Remove its largest common scalar monomial and let the remaining actual
degree be `M`. Put `n=q^2`. Both cut matrices are invertible at one.

Let `mu` be the increasing Smith orders of the row-cut reversal
`w^M K_row(w^-1)`, and let `nu` be those of the column cut at zero.
When `M>0`, the degree and common-minor arguments give

```math
\mu_1=\nu_1=0,\qquad
\sum_i\mu_i,\sum_i\nu_i\le nM,\qquad
\sum_{i=1}^q(\mu_i+\nu_i)\le qM.                \tag{2}
```

At scale `s>1`, the restored pressure is

```math
B(s)=s^M\left(\sum_i s^{-\mu_i}\right)
             \left(\sum_i s^{-\nu_i}\right).     \tag{3}
```

For `M=0`, the pressure is `q^4` and the required conclusion is
immediate.

## 2. Adjugation preserves a rank-one reshuffle

View the column cut as an operator on `Mat_q`, with the fixed column
permutation chosen so that its endpoint is the identity. Its reshuffling
back to the row cut will be denoted by `Phi`. Thus

```math
\Phi(T)_{xyuv}=T(E_{vy})_{xu},\qquad \Phi(I)=E.
```

A rank-one row reshuffling means precisely that the operator has the
form

```math
R(X)=A X^T B.                                    \tag{4}
```

The adjugate of such an operator again has rank-one row reshuffling.
Indeed, as a polynomial identity over every field,

```math
\mathrm{adj}(R)(Y)
=\varepsilon(\det A\det B)^{q-1}
  \mathrm{adj}(B)^T Y^T\mathrm{adj}(A)^T,
                                                        \tag{5}
```

where `epsilon` is the sign of transposition on `Mat_q`. One may
verify (5) for invertible matrices and then use the polynomial identity;
the formula remains valid when either matrix is singular. In particular,
no inverse at the local parameter zero is required.

Now put `T(w)=w^M K_col(w^-1)`, with the fixed permutation understood,
and set `k=mu_2`. If `k>0`, a unit pivot of `Phi(T)` produces an exact
rank-one matrix lift modulo `w^k`. Its associated operator `R` has
the form (4), and

```math
T=R\pmod{w^k}.
```

Adjugation is polynomial in matrix entries. Formula (5) therefore gives

```math
\Phi(\mathrm{adj}T)
 =\text{a rank-at-most-one matrix}\pmod{w^k}.     \tag{6}
```

Every Smith order of this reshuffling except possibly the first is
at least `k`. If `k=0`, the same lower bound is automatic. This
argument allows singular leading column matrices and nonmonomial
determinants.

## 3. The adjugate degree and its opposite profile are exact

Write

```math
c=\sum_{i=1}^{n-1}\nu_i,\qquad
d=(n-1)M-c.                                     \tag{7}
```

The actual degree of `adj T` is exactly `d`. Each of its entries is
a reversed `(n-1)`-minor of the column matrix; the minimum valuation
of those minors is `c`. Thus the largest reversed degree is (7).

The polynomial crossing

```math
A(w)=\Phi(\mathrm{adj}T(w))
```

still has endpoint `E`. Let `eta` be its row-zero profile and `theta`
its column-infinity profile, reversed at the actual degree `d`.
The preceding section gives

```math
\eta_1\ge0,\qquad \eta_i\ge\mu_2\quad(i\ge2).
                                                        \tag{8}
```

The opposite profile can be computed exactly. Up to constant
permutations, reversing the adjugate gives

```math
z^d\mathrm{adj}T(z^{-1})
=z^{-c}\mathrm{adj}K_{\rm col}(z).
```

The Smith orders of an adjugate are the determinant order minus
the original orders, in reverse order. Consequently

```math
\boxed{\theta_i=\nu_n-\nu_{n+1-i}.}              \tag{9}
```

Neither (7) nor (9) assumes that the determinant is a monomial.

## 4. Two additional all-degree constraints

The determinant of the row cut of `A` is a nonzero polynomial of
degree at most `nd`. Its order at zero is at least `(n-1)mu_2`
by (8). Hence

```math
\boxed{(n-1)\mu_2+n\sum_{i=1}^{n-1}\nu_i
       \le n(n-1)M.}                            \tag{10}
```

Next apply the ordinary shared `q`-minor theorem to `A` itself.
It is valid even when `A` has a common scalar monomial: such a
factor adds opposite charges to the two endpoint profiles. Equations
(8)--(9) give

```math
(q-1)\mu_2+q\nu_n-\sum_{i=n-q+1}^n\nu_i\le qd.
```

Substituting (7) yields

```math
\boxed{(q-1)\mu_2+q\sum_{i=1}^n\nu_i
       -\sum_{i=n-q+1}^n\nu_i\le q(n-1)M.}      \tag{11}
```

Both inequalities also hold with `mu` and `nu` exchanged. To see
the symmetry directly, exchange the middle physical legs and replace
`K(z)` by `z^M K(z^-1)`. This preserves the swap endpoint and
interchanges the two profiles.

These are scalar profile constraints. Their pressure consequence
for auxiliary cuts is transferred separately in Section 6.

## 5. A rational entropy witness at seven letters

Set `q=7`, `n=49`, and define

```math
z_i=\frac{\mu_i+\nu_i}{2M},\qquad t=M\log s.
```

Cauchy--Schwarz in (3) gives

```math
\sqrt{B(s)}\ge\sum_{i=1}^{49}\exp(t/2-tz_i).     \tag{12}
```

Average the two orientations of (10)--(11), and use (2). The three
constraints needed below are

```math
\begin{aligned}
\sum_{i=1}^7z_i&\le7/2,\\
48z_2+49\sum_{i=1}^{48}z_i&\le2352,\\
6z_2+7\sum_{i=1}^{49}z_i-\sum_{i=43}^{49}z_i
 &\le336.
\end{aligned}                                    \tag{13}
```

Take their respective nonnegative multipliers

```math
\frac{34}{875},\qquad\frac1{28000},\qquad\frac1{1200}.
                                                        \tag{14}
```

The weighted right sides sum exactly to `1/2`. Since `z_1=0`,
complete the weighted left coefficients to a probability vector by
assigning the remaining mass to its first coordinate. The resulting
probabilities, in order, are

```math
p=\frac1{84000}
 (33913,\ 4465,\ \underbrace{3901,\ldots,3901}_{5},
 \underbrace{637,\ldots,637}_{35},
 \underbrace{567,\ldots,567}_{6},\ 420).         \tag{15}
```

Thus `sum_i p_i z_i<=1/2`. Weighted arithmetic--geometric mean
applied to (12) gives

```math
\sqrt{B(s)}
\ge\prod_i\left(\frac{\exp(t/2-tz_i)}{p_i}\right)^{p_i}
\ge\exp(H(p)).                                  \tag{16}
```

The entropy comparison needed for (1) is a single integer inequality:

```math
84000^{84000}>
7^{119112}\,
33913^{33913}\,
4465^{4465}\,
3901^{5\cdot3901}\,
637^{35\cdot637}\,
567^{6\cdot567}\,
420^{420}.                                      \tag{17}
```

The checker verifies (17) exactly. Since
`119112/84000=709/500`, it says
`H(p)>(709/500)log 7`. Equations (12)--(16) therefore prove

```math
B(s)>7^{709/250},\qquad
\rho_{\rm bare}>709/1000.
```

Finally `2*709^2>1000^2`, proving the strict half-rate comparison.
No numerical minimizer, active-set assumption, or logarithm estimate
is needed for this certificate.

## 6. Transfer to fixed auxiliary bare profiles

Suppose the polynomial crossing has fixed auxiliary indices, and
a constant contraction gives a scalar crossing `K`. In each cut,
the scalar matrix is a constant left or right projection of the
rectangular auxiliary matrix.

Such a projection cannot increase Smith pressure. For a regular
matrix with orders `a_i`, its truncated multiplication matrix at
precision `P` has rank `sum_i(P-a_i)_+`. Regular projections cannot
increase these ranks. Summing them with the positive geometric weights
gives

```math
\Phi_s(M)=(1-s^{-1})^2
 \sum_{P\ge1}s^{-(P-1)}\mathrm{rank}(M\bmod\tau^P).
                                                        \tag{18}
```

Here the rank is over the coefficient field of the truncated
multiplication, or Toeplitz, operator. A common shift extends the
comparison to Laurent matrices. Therefore
each scalar endpoint pressure is at most its auxiliary counterpart.
If contraction lowers the actual polynomial degree, using the larger
degree merely adds a common regularizing shift which is restored in
the infinity pressure. Removing a scalar common monomial likewise
adds opposite factors to the two pressures.

The scalar product floor (1) thus lower-bounds the auxiliary bare
product as well. This does not assert (10)--(11) directly for
rectangular profiles, and it does not exclude adaptive cancellation
or state-dependent word factorizations.

## 7. Exact replay and the remaining boundary

Run:

```sh
python3 tools/check_adjugate_seven_letter_floor.py
```

The [checker](../tools/check_adjugate_seven_letter_floor.py) independently
constructs polynomial adjugates in small controls, checks the transpose
channel identity, verifies both profile directions and the exact degree
formula, and checks (9)--(11). A control with a singular leading column
matrix and a nonmonomial determinant exercises the scope that the
earlier inverse theorem did not cover. It also reconstructs (15) from
(13)--(14) and verifies the integer entropy comparison (17).

Evidence is saved in
[adjugate-seven-letter-floor.json](../evidence/adjugate-seven-letter-floor.json).
Eight letters remain outside the present all-degree bare exclusion.
The new constraints apply there too, but the displayed entropy witness
is specific to seven letters.

## 8. A stronger universal second-jet constraint

The adjugate retains more information than just the second order.
With the notation above, it satisfies

```math
\boxed{\eta_i\ge\min(\mu_i,2\mu_2)\quad(1\le i\le n).}
                                                        \tag{19}
```

This statement holds over every field, including singular leading
channel factors. It requires no determinant or pole condition.
The seven-letter certificate used only the weaker consequence (8).

To prove (19), choose the unit pivot used in Section 2. Its exact
rank-one lift `R` can be chosen so that `Phi(E)`, for `E=T-R`, has
zero pivot row and column. Its finite Smith orders are then exactly
`mu_2,...,mu_n`, and every entry is divisible by `w^k`, where
`k=mu_2`.

Write `R(X)=A X^T B`, put `D=det A det B`, and define the regular
transpose channel

```math
S(X)=\mathrm{adj}(B)^T X^T\mathrm{adj}(A)^T.
```

The constant and first-order terms of the polynomial adjugate are

```math
\mathrm{adj}(R+E)
=\varepsilon D^{q-1}S
 +\varepsilon D^{q-2}\bigl(\mathrm{tr}(SE)S-SES\bigr)
 +O(E^2).                                        \tag{20}
```

This follows from the derivative of the adjugate on the invertible
locus and extends as a polynomial identity. All factors in (20)
are regular because `q>=2`. Each term in the remainder contains
at least two entries of `E`, so it is divisible by `w^(2k)`.

After reshuffling, the two multiples of `S` form a single regular
rank-one matrix. The reshuffling of `SES` is a regular left and right
image of `Phi(E)^T`. Such regular maps cannot lower its individual
Smith orders. Adding one rank-one matrix can shift the increasing
Smith list down by at most one position. It follows that, modulo
`w^(2k)`, the `i`th order is at least `mu_i`, which proves (19).

Here is an elementary justification of the interlacing assertion.
Modulo `w^h`, the image of `X+uv^T` is contained in
`im(X)+O*u`. Adjoining `u` adds at most one cyclic generator. A
submodule of a module over a truncated discrete valuation ring
cannot require more cyclic generators than the ambient module,
as is seen from their socles. The number of cyclic generators in
the image of a matrix modulo `w^h` is precisely the number of its
Smith orders below `h`. Thus a rank-one addition changes that count
by at most one in the required direction. The same argument shows
that regular left or right maps cannot lower individual orders.
Use `h<=2k` to remove the remainder in (20).

Combining (19) with the determinant and shared-minor budgets gives
the additional necessary conditions

```math
\sum_{i=1}^n\min(\mu_i,2\mu_2)
 +n\sum_{i=1}^{n-1}\nu_i\le n(n-1)M,             \tag{21}
```

and

```math
\sum_{i=1}^q\min(\mu_i,2\mu_2)
 +q\sum_i\nu_i-\sum_{i=n-q+1}^n\nu_i
 \le q(n-1)M.                                    \tag{22}
```

They also hold with the two profiles exchanged. Unlike the linear
constraints used for (1), these capped-profile conditions need not
be preserved by averaging two feasible profile pairs. A symmetric
numerical search alone therefore does not certify their global
pressure floor.

The replay checks (19)--(22) in every small control. It includes a
crossing with `mu=(0,1,2,7)` and `eta=(0,1,2,6)`. Thus removing
the cap and claiming `eta_i>=mu_i` for every order would be false.

## 9. The strengthened relaxation still permits eight letters

The inequalities above do not by themselves close the eight-letter
interface. This limitation has an exact integer witness. Take
`q=8`, `n=64`, `M=140`, and the identical formal profiles

```math
\mu=\nu=(0,80^{\times7},147^{\times55},224).
                                                        \tag{23}
```

These profiles satisfy (2), (10)--(11), (21)--(22), their exchanged
versions, and every shared-minor prefix with size from two through
eight. For example, their total is `8869<=8960`, and the profile
capped at `2*mu_2=160` has total `8805`. The left sides of the
strengthened determinant and eight-minor constraints are respectively
`562085<=564480` and `70259<=70560`.

At the rational scale `s=36/35`, define the integer

```math
P=36^{224}+7\,35^{80}36^{144}
  +55\,35^{147}36^{77}+35^{224}.
```

The formal pressure of (23) is

```math
B(36/35)=\frac{P^2}{36^{308}35^{140}}.
```

Direct integer arithmetic verifies

```math
P^{500}<8^{707}36^{77000}35^{35000}.
                                                        \tag{24}
```

Thus the formal pressure exponent satisfies
`rho<707/1000<1/sqrt(2)`. The replay checks all the constraints
and (24) without floating arithmetic. It does not construct a crossing
with (23), or assert that those profiles are realizable. Further
compatibility constraints or a construction remain necessary.
