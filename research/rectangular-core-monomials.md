# Endpoint-vanishing core monomials and opposite-charge cubic pairs

## Result and scope

The corrected rectangular construction and the full-alphabet construction
admit arbitrary extreme-mode monomials whose ramified energy is strictly
positive. This includes odd degrees in characteristic two. The existing
bounded active spaces still contain all their creation and annihilation
operators, so their profiles can be calculated without growing the matrix
order with the interval length.

For the full physical alphabet, perturbations whose particle-count changes
all have the same strict sign cannot improve the profile bound. Opposite
signs are therefore the useful cubic search family. There are two negative
and eight positive cubic patterns with positive physical index increase.
All sixteen pairs change a two-copy fixed-type tensor in exact small
controls. This establishes an admissible family to search; it does not
establish a pressure improvement or a new factoring exponent.

## 1. A uniform endpoint criterion

Let the extreme modes be `0,T`, with `T=M-1>0`, and put `c=1+w`.
Choose a nonempty monomial `Q` containing at most one variable on each of
the four legs `X,Y,U,V`. An absent leg contributes neither a variable nor
an index. Define

```math
\Delta(Q)=\sum_{y_b\mid Q}b+\sum_{v_d\mid Q}d
             -\sum_{x_a\mid Q}a-\sum_{u_e\mid Q}e,
\qquad
q(Q)=\#\{Y,V\text{ variables in }Q\}
          -\#\{X,U\text{ variables in }Q\}.
\tag{1}
```

Start with the Gaussian `G` and the existing canonical quartic `Q0`.
The raw tensor is the explicitly affine expression

```math
G(w)\left(1+\frac{\gamma Q_0}{c^2}
                  +\sum_j\frac{P_j(w)}{c^{r_j}}Q_j\right).
\tag{2}
```

Here `gamma=1` uses the rectangular mixed codes, while
`gamma` different from one allows the
[full physical alphabet](full-alphabet-quartic-repair.md).
The expression is a sum; products between its added monomials are not
implicitly included. Each `P_j` is a Laurent polynomial, each `r_j` is a
nonnegative integer, and the condition is

```math
\Delta(Q_j)+\left(\mathrm{ord}_{w=1}P_j-r_j\right)T>0.
\tag{3}
```

Under `w=1+t^T`, the actual occupied-index gauge multiplies `Q_j` by
`t^Delta(Q_j)`. Equation (3) makes its full coefficient divisible by
`t`, while the gauged Gaussian is regular. The added terms are therefore
regular and vanish at the endpoint. The same physical codes and constant
correction retain exactly the old height-unitriangular endpoint.

For independently matching full-letter types on `(X,V)` and `(Y,U)`,
the actual index gauges cancel in every tensor-power entry. At the
endpoint, equal total height excludes every strict height decrease.
The restricted tensor is thus the crossing. After clearing a fixed
global power of `w`, its only possible finite pole is at `w=1`.
Regularity after ramification excludes that pole, proving polynomiality
at every length. The degree is at most a seed-dependent constant times
the length. This is an ordinary proof of the same effective crossing
interface used by the original construction.

The proof uses the characteristic-two exterior algebra, in which distinct
variables commute and each variable has square zero. It does not require
the added monomials to have even degree.

## 2. Cubics of positive index increase

For coefficients regular and nonzero at `w=1`, condition (3) reduces to
`Delta>0`. The ten cubic patterns are as follows. A dash denotes a
missing leg, and all displayed index increases are in units of `T`.

| Label | `(x,y,u,v)` | `q` | `Delta/T` |
| --- | --- | ---: | ---: |
| `minus_v` | `(0,-,0,T)` | -1 | 1 |
| `minus_y` | `(0,T,0,-)` | -1 | 1 |
| `plus_0` | `(-,0,0,T)` | 1 | 1 |
| `plus_1` | `(-,T,0,0)` | 1 | 1 |
| `plus_2` | `(-,T,0,T)` | 1 | 2 |
| `plus_3` | `(-,T,T,T)` | 1 | 1 |
| `plus_4` | `(0,0,-,T)` | 1 | 1 |
| `plus_5` | `(0,T,-,0)` | 1 | 1 |
| `plus_6` | `(0,T,-,T)` | 1 | 2 |
| `plus_7` | `(T,T,-,T)` | 1 | 1 |

The two patterns of index increase `2T` also permit a coefficient
`P/c` with `P(1)` nonzero. All other listed patterns require a regular
coefficient, unless its numerator cancels the local pole. More general
patterns with zero or negative index increase are also admissible when
the numerator vanishes to the order required by (3).

Only `plus_7` requires the additional `x_T` operator and the ten-mode
support. All other patterns lie in the seven active modes at infinity
and nine at zero. The physical alphabet can remain full in either case.

## 3. Single-sign perturbations cannot help the full alphabet

This section requires full occupation states on every physical leg.
The Gaussian, the canonical quartic, and both full-alphabet corrections
preserve total particle count. Thus the unperturbed tensor `K` has

```math
N_V+N_Y=N_X+N_U
\tag{4}
```

on its support. An added monomial of charge `q` changes the left side
minus the right side by exactly `q`, including after the correction.

If all charges have the same strict sign, any term in a tensor power
using an added monomial has a nonzero total charge. Matching the two
pairs of full-letter types forces total charge zero. Consequently the
entire fixed-type tensor power is exactly the unperturbed tensor power,
for every length, as a rational function before specialization.

There is also a direct profile obstruction, valid for every diagonal
real row and column gauge. At infinity, grade the rows by `N_V+N_Y`
and columns by `N_X+N_U`. At zero, use row grade `N_X-N_Y` and column
grade `N_V-N_U`. In both cuts the base matrix is block diagonal in
these grades. The perturbation occupies strictly one direction between
grade blocks.

Any nonzero `k` by `k` base minor has matching grade multisets on its
rows and columns. Its determinant in the perturbed matrix is unchanged:
every determinant term using a strict-sign entry would have a nonzero
sum of grade differences, whereas the fixed row and column sets have
equal total grade. In particular, the full determinant is unchanged.
Diagonal weights multiply each minor by a common scalar and preserve
these identities.

Let `a_1<=...<=a_N` and `b_1<=...<=b_N` be the weighted invariant
orders of the base and perturbed matrices. The preserved minors give

```math
\sum_{i=1}^k b_i\le\sum_{i=1}^k a_i\quad(1\le k\le N),
\qquad \sum_{i=1}^N b_i=\sum_{i=1}^N a_i.
\tag{5}
```

Thus `b` majorizes `a`. Convexity of the exponential gives, for every
positive `ell`,

```math
\sum_i e^{-\ell b_i}\ge\sum_i e^{-\ell a_i}.
\tag{6}
```

Both local profile partitions can only increase, while the physical
letter energies are unchanged. Their Gibbs pressure root can therefore
only increase at the same gauges. Optimizing all gauges does not remove
this obstruction. It applies to any collection of strictly one-sign
charges, not just one cubic.

The mixed three-state codes identify a top state with a top-plus-both
combination. They do not preserve this particle grading. This section
does not assert an obstruction for those mixed codes.

## 4. Reducing opposite-charge parameters

Consider one charge-minus-one cubic with coefficient `a*w^i` and one
charge-plus-one cubic with coefficient `b*w^j`, where `a,b` are nonzero.
On the full alphabet, scaling every output variable by `k` and every
input variable by `k^-1` fixes the base and commutes with its correction.
It changes the two coefficients to `a/k` and `bk`. Taking `k=a`
reduces the field parameters to `(1,ab)`. Nonzero field scalings are
units at both places and preserve every weighted minor valuation.

The same physical gauge with `k=w^i` changes the Laurent powers to
`(0,i+j)`. At both places this gauge adds `-i*log(s)` to the row charge
and `+i*log(s)` to the column charge. These shifts leave the paired
physical energies unchanged. Therefore the unrestricted optimized
pressure depends only on the sum of the two Laurent powers. A screen
with a bounded charge box must transport its starts or bounds explicitly;
it cannot treat this equivalence as equality at an unchanged gauge.

It suffices to search three nonzero coefficient products over `F4`, one
total Laurent power, and the sixteen pairs of modes. The two cubic
terms can still affect the typed tensor because their charges can
cancel across different tensor factors.

## 5. Exact implementation and controls

The [independent constructor](../tools/check_rectangular_core_monomials.py)
uses Gaussian minors in the original physical variables. Its interface
is

```text
(x_mode_or_None, y_mode_or_None, u_mode_or_None, v_mode_or_None,
 [(F4_coefficient, w_exponent), ...], c_pole_order)
```

`None` means that the leg is absent; the mode `0` is a present variable.
`validate_extra_monomials` combines repeated powers and checks (3)
exactly. `core_monomial_seed` supports both physical alphabets and both
base corrections, recording any global `w` power used for polynomial
numerators.

The active programs and coordinate compiler accept `extra_monomials`.
They act with the requested core CAR word on the pure Gaussian, then
apply the original physical code maps and correction. The natural
additional scalar powers, with absent terms omitted, are

```math
g_\infty=M[\,x_a\mid Q\,]-M[\,y_b\mid Q\,]
       +\sum_{u_e\mid Q}(2e-M+1)+\sum_{v_d\mid Q}(2d-M+1),
\qquad
g_0=\sum_{x_a\mid Q}(1-2M+2a)+\sum_{y_b\mid Q}(1+2b).
\tag{7}
```

A coefficient `gamma*w^ell/c^r` thus adds the local scalar
`gamma*z^(g_inf+r-ell)*(1+z)^(-r)` at infinity and
`gamma*z^(g_0+ell)*(1+z)^(-r)` at zero. The compiler checks that every
required CAR operator exists in the selected support.

The [constructor evidence](../evidence/rectangular-core-monomials.json)
contains all 24 positive-index patterns in degrees one through four,
each at three field/power choices. It checks the unchanged typed tensor
for all ten single full-alphabet cubics, and gives a nonzero two-copy
typed difference for each opposite-charge pair. It also includes a
three-mode pair and a permitted `c` pole.

The [coordinate audit](../tools/check_rectangular_monomial_coordinates.py)
compares the native CAR program with the original tensor at two points
of `F16`, in both cuts, with arbitrary quarter-weight profiles. The
[charge audit](../tools/check_rectangular_monomial_charges.py) checks
the prefix inequalities (5), field-unit gauge identities, and Laurent
charge-gauge identities independently in small full matrices.

```sh
python tools/check_rectangular_core_monomials.py \
  --output evidence/rectangular-core-monomials.json
python tools/check_rectangular_monomial_coordinates.py \
  --output evidence/rectangular-monomial-coordinates.json
python tools/check_rectangular_monomial_charges.py \
  --output evidence/rectangular-monomial-charges.json
```

The cache construction is affine in every added coefficient. One
compiled seed per unit monomial and a common base are enough to form
any finite Laurent combination exactly. The cached-matrix precision
rules still apply, and a floating screen is not an exact pressure
certificate.
