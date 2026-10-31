# Four additional extreme-mode quartics

## Result and scope

The rectangular partial repair remains valid after four further
extreme-mode quartic terms are added to the original Gaussian.
Their coefficients can be arbitrary Laurent polynomials. Each term
vanishes under the ramified energy gauge, so the existing correction,
paired-type polynomiality, and crossing endpoint remain valid.

Three terms use only the seven-/nine-mode active supports already
constructed. The fourth needs the ten-mode extension. These are
admissible search candidates, with no improved pressure asserted.

## 1. Four terms of positive energy

Write `T=M-1`, and denote a quartic by
`Q_(a,b,c,d)=x_a y_b u_c v_d`. Its physical index increase is
`Delta=b+d-a-c`. The following four terms all have `Delta=T`:

| Label | Modes `(a,b,c,d)` | Quartic |
| --- | --- | --- |
| 0 | `(0,T,0,0)` | `x0 yT u0 v0` |
| 1 | `(0,0,0,T)` | `x0 y0 u0 vT` |
| 2 | `(T,T,0,T)` | `xT yT u0 vT` |
| 3 | `(0,T,T,T)` | `x0 yT uT vT` |

Let `G` be the Gaussian and keep the original canonical quartic
`Q0=x0 yT u0 vT` with coefficient `c^-2`, where `c=1+w`.
For arbitrary Laurent polynomials `f_j(w)`, use the raw tensor

```math
G(w)\left(1+c^{-2}Q_0+\sum_{j=0}^3 f_j(w)Q_j\right).
\tag{1}
```

The sum specifies the construction directly. The displayed quartics
share physical variables, so products of distinct displayed quartic
monomials vanish in the exterior algebra as well.

After `w=1+t^T`, the full physical energy gauge multiplies a
quartic of index increase `Delta` by `t^Delta`. A Laurent
polynomial in `w` is regular at `w=1`; therefore every extra
term of (1) contains `t^T` after gauging. The Gaussian is already
regular under that gauge. Every extra term consequently vanishes
at the ramified endpoint.

The mixed code maps preserve the physical height on their support.
Applying the same constant partial correction to the energy-gauged
tensor thus leaves its old height-unitriangular endpoint unchanged.
Matching full-letter types cancel the energy factors in every
tensor power. Regularity and the strict-height endpoint argument
then prove polynomial crossing at every word length. A global
power of `w` clears any Laurent pole at zero and has value one
at the crossing point.

This is an application of the
[positive-energy endpoint argument](positive-energy-quartic-endpoints.md)
to the mixed rectangular code and its existing residual correction.

### Optional removable local denominators

The same proof permits a coefficient `P(w)/c^r` provided

```math
\Delta+\bigl(\mathrm{ord}_{w=1}P-r\bigr)T>0.
\tag{2}
```

The inequality is strict. A constant coefficient divided by `c`
has zero ramified excess for the four terms in the table. It
changes the residual endpoint and requires a different correction.
The present constructor rejects it. A numerator divisible by `c`
can cancel that denominator and satisfy (2).

## 2. Active operators and natural gauges

At infinity, the extra operator is

```math
c_{v_d}c_{y_b}\,\Gamma(\Sigma)\,a_{u_c}a_{x_a}.
```

At zero, the cut instead uses left creations on `x_a,y_b`
and right annihilations on `u_c,v_d`, in the same order as
the original canonical perturbation.

For the natural physical gauges, the quartic's additional scalar
exponents are

```math
g_\infty=2(c+d-(M-1)),\qquad
g_0=2(a+b-(M-1)).
```

Thus a coefficient `gamma*w^ell/c^r` contributes
`gamma*z^(g_inf+r-ell)*(1+z)^(-r)` at infinity and
`gamma*z^(g_0+ell)*(1+z)^(-r)` at zero. These formulas follow
by adding the four physical mode weights and substituting
`w=z^-1` or `w=z`, respectively.

Labels `0,1,3` use only the previously retained `x0` operator
and the complete `u,y,v` extreme operators. They remain in the
seven-/nine-mode family. Label `2` uses the additional `xT`
operator and therefore uses the ten-mode construction. The
implementation checks the actual required CAR names in the
supplied support; it does not silently omit a missing mode.

All scalar spectators are unchanged. Once the corrected tensor
has its crossing endpoint, the usual all-length argument also
supplies its full generic active rank.

## 3. Independent constructor and cached profiles

The [original-minor constructor](../tools/check_rectangular_extra_quartics.py)
accepts each extra as

```text
(x_mode, y_mode, u_mode, v_mode, [(F4_coefficient, w_exponent), ...], c_pole_order)
```

`validate_extra_quartics` combines repeated Laurent powers and
checks (2) exactly. `extra_quartic_seed` expands the Gaussian
minors, applies the physical mixed codes, and then applies the
original partial correction to the polynomial numerators. It
records the common `w` power used to clear Laurent terms.

The active Fock builder and coordinate builder accept the same
list through `extra_quartics`. The cache includes that list and
the constructor's source fingerprint in its key.

For a single fixed quartic `Q_j`, let `K_base` be the original
seed and `K_j` the seed with coefficient one on `Q_j`. Then

```math
K_{f_j}=K_{\rm base}+f_j(w)(K_j+K_{\rm base}).
```

Several extras are obtained by adding these differences. This
allows all coefficient and Laurent-power screens to reuse one
compiled basis matrix per quartic. Common unit denominators and
coefficient precision remain governed by the exact cached-matrix
rules, independently of the endpoint proof.

## 4. Exact evidence

The [constructor evidence](../evidence/rectangular-extra-quartics.json)
contains `36` single-extra cases, covering all four terms,
three nonzero `F4` coefficients, and three Laurent powers.
It also includes mixed three- and four-mode cases. All
`159,936` two-copy paired-type polynomial checks and all local
ramified endpoint checks pass. A residual-changing marginal
coefficient is rejected explicitly.

```sh
python tools/check_rectangular_extra_quartics.py \
  --output evidence/rectangular-extra-quartics.json
python tools/check_rectangular_extra_coordinates.py \
  --output evidence/rectangular-extra-coordinates.json
```

The independent coordinate audit compares the active CAR program
with the complete original four-mode tensor. It includes all
four extra terms, negative Laurent powers, their simultaneous
sum, and an explicitly canceled `c` denominator. The
[coordinate evidence](../evidence/rectangular-extra-coordinates.json)
contains `571,392` field-entry comparisons and twenty complete
weighted profiles, each with `192` orders. All comparisons agree.
These finite controls audit the constructor and natural gauge
formulas; Sections 1 and 2 provide the uniform proof.
