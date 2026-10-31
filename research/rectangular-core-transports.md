# Constant core transports with a certified endpoint

## Result and scope

The rectangular construction admits arbitrary invertible constant core
changes when the left and right maps are coupled to preserve a
certified ramified endpoint. The individual dressed maps can have
poles. Matching paired types still produce polynomial crossings at
every word length.

This family contains height-raising core changes that are outside the
[independent secondary-flag family](rectangular-refined-repairs.md).
It remains within the [ten-mode full-core reduction](rectangular-full-core-factorization.md).
It supplies additional admissible candidates, with no improved
factoring exponent asserted.

## 1. The endpoint criterion

Let `K_raw` be the mixed rectangular restriction of the Gaussian
quartic kernel, with constants `alpha,beta`. Its core residual matrix
is `F=F(alpha,beta)`. The
[residual-repair proof](rectangular-residual-repair.md) gives an
energy-gauged tensor `Khat(t)` that is regular after
`w=1+t^T`, where `T` is the largest selected mode. At `t=0`, its
core endpoint is `F`, and its bulk endpoint is the crossing.

The core height is `h(v,y)=h_X(v)+h_Y(y)`, with
`h_X=(0,0,1,1)` and `h_Y=(0,0,1)`. Write `H=diag(c^h)`, where
`c=1+w`. Choose constant matrices `C,W` satisfying

```math
C F W=U,
\tag{1}
```

where `U` has diagonal one and every nonzero off-diagonal entry
strictly lowers the original additive height. The identity matrix
is allowed. Define

```math
K^{C,W}=(H^{-1}CH)\ K_{\rm raw}\ (H^{-1}WH).
\tag{2}
```

The input index order is `(x,u)` and the output order is `(v,y)`.
Constant Fock permutations may interchange the order within a pair.
All maps act only on the retained twelve-dimensional core.

The complete energy gauge includes bulk index factors. Those factors
commute with `C,W`, since neither map changes a bulk occupation bit.
The energy-gauged form of (2) is therefore exactly

```math
C\,\widehat K(t)\,W.
\tag{3}
```

It is regular, and its core endpoint is (1).

In any tensor power restricted to matching full-letter types on
`x,v` and on `u,y`, the energy gauge factors multiply to one.
Equation (3) proves regularity of every restricted entry. Before
the ramification the only possible finite pole is at `w=1`.
Regularity after `w=1+t^T` therefore proves polynomiality in `w`
itself. Every off-diagonal endpoint contribution strictly lowers
an additive height, while the paired types force the total height
drop to be zero. Only diagonal terms survive, giving the crossing.

This is the complete admissibility criterion. It does not require
`C` or `W` individually to be triangular. Their endpoint product
must be checked exactly against an admissible `U`.

### Pole and complexity bounds

The raw quartic contributes at most `c^-2`. Each of the two
dressings in (2) has exponents between `-2` and `2`. Thus a
uniform local denominator `c^6` is sufficient. No minimal local
denominator is needed for the all-length theorem. Paired-type
restriction removes the pole as just proved.

The same statement allows `P(w)/(1+w)^2` with `P(1)=1`; a
global power of `w` clears any Laurent terms in `P` and has value
one at the endpoint. Degree and coefficient-encoding costs remain
polynomial in the interval length and in the supplied polynomial.
All profile calculations must restore the actual rational orders
once and must retain their independent precision conditions.

## 2. Two constructive parameterizations

### Identity endpoint from a full inverse

The residual `F` has determinant one for every `alpha,beta` in
characteristic two, including `alpha*beta=1`. For any invertible
constant matrix `V`, choose

```math
C=V F^{-1},\qquad W=V^{-1}.
\tag{4}
```

Then `C F W=I`. Equivalently, this conjugates the fully repaired
tensor by `H^-1 V H` and its inverse.

There are `132` ordered off-diagonal positions in a twelve-state
core. Thus the elementary choices `V=I+aE_rs`, with `a` nonzero
in `F4`, give `396` coupled transvection candidates. They include
actual-height-raising directions. In characteristic two each such
`V` is its own inverse. General invertible matrices, including
nonunit diagonal coefficients, are permitted too.

An input-only full correction is the particularly simple case
`C=I,W=F^-1`. It satisfies the same criterion.

### Preserve the partial-repair baseline

For `1+alpha*beta != 0`, let `C0` be the partial correction and
`U=C0 F` its already proved height-unitriangular endpoint. Set

```math
C=V C_0,\qquad
W=U^{-1}V^{-1}U.
\tag{5}
```

Then

```math
C F W=V U U^{-1}V^{-1}U=U.
```

The choice `V=I` gives the original partial repair exactly. This
parameterization therefore explores a larger family around that
baseline, without first replacing its correction by a full inverse.
Since the partial seed has at most a simple local pole, its two
additional dressings give the also valid bound `c^-5` for (5).
The general bound `c^-6` remains available uniformly for (2).

These coupled formulas do not remove the common-flag requirement
for arbitrary independently selected `V,W` in the earlier family.
Their replacement hypothesis is the explicit product identity (1).

## 3. Active supports and coordinate updates

Every entry of either correction is a matrix unit on extreme
physical core modes. The ten-mode support proof already includes
all creation and annihilation operators needed for these units.
Equations (2), (4), and (5) therefore preserve that fixed support
bound and the scalar spectators. The all-length crossing theorem
again supplies the full generic rank of the retained active matrix.

The regular physical core coordinate bases permit these corrections
to be applied without rebuilding the CAR program. Their exact
intertwining identities are

```math
L c_i=c_i^{\rm coord}L,\quad L a_i=a_i^{\rm coord}L,
\qquad
c_i B=B c_i^{\rm coord},\quad a_i B=B a_i^{\rm coord}.
```

They extend to every core matrix unit. The vacuum minor in
`LB=I_core tensor P[I,J]` acts only on the multiplicity coordinate,
so it commutes with all coordinate core maps and causes no extra
factor. A core map changes only the indicated physical core labels
and keeps the vacuum multiplicity index fixed.

At infinity the output map acts on row core `(v,y)` and the input
map on column core `(x,u)`. At zero, output correction changes
row `y` and column `v`, while input correction changes row `x`
and column `u`. The dressed scalar coefficients retain the actual
height and physical gauge differences. This is an exact algebraic
identity; numerical implementations must separately track the
precision of every cached Laurent matrix.

### Exact cached updates and precision

The [coordinate update module](../tools/rectangular_coordinate_transforms.py)
implements this identity as sparse changes of the retained core labels.
It decomposes each supplied map as `I+N`, applies the nonzero entries
of `N`, and preserves the vacuum multiplicity index. A dressed entry
has coefficient

```math
a(1+z)^d z^e.
```

Negative `d` is cleared by a common power of the valuation unit
`1+z`. The returned unit denominator records that factor. This
does not affect a Smith order, including after row and column
weights are added.

An exact input numerator produces an exact output numerator. For a
truncated input, a recorded per-cell window is not a uniform absolute
cutoff. The updater first uses

```math
K=\min\bigl(K_{\rm recorded},\min_{A_{ij}\ne0}b_{ij}+p\bigr),
```

where `p` is the retained coefficient count and `b_ij` is that
cell's Laurent shift. The recorded cutoff also covers zero cells.
It discards all coefficients at or above `K`. A stage whose least
Laurent exponent is `e_min` lowers this common cutoff to
`K+min(0,e_min)`. Multiplication by the unit numerator does not
lower it. Successive output and input stages use this rule in
order, so unknown high coefficients cannot contaminate a claimed
low coefficient.

The caller must still certify the endpoint and the full generic
rank. For an inexact matrix, its largest requested weighted order
must be below the returned cutoff after the row and column lower
bounds are included. Cached algebra does not supply either promise.

The [independent update control](../tools/check_rectangular_coordinate_transforms.py)
compares complete rational-polynomial entries with fresh CAR
programs. Its [evidence](../evidence/rectangular-coordinate-transforms.json)
contains `147,456` exact entries: both valuation places, a common-flag
repair, and a coupled identity repair containing both a height-raising
entry and a nonunit diagonal coefficient. It also truncates source
windows deliberately and checks every returned coefficient below
the conservative cutoff.

### A two-matrix cache for scalar quartic changes

Fix the core correction and coordinate bases. The tensor and all its
linear coordinate changes depend affinely on the quartic numerator
`P(w)`. Two admissible matrices, with numerators `1` and `w`, therefore
give every Laurent numerator satisfying `P(1)=1`:

```math
K(P)=K(1)+Q(w)\bigl(K(w)+K(1)\bigr),
\qquad Q(w)=\frac{P(w)+1}{w+1}.
```

Here `Q` is Laurent polynomial. In particular,
`P=1+gamma*w^ell*(1+w)` needs only the monomial
`Q=gamma*w^ell`. At infinity the local coordinate is `z=w^-1`,
so local scalar exponents have the opposite sign from `w` exponents.

`linear_combination_coordinate_matrices` clears a common valuation-unit
denominator, performs the exact Laurent sums, and propagates the
same uniform precision cutoffs. The
[scalar cache control](../tools/check_rectangular_scalar_cache.py)
checks a negative Laurent exponent and the numerator `w^3` at both
places. Its [evidence](../evidence/rectangular-scalar-cache.json)
contains `92,160` exact rational entry comparisons with independently
rebuilt programs, together with deliberately truncated controls.

## 4. Exact API and controls

The [transport module](../tools/check_rectangular_identity_transports.py)
provides the following operations:

- `validate_identity_transport(alpha,beta,C,W)` checks `C F W=I`.
- `validate_typed_transport(alpha,beta,C,W,U)` first checks the
  target's height-unitriangular condition, then checks `C F W=U`.
- `conjugation_spec(alpha,beta,V)` returns the pair in (4).
- `typed_conjugation_spec(alpha,beta,V)` returns the pair and
  target in (5).
- `build_identity_transport(...)` and `build_typed_transport(...)`
  construct polynomial numerators directly from the original
  tensor, independently of the active Fock program.

Run the exact controls with:

```sh
python tools/check_rectangular_identity_transports.py \
  --output evidence/rectangular-identity-transports.json
```

The [recorded evidence](../evidence/rectangular-identity-transports.json)
checks all `396` coupled single transvections at the two-mode
core, `132` directions preserving the mixed partial endpoint,
and three dense invertible transports with a bulk mode. It
contains `245,520` two-copy paired-type polynomial checks and
all local ramified endpoint checks. Incorrect endpoint products
are rejected. The finite controls audit the constructor; equations
(1)--(3) supply the theorem for every word length.
