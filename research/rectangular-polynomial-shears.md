# Polynomial core shears that fix the endpoint

## Result and scope

The mixed rectangular construction admits a Laurent polynomial core
change that equals the identity at `w=1`. Such a change preserves
the crossing endpoint at every word length, even when an individual
matrix entry raises the original physical height.

This supplies a new search family around the partial-repair baseline.
The construction and exact controls below establish admissibility.
They do not assert a smaller pressure or a stronger factoring exponent.

## 1. Endpoint preservation

Write the partial seed as `K0`, with energy-gauged form `Khat0` and
height-unitriangular endpoint `U`. The
[constant transport proof](rectangular-core-transports.md) establishes
that `Khat0` is regular after the ramification `w=1+t^T`.
Let `H=diag((1+w)^h)` use the same physical core height.

Choose a matrix `V(w)` whose entries are regular at `w=1` and satisfy
`V(1)=I`. Define the output change by

```math
K_V=(H^{-1}V(w)H)K_0.
```

Its complete energy-gauged form is exactly
`V(1+t^T) Khat0(t)`. It is regular and has endpoint `U`.
The bulk energy factors commute with `V`, since the latter changes
only the retained extreme core states. Matching full-letter types
on the paired legs cancel all energy factors in a tensor power.
The same strict-height argument as for `K0` then removes every
off-diagonal endpoint contribution.

If `V` is Laurent polynomial in `w`, the only possible finite poles
of `K_V` occur at `w=0,1`. A global power of `w` clears the first.
Regularity after ramification at `w=1` removes the second after
paired-type restriction. The resulting tensor is therefore a
polynomial crossing at every word length. The global power has
value one at the crossing point.

An input change `K0(H^-1 W(w)H)` has the same proof whenever
`W(1)=I`. Output and input changes may also be combined. They
need no secondary height flag, since their endpoint maps are
identities individually.

## 2. Elementary family and uniform poles

For distinct core indices `r,s`, take

```math
V(w)=I+P(w)E_{rs},\qquad P(1)=0.
```

In characteristic two this matrix is its own inverse. Its endpoint
is the identity, for every ordered pair of distinct indices.
There are `132` such directions in the twelve-state core.

A useful elementary choice is

```math
P(w)=\gamma w^\ell(1+w),
\qquad \gamma\in\mathbb F_4^\times,\quad\ell\in\mathbb Z.
```

The dressed output coefficient is
`gamma*w^ell*(1+w)^(1+h(s)-h(r))`. Heights range from zero to
two, so this coefficient has at most a simple pole at `w=1`.
The partial seed also has at most a simple pole. Thus `c^2`,
where `c=1+w`, is a sufficient local denominator for every single
such shear. The paired-type theorem removes that pole at all lengths.

More generally, `P(1)=0` means that `P` is divisible by `1+w`
in the Laurent polynomial ring, so the same local bound applies.
A finite product of these changes stays regular and equal to the
identity at the endpoint. Its degree and coefficient size grow
polynomially in the supplied list of changes.

All operations remain supported on the same active core as their
constant versions. Output changes use the seven-/nine-mode bases.
Input changes that involve the `x` core use the ten-mode extension.
The scalar spectators remain unchanged.

## 3. Exact cached implementation

Let `T_rs(K)` denote the cached result of applying the dressed
constant map `I+E_rs`. Linearity gives

```math
K_{P}=K+P(w)\bigl(T_{rs}(K)+K\bigr).
```

The intermediate constant map need not have an admissible endpoint.
It is used only in this algebraic sum. The final map has the
proved identity endpoint because `P(1)=0`.

The [polynomial shear module](../tools/check_rectangular_polynomial_shears.py)
provides `polynomial_coordinate_shear` for this construction. It
accepts a previously cached `T_rs(K)` so several Laurent powers
can reuse the same sparse core update. The underlying
[Laurent matrix helpers](../tools/rectangular_coordinate_transforms.py)
clear unit denominators and propagate conservative precision cutoffs.

The independent function `polynomial_shear_seed` instead changes
the original polynomial numerator tensor directly. It clears any
negative `w` powers and records that global monomial explicitly.
When comparing this numerator with an unscaled cached matrix,
the recorded global power must also be applied to the cached matrix.
At infinity its local exponent has the opposite sign.

## 4. Exact controls

Run:

```sh
python tools/check_rectangular_polynomial_shears.py --dense \
  --output evidence/rectangular-polynomial-shears.json
```

The [recorded evidence](../evidence/rectangular-polynomial-shears.json)
checks all `132` output directions at three coefficient/power choices.
These give `396` endpoint controls and `166,320` two-copy paired-type
polynomial checks. It also checks three input directions and applies
every chosen shear twice to recover the original rational tensor,
up to the recorded global monomial.

The independent four-mode controls include a height-raising direction,
an equal-height direction with a negative Laurent power, and a
height-lowering direction that changes the old Gaussian background.
At both valuation places, they compare `276,480` entries over `F16`
with the original dense tensor. Twelve complete weighted profiles
also agree, each with `192` orders and arbitrary quarter-integer
core and bulk-charge weights. The finite controls audit the
constructor; the endpoint argument in Section 1 proves the
all-length statement.
