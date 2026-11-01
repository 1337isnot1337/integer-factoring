# A parity obstruction at the saturated zero-energy core

**Stronger update:** the [hereditary adjugate cap](zero-energy-hereditary-adjugate.md)
now excludes the retained target without any parity assumption. The
elementary argument below remains valid. Its listed escapes avoid only
the parity hypothesis, not the stronger target obstruction.

The current compatible five-letter target starts with normalized orders
`0,6,9,9` at both cuts, has deficit six, and gives two physical labels
energy zero on each paired alphabet. If these two labels have opposite
parity, a **scalar tensor preserving parity on that core cannot realize the target**
over any coefficient field. The same conclusion holds if both actual
profile prefixes separately dominate the target prefixes.

This is an obstruction to this target sector, not a lower bound for all
parity-preserving crossings or for factoring. It leaves other profile
shapes, parity mixing within the core and same-parity zero-energy labels open.
The [six-to-five rational interface](joint-six-to-five-retractions.md)
explicitly permits parity mixing, so its whole family is not excluded.

## 1. The saturated core profiles are forced

Use the normalization and support bounds of the
[graded adjugate target](https://github.com/1337isnot1337/integer-factoring/blob/3670b7df8db2bca84bc1dc1ea492c7044fa93e60/research/graded-adjugate-caps.md).
Write `R(z)` for the polynomial zero tensor and `S(z)` for its graded
reciprocal at infinity. Restrict every physical index to the two
zero-energy symbols. This gives a two-letter strict tensor with 4-by-4
cost matrices. Let their complete increasing profiles be `a,b`.

The individual leg charges can differ, but their paired sums vanish on
these symbols. Every determinant term of the core therefore has degree
at most `4d=24`, and every term of either parity-block determinant has
degree at most `2d=12`.

Core minor valuations dominate the corresponding full-tensor prefixes.
The target's fourth prefix is 24; the cell cap bounds each core total by
24. Thus `sum a=sum b=24`.

Let `b0` be the first infinity core order, which is the minimum order
of its entries. Reciprocal support then reduces every zero-entry degree
bound by `b0`, hence reduces the core determinant degree bound by `4b0`.
Consequently `sum a+4b0<=24`. Its known total forces `b0=0`. The reverse
argument forces `a0=0`.

The shared two-by-two minor on the zero-energy symbols gives
`(a0+a1)+(b0+b1)<=12`. Each two-prefix is at least six, so both are six.
The third prefix is at least 15, making `a2,b2>=9`; the totals and sorting
force the last two orders to be nine. Both core profiles are therefore
exactly

```
0,6,9,9.
```

## 2. Parity blocks contradict the forced profile

Give the two zero-energy symbols parities zero and one. Parity
preservation means that the core row matrix splits into a two-by-two
even block on states `00,11` and a two-by-two odd block on `01,10`.
The opposite cut has the same two-block structure. Diagonal physical
height factors preserve these zero entries.

Each block determinant has degree at most 12, including independently
split leg charges: the total charge of each block determinant term is
zero. The zero core determinant has valuation 24. Its factorization
into the two block determinants therefore forces each block determinant
to have valuation 12.

The core's first order is zero, so at least one block has a unit entry.
That block has Smith orders `(0,12)`. The full core must then have an
order at least 12, contradicting its forced largest order nine. This
argument uses only scalar determinant degrees and valuations, and works
over every field.

More generally, for deficit `d>0` the same proof excludes the separate
prefix bounds `S2>=d`, `S3>=5d/2`, `S4>=4d` at both cuts when the two
zero-energy symbols have opposite parity. The forced core profile would
be `(0,d,3d/2,3d/2)`, whereas a parity block supplies order `2d`.

The same block argument also applies to the new
[six-letter target](hereditary-adjugate-caps.md). Its first four orders
are `0,6,7,11`. Saturation of the fourth prefix still forces core totals
24 and zero minima, and the shared second prefix is still six. The core
third order must be at least seven, so its only possible profiles are
`(0,6,7,11)`, `(0,6,8,10)` and `(0,6,9,9)`. None contains the unavoidable
parity-block order 12. Thus a realization of this symmetric target needs
parity mixing in an opposite-parity zero-energy core, or a same-parity
assignment of the two zero energies. This does not contradict its
survival of the stronger universal hereditary inequalities.

## 3. Exact controls and constructive consequence

The [checker](../tools/check_zero_energy_parity_core.py) reconstructs the
forced integer core profile and checks all 49 parity-block allocations
at deficit six. It checks the formal charge cancellation at both cuts,
and independently compares block profiles with fraction-free elimination
on 24 literal F4 two-letter tensors: 48 complete profiles, 192 orders,
96 block determinants and 768 endpoint entries.

It also checks all 1,296 source entries of the six-letter ambient
operator for parity conservation. In the retained five-letter ordering,
the first two labels are `0` and `8`, of opposite parity. A permitted
rational correction, row zero and bounded-space basis vector 57 with
coefficient `omega`, produces six nonzero entries between opposite joint
parities while retaining the strict endpoint. Entry `(0,1)` lies inside
the zero-energy core, so this escapes the core assumption itself.
This verifies a concrete
escape from the parity assumption; it is not a favorable profile.

Run `make verify-zero-energy-parity-core`. The
[receipt](../evidence/zero-energy-parity-core.json) retains the controls,
source hashes and explicit parity-mixing assignment. The ordinary proof
above supplies the all-field statement; finite F4 examples test the
implementation. No new Lean work, independent review or stronger
factoring exponent is claimed.
