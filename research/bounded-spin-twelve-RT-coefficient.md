# R and T inside the twelve-site coefficient algebra

## Result

In the exact ten-dimensional twelve-site coefficient table, the coefficient
unit `C_0` (the replica-swap function in the cut normalization) and the
compressed site-transposition class sum `T` generate a **five-dimensional**
commutative subalgebra. In the basis

```text
0, 2, 3, 4, 5, 6, 7, 8, 9, 10
```

the two elements are

```text
R = C_0 = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
T       = (3, 3, 1, 2, 1, 3, 3, 1, 2, 1).
```

The powers of `T` stop growing after `T^4`, and the exact relation is

```math
T^5 = 4T + 2T^3,
```

so the minimal polynomial is

```math
x(x-1)^2(x+1)^2.
```

The eight primitive central idempotents of the full coefficient algebra see
`T` with eigenvalues `0,1,4`. The two dual-number blocks are exactly the
blocks where the action has a nonzero nilpotent part. The scalar-block
eigenvalues are `0` on idempotents `0,6`, `1` on `1,5,7`, and `4` on `3`.
The two dual blocks have characteristic polynomials `(x-1)^2` and
`(x+1)^2`.

This is a coefficient-algebra computation. It is not the eight-dimensional
generic `R,T` **orbit on the 7921-state cut module** recorded in
[the twelve-site orbit note](bounded-spin-twelve-coefficient-dimension.md#the-rt-subalgebra-is-at-least-eight-dimensional); the two
dimensions measure different objects. No cut-module rank, two-cut pairing,
rank bound, or factoring exponent changes.

## Method and controls

The [table certificate](../evidence/bounded-spin-twelve-coefficient-table.json)
supplies the ten entry witnesses, the entry matrix, the multiplication
table, and the regular representation. The [idempotent
certificate](../evidence/bounded-spin-twelve-coefficient-idempotents.json)
supplies the eight central idempotents and their block bases; its
`table_sha256` pins the table used here.

`T` is evaluated as a coefficient function on the ten witness quadruples,
using the row-cut orientation

```text
T((p1,p2),(p3,p4)) = sum_{i<j} M_ij[p1,p4] M_ij[p2,p3],
```

where `M_ij` are the normalized single-transposition matrices at length
twelve. Solving the ten-by-ten entry system gives the coordinates above.
The orientation is fixed by the ten-site control: the same method recovers
the published ten-site coordinates

```text
T = (0, 2, 0, 1, 0, 2, 2, 0)
```

in the basis `0,2,3,4,5,6,7,10`, and the ten-site subalgebra has dimension
`6`, matching the published model.

At twelve sites the checker additionally verifies the `T` coordinates on
24 nonzero and 25 zero structured class-sum entries outside the witness.
It then checks the power relation, the five-dimensional span, all
one-direction extensions, and the three minimal completing pairs

```text
(coefficients 2,6), (3,5), (3,6).
```

Each single coefficient direction outside the subalgebra raises the
dimension only to `6`, `7`, or `8`; no one direction completes the full
ten-dimensional table.

## What changed

The table note left the identification of the `R,T` subalgebra inside the
ten-dimensional algebra open, and the earlier twelve-site orbit note gave
only a lower bound of eight for the orbit on the cut module. This note
supplies the exact coefficient-algebra answer: dimension five, with the
explicit minimal polynomial and block action. It also identifies the
coefficient directions that complete the algebra in pairs, a finite
structural input for the still-open cut-module rank/pairing question.

The two dual blocks are important for the next step. They are the only
blocks where `T` has a nilpotent part, so any decomposition that tries to
use the coefficient algebra to separate the eight primaries must handle
both the scalar multiplicities and the two square-zero directions. The
remaining target is still the exact row and column ranks of these
idempotents on the 7921-state cut modules, followed by the two-cut pairing.

## The separator direction at the algebra level

The completing pair `(3,5)` is the separator combination
`S=C_3+C_5` carried over from the ten-site lane. In the basis above it is

```text
S = (0, 0, 1, 0, 1, 0, 0, 0, 0, 0),
```

and it raises the five-dimensional `<T>` subalgebra to the full
ten-dimensional coefficient algebra. Its actions on the eight primitive
central idempotents are

```text
idempotent   0   1   2          3   4          5   6   7
S action     2   0   [4 4;4 0]  3   [2 4;4 3]  3   4   1
```

with the two-dimensional matrices acting on the dual blocks. The joint
`T,S` character pairs are

```text
idempotent   0     1     2     3     4     5     6     7
(T,S)       (0,2) (1,0) (1,2) (4,3) (4,0) (1,3) (0,4) (1,1)
```

so `S` separates the primaries inside every `T` block: the `T=0` block
splits as `2,4`, the `T=1` block as `0,2,1,3`, and the `T=4` block as
`3,0`. In particular, `<T,S>` already has one character for every one of
the eight primaries.

This is still a coefficient-algebra statement. It does not prove that
the same coefficient word separates the twelve-site **cut-module** coarse
blocks, and it supplies no exact rank or two-cut pairing. It does remove
the character-level ambiguity in the separator direction and gives the
next computation a fixed pair of algebra directions to test.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_RT_coefficient.py
make verify-bounded-spin-twelve-RT-coefficient
```

The checker replays the ten-site control, the twelve-site entry solve,
the structured holdouts, the minimal polynomial, the one-direction
extensions, the completing pairs, the eight `T` block actions, and the
`T,S` joint characters. The complete factoring score remains `1/2`.
