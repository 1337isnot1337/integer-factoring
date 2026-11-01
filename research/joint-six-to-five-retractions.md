# A regraded six-to-five Laurent retraction space

**Target update:** the [hereditary two-letter adjugate cap](zero-energy-hereditary-adjugate.md)
excludes the formerly compatible target used to motivate this interface,
regardless of parity. The construction below remains valid for other
profile choices, which must also satisfy the new cap.

This constructive interface compresses a 36-state joint reflection source
onto 25 states, then assigns the five physical labels the paired energies
`(0,0,12,16,18)` with declared deficit six. These differ from the natural
five-letter energy shape that is now excluded for the balanced criterion.
The aim is the newer compatible target in the
[graded adjugate note](https://github.com/1337isnot1337/integer-factoring/blob/3670b7df8db2bca84bc1dc1ea492c7044fa93e60/research/graded-adjugate-caps.md),
not the older four-letter target already excluded by cell caps.

The interface has not produced a favorable crossing. Its value is an exact
smaller search space with arbitrary internal rational poles and actual
Laurent outputs. The eight validation probes retain complete profiles;
none establishes a factoring improvement. The score remains `1/2`.

## Source and complete rational output space

Use the interval-eight reflection on physical labels
`(0,8,4,12,2,10)`. Removing its common `w^8` and deflating `z=w^2` gives
a polynomial 36-by-36 operator `A(z)` with `A(1)=I`. Retain the first five
labels on both joint directions. Let `P0,L0` be these coordinate maps,
and let `Q` insert the eleven discarded coordinates. Put

```
C0=L0 A P0,       B=Q^T A P0.
```

The [implementation](../tools/joint_six_to_five_retractions.py) constructs
only the requested Pfaffian coefficients using the existing implicit
reflection source. The [row-space routines](../tools/polynomial_retraction_spaces.py)
generalize and credit the existing rational retraction algorithm. This
source differs from both the earlier five-to-four and nine-to-eight
compression families.

The eleven rows of `B` have rational rank nine. The
[receipt](../evidence/joint-six-to-five-retractions.json) supplies a
9-by-25 polynomial basis `S`, a degree-five polynomial right inverse `R`,
and an 11-by-9 matrix `U` satisfying

```
S R=I_9,       B=U S,       rank U=9.
```

The basis row degrees are `0,7,4,7,4,9,6,9,16`. A nine-row minor of `U`
is invertible at `z=omega` in F4, certifying full rational column rank
over every characteristic-two field. This supplies a rational left
inverse `W` with `W U=I_9`.

Every rational left inverse of `P0` is `L=L0+T Q^T`, and its output is
`C=C0+T B`. If `C` is Laurent polynomial, multiplying its correction by
`R` shows that its coefficients in `S` are Laurent polynomial. Strictness
at one makes each coefficient divisible by `z+1` in the Laurent ring.
Consequently the complete strict Laurent-output family is

```
C=C0+(z+1)G(z)S(z),       G a Laurent-polynomial 25-by-9 matrix.
```

Conversely `T=(z+1)G W` realizes every member. This removes arbitrary
internal pole locations from the output search without assuming the
output itself has no pole at zero.

## Bounds at the new energies

Assign physical heights `h=(0,0,6,8,9)` and paired heights
`H_(x,y)=h_x+h_y`. Set `mu_zero=10` and `mu_infinity=-16`. Require

```
val C[i,j] >= 10-H_i,       deg C[i,j] <= 16+H_j.
```

Both normalized cuts are regular, but the raw output may have poles.
The declared deficit is `-mu_zero-mu_infinity=6`. The reference `C0`
obeys these bounds. A correction fixing the first `k` raw matrix jets
at both cuts has lower bound `10+k-H_i` and upper bounds `16+H_j-k`.

The constant rows of `S(0)` and its leading rows are independently
linearly independent. They give predictable minimum valuation and
maximum degree of a Laurent row combination. Therefore the coefficient
shifts in each row range from its lower bound through
`max(upper)-degree((z+1)S_i)`. An exact binary kernel imposes all column
bounds. This is a complete bounded Laurent space, rather than a sampled
pole-order window.

For `k=0`, the dimensions at the ten occurring paired heights are:

| Paired height | Row dimension | Common row clearing |
| --- | --- | --- |
| 0 | 58 | 0 |
| 6 | 110 | 0 |
| 8 | 128 | 0 |
| 9 | 137 | 0 |
| 12 | 164 | 2 |
| 14 | 182 | 4 |
| 15 | 191 | 5 |
| 16 | 200 | 6 |
| 17 | 209 | 7 |
| 18 | 218 | 8 |

## Full profiles and cost restoration

The native interface handles each entry's actual Laurent shift at both
cuts, including the opposite reshuffle. Its reported raw zero and
infinity orders can be used directly with paired energies `2h`.
If normalized orders are used instead, the full balanced critical
pressure is

```
6t + log sum exp(-t r) + log sum exp(-t s)
   - sqrt(2) (log sum exp(-t X/sqrt(2))
              + log sum exp(-t Y/sqrt(2))).
```

The `6t` term is required. Each validation case independently clears
all rows with one global monomial, computes both profiles by Python
fraction-free polynomial elimination and restores the opposite shifts.
This agrees with the native entrywise-Laurent calculation on all 400
complete orders. Five thousand endpoint entries are also checked.

The [checker](../tools/check_joint_six_to_five_retractions.py) additionally
checks 1,296 source endpoints, 625 retained entries against the earlier
five-letter formula, the polynomial inverse and source containment,
70 complete bounded spaces, and 63 space regressions against the
separately implemented earlier polynomial module.

Run `make verify-joint-six-to-five`. The focused replay postdates the
clean full-suite source `f2a16c0`. No new Lean work or independent
mathematical acceptance is claimed. Choosing a favorable combination of
the output rows, or changing the embedding as well, remains open.
