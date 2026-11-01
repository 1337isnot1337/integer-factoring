# Exact pressure certificates without large integer powers

## Result

The complete ninety-six-letter physical code has the exact uniform bound

```math
\boxed{\rho<\frac{710005874881191}{10^{15}}
              =0.710005874881191.}
```

The terminating decimal is the rational bound itself. It is not a
rounded numerical optimization result. The proof evaluates rational
polynomials and logarithms with exact outward intervals. No floating
comparison decides whether the certificate passes.

This improves the separate [integer-power certificate](hot96-rational-pressure.md)
`rho < 49/69` for the same code. The new bound still exceeds
`1/sqrt(2)`, so it does not improve the established factoring score.

The complete 9,216-channel profile and its construction are recorded in
[hot96-charge-profile.json](../evidence/hot96-charge-profile.json).
This note proves and checks the pressure arithmetic on that frozen
profile. It does not claim a second independent construction of the
large polynomial matrix.

## 1. A rational inequality with the full physical type

Let sigma_i be the complete nonnegative natural orders after cancelling
the opposite scalar shifts of the two cuts. Let E_a be the physical
energies of the identical paired code. For positive integers A,B and a
rational number `0<z<1`, put

```math
R=A/B,\qquad
P=\sum_i z^{A\sigma_i},\qquad
Q=\sum_a z^{BE_a}.
```

Both lists contain zero, so `P,Q >= 1`. It suffices to prove

```math
B\log P<A\log Q.                                      \tag{1}
```

This is equivalent to `P^B<Q^A`, but computing those large rational
powers is unnecessary. Take the exact positive rational type

```math
p_a=z^{BE_a}/Q,\qquad s=A\log(1/z).
```

Its entropy is `H=log Q+(s/R) E_p`. The complete natural numerator is
`2 log P+2s E_p`, including both physical means. Consequently (1)
gives

```math
\frac{2\log P+2sE_p}{4H}<\frac R2.
```

All types, the code, and the coefficient field are fixed. Clearing the
rational type denominators gives admissible finite blocks. The strict
gap absorbs the polynomial coefficient and reconstruction factors.
The [finite amplification theorem](finite-block-crossing-amplification.md)
therefore gives the same strict bound at every external label width.
This is the same normalization used in the earlier exact integer
certificates, with a more efficient proof of their scalar inequality.

## 2. Exact polynomial enclosures

Fix an integer precision b and write `S=2^b`. An interval `[l,u]`
represents the rational interval `[l/S,u/S]`. For nonnegative intervals,
outward multiplication is

```math
[l_1,u_1]\,[l_2,u_2]
\subseteq
\left[\left\lfloor\frac{l_1l_2}{S}\right\rfloor,
      \left\lceil\frac{u_1u_2}{S}\right\rceil\right].   \tag{2}
```

Start from `floor(Sz),ceil(Sz)` and apply (2) in binary exponentiation.
The exponent-zero result is exactly `[S,S]`. Induction over the
multiplications and squares proves enclosure of every required power.
Adding the intervals with their nonnegative multiplicities encloses P
and Q. Equal exponents are computed once.

Since `0<z<1`, the scaled power endpoints remain between zero and S.
An extremely small power may have lower endpoint zero; its upper
endpoint is retained. The zero-order term guarantees that the complete
polynomials still have strictly positive lower endpoints.

The cost depends logarithmically on the largest exponent. It does not
construct a denominator whose bit length is proportional to that
exponent. The large final powers in (1) are also avoided.

## 3. Exact logarithm enclosures

For any positive rational x, choose an integer k by exact comparisons
such that

```math
1\le y=x/2^k<2,\qquad
t=\frac{y-1}{y+1}\in[0,1/3).
```

The ordinary power series for the logarithm gives

```math
\log y=2\sum_{j=0}^{J-1}\frac{t^{2j+1}}{2j+1}+\mathcal R_J,
\qquad
0\le\mathcal R_J\le
\frac{2t^{2J+1}}{(2J+1)(1-t^2)}.                       \tag{3}
```

The tail bound follows by replacing every later denominator by
`2J+1` and summing the geometric series. All terms and the tail bound
are rational. The same formula at `t=1/3` encloses `log 2`.
Now use `log x=k log 2+log y`, reversing the endpoints of the
`k log 2` interval when k is negative. Round the resulting rational
bounds outward to multiples of `1/S`.

Normalize each rational endpoint separately. Thus an input interval
that crosses a power of two needs no special approximation rule.
Monotonicity gives `log(P_lower) <= log P <= log(P_upper)` and the
corresponding Q inequalities.

If L_Q is the scaled lower bound for `log Q` and U_P is the scaled
upper bound for `log P`, the certificate accepts only when the integer

```math
A L_Q-B U_P
```

is strictly positive. Every rounding direction therefore weakens the
desired strict inequality. Equality or an unresolved interval fails
the gate.

## 4. The ninety-six-letter certificate

The exact parameters are

```math
A=710005874881191,\quad B=500000000000000,\quad
z=\frac{5571156672673465}{5571156672673466},\quad
b=192,\quad J=96.
```

The saved lower gap numerator is

```text
11169158249146545798852467616219265398210622320232377963976
```

over denominator `2^192`, so (1) holds strictly. The
[certificate record](../evidence/hot96-interval-pressure.json) stores
every polynomial and logarithm interval, the exact parameters, and
the source and complete-profile hashes. The arithmetic replay took
less than one second in the recorded environment. The separate full
profile construction took about 25 minutes; its cost is not included
in that arithmetic timing.

The rational parameters were proposed from a numerical search. Their
acceptance and the stated bound depend only on the exact test above.
No global minimum of the physical-code family is asserted.

```sh
python3 tools/certify_natural_pressure_interval.py \
  evidence/hot96-charge-profile.json \
  --rho 710005874881191 1000000000000000 \
  --base 5571156672673465 5571156672673466 \
  --output /tmp/hot96-interval-pressure.json
```

## 5. Independent controls and input scope

The [implementation](../tools/certify_natural_pressure_interval.py)
validates the identical paired downsets, their physical energies, the
total channel count, and the agreement between the saved block orders
and the flattened profile. It binds the complete input file by hash.
It does not infer a valid polynomial crossing from an arbitrary list
of integers; that input obligation is discharged by the separately
documented charge-block construction.

The [independent checker](../tools/check_natural_pressure_intervals.py)
uses exact rational powers and a different logarithm expansion. Its
alternating series for `log(1+u)` at `u<=1/2` supplies rigorous bounds
for products and reciprocals of small rationals. These exercise
negative k, normalization boundaries, and logarithm addition without
using the certificate's atanh-series recurrence.

The [control record](../evidence/natural-pressure-interval-controls.json)
contains 911 exact power enclosures, 125 independent logarithm
enclosures, 24 direct integer-power pressure sign comparisons, two
previous large-integer certificate controls, and seven rejected invalid
inputs. An [independent mathematical and executable audit](../evidence/hot96-interval-independent-audit.json)
also checked the tail, rounding, type normalization, and fresh fine
certificate replay. Its high-precision Decimal comparisons are controls;
the acceptance proof is the exact outward interval calculation.

```sh
python3 tools/check_natural_pressure_intervals.py \
  --output /tmp/natural-pressure-interval-controls.json
```

The earlier `49/69` and `27/38` certificates remain separate controls
using full denominator clearing. Their provenance and a command that
reconstructs all 9,216 polynomial channels appear in
[hot96-rational-pressure.md](hot96-rational-pressure.md).
