# A seven-letter asymmetric profile survives inverse-degree compatibility

## Result and scope

The [shared-minor relaxation](shared-minor-smith-floor.md) identifies
`M=12` and `mu=nu=(0,7^6,13^42)` as the formal `q=7` profile optimizer.
The [inverse-degree theorem](inverse-degree-crossing-compatibility.md)
excludes that exact profile in the quadratic-inverse family at every
integer ramification, because `mu_2=7>2(13-12)=2`.

A nearby asymmetric profile satisfies every published necessary
condition, including the inverse-degree bound, and still crosses the
profile threshold.  Take

```math
M=12,\qquad
\mu=(0,7^6,13^{42}),\qquad
\nu=(0,7^6,12^3,13^{38},16).
```

Both lists have 49 entries and are nondecreasing.  Their sums are
`49*12=588`; the first seven entries of each sum to 42, so the shared
minor prefix saturates `7*12=84`.  The inverse-degree condition is
satisfied with one unit of slack:

```math
\mu_2=7\le 2(\nu_{49}-M)=2(16-12)=8.
```

At `s=34/25`, the profile pressure is

```math
B(s)=\left(\sum_i s^{M-\mu_i}\right)
     \left(\sum_i s^{-\nu_i}\right)
   =244.848867\ldots .
```

The exact certificate is

```math
B(34/25)^{169}<7^{478},
\qquad
478^2<8\cdot169^2,
```

so `B(34/25)<7^{478/169}<7^{2\sqrt2}`.  A finer sample near
`s=1.3618` gives `B=244.835361...`, below the threshold
`7^{2\sqrt2}=245.639508...`.  Thus this profile has
`rho=log B/(4 log 7)≈0.70669<1/sqrt2`.

This is a necessary-condition result, not an algebraic crossing.  No
tensor, auxiliary matrix, or polynomial endpoint is constructed, and no
factoring score is claimed.  It does show that the inverse-degree
exclusion is not a blanket obstruction to the seven-letter window: the
remaining question is realization of an asymmetric profile, not the
profile inequalities alone.

## What remains

The exact optimizer was excluded by a one-unit inverse-degree violation.
The present profile uses the allowed slack by lowering three zero
orders from 13 to 12 and raising the top zero order from 13 to 16.  A
construction must realize this asymmetry, or another profile in the same
feasible region, with an actual polynomial crossing.  The quadratic
inverse family with the displayed `q=7` transpose channel has degree 24
and does not directly supply `M=12`; a shorter kernel flag or a
higher-degree inverse with a rank-one first jet remains to be built.

## Verification

`python3 tools/check_seven_letter_asymmetric_profile.py` checks the
monotonicity, sums, shared prefix, inverse-degree inequality, and the
exact `B^169<7^478` comparison.  The [saved record](../evidence/seven-letter-asymmetric-profile.json)
reconstructs the same certificate.
