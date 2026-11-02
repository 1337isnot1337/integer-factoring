# A historical seven-letter asymmetric profile passes the older inverse bound

**Current status: excluded, not an active construction target.** The later
[adjugate compatibility theorem](adjugate-seven-letter-pressure-floor.md#4-two-additional-all-degree-constraints)
excludes this profile at every positive integer ramification. The stronger
[box-minor theorem](bare-seven-to-twelve-letters.md) excludes the entire bare
criterion through twelve physical letters, over every field and degree.
The old inverse-bound pass and exact formal pressure below remain valid
as historical relaxation calculations.

## Result and scope

The [shared-minor relaxation](shared-minor-smith-floor.md) identifies
`M=12` and `mu=nu=(0,7^6,13^42)` as the formal `q=7` profile optimizer.
The [inverse-degree theorem](inverse-degree-crossing-compatibility.md)
excludes that exact profile at every integer ramification, because `mu_2=7>2(13-12)=2`.

The nearby asymmetric profile below satisfies the older shared-minor and
inverse-degree conditions and crosses their relaxed profile threshold.
It fails the later adjugate condition.  Take

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

This is a historical relaxation point, not an algebraic crossing. No
stronger factoring score is established.

## The later constraint excludes every ramification

The [adjugate determinant inequality](adjugate-seven-letter-pressure-floor.md#4-two-additional-all-degree-constraints)
is necessary for an actual polynomial crossing:

```math
(n-1)\mu_2+n\sum_{i=1}^{n-1}\nu_i\le n(n-1)M.
```

For the displayed profile its left side is `28364`, its right side is
`28224`, and the excess is `140`. Multiplying every order and the degree
by any positive integer `r` multiplies that violation by `r`. Thus this
profile is impossible at **every positive integer ramification**, over
every field. Passing the older `mu_2<=2(nu_n-M)` test is insufficient.

The broader bare seven-through-twelve-letter exclusion is already proved;
changing the degree or realizing another bare profile in that range
cannot improve the score. Physical gauges, collective recombinations and
other transfer criteria require their own analysis.

## Verification

`python3 tools/check_seven_letter_asymmetric_profile.py` checks the
monotonicity, sums, shared prefix, older inverse-degree inequality, exact
`B^169<7^478` comparison, and the later adjugate violation of 140.  The [saved record](../evidence/seven-letter-asymmetric-profile.json)
reconstructs the same certificate.
