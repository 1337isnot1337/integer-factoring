# Degenerate asymmetric scaling of complete rectangular two-cut profiles

## Result and scope

The [constant-companion search](constant-companion-reflection-search.md)
recorded `1.001436002791` as the best asymmetric support value on the
complete q4 rank-one frame catalogue.  That search constrained the
supporting angle to `0.03 <= theta <= pi/2 - 0.03` and the gauge scale to
`0.0005 <= log(s) <= 2`.

Replaying the winner with those bounds relaxed shows that the recorded
value was bound-constrained.  Every rectangular `A,B,B,A` complete profile
has the exact zero-scale surface

```math
R(\log s\to0,\theta)=\cos\theta+\sin\theta .
```

The q4 finite-scale values approach the corner
`(theta,log(s))=(0,0)` from above.  The smallest tested value is
`1+1.37e-9`; no tested point is below one.  This note records a numerical
diagnostic.  It does not prove a uniform floor, and it does not change the
certified balanced pressures.

The concurrent [q96 unequal-scale screen](q96-asymmetric-scale-screen.md)
reports the same boundary layer for six complete 96-letter profiles: its
best interior value is `1.0037692355`, and its boundary values decrease
toward one as the supporting angle tends to zero.  That screen treats the
near-one values as a one-cut boundary artifact.  The present note adds the
exact zero-scale identity for the rectangular family, the q4 rank-one
catalogue replay, and a rank-two square-zero screen.

## 1. The recorded asymmetric winner

The best recorded case is catalogue generator 4035, constant one, pole
four.  Its stored parameters have `theta = 0.030000000000000023`, exactly
the lower angle bound, and `log(s) = 0.002712151136578734`.  The stored
support value is `1.0014360027905598`.  The second recorded case,
generator 63, constant one, pole eight, has `1.0014407204155154`.

Relaxing only the angle, with the stored gauges and scale, gives
`1.0009658301049553` at `theta = 0.01543234352330092`.  A full
gauge-and-type re-optimization with `theta` allowed down to `1e-6` gives
`1.000003353164649` at `theta = 4.765350789642891e-6`.  Continuing with
`theta = 1e-8` and then scaling the whole gauge down gives the values in
Section 3.  The prior lower angle bound, not the profile, set the scale of
the recorded `1.0014` value.

## 2. The zero-scale identity

Let the seed have alphabets `A,B,B,A`, so each cut has `n=AB` Smith
channels.  The letter alphabets contain `A+B` letters and `AB` ordered
pairs.  Send every physical height to zero, i.e. `log(s) -> 0` with the
relative gauge ratios fixed.  Every Smith order and every letter energy
tends to zero.  The support-value equation is then

```math
a\log n+b\log n
 = R\bigl(\log A+\log B\bigr),
\qquad a=\cos\theta,\ b=\sin\theta .
```

Because `\log A+\log B=\log n`, this reduces to

```math
R=a+b=\cos\theta+\sin\theta .
```

The identity is independent of the seed tensor, the field, the pole order,
and the choice of cut.  It uses only the rectangular alphabet product.
Since `cos(theta)+sin(theta) >= 1`, with equality only at `theta=0`, the
zero-scale surface is a universal lower boundary at one.  The physical
pressure is the finite-scale correction below this surface.  The correction
can be negative: at `theta=0.03`, `log(s)=0.0027`, the recorded winner is
`1.0014`, below `cos(0.03)+sin(0.03)=1.0295`.  The open question is whether
the correction can drive the value below one.

## 3. Numerical collapse in the q4 family

The cleanest tested case is generator 5056, constant one, pole two.  The
gauge ratios are fixed from the limit optimization; `theta=1e-10`.  The
support value satisfies `R-1 = delta*log(s)^2` with `delta = 8.786105`.

| `log(s)` | `R` | `(R-1)/log(s)^2` |
|---:|---:|---:|
| `2e-4` | `1.0000003514430602` | `8.786` |
| `1e-4` | `1.0000000878650024` | `8.787` |
| `5e-5` | `1.0000000219678276` | `8.787` |
| `2.5e-5` | `1.0000000054926776` | `8.788` |
| `1.25e-5` | `1.0000000013735213` | `8.791` |

The quadratic coefficient is positive and stable.  The extrapolated
intercept is one; the finite-scale value is above one at every tested
point.

Generator 4035, constant one, pole four, gives the same qualitative
behaviour but a less clean quadratic fit at the smallest scales:

| `log(s)` | `R` |
|---:|---:|
| `2e-4` | `1.0000000244048803` |
| `1e-4` | `1.0000000076146291` |
| `5e-5` | `1.0000000026597269` |
| `2.5e-5` | `1.0000000010431750` |
| `1.25e-5` | `1.0000000004498815` |

The first three values fit `delta = 0.620654`; the later values are at the
level of the oracle's numerical resolution.  All are above one.

The zero-scale identity was checked independently at `log(s)=1e-8` for
both cases and for angles `1e-10`, `1e-4`, `1e-3`, `1e-2`, `0.03`, `0.1`,
and `pi/4`.  The maximum deviation from `cos(theta)+sin(theta)` was
`5.77e-8`.

A bounded re-optimization of the top twenty recorded asymmetric
refinements at `theta=1e-8` and `log(s)>=0.0005` found limits between
`1.000001088287796` and `1.0001074779055672`.  Every limit was above one.
This scan is exploratory and is not replayed by the checker.

## 4. Rank-two square-zero frames

The rank-one catalogue is the complete `G^2=0`, `rank G=1` family.  The
next square-zero class over `F4^4` is the nilpotent type `(2,2)`, with
`rank G=2`.  A random screen of forty such frames, conjugated by random
invertible `F4` matrices, gave natural profile values in `[0.8597,0.9411]`
and full-gauge asymmetric support values above `1.11`.  The best natural
value was `0.8597077729600762`; the best asymmetric replay was
`1.1132973060388962`.  Both are worse than the rank-one family, whose best
balanced value is about `0.7311`.

This is a bounded random screen, not a classification.  It shows that
moving from rank one to the full square-zero class does not by itself
produce a better asymmetric frame.  The rank-one family remains the best
sampled class.

## 5. Consequences and next test

The recorded `1.001436` should be read as a bound-constrained search
observation.  The new evidence puts the infimum of the tested q4
asymmetric surface at one, approached from above in the degenerate corner.
There is no strict crossing in this family.

The zero-scale identity uses only the rectangular product `n=AB`.  A
strict crossing requires changing that counting, not rescaling the same
profile.  The [q96 unequal-scale screen](q96-asymmetric-scale-screen.md)
varies the two cut scales `(a,b)` but keeps the rectangular `96,96`
channel product; the identity applies and its near-one values are the
same corner.  Unequal scales are not unequal channel counts.  Concrete
escapes are:

- a cut with `n0 != A B`, i.e. unequal channel counts or a non-rectangular
  leg partition;
- a nonzero order offset that survives the zero-scale limit, rather than
  a correction proportional to `log(s)`;
- a different protected core or a joint boundary compression that changes
  the channel algebra.

The existing rectangular oracle already exposes both pairings as `cut=0`
and `cut=1`; both still have `AB` channels, so the cross-pairing alone does
not change the zero-scale count.  A genuine escape needs a `3+1` or `1+3`
partition whose Smith matrix is rectangular.

The next experiment should construct one explicit unequal-channel or
non-rectangular two-cut profile and evaluate its support value.  A tested
value below one would be a genuine crossing candidate; a value above one
would sharpen the rectangular degenerate floor into a statement about the
cut algebra.

## Reproduction

```sh
python3 tools/check_asymmetric_scaling_floor.py \
  --output evidence/asymmetric-scaling-floor-check.json
```

The checker reconstructs both frames, verifies their tensor hashes,
replays all scaled rows, checks the zero-scale identity within `2e-6`, and
checks the positive quadratic coefficient.  The evidence is
[asymmetric-scaling-floor.json](../evidence/asymmetric-scaling-floor.json).
The source results are
[constant-companion-reflection-search.md](constant-companion-reflection-search.md)
and [dense-rank-one-coupled-frames.md](dense-rank-one-coupled-frames.md).
The asymmetric transfer criterion is
[asymmetric-support-line-crossings.md](asymmetric-support-line-crossings.md)
and [asymmetric-modular-transfer.md](asymmetric-modular-transfer.md).
