# A complete q96 hot-code profile and a rational pressure certificate

The first 96 hot masks in the code specified below give the exact bound

\[
\rho<\frac{49}{69}=0.710144927536\ldots .
\]

This improves the earlier [q64 rational bound](hot64-rational-pressure.md)
\(32/45\). It remains above \(1/\sqrt2\), so this result does not change
the proved factoring exponent or the competition score. Numerical
optimization of the same complete profile gives
\(\rho\approx0.7100058748811902\); that decimal is not the assertion of
this rational certificate.

## 1. The fixed physical code

For a twelve-bit mask \(a\), set

\[
E(a)=\sum_{j:a_j=1}(8j+4).
\]

Order the masks in \(\{0,\ldots,2^{12}-1\}\) by \((E(a),a)\), and retain
the first 96. Every bit deletion strictly lowers energy, so this is a
downset. The maximum retained energy is 100. The fixed reflection interval
\(M=32\) has enough modes and satisfies \(100<4M\), so the
[stabilized realization](stable-hot-reflection-codes.md) applies.

The [charge-block constructor](hot-code-charge-profile-constructor.md)
recovers all \(96^2=9,216\) common natural orders. The block sizes are

\[
1,24,216,932,2136,2738,2044,888,212,24,1.
\]

Their total is 9,216; the largest block is 2,738. The complete nonnegative
orders have minimum 0 and maximum 240. The full construction took
approximately 1,517.281 seconds in the recorded run. No unresolved or
partially computed block is included.

The frozen complete profile is
[`hot96-charge-profile.json`](../evidence/hot96-charge-profile.json), with
file SHA256
`0623bf18f2dab0766c9137ccf05b41c8cd54a9c7be8118a74c3d918b023865d4`.
It contains the labels, energies, every complete order, and every charge
block. The profile construction uses the exact field and polynomial
operations of the already controlled charge engine; no real-number pivot
or numerical rank tolerance is used.

## 2. Exact integer comparison

Write

\[
P(x)=\sum_{i=1}^{9216}x^{\sigma_i},\qquad
Q(x)=\sum_{a\in C}x^{E(a)}.
\]

Choose

\[
z=\frac{768}{769},\qquad R=\frac{98}{69}.
\]

The checker proves, after exact denominator clearing,

\[
P(z^{98})^{138}<Q(z^{69})^{196}.
\tag{1}
\]

The two compared positive integers each have 31,116,901 bits. The
certificate stores their bit lengths, SHA256 digests, and the digest of
their strict positive difference. It also stores the positive integer
type-count formula and its exact normalizing sum. The large decimal
integers are not needed in the repository because the checker reconstructs
them from the small rational parameters and the complete order list.

For completeness, an independent smaller comparison with
\(z=211/212\) and \(R=27/19\) also passes. It uses two 1,903,011-bit
integers and proves the weaker bound \(\rho<27/38\). Both comparisons are
in [`hot96-charge-rational-pressure.json`](../evidence/hot96-charge-rational-pressure.json).

## 3. Entropy and uniform transfer

Use the positive rational physical type

\[
p(a)=\frac{z^{69E(a)}}{Q(z^{69})}
\]

on each paired alphabet, and set \(t=98\log(1/z)\). Then

\[
H(p)=\log Q(z^{69})+\frac{t}{R}\mathbb E_p E.
\]

Taking logarithms of (1) gives
\(\log P(z^{98})<R\log Q(z^{69})\). The full natural rank numerator is
\(2\log P(z^{98})+2t\mathbb E_pE\), and its denominator is \(4H(p)\).
The two physical means therefore give the exact bound \(\rho<R/2=49/69\).

The code, field, and ambient interval are fixed constants. Rational type
denominators can be cleared at admissible block lengths. The strict
inequality absorbs the polynomial coefficient and reconstruction factors,
as in the [q64 proof](hot64-rational-pressure.md) or the
[finite-block amplification theorem](finite-block-crossing-amplification.md).
This uses the full physical entropy and both complete natural profiles.

## 4. Reproduction and provenance

The full independent replay command reconstructs every charge block before
checking the integer comparison:

```sh
OPENBLAS_NUM_THREADS=1 python tools/certify_hot_charge_pressure.py \
  --size 96 --modes 12 --base 768 769 --support 98 69 \
  --output /tmp/hot96-full-replay.json
```

The [certificate engine](../tools/certify_hot_charge_pressure.py) is unchanged
from the q64 checkpoint. The recorded q96 certificate reused the complete
profile from the earlier exact construction and performed a fresh integer
comparison. Its provenance reports the original profile runtime and hash;
it does not describe that reuse as a second independent large elimination.
The final search archive also preserves the original constructor record.

The floating optimum remains useful for choosing rational parameters.
Only the strict integer inequalities above are asserted here.
