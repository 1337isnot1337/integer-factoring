# Single Laurent updates of a joint rectangular retraction

Every declared one-monomial update in this note has balanced natural
`rho >= 3/4` at all scales. The complete finite family consists of 14,400
row or column updates over F4, with Laurent powers from -12 through 12.
This changes the joint embedding or its retraction and permits unequal
energy cells. It is a larger family than constant homogeneous moves, but
the result does not cover combinations, arbitrary Laurent powers, other
coefficients or independently varied boundary charges. No stronger
factoring exponent is obtained.

## The maps and their pole accounting

Let `A(z)` be the existing literal 25-by-25 five-letter source operator,
with `A(1)=I`. Let `P0,L0` retain the four-letter Cartesian grid of physical
heights `(0,1,3,4)`, so `L0 P0=I_16`. There are nine discarded ambient
coordinates. For a retained target coordinate `i`, discarded coordinate
`a` and `f(z)=c z^k`, either use

```
P = P0 + f e_a e_i^T,       L = L0,
```

or use

```
P = P0,                    L = L0 + f e_i e_a^T.
```

Both satisfy `LP=I` identically, since the discarded coordinate is absent
from the retained maps. They need not factor into separate physical-leg
maps. The strict output is `C=LAP`, and `C(1)=I` because the original source
is identity there. In the first case only column `i` changes; in the second
only row `i` changes.

The [constructor](../tools/joint_moving_retractions.py) explicitly clears
`z^-m`, where `m=max(0,-k)`, and keeps the clearing exponent. It constructs
both actual flattenings of the polynomial numerator. For the original
Laurent output, its zero orders are the numerator orders minus `8+m`,
and its infinity orders are the numerator orders plus `8+m`. Thus the
common clearing has opposite, cancelling effects in the balanced product.
Poles are retained in the profiles; they are not discarded as free gains.
The endpoint supplies the generic full rank for both 16-by-16 cuts.

## Complete finite family and exact scale certificates

The [screen](../tools/search_joint_moving_retractions.py) covers every
combination of:

- 16 target coordinates and nine discarded source coordinates;
- either the row or column update;
- every integer power from -12 through 12;
- scalar `1` or `omega` in F4.

This gives `16*9*2*25*2=14,400` cases. The remaining nonzero scalar
`omega+1` is the coefficientwise Frobenius image of `omega`; the source is
binary, so its full Smith profiles coincide. The
[complete receipt](../evidence/joint-moving-monomials.json) contains every
case, all 1,745 distinct complete profile pairs, source hashes and the
finite numerical scale grid. The sampled minimum rho is approximately
0.750861; this sample value is not a global minimum certificate.

For each profile pair, define `Z0(x)=sum x^nu` and `Zinf(x)=sum x^mu`,
including negative exponents when present. Let `F(x)=1+x+x^3+x^4`.
The balanced natural floor `rho >= 3/4` follows from

```
Z0(x^3) Zinf(x^3) - F(x^4)^3 > 0,        0 < x < 1.
```

The [checker](../tools/check_joint_moving_retractions.py) removes the
leading Laurent monomial and proves positivity by integer Bernstein
coefficients for every class. The largest required Bernstein degree is
451; the smallest final integer coefficient over the family is one.
These are exact all-scale inequalities for the declared profiles.

## Literal reconstruction and verification limits

A separate construction starts with the unreshuffled `five_letter_source`
array and explicit sparse polynomial maps `L,P`. It checks `LP=I`, computes
every entry by polynomial multiplication and restores the monomial clearing.
All 14,400 cases are reconstructed this way, giving 3,686,400 endpoint
entries and 460,800 complete native weighted orders. Every case's actual
profile is checked against its class assignment.

For every one of the 1,745 classes, separate Python fraction-free
polynomial elimination reproduces both full profiles: 55,840 independently
computed orders and 446,720 tensor entries. Thirty extra map controls cover
negative, zero and positive Laurent powers and all three nonzero F4
scalars; they also check the Galois relation directly. The
[replay receipt](../evidence/joint-moving-replay.json) records these counts
and the common hash of the regenerated integer certificates.

Run the focused exact replay with:

```
make verify-joint-moving-retractions
```

Regenerate the finite discovery receipt with:

```
python tools/search_joint_moving_retractions.py --output /tmp/joint-moving.json
```

The focused replay proves neither optimality over general retractions nor
a bound for asymmetric charges. It postdates the earlier full-suite source
pin at `decce64`; no new full-suite pass, Lean result or independent
mathematical review is claimed. The factoring score remains `1/2`.

Combinations of discarded rows are a distinct next search: several terms
may cancel poles or low-order minors that a single monomial cannot cancel.
The existing [rational row-module method](joint-hot-rational-retraction-module.md)
provides a way to enumerate polynomial outputs of such rational maps,
without imposing an arbitrary pole-location cutoff.
