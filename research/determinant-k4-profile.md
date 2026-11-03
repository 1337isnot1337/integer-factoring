# A complete k=4 determinant-code profile stays above the critical target

## Result and scope

The diagonal retraction of the k=4 determinant-minor code has both complete
natural Smith profiles resolved: 4,900 orders in each cut, zero unbounded
tails at cutoff 511, and largest orders 450 (zero cut) and 186 (infinity
cut). The numerical natural-scale fit is

```math
\rho \approx 0.7394802516174483
```

at scale `0.7413880836`. This is a complete finite diagnostic, not a global
scale optimum, an all-scale floor, a crossing certificate, or a factoring
improvement. It remains above `1/sqrt(2)`, so the complete exponent remains
**1/2**.

The value extends the registered k=2,3 panel (`0.73964`, `0.73954`) by one
exact case. The three sampled values are nearly flat, so this does not
suggest that raising k alone approaches the threshold. The k=4 case is also
the largest sparse map the prescribed F8 Cauchy construction supports.

## 1. What was computed

Run `make verify-determinant-inverse-k4` for the source-pinned replay and
selected entry controls; run `make verify-determinant-inverse-k4-full` to
reconstruct all 9,800 orders again (about 24 minutes in the recorded
environment). The [complete profile receipt](../evidence/determinant-inverse-k4-diagonal-profile.json)
retains:

- both complete order lists, all 70 labels and their 16 physical heights;
- all charge blocks, their sizes, multiplicities and retained orders;
- the numerical fit and its physical probabilities;
- 130 selected complete polynomial entries from every charge and both cuts,
  checked against determinants of the original one-particle matrix;
- source fingerprints for the constructor, native F8 engine and entry oracle.

The selected entries are not a second reconstruction of all 4,900 matrix
entries. The complete profile was produced by the same native engine used
for the registered panel; the new receipt records its exact source pin.

## 2. Relation to the existing panel

The [determinant-code note](determinant-inverse-codes.md) registered the
k=2,3 cases and explicitly declined a bare k walk. This checkpoint adds the
last F8 case because the sparse map constructor stops at k=4 and because a
complete profile, not another near-critical fit, is needed to see whether
the sampled pressure is moving. It is not moving: `0.73964`, `0.73954`,
`0.73948`.

Changing the retraction at k=3 moved the sampled value by less than
`3e-4`; the k=4 diagonal case is therefore representative evidence for the
family, not a global optimum over retractions. Larger k needs a larger
field and a different acquisition/height regime, and remains outside this
finite statement.

## 3. Consequence for the search

The determinant-minor family is not a promising route to a crossing below
`1/sqrt(2)` at fixed small k. A continuation should change the physical
map space, the kernel, or the protected low structure, rather than enlarge
the k panel again. This checkpoint does not prove that the family is
excluded asymptotically; it records the exact finite evidence and stops
the local k walk.

Credit: classical determinant derivative/minor spaces; the existing
coherent inverse-system, dual-space and F8 compound constructors. This
checkpoint was prepared with OpenAI assistant support. No new Lean
formalization or independent mathematical acceptance is claimed.
