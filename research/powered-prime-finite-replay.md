# Finite replay of powered-prime norm witnesses

The [powered-prime theorem](powered-prime-smooth-norms.md) has a
separate finite replay that produces actual primitive norm witnesses
and then recovers their coordinates using the independent
[orientation decoder](factored-norm-orientation-search.md).

The checker deliberately enumerates the full class groups at
discriminants `-41612` and `-300004`. This is a small-instance
verification cost, not the proposed norm-selection algorithm.
Its prime cutoff is `32`. The smooth sets use one or two distinct
prime factors, respectively, and retain only products below a
radius `R` satisfying `4R^2 < abs(Delta)`.

It constructs the exact rational prime distribution and convolves it
until its exact `l1` distance from uniform is at most one quarter
of the smooth set's density. This directly bounds every nonprincipal
character bias without floating-point Fourier calculations. The
positive correlation criterion then applies to every target class.

For each target, the checker records a product-walk decomposition
into two smooth classes and a prime walk, cancels opposite prime
orientations, constructs the corresponding norm, and supplies only
that factored integer to the orientation decoder. The recovered
coordinates must satisfy the original target form exactly and have
gcd one. It also checks coprimality to the discriminant, the prime
cutoff, the multiplicity budget, and both norm-height bounds.

The first case includes target classes missed by the one-prime
version at the same finite cutoff. Powering repairs that coverage.
The second case checks the same construction with two-factor
squarefree starting norms.

Run:

```sh
python3 tools/check_powered_prime_norm_witnesses.py
```

The [checker](../tools/check_powered_prime_norm_witnesses.py) and
[evidence](../evidence/powered-prime-norm-witnesses.json) record
the exact cases, walk lengths, coverage counts, largest norm sizes,
and representative primitive coordinates. These finite cutoffs are
certified by complete class enumeration. They are not substituted
for the much larger uniform analytic cutoff in the theorem.

This note used substantial OpenAI assistant support. The replay
establishes no stronger complete factoring score.
