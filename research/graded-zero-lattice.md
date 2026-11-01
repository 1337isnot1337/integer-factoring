# A constructive graded lattice with zero profile beginning 8,10,10

The retained five-letter crossing admits an explicit polynomial change
whose normalized zero profile begins **8,10,10**, improving the original
8,9,11 prefix. Its actual opposite profile becomes substantially worse.
The complete paired pressure is approximately **0.9964095025**, so this
construction does **not** improve the factoring score, which remains 1/2.

This is a finite explicit construction beyond the simple-pole pencil,
not a contradiction to that pencil's polynomiality restrictions. The
source crossing, joint grading and matrix layout are inherited from
[joint energy conjugations](joint-energy-conjugations.md). The
[graded Newton ceiling](joint-graded-newton-ceiling.md) remains the
relevant upper target; the new construction does not attain its whole
profile.

## Construction and complete costs

Write the source joint operator as `A(t)` on rows `x+5*y` and columns
`v+5*u`, with joint energy `H`. The leading coefficient of `t^H A(t)`
at order eight is a rank-one idempotent `Pi`. Put `Q=I+Pi`.
The next compressed coefficient `Q R_9 Q` is rank one and square zero.
Its image and covector coincide, with two nonzero coordinates in each
of the energy-one, energy-five and energy-nine blocks.

In each affected grade choose a constant basis whose first column is
that image vector, second column is an original coordinate pairing to
one with the covector, and remaining columns complete the basis.
Use the identity on the other grades; call the resulting binary matrix
`W`. Give its first two coordinates potentials `+1,-1` in each changed
grade, and all other coordinates potential zero. The constructor
checks invertibility, grade compatibility and every endpoint entry.

With these potentials `g`, define

```math
K(t)_{r,c}=t^{2+g_r-g_c}(W^{-1}A(t^2)W)_{r,c}.
```

The outer exponents are all nonnegative. Thus `K` is polynomial, of
actual degree 38, and `K(1)=I`. Ramification doubles physical energies;
the common monomial shifts the two cut orders by opposite constants
and cancels from the paired pressure. Dividing out these bookkeeping
changes gives the complete zero profile

```text
8,10,10,12,13,14,14,16,16,17,18,18,18,19,19,
20,20,20,21,22,22,22,23,24,28.
```

The complete raw infinity profile, before undoing ramification and
clearing, is

```text
-20,-20,-20,-19,-19,-19,-19,-19,-19,-18,-18,-18,-17,
-17,-17,-17,-16,-16,-15,-15,-11,-9,-2,4,6.
```

Its normalized orders are `(raw+2)/2`. The exact raw zero orders are
`2*normalized+2`, and their sum is `2*444+25*2=938`.
Both complete profiles are independently computed by untruncated
packed-F4 fraction-free polynomial elimination. Binary coefficients
embed in every characteristic-two extension without changing these
minor valuations.

Substituting the *old* infinity profile would yield a numerical value
near `0.6944141674`. That profile belongs to a different kernel and is
an unrealized target; it cannot be used in the construction's cost.
The actual paired value is the poorer `0.9964095025` stated above.

## Polynomial units preserving the improved zero profile

Let `N` be any energy-preserving square-zero matrix over F4, and `k>=0`.
Then `U=I+t^k N` and its inverse are polynomial units. Conjugating `K`
by `U` preserves the entire weighted zero Smith profile: `U` commutes
with the row grading, and both left and right factors are units at zero.
The strict endpoint identity is also preserved. The opposite profile
can change. The [constructor](../tools/graded_zero_lattice.py) supplies
this family with exact energy, square-zero and endpoint checks.

A registered finite screen tested every triple of square-zero two-by-two
F4 matrices on the three declared subspaces, except the all-zero triple,
at powers zero through four: `(16^3-1)*5=20,475` candidates. Every complete
zero profile matched the target. The screen found 523 complete opposite
profiles; its best numerical paired value was `0.9930831371`. This is
not an all-field or all-unit obstruction. The
[packet](../evidence/graded-zero-positive-dressings.json.xz) retains the
registration, every result, driver and independent complete replay of
both profiles of the best candidate.

An additional [physical retraction packet](../evidence/graded-zero-physical-retractions.json.xz)
records all 496 binary rank-four split retractions of the five-letter
space, applied to each of the reference and new kernels. All 992
declared cases completed. Their 975 distinct cost profiles gave no
improvement; the best numerical value was the reference `0.7310809360`.
Inherited physical heights and only this finite binary family were tested.

## Reproduction and limits

Run `make verify-graded-zero-lattice`. The
[checker](../tools/check_graded_zero_lattice.py) verifies 625 strict
endpoint entries, both complete profiles and two nonbinary zero-preserving
dressings, totaling 100 independently computed orders. The retained
[receipt](../evidence/graded-zero-lattice.json) also includes the optional
numerical pressure comparison (`--numerical`).

A full default pass at [2528317](../evidence/verify-2528317.json) predates
this construction and is retained with compressed stdout. It does not
verify later source changes. No new Lean formalization, globally optimal
pressure or improved uniform factoring exponent is claimed. Substantial
OpenAI assistant support was used for the construction, code and text.
