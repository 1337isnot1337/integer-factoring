# Recombination before orientation leaves the Gaussian width regime

## Status

The adaptive global-pairing step is rejected in the
[global-pairing audit](two-axis-global-pairing-repair.md).  The
[favorable-reweighting diagnostic](cube-favorable-reweighting.md) excludes
all literal favorable orientations of the fixed `m=5`, weight-1 Gaussian
coefficient family.  A remaining question is whether **new recombinations
before orientation** can restore the endpoint while keeping a
maximum-width-hiding orientation.

For the same fixed single-block family, the answer is: the endpoint can be
restored, but only by a recombination that leaves the `5/8` width regime.
The complete factoring score remains `1/2`.

## Exact result

Use the 22 coefficient tensors `A_i`, the endpoint `E4`, and the three
fixed parent pairings `P_c` from the audit.  Each `P_c` is an invertible
permutation of the eight leg labels.  Hence for each `c` the equation

~~~
P_c(R)=E4
~~~

has the unique preimage `R_c=P_c^{-1}(E4)`.  The checker verifies that the
three preimages coincide at one recombination `R`, and that `R` is the XOR
of 11 of the 22 original coefficient tensors.  The receipt records the
explicit exponent subset.

The exact cut ranks of `R` are

~~~
(r_0,r_1,r_2)=(16,16,16),
~~~

so the exact channel widths are `(4,4,4)` and

~~~
S_R = 4+4+4 = 12.
~~~

At `lambda'=2`, the inherited Gaussian budget is

~~~
6 rho lambda' = 6*(5/8)*2 = 15/2 = 7.5,
~~~

so `S_R>15/2`.  The per-cut cap `2 lambda'=4` is tight but not violated.
The effective one-piece ratio is `S_R/(6 lambda')=1`, not `5/8`.

The favorable two-entry dual `L(T)=T[240]+T[255]` has value one on both
`E4` and `R`; it does not separate this recombination.  The width sum is
the obstruction, not the endpoint identity.

## Scope

- **Shown:** in the fixed 22-piece family, every recombination that one
  fixed parent pairing sends to `E4` is this full-rank recombination `R`.
  The endpoint is restored, but the Gaussian pressure gain is lost.
- **Not shown:** larger native alphabets, other base families, coupled
  blocks, or a new adaptive identity are impossible.  This is a finite
  single-block comparison, not an asymptotic theorem.
- No stronger complete factoring exponent, Lean result, or independent
  mathematical acceptance is claimed.

## Checks

Run

~~~
make verify-recombination-before-orientation
~~~

The [checker](../tools/check_recombination_before_orientation.py)
reconstructs the actual coefficient family through the rejected-pairing
audit, verifies the three fixed pairings of the original pieces sum to
`E4`, that the three preimages coincide, that `R` lies in the span of the
original pieces, and the exact ranks, widths and budget comparison.  The
[receipt](../evidence/recombination-before-orientation.json) records the
recombination subset and all controls.
