# Integer targets of degrees 22 and 24 excluded by box minors

**Status correction (2026-11-01).** Both profiles below are impossible for
bare polynomial crossings. They pass the newer cofactor constraints listed
here, but each violates sixteen of the earlier
[box-minor budgets](bare-seven-to-twelve-letters.md). The original package
omitted those checks. The exact formal rate calculations remain valid for
the proposed lists, and are retained as historical relaxation calculations.
They are not open polynomial construction targets. The complete factoring
score remains `1/2`.

There is an integer twelve-letter profile point satisfying every declared
cell/shared, inherited square adjugate, rectangular first-jet, multichannel,
second-jet and [higher-word constraint](higher-jet-word-caps.md). Its formal
bare rate at the rational scale `s=17/14` is **less than `707/1000`**, hence
less than `1/sqrt(2)`.

A smaller degree-22 point passes the same constraints and has formal rate
below `7071/10000`, also less than `1/sqrt(2)`. The degree-24 point has
the stronger certified rate margin; neither is a polynomial tensor.

These are points of an incomplete spectral relaxation. **No bare polynomial
tensor can have either profile pair**, by the box budgets below. They are
different from the preceding twelve-letter target excluded by the word
theorem, but changing that target did not evade the older box obstruction.

## 1. The exact integer point

Use degree/deficit `M=24`, zero paired energies on both twelve-letter
alphabets, and identical proposed 144-order profiles. In run notation:

```
0,6,9,11,13,14,15^5,16,19^11,20,22^11,23^7,24^17,
25^30,26^13,27^11,28^29,29,39.
```

The superscripts denote multiplicities. The full lists are retained in the
[point](../evidence/q12-integer-word-target.json). Their first twelve sum to
144, and their total is `3449<=144*24=3456`.

The [checker](../tools/check_q12_integer_word_target.py) verifies the original
cell/shared caps, 1,452 oriented square caps, 6,578 rectangular first-jet
caps, 6,578 multichannel caps, 1,452 second-jet caps and 1,452 higher-word
caps. Both rectangular cut orientations are identical for this equal pair.
Their minimum integer slacks are respectively `22,0,1,22,4` for the five
compiled families. The nonlinear word and second-jet functions are
evaluated separately on the two complete proposed profiles.

## 2. Rate certificate using integers only

Let `A=39` be the largest order and define

```
P = sum_i 17^(A-r_i) 14^(r_i).
```

For the restored formal bare product,

```
B(17/14) = (17/14)^24 (sum_i (17/14)^(-r_i))^2
         = P^2 / (17^(2A-24) 14^24).
```

The checker verifies the strict integer inequality

```
P^500 < 12^707 * 17^(250*(2A-24)) * 14^(250*24).
```

Therefore `B(17/14)<12^(707/250)` and
`rho_formal=log(B)/(4 log12)<707/1000`. The threshold comparison is exact:
`2*707^2=999698<1000000`. The compared integers have 80,559 and 80,560 bits; their
hashes are retained in the [receipt](../evidence/q12-integer-word-target-controls.json).
No floating arithmetic or optimization is used in this certificate.

## 3. Discovery and the missing compatibility constraint

The point is an incumbent from a bounded integer search using an inner
linearization of the word-prefix functions. A selected set of `k` words
gives a linear upper bound on the minimum `k`-word sum; satisfying the
selected linear constraints is sufficient for those nonlinear conditions.
The solver reached its 60-second time limit. Neither integer optimality
nor a minimum possible degree is claimed. The point is verified independently
of the solver and of this discovery heuristic.

The earlier box theorem already excludes the degree-24 point. For a
`2 by 12` box, let `N=24` and `r=4`. Its coefficient vector on one
profile is `(4,4,4,4,4,4,0,...)`, while the other contributes its first
24 orders. Thus the required inequality is

```
4 * sum_{i=1}^6 mu_i + sum_{i=1}^24 nu_i <= 24 M.
```

For this point the left side is `585` and the right side is `576`.
Both cut orientations fail. The boxes `(alpha,beta)=(a,12)` for
`a=2,3,4,5,7,8,9,10` fail in both orientations, giving sixteen exact
violations. The checker now evaluates all 156 oriented box budgets and
records every violation, in addition to replaying the original cofactor
and integer-rate checks. The zero paired energies and degree convention
here agree with those of the box theorem; no gauge conversion is needed.

More generally, the older entropy-dual certificates exclude the desired
bare rate for **every** twelve-letter profile, at every degree and field
under that theorem's assumptions. Enlarging the degree or finding another
point of the incomplete relaxation cannot evade it. The thirteen-letter
relaxation in the older note is not covered by that exclusion; actual
realization and any additional compatibility constraints remain separate
obligations there.

Run `make verify-q12-integer-word-target`. This point package has a separate
focused replay; earlier full-suite source pins do not establish its checks.
No new Lean work or independent mathematical review is claimed.
Substantial OpenAI assistant support was used.

## 4. A smaller degree-22 alternative

The [smaller point](../evidence/q12-degree22-word-target.json) has `M=22`,
the same zero energies, and equal profiles

```
0,6,9,10,11,12,14^6,17^5,18^9,19^7,20^3,21^12,22^12,
23^28,24^14,25^12,26^25,27^4,36.
```

Its first twelve sum to 132 and its total is `3162<=3168`. All originally
declared cofactor constraints pass: the square/first-jet/multichannel/second-jet/word minimum
slacks are `20,0,0,14,4`. At `s=57/46`, put
`P=sum_i 57^(36-r_i)46^(r_i)`. The independent checker verifies

```
P^5000 < 12^7071 * 57^(2500*(72-22)) * 46^(2500*22).
```

Both sides have 1,058,257 bits. This proves
`rho_formal<7071/10000<1/sqrt(2)`; the threshold check is
`2*7071^2=99998082<100000000`. The combined receipt retains both integer
comparisons, complete profiles, all compiler hashes and input hashes.
Both certificates run under the same verification target.

The bounded degree-22 discovery also timed out with an incumbent; no
optimality or minimum-degree result is claimed. Its `2 by 12` box instead
requires `535<=528`, which is false. The same eight box shapes fail in
both orientations, with maximum excess eight. Both original point files
and their discovery provenance are retained unchanged. The combined
receipt now proves exclusion as well as the historical formal inequalities;
neither point is a remaining polynomial realization problem.
