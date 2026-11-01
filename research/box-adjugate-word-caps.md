# Apply box budgets to the word-controlled adjugate

The [graded extension](graded-box-minor-caps.md) retains independent
physical energies through the original box and adjugate budgets. Its
necessary inequalities do not supply a uniform pressure floor.

## Question and scope

Base: `4c3d135`. The box budgets and higher-word adjugate bounds have so
far been imposed separately: the newer adjugate checks apply cell/shared
budgets to the derived crossing. Applying **every box budget to that same
crossing** gives another necessary coupling. This changes the compatibility
theorem, rather than extending a bounded profile search.

This note concerns scalar polynomial strict crossings over every field,
bare profiles and zero physical energies. Auxiliary matrix entries and
general graded physical energies are outside the statement. The intended
next checks are complete literal adjugate controls and an exact audit of
the thirteen-letter point that passes the previously combined families.
The complete factoring score remains `1/2`; no universal thirteen-letter
pressure exclusion follows merely from rejecting that point.

## Composed inequalities

Let `F(z)` be a polynomial crossing with entry degree at most `M`, endpoint
`F(1;x,y,u,v)=[x=v][y=u]`, and `q` letters. Use one ordinary row cut at
zero, with complete nonnegative profile `r`, and the opposite cut of
`z^M F(1/z)` at zero, with profile `s`. The two cuts may be exchanged.

For `2<=m<=q`, put `n=m^2` and define

```
L_i(r;m) = sum_{j=1}^n r_j - r_(n+1-i),   1<=i<=n,
H_i(s;m) = the m-letter ordered-word profile of (s_1,...,s_n),
C_m = (n-1) M.
```

The word profile retains the empty word, all ordered words of length at
most `m-1` in the costs `s_2,...,s_n`, and arbitrarily many copies of
`m*s_2`, taking the first `n` costs with multiplicity. For every box
budget `(a,b,N)` from the
[box theorem](bare-seven-to-twelve-letters.md), with coefficient vectors
`a,b` on `n` channels and `N=alpha*beta`, the necessary inequality is

```
sum_i a_i L_i(r;m) + sum_i b_i H_i(s;m) <= N C_m.       (1)
```

Include both box orientations and both source orientations `(r,s)` and
`(s,r)`. These constraints include the former shared-word conditions
at `alpha=1`; boxes with larger first side provide additional conditions.

## Proof at the full alphabet

Regard the row cut as a matrix with columns `(v,u)`, so its endpoint is
the identity. Form the scalar polynomial crossing

```
B(z) = adj( row(z^M F(1/z)) ),        C=(q^2-1) M.
```

Its endpoint is again the identity. Reversing its parameter and restoring
degree `C` gives exactly `adj(row(F(z)))`. Smith diagonalization therefore
gives its normalized infinity row orders exactly as
`L_i(r;q)=sum r-r_(n+1-i)`. The higher-word theorem gives the opposite
zero orders `eta_i>=H_i(s;q)`, with the positive-minimum restoration
justified below when `s_1>0`.

The box budgets apply to the actual two profiles of `B` with homogenizing
degree `C`, even when their first orders are positive. To see the latter
point directly, let these first orders be `u,v`. Divide the common zero
monomial `z^v` from `B`; its actual degree is `C-u-v`. The normalized
profiles subtract `u,v`. Each box coefficient vector has sum `N`, so
restoring both shifts gives precisely the original right side `NC`.
No monomial determinant assumption is used.

All box coefficients on `eta` are nonnegative. Substituting its word
lower bound proves (1) for `m=q`. Reversing the original parameter and
exchanging its middle physical legs gives the other source orientation.

## Restriction and positive minima

Restrict to any `m` letters on each matched pair. Determinantal
interlacing gives its actual source orders `r'_i>=r_i`, `s'_i>=s_i`
for `i<=m^2`. Every `L_i` is the sum of all but one source order,
so is coordinatewise nondecreasing.

If the restricted opposite profile has minimum `b>=0`, remove its common
monomial before applying the word theorem. The actual adjugate orders
are bounded below by

```
(n-1)b + H_i(s'-b;m) >= H_i(s';m) >= H_i(s;m).
```

For each word the first inequality follows because its length is at most
`m-1<=n-1`; for the formal cap it follows from `m<=n-1`. Sorting the
candidate costs preserves the inequality. The word-cost functions are
coordinatewise nondecreasing, giving the second comparison. Thus the
restricted crossing proves (1) using just the first `m^2` source orders.

The proof uses the earlier box theorem, exact adjugate reversal identity
and higher-word theorem, retaining their attribution. Substantial OpenAI
assistant support; no Lean or independent mathematical review is claimed.

## Exact controls and target rejection

The [checker](../tools/check_box_adjugate_word_caps.py) constructs twelve
literal polynomial crossings including restricted controls, over F2, F3,
F5 and F7. It expands all 322 complete cofactor polynomials, retaining
every coefficient through the full degree bound, and verifies
`T adj(T)=det(T) I` coefficientwise. It independently reconstructs 232
complete source/adjugate orders, including nontrivial qutrit word jets,
the opposite source orientation, physical restrictions and positive minima
at both places. All 192 composed inequalities and 84 literal adjugate-box
checks pass. The [receipt](../evidence/box-adjugate-word-controls.json)
retains the source polynomials, full adjugates, all profiles and hashes.
The universal claim follows from the composition proof, not those samples.

The [degree-2000 thirteen-letter point](joint-bare-profile-feasibility.md)
passes its original seven families but violates **52 oriented conditions**
from (1). For `m=13`, `alpha=7`, `beta=13`, use `N=91` and
`a=(49,42,0,...)`, `b=(1^91,0,...)`. The exact required inequality is

```
49 L_1(r;13) + 42 L_2(r;13) + sum_{i=1}^91 H_i(s;13)
    <= 91*168*2000 = 30576000.
```

Here `L_1=334136`, `L_2=335027` and the word prefix is `145158`.
The left side is `30588956`, an excess of `12956`. The original point and
its exact formal-rate calculation are preserved, but polynomial realization
is excluded. These are newly composed necessary conditions; passing the
original smaller family was a correct statement about that relaxation.
This rejects a point, not every thirteen-letter profile or a factoring
algorithm. The complete score remains `1/2`.

Run `make verify-box-adjugate-words`. A weaker linear compiler also uses
`H_i(s)>=s_2` for `i>=2`; its coefficients are nonnegative, and 1816 exact
comparisons verify that it is weaker on the recorded point. No global
entropy optimization or new alphabet sweep is asserted here.
