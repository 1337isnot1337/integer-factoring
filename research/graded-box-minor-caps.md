# Box-minor budgets with independent physical energies

## Registered extension

Base: `37dc71d`. Extend the bare box theorem to the actual graded setting
used by weighted physical types, then apply it to the graded adjugate.
The mechanism is the determinant's exact physical charge: each selected
X letter occurs `beta` times and each selected Y letter occurs `alpha`
times. The next test is literal Laurent tensors with independent heights
at zero and infinity, complete orders, direct box determinants and the
full adjugate cofactor identity. No numerical profile sweep is registered.
The complete factoring score remains `1/2`.

Use the [cell-cap normalization](cell-minor-caps.md): nonnegative sorted
profiles `r,s` with first order zero, normalized deficit `d>=0`, and
nonnegative paired energy lists `X,Y` with minimum zero. Alphabets have
the same size `q`; individual leg charges can have either sign. The claims
below concern scalar Laurent polynomial crossings over every field.

## Original crossing

For `1<=alpha<=beta<=q`, let `N=alpha*beta`, `h=alpha^2`, and let
`a,b` be either oriented box coefficient pair from the bare theorem.
Write `Low_k(X)` for the sum of the k smallest energies. Then

```
a dot r + b dot s <= N*d + min(
    beta*Low_alpha(X) + alpha*Low_beta(Y),
    alpha*Low_beta(X) + beta*Low_alpha(Y)).              (1)
```

At zero energies this is exactly the earlier box budget. At `alpha=1`
it recovers the shared-prefix bound; it is not an independent replacement
for the arbitrary cell-set bounds.

To prove it before normalization, restrict a determinant box to alpha
matched X letters and beta matched Y letters. The ordinary-cut determinant
order is bounded below by its source prefix. The opposite cut's integral
Smith layers realign to rank at most `min(alpha,beta)^2`, giving the other
weighted prefix by the bare argument. A common Laurent shift restores
the same inequality if either profile begins at a negative order, since
each coefficient vector has sum N.

The box determinant is nonzero at the strict endpoint. Its raw zero and
infinity orders therefore have sum at most zero before gauges. The gauges
add exactly `beta*sum(G_X)+alpha*sum(G_Y)`. Subtract the two profile minima
and the two physical energy minima; this gives `N*d` plus the normalized
selected energy sums. Choose the lowest energies, and then the stronger
of the two physical placements, to obtain (1). Both cut orientations hold.

## Word-controlled adjugate

For a restriction to the lowest m energies on each alphabet, put

```
n=m^2,
C_m=(n-1)*d + m*(sum(X[:m])+sum(Y[:m])),
L_i(r;m)=sum_{j=1}^n r_j-r_(n+1-i).
```

Let `H_i(s;m)` be the higher-word profile used in the
[bare composition](box-adjugate-word-caps.md), and let `Top_k` sum the k
largest energies within these m-letter cores. Then

```
a dot L(r;m) + b dot H(s;m) <= N*C_m - max(
    beta*Top_alpha(X[:m]) + alpha*Top_beta(Y[:m]),
    alpha*Top_beta(X[:m]) + beta*Top_alpha(Y[:m])).       (2)
```

Again include both box and source orientations, at every `2<=m<=q`.

The [graded cofactor identity](graded-adjugate-caps.md) gives the actual
adjugate infinity profile `L(r';m)-C_m` and paired energies `-X,-Y` on
the restricted core. Apply the raw form of (1) to that crossing. Selecting
the lowest negative energies selects the highest original energies;
moving the `-C_m` shift, whose coefficient sum is N, gives the right
side of (2). Interlacing makes `L(r';m)>=L(r;m)`. The word theorem and
its positive-minimum restoration give the other profile at least `H(s;m)`.
Every coefficient is nonnegative, proving (2).

For `alpha=1`, (2) recovers the inherited shared-word inequalities. At
zero energies it is exactly the new bare box-adjugate family. Integer
heights give the polynomial proof; rational heights follow by ramification
and real heights by continuity of the finitely many determinantal-prefix
functions. No monomial determinant assumption or uncharged pole is used.

This is a necessary compatibility theorem, not a global entropy floor,
a tensor construction or a factoring improvement. Previous box, graded
cofactor and word proofs retain credit. Substantial OpenAI assistant
support; no Lean or independent mathematical review is claimed.

## Literal verification

The [checker](../tools/check_graded_box_minor_caps.py) uses twelve actual
scalar polynomial sources from the preceding adjugate controls, including
restrictions, a nontrivial qutrit cell and its opposite orientation. For
each source it assigns zero heights and two independent signed integer
height patterns on all eight physical legs. These **36 Laurent cases**
span F2, F3, F5 and F7, including positive raw minima and nonuniform
energies with unrelated signs on the individual legs.

It reconstructs 696 complete source and derived-adjugate orders. For
every original box it expands the complete endpoint determinant before
gauges, then restores the two raw orders and their physical charges
separately. All **348 literal box determinants**, **252 original caps**
and **576 composed caps** pass. For the actual graded adjugate it acquires
every cofactor using exact monomial row/column charges, verifies
`S adj(S)=det(S) I` coefficientwise, and checks both the exact infinity
complement profile and the opposite word lower bound. It then applies
the original graded boxes to that actual adjugate as a separate check.

The [receipt](../evidence/graded-box-minor-controls.json) retains all signed
heights, normalization constants, original and adjugate polynomials/profiles,
literal determinant polynomials, physical budgets and source hashes. At
zero energies the composed records agree exactly with the bare compiler.
The box prefix is evaluated from cumulative sums, avoiding dense coefficient
arrays; no numerical optimization or truncated tail is used.

Run `make verify-graded-box-minors`. These finite controls validate the
implementation and charge accounting; the preceding ordinary proof gives
the uniform result. No global weighted pressure exclusion is inferred.
