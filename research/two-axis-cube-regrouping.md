# Audit of the rejected two-axis cube regrouping

## Result

The proposed regrouping in `ddb303e` / `00e00c9` does **not** establish
an exponent above 1/2. Its one-cut transfer is invalid. The
[original proposal](two-axis-cube-regrouping-proposal-rejected.md) and
[original receipt](../evidence/two-axis-cube-regrouping-proposal-rejected.json)
are retained as rejected historical work, including the original claims.
The corrected [audit receipt](../evidence/two-axis-cube-regrouping.json)
checks the exact counterexamples below. The complete factoring score
remains **1/2**.

The valid Gaussian cube, its three-cut profiles and its rational pressure
certificates are unchanged. The error is in changing the global grid
transfer, not in those local computations. This is an ordinary source
audit with exact controls and OpenAI assistant support; it is not a new
Lean theorem or independent mathematical acceptance.

## 1. The established transfer has one global plane orientation

The [multiaxis transfer](multiaxis-cube-transfer.md) fixes the actual cube
legs at the parent vertices and, for a complete piece assignment, chooses
one axis a minimizing

~~~
min_a sum_cubes r_(piece(cube))^a.
~~~

All planes are then formed along that same axis. A tensor factor is local
to one of its two coordinate planes because the selected cut is exactly
that plane partition. The source's average-width estimate bounds this
minimum of global sums.

The rejected proposal instead uses

~~~
sum_cubes min_a r_(piece(cube))^a.
~~~

These quantities need not coincide. For widths (0,3,4) and (3,0,4), the
second quantity is zero and the first is three. The inequality between
them points in the direction that makes independent minima smaller; it
does not provide an implementable common plane cut. The numerical width
and child-volume calculations in the original checker assume this missing
geometric step. Their arithmetic remains true as arithmetic, without
proving the step.

## 2. A hidden-axis cut is not a cut on the fixed parent labels

Fix cube coordinates 0 and 1 as the two parent axes. Each of the four
parent superlegs bundles the two values of coordinate 2. The cuts along
axes 0 and 1 partition whole superlegs. The cut along axis 2 splits **all
four** superlegs into their internal halves.

Thus the third cut cannot be inserted as a factor between two existing
parent planes. Changing which axis is hidden changes which actual legs
belong to a parent label. Such a change requires a new wiring argument;
endpoint invariance by itself does not relocate vertex-table arguments.
Even the two available parent orientations must respect the common global
plane partition used by the cited transfer.

The checker verifies this incidence directly on the eight physical legs.
All three pairings still give the correct endpoint when used consistently.
Those valid pairing identities do not authorize independently changing the
pairing after each summand is selected.

## 3. Endpoint symmetry does not permit separate summand permutations

If E=sum_i T_i and a common permutation P preserves E, then

~~~
E = sum_i T_i composed with P.
~~~

For permutations P_i depending on i, the corresponding identity does not
follow. An exact characteristic-two example already occurs on eight bits.
Let E_8 be the product of the four antipodal equality indicators. Set

~~~
A(x)=x_001,     T_0=E_8+A,     T_1=A.
~~~

Then T_0+T_1=E_8. Swap cube coordinates 0 and 2 in T_1 only. This
permutation preserves antipodes, but the new sum is

~~~
E_8 + x_001 + x_100.
~~~

For the assignment x_001=1 and every other bit zero, E_8=0 and the
new sum is one. Of the 256 bit assignments, exactly 128 disagree. The
checker retains this complete finite counterexample.

This does not claim that every piece-specific reindexing is impossible.
A special family could have additional identities. Such an identity must
be proved for the actual pieces and wiring; the rejected proof supplies
neither. Its claim that a coordinate permutation only relabels the same
piece omits the change of that piece as a function of fixed leg arguments.

## 4. The claimed uniform fixed-cut conclusion contradicts rank

There is a separate asymptotic check. At a fixed row cut,

~~~
E_4(x,y,u,v) = [x=v][y=u]
~~~

is a permutation matrix on row labels (x,y) and column labels (u,v).
With lambda'-bit letters its rank over every field is exactly 2^(2lambda').

Suppose polynomially many pieces, at most (lambda'+1)^C, summed to E_4
and **every piece at that same fixed cut** had channel width at most
2rho*lambda'+C'log_2(lambda'+2). Rank subadditivity would give

~~~
2^(2lambda')
 <= (lambda'+1)^C * 2^(2rho*lambda'+C'log_2(lambda'+2)).
~~~

For any fixed rho<1 and constants C,C', this fails for sufficiently large
lambda'. Polynomial factors cannot compensate for the linear exponent gap.
The claimed one-cut construction purports to give precisely such a uniform
fixed-row decomposition by choosing each piece's minimum width. Hence
preservation of the endpoint, the asserted row width and polynomial piece
count cannot all hold.

For a finite arithmetic illustration with zero additive error,
lambda'=32, C=4 and rho=5/8 already give

~~~
33^4 * 2^40 < 2^64.
~~~

The checker also builds the endpoint permutation matrix at a four-element
alphabet and checks its full rank 16. These finite cases illustrate the
all-parameter proof; they do not replace its asymptotic quantifiers.

This rank argument concerns the proposed **uniform fixed-cut conclusion**.
It does not invalidate the original multiple-cut crossing theorem, in
which different pieces can have large rank at different cuts, and the
single orientation is chosen from total assignment widths.

## 5. Surviving results and reopening conditions

The following original controls remain valid and are retained:

- Three consistent pairings of the antipodal endpoint into four superlegs.
- The elementary minimum/average inequalities and formal volume formulas.
- The Gaussian one-particle profiles (2,4,8,10) and their source certificate.
- The numerical facts 2(5/8)^2<1 and the sharper rational comparison.

They do not establish the missing transfer. In particular, substituting
the numeric pressure into log_2(1/rho) is a formal rate calculation until
a valid growth-base-two evaluator has been supplied. The existing
three-axis exponent remains log_3(1/rho), below 1/2 for this family.

A replacement must give a consistent exact tensor identity on the actual
parent labels and a global child-network construction with charged table
queries, edge widths and regrouping costs. Merely repeating independent
minima, changing local labels or increasing the arithmetic test panel does
not resolve the failure. Other tensor families, structured contractions
and different global networks remain open.

Run make verify-two-axis-cube-regrouping. Its success means the audit and
counterexamples pass, and its explicit factoring_exponent_improved field
is false. The rejected historical JSON is not a current validation receipt.
