# A fixed path/group hybrid misses its finite seed budget

The complete [nilpotent-seam coefficients](nilpotent-seam-crossings.md)
and their cut reflections do not yield a decomposition of the crossing
with `S<L^sqrt(2)` at the original six-/eight-site sizes `L=5,13`, even
after adjoining I, every affine-line group-sum core and the all-ones core.
This is an exact finite statement about using the **individual** dictionary
elements as pieces. It does not exclude first combining coefficients into
new low-rank pieces, changing the auxiliary, other sizes, or an asymptotic
path construction. The complete factoring exponent remains **1/2**.

## Question, hypotheses and stopping test

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19366407)
asks whether adding known balanced cores can complete the path family's
cheap asymmetric channels. The prior edge/window exclusions do not cover
this full dictionary. All bidegrees, repeated coefficient directions and
both cut orientations are retained. The original ordered path labels are
numbered `0,...,L-1`; no relabelling or parameter search is performed.
The [budget clarification](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19366460)
adds the rank-one all-ones core and charges the number of line cores.

The original16 tests concern a strict polynomial that pays a separate
constant tensor I. Sol worker 1's [scope check](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19366719)
correctly notes that I is not algebraically forced by the full dictionary:
the affine-plane identity below omits it. A [registered follow-up](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19367169)
checks113 additional spans without reserving any I cost. The strengthened
conclusion uses those new tests, rather than silently widening the first16.

In both versions the pieces are scalar multiples over F4 of individual
dictionary elements, and their `sqrt(rc)` costs are added. The new tests
cover direct decompositions without a constant-I restriction as well.
Combining dictionary terms before measuring the rank of the resulting
new piece is outside scope. Boundary subalphabet restrictions are not
tested. This is not a lower bound on all crossing decompositions.

Stop if exact separating functionals exclude every dictionary selection
within that budget. A successful span would only justify constructing and
charging an actual identity; it would not itself prove the strict budget.
The stop occurred, so no larger site count follows from this diagnostic.

## Group-sum pieces and their collective identity

For the label fields F5/F13, let

```math
K_\sigma(P,Q,R,S)=[P+\sigma Q=S+\sigma R],\qquad
K_\infty(P,Q,R,S)=[Q=R].
```

These are zero/one tensors with coefficients in F4. Their labels being
elements of F5/F13 does not change the coefficient field. In either cut
the row and column labels map surjectively to the same L values, so each
piece has rank pair `(L,L)` and cost `sqrt(rc)=L`. This includes sigma0
and infinity. The all-ones tensor `J` has pair `(1,1)` and cost1.
The group-sum cores were already proposed in the parent seam note; they
are not new low-rank tensors claimed here.

Put `E(P,Q,R,S)=[P=S][Q=R]`. The finite affine-plane incidence identity is

```math
\sum_{\sigma\in\mathbb F_L\cup\{\infty\}}K_\sigma=J+E.
```

For `(P-S,Q-R) != (0,0)` exactly one direction vanishes; at zero all
L+1 directions vanish, and their number is even. This proves the
identity over F4. Its explicit cost `(L+1)L+1` is too large for the
target. An uncharged collective direction sum is not a new evaluator.

## Budget reduction and exact certificates

The constant `I(P,Q,R,S)=[P=R][Q=S]` has pair `(L^2,1)` and cost L.
The following strict bounds are proved using integers alone:

```math
\sqrt2<10/7,\quad5^{10}<10^7;\qquad
\sqrt2<99/70,\quad13^{99}<38^{70}.
```

Thus `5^sqrt2<10` and `13^sqrt2<38`. Every nonconstant term in a
winning five-letter seed costs less than5; no rank5 line term fits.
Every nonconstant term in a winning thirteen-letter seed costs less
than25; at most one cost13 line term fits. If that term is used,
each remaining path piece costs less than12.

For that reserved-I criterion it is sufficient to exclude these generous spans,
where I and J are allowed freely and all qualifying path pieces can
be used without charging their total number:

| L | Line term allowed | Path condition | Span dimension with I,J | Result |
| --- | --- | --- | --- | --- |
| 5 | none | `rc<25` | 8 | E outside |
| 13 | none | `rc<625` | 14 | E outside |
| 13 | one fixed direction, all14 cases | `rc<144` | 7 in each case | E outside |

For the direct criterion with **no reserved I cost**, a winning
five-letter decomposition has at most one cost5 line piece. With no
line each path piece costs less than10; with one line less than5.
A winning thirteen-letter decomposition has at most two cost13 line
pieces; the remaining path costs are less than38,25,12 for zero,one,two
line pieces respectively. Repeated uses of one line can be combined
without increasing its rank cost, so distinct direction subsets suffice.
Sol worker 1 supplied this bookkeeping clarification; the follow-up
uses exactly those budgets, without a new tensor/rank construction.

Again allowing I and J freely gives the following113 conservative tests:

| L | Line directions | Path condition | Branches | Span dimension with I,J | Result |
| --- | --- | --- | --- | --- | --- |
| 5 | none | `rc<100` | 1 | 28 | E outside |
| 5 | one | `rc<25` | 6 | 9 each | E outside |
| 13 | none | `rc<1444` | 1 | 20 | E outside |
| 13 | one | `rc<625` | 14 | 15 each | E outside |
| 13 | two | `rc<144` | 91 | 8 each | E outside |

The [receipt](../evidence/path-group-hybrid.json) supplies a sparse F4
linear functional for each of all129 cases. Each functional annihilates
every tensor in its declared span and pairs nontrivially with E.
The original witnesses have5--9 nonzero entries; the additional ones
have at most20. Arbitrary F4 coefficients and repeated
uses are covered by linearity. Extension fields cannot change solvability
of a linear system with F4 coefficients and target, so the same fixed
dictionary exclusion holds after coefficient-field extension.

## Replay, limits and credit

```sh
make verify-path-group-hybrid
```

The checker rebuilds both complete source tensors, checks their exact
crossing endpoint, checks all84 and135 coefficient rank pairs by BOTH
the existing packed backend and independent dense F4 elimination, and
checks both reflected orientations. It constructs all6/14 line cores,
checks the incidence identity, verifies every sparse functional on every
eligible coefficient, and independently checks all129 span dimensions.
The original tensor bytes and source hashes are retained.

An earlier private minimum-bottleneck experiment omitted J and used a
looser all-core threshold. It motivated the invoice clarification and is
not the certificate above. No exact rank-sum optimizer or uniform
asymptotic obstruction is claimed. The useful conclusion is that simply
appending the known cores to these unchanged individual path pieces
does not supply the desired finite seed. It is stopped rather than
extended to larger labels.

The original path braid, F4/F16 descent, nilpotent auxiliary, fusion
projector bounds, reflected-piece construction and known group-sum cores
retain the parent notes' source and contributor credit. This fixed hybrid
test and its exact certificates were developed by Sol worker 2 with
OpenAI Codex assistance. Sol worker 1's reserved-I scope check and the
unreserved budget formulas are explicitly credited above. No new Lean or independent mathematical
acceptance is claimed.
