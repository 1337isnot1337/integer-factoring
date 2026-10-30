# A path obstruction for universal grouping

## Result and scope

The sufficient grouping theorem in
[`grouping-improvement.md`](grouping-improvement.md) uses
`g = O(D^2 K^2 log K)` vertices. A unit-weight path gives a lower bound of
order `K^2` for fixed discrepancy `D`, even when the maximum weighted degree
is only two. The bound is exact for this path family. It does not show that
the logarithmic factor or the quadratic dependence on `D` in the sufficient
threshold is necessary.

For integers `K >= 3` and `D >= 2`, define

\[
g_*(K,D)=\left\lceil\frac{D(K^2+2)}{2(D+1)}\right\rceil.
\]

Consider the grouping conditions from the grouping theorem, with `g` vertices,
maximum weighted degree `P`, ordered edge total `U`, class sizes `s_a`, and
cross-class edge weights `Y_ab`:

\[
K s_a\le 2g,\qquad D K^2Y_{ab}\le DU+gP\quad(a\ne b).
\]

**Theorem.** For the unit-weight path on `g >= 3` vertices, a valid grouping
exists exactly when `g >= g_*(K,D)`. Consequently, any universal sufficient
vertex threshold for these conditions must be at least `g_*(K,D)`. For fixed
`D`,

\[
g_*(K,D)=\frac{D}{2(D+1)}K^2+O(1),
\]

so a universal threshold that is `o(K^2)` cannot hold at fixed discrepancy.

The necessary inequality is also formalized in
[`formal/GroupingObstruction.lean`](../formal/GroupingObstruction.lean), using
the pinned upstream's actual `grpSize`, `grpY`, and `totU` definitions. The
formal file separately proves symmetry, a zero diagonal, maximum weighted
degree at most two, and ordered total `2(g-1)` for the path. The contiguous
construction proving sufficiency below is handwritten and tested in Python;
it is not part of the Lean theorem.

The source reference is
`IsaacRudich/SubexponentialDeterministicFactoring` revision
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`; the grouping definitions and
factoring engine are due to Isaac Rudich and Louis-Martin Rousseau. The proof
and checker in this note were developed with substantial OpenAI Codex assistance.

## Proof

For a path on `g >= 3` vertices, every edge has weight one, the maximum
weighted degree is `P=2`, and the ordered edge total is `U=2(g-1)`. The
right-hand side of the pair inequality is therefore

\[
DU+gP=2D(g-1)+2g=2(D+1)g-2D.
\]

Suppose

\[
g<\frac{D(K^2+2)}{2(D+1)}.
\]

Then `2(D+1)g-2D < DK^2`. If any path edge joined two distinct color
classes `a,b`, its unit weight would give `Y_ab >= 1`, contradicting the pair
inequality. Every path edge must therefore be monochromatic. Connectivity
forces all vertices into one class, which violates `K s_a <= 2g` because
`K >= 3` and `s_a=g`.

Now take `g >= g_*(K,D)` and partition the path into consecutive blocks of
size at most `ceil(g/K)`, using at most `K` blocks. Such a partition exists
because `K ceil(g/K) >= g`. Since `D/(D+1) >= 2/3` and
`(K^2+2)/3 >= K` for `K >= 3`, we have `g_*(K,D) >= K` and hence `g >= K`.
Therefore

\[
K\left\lceil\frac gK\right\rceil\le K\left(\frac gK+1\right)
=g+K\le2g,
\]

which proves the size condition. Each pair of distinct blocks has at most
one path edge between them, so `Y_ab <= 1`. By the definition of `g_*`,

\[
2(D+1)g-2D\ge DK^2,
\]

and the pair condition holds. This proves existence at every integer
`g >= g_*`. Together, the two directions give the exact threshold.

## Exact finite checks

[`tools/check_grouping_obstruction.py`](../tools/check_grouping_obstruction.py)
computes the threshold using integer ceiling arithmetic. It exhausts every
coloring at `g_*-1` for four small `(K,D)` pairs and checks a contiguous-block
witness at `g_*` with the independent grouping verifier. The accompanying
unit tests include the equality case `(K,D,g)=(4,2,6)`, input validation,
and additional thresholds.

The four nonexistence checks enumerate respectively 27, 81, 1,024, and 65,536
colorings. The equality case `(K,D,g)=(4,2,6)` has a valid grouping, so the
strict inequality in the obstruction is essential. Full receipts are in
[`evidence/grouping_obstruction.json`](../evidence/grouping_obstruction.json).

Reproduce the finite checks with:

```sh
python3 -m unittest discover -s tests -p 'test_grouping_obstruction.py' -v
python3 tools/check_grouping_obstruction.py
make verify-formal
```

These finite checks illustrate the proof; the exact all-`K,D` threshold
follows from the argument above, not from exhaustive search. The result is a
lower bound on a universal grouping lemma only. It is not a lower bound on
the integer-factoring algorithm or on its achieved time exponent.

## Consequence for the parameter-growth route

For `D >= 2`, the necessary condition gives `g >= K^2/3`. Thus a schedule
that insists on a grouping guarantee for **every** weighted graph and uses
the current grid relation `K_next = g^2` has

\[
K_{d+1}\ge K_d^4/9.
\]

Writing `L_d = log K_d - (log 9)/3` gives `L_(d+1) >= 4 L_d`.
For the current `K_0=48`, `L_0>0`. Such a schedule therefore has at most
`log_4 log P + O(1)` recursive levels before reaching group count `P`.
Sharpening the universal grouping theorem alone cannot replace that base
four by a smaller leading growth base. Removing `log K` or improving the
dependence on `D` remains compatible with this obstruction.

The actual row and column networks form a restricted graph family. This
path example does not prove that those networks require the universal
threshold. Exploiting their special structure, changing the grid relation,
or changing the crossing decomposition lies outside this lower bound.

The next concrete question is whether the grid networks admit a stronger
grouping guarantee using properties that exclude these sparse paths. A
universal theorem with a subquadratic vertex threshold is closed by this
obstruction and should be reopened only with different hypotheses.
