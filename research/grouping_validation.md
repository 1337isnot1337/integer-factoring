# Finite validation of the grouping checker

`tools/grouping.py` implements exact integer certificate verification and a
bounded lexicographic search. The input is a loopless graph with nonnegative
integer edge weights. Repeated undirected edges are merged by addition; zero
weights are discarded. It checks the two inequalities

\[
K|V_a|\le2g,
\qquad DK^2Y_{ab}\le DU+gP,
\]

with `U` the ordered total, so every edge contributes twice to `U`.

The test command is:

```sh
python -m unittest discover -s tests -p 'test_grouping.py' -v
```

All 10 tests passed. The exhaustive test considers every labeled graph on zero
through four vertices with each undirected edge weight in `{0,1,2}`. This is
761 graphs. For each graph it checks every coloring for `(K,D)=(2,2),(3,2),
(3,4)`, giving 131,513 comparisons against a separately written specification
of the ordered-edge inequalities. It also verifies that the search returns the
lexicographically first valid coloring, or reports exhaustive nonexistence.

The other tests cover matching and star graphs, duplicate edges, zero weights,
an empty graph, malformed endpoints and weights, invalid parameters and color
assignments, a too-small degree bound, and explicit search limits. The cap is
checked before enumeration. A saturating exponent calculation avoids creating
an astronomical integer just to decide whether `K^g` exceeds the cap.

The certificate examples in `evidence/grouping_examples.json` are regenerated
by `python tools/grouping.py --demo`. They include a two-vertex, one-edge graph
at `K=3,D=2` for which all nine colorings fail. Its failure is consistent with
the lemma because the theorem's sufficient size condition is not met.

All these finite examples are below the theorem's sufficient threshold. They
test the implementation's exact inequalities and behavior; they do not prove
the concentration theorem or the proposed factoring runtime. The mathematical
proof of the general statement is separate.
