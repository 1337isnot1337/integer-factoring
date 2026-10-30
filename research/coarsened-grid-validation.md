# Finite coarsening and clique-table check

The integration proposal requires a row or column network to be grouped and
materialized as a clique with a common edge width. The existing grid check
stops at the network contraction. `tools/coarsened_network.py` provides a
small exhaustive reference for the next step; it is not a factoring engine
or a proof of the uniform machine cost.

For each group, the reference sums over its original node settings. Internal
edges are contracted inside that sum. For each pair of groups it packs the
original crossing channels in edge order, selecting the left or right factor
according to the original endpoint, and pads to the chosen common width.
Unused high bits are forced to zero. The resulting clique tables are then
contracted over all labels. This mirrors the semantics of `Net.coarsen` in
the upstream pinned `FactoringAux/Net.lean`; it does not import or formally
verify that implementation.

`tests/test_coarsened_network.py` checks:

* A signed four-node network with internal edges, reversed group orientation,
  multiple crossing factors on one group pair, a zero-bit factor, an empty
  group and padding. The checker compares the grouped clique against a
  separately enumerated original network, over integers and modulo 4, 6, 9,
  and 15; an insufficient common width is rejected.
  Its value is `-720`; reversing one original edge's factors changes it to
  `1440`, so this fixture can detect an orientation mistake.
* Each signed width-one row-network piece from the existing 2-by-2 grid
  certificate, with both split and merged groups. It compares the row
  contraction, the original network value, and the grouped clique value,
  over integers and modulo 4 and 9. Summing the two coarsened pieces yields
  `-1084 + 700 = -384`, eight times the original clique contraction `-48`,
  as required by the multiplier.

The tested groupings exercise value preservation independently of whether
they satisfy the grouping theorem's width condition. The toy parameters do not
meet the theorem's sufficient size hypothesis. These tests do not establish
the general grouping theorem, a multi-level evaluator, the modular positivity
test, or any asymptotic or machine-level factoring claim.

Run `python3 -m unittest discover -s tests -p 'test_coarsened_network.py' -v`.
