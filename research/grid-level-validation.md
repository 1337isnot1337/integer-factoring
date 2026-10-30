# Exact certificate for one complete grid step

`tools/check_grid_level.py` checks the interaction of the new crossing
pieces with actual vertex tables on a `2 by 2` grid. This is a finite
integration certificate. It does not construct or certify the final
factoring machine.

For each case it computes four quantities independently:

1. The original four-node clique contraction, summing its six undirected
   edge labels.
2. Each modified piece contraction, freeing the two diagonal edge copies
   and summing eight independent labels.
3. The row network contraction through explicit left/right factor
   channels and the two remaining vertical equality channels.
4. The column network contraction through its factor channels and the
   two remaining horizontal equality channels.

The direct, row, and column values agree for every piece. Their sum is
`8^ceil(lambda/3)` times the original contraction. There is one rectangle
on this grid; the test enumerates its piece choices directly, without
the source program's redundant non-rectangle choice indices.

The exact results are:

| Label width | Tables | Original value | Modified piece values | Sum |
| --- | --- | ---: | --- | ---: |
| 1 | All ones | 64 | 256, 256 | 512 |
| 1 | Signed seed 100 | -48 | -1084, 700 | -384 |
| 1 | Signed seed 101 | 87 | 408, 288 | 696 |
| 1 | Signed seed 102 | 22 | 116, 60 | 176 |
| 1 | Signed seed 103 | -194 | -212, -1340 | -1552 |
| 2 | All ones | 4096 | 16384, 16384 | 32768 |
| 2 | Signed seed 200 | -165 | -7812, 6492 | -1320 |
| 2 | Signed seed 201 | 169 | -1240, 2592 | 1352 |

All eight cases also pass actual modular row and column contraction
checks modulo `4,6,8,9,15`. These include moduli dividing the multiplier.
The checks use integer arithmetic and modular reduction; no inversion
or division modulo a composite is used. Passing a residue identity when
the multiplier is zero modulo a modulus is not used to infer positivity.

For `lambda=1,2`, the minimum label width selects the trivial `2lambda`
fallback. The tests explicitly check this branch and separately check
the compressed channels:

- At `lambda=3`, all 4096 block quadruples are checked for a compressed
  row piece and its corresponding column piece.
- At `lambda=6,t=2,i=1`, the actual compressed width is `11<12`.
  Five deterministic examples, including positive and negative values,
  check the mixed rank-8/rank-56 decoder in both orientations. Each
  decoder sums at most 2048 padded channels, without building a large
  grid or materializing all pairs of six-bit strings.

The rank-8 and rank-56 block functions are imported from
`tools/check_crossing.py`. The grid implementation independently assembles
their channels into row and column node tables. This tests their use in
the grid transformation, beyond the original base tensor identities.

Run:

```sh
python -m unittest discover -s tests -p 'test_grid_level.py' -v
python tools/check_grid_level.py --json
```

All five grid test methods pass. The full-grid verifier rejects more than
100,000 expanded states by default, before constructing the state space.
The largest tested grid has 4096 original assignments and 65,536 expanded
assignments. The zero-width and invalid-input cases are also checked.

`evidence/grid_level.json` records the full signed input tables, exact
piece values, channel-width decisions, modular results, and the five
mixed-channel certificates.
