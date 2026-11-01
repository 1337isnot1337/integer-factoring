# Mixed occupation and height in six-to-five split graphs

**Status, 2026-11-01:** exact constructor and independent finite controls;
the following numerical experiment is registered before its pressure screen.
There is no stronger factoring exponent or strict supporting-line certificate.
The complete score remains **1/2**. No new Lean formalization is attempted.

The source is the existing reflection-pair Gaussian at interval eight,
restricted to occupation words `(0,2,4,8,12,1)`. Deflate the canonical even
kernel from `w` to `z=w^2` and multiply its selected Pfaffians by `z^8`.
Every resulting entry is polynomial, the degree is 22, and all 1,296
endpoint entries equal the six-letter crossing. Restricting each leg to
the first five words reproduces all 625 entries of the degree-18 reference
in [the earlier recentered search](../tools/search_f4_recentered_crossings.py).
No implicit Gaussian formula replaces the profiles of a compressed tensor.

For `a,b` in `F4^5`, take the right columns and left rows

```
P_j = e_j + a_j e_5,       L_i = (delta_ij + b_i a_j, b_i).
```

In characteristic two, `LP=I`. Choose these independently for X/V and
Y/U and use `(Lx,Ly,Py^T,Px^T)` on the four physical legs. The existing
[physical split implementation](../tools/physical_split_crossings.py)
already proves and implements this endpoint-preserving operation. This
note introduces a finite chart to search, rather than a new split theorem.
The resulting five-letter tensor has both complete generic ranks 25.

When `a_0` is nonzero, the right image mixes the vacuum and the additional
one-particle word. It can omit the pure vacuum and fail deletion closure
and height homogeneity. Consequently the old homogeneous compound formulas
do not calculate its compressed profiles. The new implementation contracts
all entries and uses the complete native polynomial oracle. It charges
both physical types and all 32 nonreference heights at the two places.
The [fixed-gauge axis theorem](polynomial-crossing-axis-obstruction.md)
still applies. The bounded interior interval below is a search choice,
not a global exclusion of other angles or moving gauges.

## Independent controls

[The checker](../tools/check_mixed_graph_retractions.py) verifies:

- All ambient endpoints and the 625-entry reference identity.
- 1,295 independent determinant identities on the selected kernel minors,
  using permutation subset DP rather than the Pfaffian recurrence. For
  selected size `k>=2`, `det(z^8 K_I)=z^(8k-16) T_I^2`; the empty and
  singleton cases are checked separately.
- 1,875 entries from direct four-index contraction with separately packed
  F4 arithmetic, plus 1,875 coefficientwise Frobenius identities.
- 150 independently recomputed, untruncated fraction-free Smith orders
  at new integer heights, including nonbinary maps.
- 102 numerical derivative comparisons and nine comparisons between the
  native objective and a closed formula for a fixed unit profile pair.

The controls are [retained here](../evidence/mixed-graph-retractions-controls.json).
The determinant, coefficient and integer-profile comparisons are exact;
the derivative and pressure checks are floating comparisons. The native
512-coefficient search window must resolve every promised pivot. Errors
are recorded and rejected, never counted as zero tail coefficients.

## Registered campaign

Base main: `a32d9b4056314b546902f2495d2e9dbf6740706d`.
The [manifest](../evidence/mixed-graph-retractions-manifest.json) records
source hashes, the domain and budgets. Run only from a clean frozen commit;
the driver records that commit and its source hashes in every summary.

1. For each graph, allow `a` and `b` to be zero or have one nonzero
   coordinate. There are 16 choices each. Place the graph on X alone,
   Y alone, or both; remove duplicates. This gives exactly 766 choices.
   Frobenius `omega -> omega^2` has 106 fixed choices and 436 orbits.
   Evaluate every orbit representative, all 50 orders each. Because the
   ambient tensor is binary, Frobenius preserves every weighted profile
   at every gauge. This quotient is valid for later free-height fits too.
2. At the inherited natural heights, fit the balanced scalar curve and
   three interior supporting-angle starts. Cache only identical complete
   unit profile pairs. The natural physical energy is `(0,10,6,2,8)` on
   both pairs. Its [global shape floor](natural-five-letter-pressure-floor.md)
   forces balanced `rho>=5/7`; this is a negative control, not a promising
   balanced escape. Use log scale in `[-9,log(4)]` and log tangent in
   `[-2,2]`, starting at angles `0.35,pi/4,pi/2-0.35`.
3. Select the reference and the 16 distinct nonreference literal tensors
   with lowest fitted interior support. Fit all 32 heights, log scale
   and log tangent from three angle starts for each. Heights lie in
   `[-64,64]`. Each SLSQP trial has at most 100 iterations and five seconds,
   checked between complete queries; the stage has 600 seconds total.
4. Starting from that stage's best successful point, run a **7,200-second**
   coefficient/gauge walk. Change one of the 20 graph coefficients by
   a uniformly chosen nonzero F4 increment. Recompute the literal tensor,
   canonicalize by Frobenius and warm-start its free-height fit, with the
   same per-trial bounds and caps. Accept improvements, or uphill changes
   with probability `exp((old-new)/temperature)`. Temperatures cycle
   through `(.001,.005,.02,.1)` in 50-proposal blocks; restart from the
   global best every 200 proposals. Record an OS-random 128-bit seed.

The whole independent chart has `4^20` coefficient tuples. The sparse
screen is not a census of that chart, and the walk is not exhaustive.
The screen budget is 600 seconds. A stage is marked incomplete if its
total budget prevents its declared finite inventory from finishing.
Per-trial timeouts and nonconvergence retain their explicit statuses.

The driver streams every trial's coefficients, initial parameters, best
and last successful parameters, native profile details, optimizer status,
and every rejected query's coordinates and reason. It counts successful
queries without retaining each intermediate successful query. Best literal
tensors and atomic summary checkpoints are retained. A floating support
below one requires fresh independent exact profiles and a rational outward
pressure certificate before any recursion or score claim.

```
make verify-mixed-graph-retractions
python3 tools/search_mixed_graph_retractions.py --stage screen --seconds 600 --output /path/screen.json
python3 tools/search_mixed_graph_retractions.py --stage free --seconds 600 --input /path/screen.json --output /path/free.json
python3 tools/search_mixed_graph_retractions.py --stage walk --seconds 7200 --input /path/free.json --output /path/walk.json
```

## Completed initial stages at `1461d7d`

The [complete initial archive](../evidence/mixed-graph-initial-screen-complete.json.xz)
retains every screen case, every selected free-height trial and their
summaries. All 436 Frobenius representatives resolved 50 orders each,
with no rejection. They give 436 distinct literal tensors and 158 distinct
unit profile pairs. The best balanced rate is the unchanged reference
`0.7263479477044683`; the best inherited-height interior support is
`1.0183848616343762`. The screen took 12.86 seconds.

All 51 selected free-height trials were attempted, making 9,614 successful
complete queries and one rejected query in 14.27 seconds. Only three trials
reported optimizer convergence. Forty-seven reached the 100-iteration
limit, and one stopped on an unresolved native prefix; its full rejected
coordinates and reason are retained. The best successful numerical query
has support `1.0104250223304403`, above the critical value one. Thus this
archive is a completed attempt ledger, not 51 converged fits or an optimum.
The two-hour coefficient/gauge walk remains a separate registered run.

## Completed two-hour walk at `1461d7d`

The [summary](../evidence/mixed-graph-walk-complete.json.xz) and
[full proposal ledger](../evidence/mixed-graph-walk-complete.jsonl.xz)
retain every proposal, mutation, initial coordinates, best and last complete
query, optimizer outcome, acceptance decision and rejection. The declared
7,200-second budget completed with **19,059 proposals**, 9,998 accepted
moves, **3,054,386 successful complete-profile queries** and 86 rejected
queries. Of the optimizer outcomes, 4,356 reported convergence, 14,616
reached the iteration limit, 85 stopped on an unresolved native prefix,
one had no positive common-type root, and one reached the final wall-clock
cap. The rejected coordinates and reasons remain explicit.

The best successful support was `1.0097450553486638`, still above one.
This is a finite numerical campaign, not an optimum or an obstruction to
the full `4^20` chart. It produced no stronger complete factoring bound;
the score remains `1/2`. Further coefficient walks are not justified by
this result. A continuation must change a structural assumption rather
than extend the same search budget.

Credit: this uses the existing reflection Gaussian source, generic F4 split
maps, complete rectangular Smith backend, Gibbs pressure and asymmetric
chain rule. Cinna Davis / 1337isnot1337 directed the campaign; substantial
OpenAI assistant support was used. The current contribution is a literal
mixed-occupation chart, independent controls and a reproducible numerical
campaign. It is neither a new Gaussian reduction nor an improved algorithm.
