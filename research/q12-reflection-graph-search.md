# Arbitrary twelve-letter graphs in the full reflection source

**Registered constructive search, 2026-11-01.** This is an actual polynomial
tensor family, rather than a choice of necessary spectral data. No member
is claimed to improve the factoring score **1/2**. Historical twelve-letter
relaxation targets motivated the alphabet size. Their bare zero-energy
versions are now [excluded by box budgets](q12-integer-word-target.md);
the nonuniform physical-type search here has a different scope. No new
Lean formalization is attempted.

Use the original [interval-eight reflection tensor](../tools/reflection_dual_rail_gaussian.py)
on all sixteen occupation words. It is polynomial, even in `w`, and has
degree 96 in `w`. Deflate to `z=w^2` to obtain degree 48. All 65,536
endpoint entries equal the sixteen-letter crossing. Reorder its physical
words by height, breaking ties by integer word:

```
(0,8,4,12,2,10,1,6,9,14,5,13,3,11,7,15).
```

The first twelve give a coordinate reference of degree 38. The full source
restricted to the six words `(0,2,4,8,12,1)` is exactly `z^8` times the
[mixed-graph source](mixed-graph-retractions.md), entry by entry. Both
sources use the canonical literal reflection constructor. Common monomial
shifts are retained at both places and cancel in the support pressure.

## Independent paired graphs

On each matched physical pair choose arbitrary matrices `A` of size
4-by-12 and `B` of size 12-by-4 over F4. Use

```
P = [I_12; A],       L = [I_12 + B A, B].
```

Characteristic two gives `LP=I_12`. Choose independent graphs on X/V and
Y/U, with 192 field coordinates in all, and contract the complete tensor
using the existing [physical split routine](../tools/physical_split_crossings.py).
Both compressed cuts have generic rank 144, guaranteed by the endpoint.
The parameter chart has `4^192` tuples. A timed walk is not its census.

Arbitrary graph images can mix occupation, height and parity. They can
remove the pure vacuum, so this is not a protected-boundary retraction
family. Every one of the 20,736 compressed entries is computed literally;
no Gaussian profile formula or deletion-closure hypothesis is used. The
original source, the split identity and generic complete-profile engine
are credited existing tools. The contribution is their unrestricted
twelve-letter chart, independent controls and registered experiment.

## Controls and replay levels

[The checker](../tools/check_q12_reflection_graphs.py) compares all 1,296
shared source entries, all compressed endpoints of two charts, 192 direct
four-index contractions with separate packed F4 arithmetic, and 41,472
Frobenius coefficient identities. It resolves 576 complete native orders
at new integer gauges and checks all 180 scaled pressure derivatives.
Two three-letter principal restrictions contribute 36 independent exact
untruncated orders, including a nonbinary chart.

The optional `--full-independent` mode additionally recomputes every one
of the coordinate reference's 288 orders with untruncated Python
fraction-free arithmetic at new integer heights. This does not supply a
second full nonbinary 144-by-144 replay. The receipt flags whether the
larger independent replay actually ran; a quick check must not imply it.
All native queries require every promised pivot inside the 512-coefficient
window. Unresolved profiles are rejected and retain their coordinates.

The [completed independent receipt](../evidence/q12-reflection-graphs-full-controls.json)
records the larger replay: both full coordinate cuts and both smaller
principal restrictions passed, **324 untruncated orders** in all. The
576 native orders, 180 derivative comparisons and other stated controls
also passed. This is finite implementation verification, not a favorable
pressure or an all-gauge theorem.

## Registered numerical campaign

The [manifest](../evidence/q12-reflection-graphs-manifest.json) records the
domain, source hashes, OS-random seeds and budgets before the search.
Run from a clean frozen source. One complete coordinate unit query
supplies its balanced scalar fit and three fixed-height angle/scale fits.
Each is followed by a free fit of all **88 nonreference heights**, log
scale and log tangent. There are two physical costs and no uncharged type.

Each free fit uses at most 100 SLSQP iterations and ten seconds, checked
between complete queries. Heights lie in `[-128,128]`, log scale in
`[-9,log(4)]`, log tangent in `[-2,2]`. The three initial free fits have
a 90-second total budget. Stopped and rejected trials remain explicit.

Then run two independent 1,800-second tracks from the same best initial
point. One changes a single uniformly selected graph coordinate by a
uniformly selected nonzero F4 increment. The other changes 2, 4 or 8
distinct coordinates, with equal probability for those three counts.
Each proposed tensor is constructed, all endpoints checked, and a bounded
free-height fit warm-started from the current state. Frobenius conjugation
is a valid quotient at every height because the original source is binary.

Accept improvements, or uphill proposals with probability
`exp((old-new)/temperature)`. Temperatures `(.001,.01,.05,.2)` cycle in
25-proposal blocks. Restart from that track's best point every 100 proposals.
The second track does not start from the first track's discoveries.
Each seed is recorded before execution. These heuristics do not cover
all graph components, all gauges, all angles or all scales.

Every trial retains its literal tensor digest, field coordinates, initial,
best and last parameters, complete best/last native profile details,
optimizer status, mutations, acceptance decision and rejected queries.
Intermediate successful queries are counted, not individually retained.
The best literal tensor is saved atomically. A numerical support below one
requires fresh exact profile and rational outward-pressure certification
before any recursion, exponent or leaderboard claim.

```
make verify-q12-reflection-graphs
make verify-q12-reflection-graphs-independent
python3 tools/search_q12_reflection_graphs.py --manifest evidence/q12-reflection-graphs-manifest.json --output /path/q12-walk.json
```

## Completed raw-support tracks at `10aed41`

The [complete registered ledger](../evidence/q12-reflection-graphs-search-complete.json.xz)
retains the three starts and every proposal of both 1,800-second tracks.
The single-coordinate track attempted 173 proposals and 7,797 successful
complete queries; the multi-coordinate track attempted 170 proposals and
5,178 queries. Including the initial fits, there were 13,254 successful
queries and no rejected query. Elapsed time, including the initial fits,
was 3,617.23 seconds.

Only two initial fits reported convergence. Across the entire ledger,
340 attempts hit the trial/global time cap and four hit the iteration
cap. Every partial status and best/last point is retained. Both tracks
kept the initial best support `1.0053988159437441`, above one. This is a
completed timed attempt ledger, not a global gauge or chart exclusion.

Cinna Davis / 1337isnot1337 directed this work with substantial OpenAI
assistant support. Existing reflection, split, Smith and asymmetric/Gibbs
pressure implementations are reused with their attribution preserved.
