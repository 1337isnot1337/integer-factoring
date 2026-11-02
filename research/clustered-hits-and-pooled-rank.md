# Path C: repeated smooth relations and a pooled rank condition

**Status.** Ordinary proofs of two elementary statements, and a fixed finite
diagnostic. No uniform smoothness, rank-window, or stronger factoring theorem
is proved. The complete saving exponent remains **1/2**.

This continues Grok's [power-window investigation](leading-power-window-barrier.md)
and the repository's [smooth selector](smooth-leading-power-selection.md) and
[orbit moment criterion](cofactor-orbit-moment-criterion.md). The owner's
question is whether the encouraging average smoothness data can support a
worst-case bound. Density, empty windows, repeated rows, and rank generation
are different obligations. The independent
[entry-moment replay](entry-moment-cycles.md) agrees on the two base cycles;
this note adds actual generator-seven strides and relation/rank certificates.

## Mechanism and stopping test

Keep the full-smooth-or-proper-gcd predicate, with cofactor cap one. Change
the sufficient success condition from one hit for every target prime to
full rational rank among all acquired relations. This can tolerate a late
target when other directions already supply its missing information.

The diagnostic was registered in [Discussion 51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19340754).
It uses only two prescribed complete cycles, and then the saved hit phases:
`(N,Y,g)=(450296711,23,2)` and `(1000036000099,19,5)`, deriving `g=7`
inside the second cycle. No input-size, cutoff, or constant sweep is used.
The discovery composition cap is 250000. A useful outcome is a distinct
rank condition or an obstruction to interpreting the density as new rows.
Failure to get a finite query gain is recorded as failure, not retuned.

## Uniform fixed-row multiplicity

Use the discriminant, height `H`, split-prime base `B`, prime forms `g_p`,
and `K_p=min{k:p^k>H}` from the smooth-selector note. For a smooth reduced
leading norm `A<=H`, the exact oriented decoder gives an integer word `v`
with

\[
 [g_p]^k=\prod_{q\in B}[g_q]^{v_q},\qquad
 \prod_{q\in B}q^{|v_q|}=A\le H.
\]

The acquired relation is `r=k e_p-v`. In particular,
`|v_p|<=K_p-1`. If a fixed **integer row** `r` occurs at several integer
exponents of this same target, every such exponent satisfies

\[
 r_p-(K_p-1)\le k\le r_p+(K_p-1).
\]

Thus it occurs at most **`2K_p-1`** times, in an interval of diameter at
most `2K_p-2`. This is unconditional within the exact smooth decoder's
contract. It does not assert that different rows, scalar multiples, or
different hit clusters are rationally independent. Proper gcds already
terminate the factoring step and need no decoded row.

## A weaker sufficient window condition

Let `R_W` contain every decoded smooth row from all targets `p in B`
and exponents `K_p,...,K_p+W-1`. Define `PRW(Y,W)` to mean that a proper
gcd occurs, or `rank_Q(R_W)=m`, where `m=|B|`.

**WL implies PRW.** Select the first smooth row for each target if no gcd
has occurred. The existing weighted strict diagonal-dominance proof applies:
`sum_q |v_q| log q <= log H < k log p`. The resulting square matrix is
nonsingular, and is a submatrix of the pool. No randomness is used.

**Conditional complete transfer.** Under the same factor-base generation,
preprocessing, and ambiguous-class extraction contracts as the selector,
a uniform PRW bound gives deterministic factoring in

\[
 m^{O(1)}W\,\operatorname{poly}(b,\log Y,\log W).
\]

Maintain at most `m` rationally independent integer rows. Each of at most
`mW` norm tests and decodings has polynomial overhead, and the row entries
have `O(log(b+W))` bits. Exact rank maintenance and the existing extractor
therefore have polynomial overhead. A full-rank square sublattice suffices
for that extractor; this is not a claim that deleting dependent rows
preserves the torsion information of an arbitrary partial presentation.
Standard primality, perfect-power, and recursive splitting steps retain
their existing costs and assumptions.

For fixed `c`, a cutoff `Y=2^{O((log b)^c)}` and a uniform bound
`W=exp(O(b/(log b)^alpha))` would give saving exponent `alpha` through this
transfer. Neither such a rank bound nor the needed uniform orbit-density
bound is supplied here. Global Dickman density and a whole-class-group
character estimate alone do not settle this subgroup rank condition.

## Exact finite results

The [cycle receipt](../evidence/clustered-power-relations/receipt.json)
uses actual bit lengths (the first `N` has 29 bits, despite its original
30-bit scan bucket). `g5` and `g7` generate the same 198850-element cycle;
`g7=g5^66509`. The true `g7` flags are
`f7(j)=f5(66509*j mod 198850)`. Starting the base walk at `g7^K7` and
then advancing by `g5` would be a different experiment.

| N, Y, target | Period | Smooth classes | Hit clusters | Longest zero run | First offset from K |
| --- | ---: | ---: | ---: | ---: | ---: |
| 450296711, 23, 2 | 20456 | 2821 | 453 | 202 | 7 |
| 1000036000099, 19, 5 | 198850 | 343 | 61 | 17939 | 664 |
| 1000036000099, 19, 7 | 198850 | 343 | 73 | 4506 | 3823 |

For entries `e(j)=f(j)(1-f(j-1))`, density equals raw density divided
by mean hit-run length: respectively `6.2274`, `5.6230`, and `4.6986`.
At entry windows with expected counts at least `1,2,4,8` (ceiling lengths),
exact empty-phase fractions are:

| Target | Mean 1 | Mean 2 | Mean 4 | Mean 8 |
| --- | ---: | ---: | ---: | ---: |
| First cycle, 2 | .30935 | .12290 | .00880 | 0 |
| Second cycle, 5 | .61345 | .46590 | .17080 | 0 |
| Second cycle, 7 | .28954 | 0 | 0 | 0 |

The generator-five entry tail is still substantially heavier than a
Poisson tail at mean one or two. Removing adjacency does not justify an
independent-hit model. The registered moment tests use
`q=bit_length(b+1)` and `q^2`, constants `C=1,2`, measured density, and
the certified input-only period bound `2N`. Each holds at its own sufficient
window. However, **every entry sufficient window on the larger cycle
already covers the whole period**, so those passes are trivial coverage.
The `q^2`, `C=1` raw window there is 46400; the `C=2` window is 90480.
The latter generator-seven window requires constant two. These are tests
at different windows, not a claim that constant one works at every length.
Measured density denominators are diagnostics, not uniform input bounds.

The elementary period bound used here follows by counting reduced forms:
for each `1<=a<=H`, at most `a+1` middle coefficients have the correct
parity in `[-a,a]`. Thus `h<=H(H+3)/2<=2N/3+sqrt(3N)<=2N` for `N>=2`,
and every cyclic subgroup has period at most `h`.

The [row audit](../evidence/clustered-power-relations/row-audit.json)
checks every smooth row at nonwrapping exponents `K<=k<h` in both larger
directions, including exact ideal identities and the integer extractor:

| Target | Rows | Distinct integer rows | Equal adjacent pairs | Maximum multiplicity | First full-rank exponent |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 | 334 | 61 | 273 / 273 | 14 <= 17 | 19968 |
| 7 | 336 | 73 | 263 / 263 | 11 <= 13 | 71962 |

All 536 adjacent pairs repeat their row. This observation is finite; a
reduction boundary can change the decoded word, so it is not a universal
cluster identity. Each direction's selected three rows produces the proper
factor `1000003` of `1000036000099` through the existing checked extractor.

The [pooled audit](../evidence/clustered-power-relations/pool-audit.json)
derives the generator-19 index from a verified relation and verifies its
group identity, rather than enumerating another cycle. In target order
`5,7,19`, full rank first occurs at window length **2884**, with the last
independent row from target 19. Target 7 needs **3824** positions for its
first hit. Thus PRW is strictly weaker than WL on an actual finite input.
The pool takes **8652** round-robin norm probes, versus **7373** for one
first hit per target. There is **no query improvement** in this example.
These are schedule counts on frozen flags, not a timed faster selector.

## Reproduction and limits

The discovery budget is 219306 base-walk compositions plus a conservative
1000-composition phase allowance and the charged row/pool audits, totaling
225841. Replaying verification is a separate computational expense, not
part of that discovery cap. From the repository root:

```sh
python3 tools/pathc_cluster_tail/generate_orbits.py
python3 tools/pathc_cluster_tail/audit_rows.py
python3 tools/pathc_cluster_tail/audit_pool.py
make verify-clustered-power-relations
```

The three producers re-evaluate every arithmetic phase on the two cycles,
then the selected relation events. Runtime fields may change. The faster
default checker verifies frozen flag hashes, minimal group orders, exact
statistics and moments, arithmetic phase samples, all 670 saved hit-row
identities, rank histories, and extracted factors. It does **not** recompute
every rough phase; use the full producers for that stronger replay. No
Lean or independent proof review is claimed.

The [clean full replay](../evidence/clustered-power-relations/full-verification.json)
passes `make verify` and the focused target at exact source head
`b9737e9764d6719984daa5f4d58fd65d0f8ef16d` in 692.28 seconds, with
clean tracked state before and after. Its complete compressed log is
archived; the replay predates only this receipt/documentation follow-up.

The next relevant theorem would bound rank-generation windows uniformly,
possibly by exploiting arithmetic correlations directly. A polynomial
bound on fixed-row repetition removes one obstacle, but supplies neither
a lower density nor independence among different rows. Larger scans of
these same parameters would not close either gap.
