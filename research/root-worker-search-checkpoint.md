# Completed root-worker searches and resumed physical-frame runs

## Result and scope

This checkpoint preserves six interrupted physical-frame searches and
their completed continuations, together with the independent-charge,
dense-target, filtered-Laurent, and truncated-affine searches. None supplies
a new complete factoring exponent or a certified crossing improvement.

The [compact archive](../evidence/root-worker-search-archive.json.gz)
contains all 6,722 saved physical candidate scores and canonical record
hashes, full winning records, and the original source and continuation
metadata. The seven smaller target/search reports are retained in full.
The original physical optimizer traces are summarized, not reproduced
byte for byte. Their original-byte hashes and lengths remain recorded;
the local raw files are unchanged. In particular, this is not a claim to
reproduce every optimizer trajectory from its seed.

The [replay](../evidence/root-worker-search-replay.json) separates exact
coefficient, endpoint, finite-rank, and linear-system checks from numerical
pressure calculations. Native weighted-profile calculations use exact
`F4` polynomial coefficients but floating valuation comparisons. The
reported pressure values are diagnostics, not rational certificates.

## 1. The six interrupted and resumed families

| Family | Preserved prefix | Continuation | Saved total | Best recorded candidate objective |
| --- | ---: | ---: | ---: | --- |
| Four-letter physical bases | 1,458 | 135 | 1,593 | balanced `rho = 0.731080936633` |
| Four-letter asymmetric bases | 954 | 751 | 1,705 | support value `1.001432887611` |
| Sixteen-letter local bases | 20 | 58 | 78 | balanced `rho = 0.716653659728` |
| Polynomial physical frames | 927 | 1,482 | 2,409 | balanced `rho = 0.731080935973`; support `1.001438165356` |
| Rank-two square-zero frames | 226 | 411 | 637 | nonbaseline balanced `rho = 0.793662338724`; support `1.002018244245` |
| Affine coupled frames | 32 | 268 | 300 | nonbaseline balanced `rho = 0.737916287620`; support `1.001530331311` |

The square-zero and affine workers also retain their original baseline
`rho = 0.731080935973`. Their top-level best-score fields include that
baseline; the table distinguishes the best actually tested deformation.
All balanced values stay above `1/sqrt(2)`, and all asymmetric support
values stay above one.

Every prefix explicitly reports interruption. Its continuation records
the prefix byte hash, saved candidate count, elapsed checkpoint time,
remaining declared budget, original source hash, and the hash of the
instrumented initialization. The replay applies the recorded textual
initialization changes to the frozen source and verifies that hash.
The continuation uses a new declared random stream. Lost random-generator
and optimizer states were not restored, and unsaved interrupted work is
not counted as a completed candidate.

Five families reached their remaining time budgets, with small final
candidate overruns. The affine family instead exhausted its selected
300-candidate refinement queue: its recorded combined time is about
2,857 seconds of the original 3,600-second limit. It does not claim a
completed hour of computation or an exhaustive affine-space search.

## 2. Additional completed workers

| Worker | Recorded coverage | Best recorded result and stopping scope |
| --- | --- | --- |
| Independent ternary charges | All `3^6 = 729` masks; 364 distinct tensors; 80 refinements | `rho = 0.731080935973`; support `1.001104166005`; selected refinement queue exhausted |
| Dense complementary-rectangle targets | 5,610 of 5,610 systems; 162 feasible targets; 157 distinct seeds | 59 retained seeds; best retained natural diagnostic `rho = 1.160587126274` |
| Dense paired Bell targets | 18 exact systems | All 18 recorded inconsistent; these systems are preserved rather than independently rerun here |
| Affine Hasse targets | 218 of 624 planned systems | No feasible recorded target; time budget exhausted, target queue incomplete |
| Filtered sixteen-letter Laurent series | All 64 elements of the three-dimensional `F4` parameter space; 10 refinements | `rho = 0.716653659728`; refinement queue incomplete |
| Degree-16 truncated affine lattices | 16,352 walks; 2,467 distinct target vectors; 7,401 full profile queries | 5,219 profile classes; best `rho = 0.887618577263`; time budget exhausted |
| Degree-22 truncated affine lattices | 49,784 walks; 1,687 distinct target vectors; 5,061 full profile queries | 2,241 profile classes; best `rho = 0.737485307007`; time budget exhausted |

The Hasse report's `completed` flag means its timed worker returned
normally. Its `target_queue_exhausted` flag is false. The same distinction
between a finished timed run and an exhausted search space applies to
the truncated and refinement searches.

## 3. Fresh replay

Run the checker against the committed archive:

```sh
python tools/checkpoint_root_worker_searches.py
```

To rebuild the same compact archive from the preserved local original
reports before replaying it, use `--build-archive`. This reads the raw
reports and writes only the archive and replay; it does not alter any
search source or original report.

The checker rebuilds each physical family winner and checks all endpoint
entries, the saved tensor digest when available, every complete weighted
profile order, and its numerical objective. It also:

- Reconstructs all 729 independent masks, their three-character identities,
  distinct tensor hashes, and degree histogram, then replays both winners.
- Checks all 59 retained dense witnesses at both cuts and both nontrivial
  `F4` places, including the imposed zero rectangles and protected constant
  coefficient when required.
- Reconstructs the first and last saved Hasse systems from the exact
  separately archived 1,911-dimensional span and verifies their ranks and
  nonzero residuals. It does not extrapolate these two fresh checks to an
  exhaustive replay of all 218 systems.
- Reconstructs all 64 filtered Laurent tensors and all 32,768 natural
  profile orders. The 16,384 infinity orders remain unchanged.
- Rebuilds both truncated regular Smith frames, checks the best complete
  tensor against its affine constraints and endpoint, and recomputes its
  32 integer profile orders and numerical objective.

The input bindings identify the existing dense-frame catalogue and
polynomial-span archive by byte hash. This checkpoint does not duplicate
the separately owned fused-basis, Laurent-four-letter, filtered-affine,
implicit-code, restricted-place, or finite-fiber archives.

## 4. Publication manifest for the arithmetic results

The [constructor and trace manifest](../evidence/bound-search-publication-manifest.json)
lists 25 frozen arithmetic packages and their 77 proof, checker, module,
and evidence files in mathematical dependency order, followed by this
search-archive package and its source dependencies. Every listed file has a byte hash and
length, and each package has its checker command. The manifest covers
earlier published prerequisites as well as the later constructor and
trace packages; the local tracked-file field is only a snapshot for
integration, not a claim about the current remote branch.

The manifest's scope cautions retain the relative-linearized
`(s+2)^{O(J)}` operation count and the unresolved joint-contraction cost.
The constructor and local arithmetic improvements remain distinct from
the missing crossing improvement needed for a factoring score above
one half.
