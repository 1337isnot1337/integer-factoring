# Contribution review and disposition ledger

This is the maintenance checkpoint of 2026-10-30 against `main`
`e91442d3162cae7f775e008efb2ea998c4478344`. GitHub is the
[live PR queue](https://github.com/1337isnot1337/integer-factoring/pulls).
Recheck heads before acting; this inventory is not a correctness review of
research drafts. See [REVIEWING.md](REVIEWING.md) for the evidence policy.

## Completed external verification request

[Issue #28](https://github.com/1337isnot1337/integer-factoring/issues/28)
requested packaging for an external Lean rerun. The reviewer subsequently
reported that packaging was unnecessary, independently rebuilt `e91442d3`,
and explicitly invited closure. The
[verifier run](https://github.com/bigospeedrun/verifier/actions/runs/41108289491)
at verifier commit `e6f6793aad7f517fe660d49caaa54295e6350607` succeeded;
its log reports that nanoda passed. This resolves the packaging request and
supports the README's precise external-verification status. It does not
establish an exponent above `1/2`.

## Archived predecessor routes

Both PRs are by **Cinna Davis / 1337isnot1337**, with OpenAI assistant support
disclosed in their contribution records. Their heads contain unique material;
they are archived as superseded routes, not represented as incorporated.

| PR and exact head | Retained source | Disposition and scope |
| --- | --- | --- |
| [#1](https://github.com/1337isnot1337/integer-factoring/pull/1), `18d788a0d8a37e9deb0a35389b4f1ddf0d4bee02` | [archive/pr-1-exact-crossing-width](https://github.com/1337isnot1337/integer-factoring/tree/archive/pr-1-exact-crossing-width) | Archived: integer widths for an earlier conditional exponent about `0.07551`, below the complete `1/2` theorem. The unique implementation and boundary tests remain available at the tag. Their correctness was not newly reviewed in this maintenance pass. |
| [#14](https://github.com/1337isnot1337/integer-factoring/pull/14), `a064085db89d1f54d1c75dc42ae22408e5f60d22` | [archive/pr-14-sixth-rate-bridge](https://github.com/1337isnot1337/integer-factoring/tree/archive/pr-14-sixth-rate-bridge) | Archived: conditional sixth-power rate arithmetic for the earlier `1/30` route. Complete k6 and `1/2` machine theorems now exist separately. The unique Lean file and recorded evidence remain at the tag; its Lean build was not rerun here. |

Archive tags preserve the exact heads before remote branch deletion. Local
branches, local worktrees, and their uncommitted changes are retained. A future
contributor can cite an archived result or open a new focused PR from its source.

## Work retained for review

The following inspected heads contain unique work absent from `main`; they stay
open. These entries record inventory and scope only. All authors at this
checkpoint are Cinna Davis / 1337isnot1337; retain their detailed source credits
and AI disclosures when reusing the work.

| PR | Inspected head | Remaining review scope |
| --- | --- | --- |
| [#4](https://github.com/1337isnot1337/integer-factoring/pull/4) | `16b9dba1839aa671a68802309d08d9817ee8f799` | Independent replay and review of finite coarsened-network/clique adapters, including the orientation negative control. Unique tests can remain useful after their original target is superseded. |
| [#22](https://github.com/1337isnot1337/integer-factoring/pull/22) | `d94e4a99fec453768d69ae083038ba65cf31ae9e` | Fixed-characteristic transfer, Pfaffian evaluator, and focused exact checks. |
| [#23](https://github.com/1337isnot1337/integer-factoring/pull/23) | `147652574e7ce1ff68b2c378712846cbf438da34` | Half-rate space bounds, direct counts, and construction-specific obstructions. |
| [#24](https://github.com/1337isnot1337/integer-factoring/pull/24) | `05d9c1b8aed4f8494ce230f0869f9048ef879474` | Effective crossing families, finite rank optima, and all-field barriers. |
| [#26](https://github.com/1337isnot1337/integer-factoring/pull/26) | `dd28d375ea65863356c67b051f5d4f217562475d` | Factoring transfers, divisor bounds, and new crossing certificates. |
| [#27](https://github.com/1337isnot1337/integer-factoring/pull/27) | `b49b9e19480c65d139ef15beecacb92064c1a960` | Factoring constructions and conservation-space barrier. |
| [#29](https://github.com/1337isnot1337/integer-factoring/pull/29) | `61a4f4537488900a77d6cfb18e90e88fd8165bd4` | Crossing constructions, finite costs, and operator-preserving rank methods. |
| [#30](https://github.com/1337isnot1337/integer-factoring/pull/30) | `28105165edb8e9a0427468dc833f3c09bc0a7fb7` | Factoring recovery and algebraic crossing interfaces. |
| [#31](https://github.com/1337isnot1337/integer-factoring/pull/31) | `32e758239d1b51210508adb6179d6ba56e8aed5d` | Structured frontier and preprocessing bounds; stacked on #29, whose branch must be preserved as its base. |

## 31 October frontier integration: PRs 39 and 43

This is an integration record, not independent mathematical review. Both
source contributions are by Cinna Davis / 1337isnot1337 with substantial
OpenAI Codex assistance; their detailed predecessor credits are retained.
The complete factoring exponent remains `1/2`, and no Lean source changed.

| Source | Exact inspected head | Integration scope |
| --- | --- | --- |
| [#39](https://github.com/1337isnot1337/integer-factoring/pull/39) | `a3c5bf15fcc80ea8ab9af7b4853a1a70ae6a88b0` | Retained the constructive typed/filter families, scoped ordinary obstructions, checkers, and finite search receipts. Resolved additive Makefile/index conflicts against main `a5e05feb756b4ad35ae1877d54a5f2e4eb387c19` by keeping both sets. Integrated in `dd386d8`; the broad local verification log completed. Large optional numerical searches and independent proof review were not repeated. |
| [#43](https://github.com/1337isnot1337/integer-factoring/pull/43) | `307d752cebb528dd3ee58c68abaf3e9747727645` | Retained cyclic/register and trace constructions, convex fixed-type Smith pressure, matrix-vacuum filters, Steiner and p-group cubics, and the actual-rank Steiner coefficient obstructions. Resolved additive README/Makefile conflicts after the prerequisite. Clarified that the fully multivariate claim applies at `q>=7`, while the three-letter middle-coefficient proof applies to the univariate span and a sublinear number of parameter classes. |

The integrity check found that main's PR #38 merge `a5e05fe` had omitted
141 new tool files present in source parent `c2e951f`; none existed in that
merge's common base. This caused 138 broken local documentation links.
The files were restored byte for byte from `c2e951f` in `dd386d8`, with
the source history retained. The analogous audits of the PR #36 and #37
merges found no omitted new source files.

On the combined tree, `make verify` passed with exit code zero. All ten
focused targets also passed: cyclic algebra, shift registers, copied-wire
barrier, quadratic registers, quantum Pascal, matrix vacua, Steiner cubic,
group shears, Steiner coefficients, and Smith capacity. These checks
replay their documented finite scopes; the optional full large catalogues,
numerical optimization searches, and Lean gates were not rerun in this
integration. Integration and passing checks do not establish a stronger
factoring machine or independent acceptance of the author proofs.

Concurrent main reached `6646d5a` during this integration, including the
GitHub merges of #39 (`e19d077`) and #42 (`65823ce`). It was reconciled
without dropping either set of research targets, and the combined broad
gate passed again. A later eight-file arithmetic addition at `601735e`
merged without conflicts; its smooth-norm relation and gap checkers passed,
as did the final documentation gate. Those later files did not modify the
frontier sources covered by the ten focused replays.

## 2026-10-31 integrated Bell and arithmetic checkpoints

These are author integration records for work by Cinna Davis with substantial
OpenAI Codex assistance. They preserve the ordinary mathematical and finite
evidence scope; they are not independent mathematical acceptance.

| PR and inspected head | Integration | Evidence and remaining scope |
| --- | --- | --- |
| [#41](https://github.com/1337isnot1337/integer-factoring/pull/41), `8c95d717cc30d1826f64ea2bfa93fc87a924ee50`; mathematical checkpoint `dbbd7b14d52f4ded9a85e694ce2f10573d9f3564` | Merged at `3a1cd6758393c70ee2e844ec343ecb077c16cf6f` | Bell transforms, folded/code constructions, scoped rank and pressure floors, rigidity interface, and completed finite searches. All 40 research/tool/evidence files from the checked mathematical checkpoint were compared with integrated main and retained byte for byte. `make verify`, focused Bell checks, and independent profile readbacks passed at the recorded source. The lane is released; no stronger factoring bound is claimed. |
| [#47](https://github.com/1337isnot1337/integer-factoring/pull/47), `80d3af0602476648d13601f28ff721822aec1187`; mathematical checkpoint `2aef203c88ae8de0fa19d103bc1d70f3723c402c` | Merged at `a24247394438492dc8d1689dbfb15cdf295774e8` | Arithmetic inverse sections, complete-orbit aggregation, global modular product circuits, and valuation-separated local jets. `make verify` and the class-group integration gate passed before merge; no new Lean source or local Lean build. Global rational-factor extraction remains unresolved. The later exact rational-part identity, specialization caveat, and 192 pole controls are retained at `4f752e4c2301c899d9af51415cc2e7a844e880e6`. |

The combined local default gate passed on the rational-part source before
its final rebase; that rebase also incorporated concurrent code changes.
The final documentation gate passed on `4f752e4` across 289 Markdown files.
After those integrations, `make verify` passed again at the exact main head
`38030d3f7f0d56c609a729ee69c294988edfbaaf`. These checks are local replay,
not independent mathematical review. The intervening missing-tool restoration
was coordinated with its source owner in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51).
The active global elimination question is
[Discussion #53](https://github.com/1337isnot1337/integer-factoring/discussions/53).
The strongest complete factoring score remains `1/2`.

## Maintaining this record

### 31 October actual-rank and supervisor continuation: PR 42

PR [#42](https://github.com/1337isnot1337/integer-factoring/pull/42) merged
at `1aa090a560998ccf4bdca86de792def642700eb5`, with exact source head
`236fd85832c8d7c8cc9f705e9fea4500cab33cb5`. Its earlier mathematical
checkpoint `6440212` had already been integrated through `65823ce`.
The incremental work retains the paired-downset normal form, the source-
pinned updated receipts, and the owner-requested Astra supervisor prompt
with current Discussions/main-first instructions. Cinna Davis /
1337isnot1337 and substantial Codex AI assistance retain their source
credits. This is an integration record, not independent proof review.

The broad `make verify` gate passed on integrated candidate `39deedf`:
50 unit tests, five crossing tests, default exact certificates, and links
in 291 Markdown files. The binary-site-symmetry, downset and regular-code
focused gates passed on the unchanged research sources. The new paired
receipt covers all 361 ordered three-bit downset pairs, 160,000 literal
polynomial entries and 4,670 actual rank queries. Later reconciliation
retained the incoming new standalone research files and all incoming
README/index link targets; the final local link gate checked 296 files.
The later optional arithmetic, flag and matching searches were not rerun
as part of this integration. No Lean source changed or Lean gate ran.

The integration retained the other frontier agent's byte-for-byte tool
restoration and full research-note/index restoration. The resulting diff
against the reconciled main source contained only the intended nine-file
PR delta. A Markdown coefficient-extraction ambiguity was reconciled to
the restored source's `operatorname` notation. A merge or passing gate
does not establish a new factoring bound; the complete exponent is `1/2`.

For each later decision, record the exact reviewed head, changed-file coverage,
checks actually run or reused, accepted scope, exclusions, source attribution,
retaining commit or archive tag, and actual GitHub disposition. A snapshot or
successful finite replay alone must not be recorded as mathematical acceptance.
When a PR changes, retain the historical checkpoint and review the delta.

## 2026-10-31: PR #50 integration by research director

- Inspected PR head `e5ba011e5264c6e78086f9b57da7a5be60048dc6` against
  main `38030d3` in an isolated checkout. Contributor: orda-aska / Claude Code.
- Preserved main's already-restored tools and subsequent research. Resolved
  additive README/index conflicts by retaining current content and restoring
  missing predecessor paragraphs/entries. Retained the 30 missing Make targets.
- The incremental code changes are lazy SciPy imports and a source-hash-checked
  embedded-reference fallback for saved walk records. No research claims changed.
- Passed `make verify`, the five exact binary/flag/target checks with SciPy
  available and with all SciPy imports deliberately blocked, and
  `make verify-coefficient-height-walks` (both embedded-reference records).
- This is source inspection and integration testing, not independent acceptance
  of the mathematical claims. No Lean rebuild or PR approval was performed.

## 2026-10-31: PR #49 integration by research director

Exact head 78b3cca917583eff2324671d6fd241c28ebd7019; contributor orda-aska / Claude Code. Inspected incremental source, preserved current main, and replayed make verify-cell-minor-caps byte-identically against evidence. Full make verify passed on the integration tree (315 docs, 50+5 tests). Finite floating findings remain distinct from global positivity. Integration is not independent mathematical acceptance; no new factoring bound, Lean rebuild, or PR approval.

## 2026-10-31: PR #54 integration by research director

Exact head fb5d2fe333a2c21b41c7f092ee25ef87d4801921; contributor orda-aska / Claude Code. Inspected scope and checker, preserved current main, replayed make verify-auxiliary-binary-floor byte-identically:230 cases,3320 rectangular profiles,128 all-minor controls,68 larger constant-contraction pressure cases. Full make verify passed on integration tree (322docs,50+5tests). The result concerns the stated local Smith criterion, not actual coefficient ranks or all auxiliary mechanisms. Integration/source inspection is not independent mathematical acceptance. No Lean rebuild or PR approval.

## 2026-10-31: PR #55 integration by research director

Exact head ba5974a05a998572d14763c3c40997aabf65f217 after its #40 prerequisite. Contributor orda-aska / Claude Code. Inspected checker and scope, retained optional generator without rerunning its floating search. Focused exact replay passed: six entropy certificates, q13 feasible profile, 140 literal crossings,1200 boxes,280 adjugate controls and mutation rejections. Full make verify passed on reconciled local tree (327docs,50+5tests); restored focused/default wiring while preserving current targets. The q13 profile is not a realizing tensor, and these are bare-criterion obstructions only. Integration is not independent mathematical acceptance; no Lean rebuild or PR approval.

## 2026-11-01: PR #44 glossary integration

Exact source head `b89916cde2b30c9d1c5bfc9500482856999f6b4e`, by Hao Yan / FanPUYen with disclosed OpenAI Codex assistance. Retained the seven-term glossary and quickstart link. Preserved main's newer contribution policy, discussion links, and precise distinction between the existing Lean endpoint and the ordinary general transfer. All glossary fragment targets resolved. `make verify-docs` and the full `make verify` passed on the integration tree (339 Markdown files, 50 unit tests and five crossing tests). No theorem change, new Lean work, PR approval, or independent mathematical acceptance.

## 2026-11-01: PR #4 coarsening checker integration

Exact head `16b9dba1839aa671a68802309d08d9817ee8f799`, by Cinna Davis / 1337isnot1337 with disclosed AI assistance. Retained the exact finite network-to-padded-clique reference, four unit tests, orientation negative control, and scoped note; added current README/index links without reverting later work. All four focused tests passed. Full `make verify` passed on the integration tree (340 Markdown files, 54 unit tests and five crossing tests). These are finite integration checks, not a uniform machine proof or a new exponent. No Lean rebuild, PR approval, or independent mathematical acceptance.

## 2026-11-01: PR #22 fixed-characteristic checkpoint integration

Exact head `d94e4a99fec453768d69ae083038ba65cf31ae9e`. Retained the original source attribution to Claude / Anthropic Appendix B and Isaac Rudich / Louis-Martin Rousseau, and the contribution's OpenAI assistant disclosure. Inspected the incremental source and preserved all current main work. Restored the six notes, seven implementation/checker files, six evidence records, default verification wiring, and extended `verify-modular` target. The weighted transfer is conditional on its stated evaluator; the Pfaffian algorithm evaluates profiles, and the finite rank/Schur/grid controls do not prove a stronger exponent.

`make verify-modular` passed, and all six JSON results matched the saved evidence exactly. The extended Pfaffian replay includes 45,122 source piece entries, 80 independent perfect-matching comparisons, 992 product/division controls, and 768 folded-block controls. Full `make verify` passed on the integration tree (347 Markdown files, 54 unit tests and five crossing tests). This is source inspection and integration testing, not independent mathematical acceptance. No new Lean work or PR approval review was performed.

## 2026-11-01: PR #23 time-space checkpoint integration

Exact source head `147652574e7ce1ff68b2c378712846cbf438da34`, with original author and source credits retained. Inspected the delta and reconciled it with current main, retaining 200 added research, checker, implementation and evidence files. The 60 notes are indexed separately; the [exact source README](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/README.md) remains available for its original summary. The time-space claims are ordinary proofs, while the rank catalogues are finite or explicitly scoped obstructions. Neither is promoted by this merge to independent mathematical acceptance or a stronger complete time exponent.

Preserved main's existing verification snapshot and saved the source's historical snapshot as `evidence/pr23-beyond-half-verification.json`, updating the two progress-note links. Its recorded hashes describe its historical source, not all of current main. Four false Markdown links in three formula tables were repaired with whitespace only. Preserved current README, Makefile dependencies, index and subsequent modular-crossing work; retained the regular-fusion cross-reference.

Full `make verify` passed on the reconciled tree: 424 Markdown files, 59 unit tests, five crossing tests, the 61 newly restored source checkpoint commands and the existing default gates. The dedicated `verify-pr23` target exposes those same 61 commands; no second redundant replay was claimed. No new Lean rebuild or PR approval review was performed.

### PR #56 integration record

Exact source `9c35bd30544a43291a1b63293c2c2656e97f098e` by Hao Yan / FanPUYen, merged at `db73b92725469b2e5ab4f49ac01f9c43f04d2f7d`. All eight focused offline overlap-search tests and the default gate passed (362 Markdown files, 59 unit tests and five crossing tests). A live read-only test of source blob `28d396324513a2de64a4264a704c002de3701763` scanned two discussion posts and 36 comments/replies and found a known nested reply. GitHub subsequently confirmed the PR merged. This repairs research coverage; it is not a novelty guarantee or mathematical result.

## 2026-11-01: PR #24 effective-family and finite-rank integration

Exact source head `05d9c1b8aed4f8494ce230f0869f9048ef879474`, with all source attribution and AI disclosures retained. Restored 85 missing files and indexed 19 notes. Preserved main's current half-rate evidence/leaderboard wording rather than restoring the source's outdated pending-acceptance statement. Retained the quantitative modular endpoint cross-reference and all existing default Makefile commands.

All 23 focused checkpoint commands passed under `make verify-pr24`. This includes compilation of the retained C verifiers, independent RUP verification of 8,532 additions against 9,802 initial clauses with four negative proof controls, and the independent classification of 47,677 rank-one outer products in the declared three-label problem. These certify the stated finite problems, not a uniform new rate. Full `make verify` passed on the integration tree: 448 Markdown files, 59 unit tests, five crossing tests and all default gates. The source's conditional rational-family specialization and exact modular half-endpoint remain ordinary proofs. No new Lean work, PR approval review, independent mathematical acceptance, or exponent above `1/2` is claimed.

## 2026-11-01: PR #26 effective-transfer checkpoint integration

Exact source head `dd28d375ea65863356c67b051f5d4f217562475d`, with original author/source attribution and AI disclosures retained. Restored 260 missing files and indexed 66 notes. Kept current main's precise half-rate and general-transfer evidence wording. Restored `verify-further` as a default dependency and the source's four optional extended/native targets without removing any current checks. The universal-search observation is explicitly not a new numerical factoring bound; conditional interfaces and construction-specific barriers retain their original scope.

Full `make verify` passed on the reconciled tree: 514 Markdown files, 59 unit tests, five crossing tests, the newly restored further checks and the existing default gates. Both 300-second tagged/bivariate checks actually completed with `status: passed`, rather than timing out. The bivariate replay checked all 2,048 tensors with both block and independent dense ranks and recorded the finite exact optimum. All five generated two-patch/shear verification records reported passed with their finite scope retained. The polynomial-shear diagnostic does not turn its retained partial search into an exhaustive result.

Not newly run: optional `verify-energy-full`, `verify-search-native`, `verify-weyl-orbit`, `verify-witt-full`, or Lean targets. Those source implementations/evidence are retained without a fresh full-run claim. Integration and finite checks are not independent mathematical acceptance, a PR approval review, or an exponent above `1/2`.

## 2026-11-01: Second spin coefficient checkpoint and session handoff

Exact published mathematical source `c31de0f66d3680311482e5cd8785bac2849f8537` adds the ordinary proof of `3C_2=sum T_ij tensor T_ij+k(k+1)I`, independent literal-transposition entry controls, and an interpretation of the existing eight-site row model. The constructor and original rank certificate remain credited to the bounded-spin lane. The checker verifies 198 single-copy class-sum entries and 104 second-coefficient controls, then reads the published polynomial model to label its blocks; that step does not rebuild the dense rank certificate.

`make verify-bounded-spin-second`, whitespace checks, and the full `make verify` passed before integration, with this additive patch on base `e8787861611aea3a43c1238fa0276b81eb964f41` (521 Markdown files, 59 unit tests and five crossing tests). Subsequent integration replay is limited to the focused checker and documentation gate. The latter initially found a newly introduced reference to the absent `smooth-set-fourier-norm.md`; removed the broken reference without claiming to have reviewed the pre-existing class-number bound. The final focused and documentation replays passed. No new Lean build, independent mathematical acceptance, favorable all-length cut-rank theorem, or factoring exponent above `1/2` is claimed. The owner requested wrap-up; this session's arithmetic and additive spin-lemma scopes are released in Discussions #53 and #51.

## 2026-11-01: Joint-energy conjugations and polynomial dressing reduction

Exact mathematical sources `795825c`, `c5701ff` and `74dde47` supply a one-place-isospectral joint-pair construction, its complete 1,296-member F4 leading-coefficient-preserving orbit at the five-letter reference's natural heights, a formal leading-rank-one factorization, and an all-field/all-degree reduction for a square-zero polynomial sector. A unipotence argument extends the latter to every rank-one scalar-polynomial energy generator under the stated pole and leading-rank conditions. The Gaussian reference and its authorship remain credited. These are scoped construction/rigidity results; the complete factoring exponent stays `1/2`.

The full orbit replay recomputed 1,298 profiles with fraction-free polynomial arithmetic, checked every endpoint and all integer prefix inequalities, regenerated a 16-dimensional bimodule orbit, and proved completeness by the exact unit counts 186,624/144. The fast companion checkers independently computed the 16-dimensional tangent identity, all 180 quadratic identities and both symbolic leading minors. All focused checks passed. Full `make verify` passed on integration tree `8f7984459fdf01c9a1e21f1a6ed3e75a9822cdb0`; subsequent concurrent integrations retained the exact research source fingerprints and used the documentation and focused gates. This is author computation and an ordinary proof, without independent mathematical review or new Lean work.

Integration repairs: `645fad8` replaces a then-absent Fourier-note citation in `smooth-band-density.md` with the published quadratic-order class-number bound of Adzaga et al., Lemma 4.3, and the explicit comparison to the stated bit-length estimate. Two 124-bit replay files missing from main were restored byte-for-byte from its ancestor `11fc0aa`; the restored checker replayed all 942 supplied norms and recovered factor `2305843009213693951`, with unchanged evidence bytes. This restoration adds no new arithmetic theorem or search-cost claim.

## 2026-11-01: PR #27 conservation checkpoint integration

Exact source head `b49b9e19480c65d139ef15beecacb92064c1a960`; original source attribution and AI disclosures retained. Integrated the missing implementations, evidence and 28 research notes, with current main's half-rate/general-transfer evidence wording preserved. Restored all 30 source default commands, a focused `verify-pr27` target, and the optional native width-five target.

The full-conservation theorem permits every binary conservation kernel in characteristic two on the full cube. The elementary code-span theorem permits arbitrary code pairs and coefficient fields but a narrower kernel space. A related-scope section now records this distinction and the earlier PR #27 full-cube result; neither theorem is treated as superseding the other or as a factoring lower bound.

Full `make verify` passed on the integration tree: 552 Markdown files, 59 unit tests, five crossing tests, all 30 restored source commands and all other default gates. The conservation replay includes 65,812 exhaustive binary function-degree controls and 336 extension-field controls. Optional `verify-code-width5` regeneration and Lean targets were not newly run. Integration is not independent mathematical acceptance, a PR approval review, or an improved complete factoring exponent.

## 2026-11-01: All-rank dressings, formal polynomial rigidity and discrete lattices

Mathematical source `22b977f` removes the rank assumption on scalar-polynomial square-zero energy generators. Its quadratic/Frobenius elimination starts from all 69 variables, treats degrees 1 through 17 exactly, and has a separate argument for every higher degree over every characteristic-two field. The checker regenerates every linear and pure-square consequence and independently controls 58,125 entries with dense F4 products. It does not specialize squares to binary coordinates or depend on the exploratory SAT solver.

Source `7ea5f7d` proves two direct-sum identities, with rank triples `(59,6,65)` and `(53,6,59)`. They classify polynomial tangents of every degree and give the ordinary nonlinear lifting proof for every identity-based formal polynomial gauge retaining the reference pole and leading rank. Such a gauge changes the tensor only by a constant conjugation. Separate components and transformations relaxing the hypotheses remain open.

The discrete-lattice companion classifies every integer coordinate charge preserving polynomiality and the original pole, modulo two neutral support-component shifts. It replays all 12 effective tensors and both complete profiles, 26 fraction-free profiles including the reference. Both profile prefix comparisons hold for every candidate; no numerical optimization is needed for the all-scale conclusion. The noncoordinate idempotent interface is supplied as an open construction target, not an existence result.

Full `make verify` passed for source `08bb5b71c7b6cfc987ed6c773ca4fe8b6cd2b6d5`: 598 Markdown files, 63 unit tests, five crossing tests, and the default research gates. This includes repair of a then-unavailable smooth-power link in incoming arithmetic work. Later main integrations substantially expanded the default suite; that expanded whole suite was not newly rerun. The research dependencies retained their exact fingerprints and the new focused Riccati, formal-tangent, monomial-lattice and documentation checks passed after integration (686 Markdown files). These are author calculations and ordinary proofs, without new Lean, independent review, or a factoring exponent above `1/2`.

The subsequent projector checkpoint closes arbitrary energy-block idempotents and every positive lattice depth. The one-step pole/rank system has exact all-field dimension chain `69 -> 13 -> 2`; for depth at least two the first two cold/infinity jets alone give `69 -> 33 -> 10 -> 3 -> 2`, without a leading-rank condition. The terminal space consists of the two neutral support-component projections, so idempotence gives exactly four solutions and every transformed tensor is unchanged. The standard-library checker regenerates all formal polynomial row spaces and Frobenius-root consequences and controls 20,000 entries against independent dense F4 multiplication. This is a new focused replay and ordinary proof, not a rerun of the expanded whole suite, a finite-field census, or a factoring lower bound. Other finite places and interacting projection flags remain open.

The stronger one-step replay removes the cold-place condition entirely and gives `69 -> 14 -> 2`. It proves the same four-projector conclusion for every scalar rational multiplier `f(t)=t+O(1)` at infinity, including `t-a` at every finite point over every characteristic-two field. Rational multipliers of degree zero at infinity remain a distinct target; the note gives their exact two-finite-place polynomiality interface. The focused projector and documentation checks passed after this extension.

The single-evaluation algebra checkpoint closes that scalar-rational three-place target uniformly. Polynomial Euclidean module generation produces 169 members with a triangular minor exactly `t^1452(t+1)^350`. Thus every finite evaluation away from zero and one, over every characteristic-two field, generates the entire constant coefficient algebra. A graded invariant-subspace argument shows that every nonconstant scalar-rational projector gauge regular and invertible at zero, one and infinity which yields a polynomial seed is trivial. The checker retains the complete polynomial basis, processes all 632 generated members, replays 942 final closure controls and 3,620 coefficient-algebra memberships, and independently generates the F4 point algebras of dimensions `19,11,169,169`. The all-field result follows from the symbolic minor; finite point dimensions are controls. No new full-suite, Lean, independent-review or stronger factoring-bound claim is attached to this focused checkpoint.

The same determinant also proves generation over `K[t,1/t,1/(t+1)]`. This yields global rational rigidity for arbitrary energy-block matrix gauges regular and invertible at zero, one and infinity, with no projector or formal-family restriction: conjugation sends each constant coefficient-algebra member into that localized ring, and regularity at the three distinguished places makes it constant. Every such gauge factors as its constant value at zero times a full-commutant rational unit. The checker additionally verifies full rank of the selected constant-algebra coordinate map. Gauges singular at one remain outside the theorem even if their transformed crossing endpoint survives. Full `make verify` subsequently passed on isolated exact source `e93fb4732c5554da0a06ef600a5e35ed27080c4e`: 689 Markdown files, 63 unit tests, five crossing tests and all default research checks. The later global-rigidity extension `3657fa7` retained its separate focused replay and documentation check; the full gate is pinned to the stated earlier source.

The singular-endpoint checkpoint constructs genuine new strict crossings outside the regularity theorem. The first endpoint-jet commutant has dimension five and a three-dimensional square-zero sector; modulo two tensor-preserving directions every regular-at-zero, identity-at-infinity rational gauge in that sector reduces to three principal parts at one. The all-field proof uses 114 quadratic identities, commutator order exactly four, a two-square minor `z(1+z)`, and a nine-square minor with formal-parameter coefficient valuations `52,34,18,4`. A degree-69 pressure difference factors as `x^3 Q(x)` with all 67 scaled degree-66 Bernstein coefficients positive. This proves typed Smith parameter at least `3/4` for every nontrivial member over every characteristic-two field at the inherited heights. A separate complete F4 replay checks all 64 triples and 130 fraction-free profiles; all are prefix-dominated by the reference. The focused checker and documentation gate are the new checks, without a new whole-suite, Lean, independent-review or factoring-score claim.

The all-rank rational extension removes both the initial endpoint-commutant assumption and the pole-location restriction for a fixed constant square-zero generator. At any finite point other than zero and one, the top Riccati equation makes its graded kernel and image invariant; single-evaluation generation then forces the generator into the full commutant and the tensor is unchanged. At one, seventeen exact quadratic/Frobenius reductions and the separate degree-at-least-18 argument force every rank and pole order into the same endpoint sector. The checker retains every stage constraint and terminal basis and independently controls 58,125 new dense F4 entries. Thus every nontrivial normalized scalar-rational square-zero energy gauge has the same three-parameter form and the all-field `rho>=3/4` floor. Higher-nilpotency generators and general rational matrices singular at a distinguished place remain open.

The higher-nilpotency monomial checkpoint uses the maximal energy-block dimension four to bound the nilpotence index. Exact expansion makes polynomiality and the strict endpoint equivalent to `P^4=0` and the displayed degree-at-most-four equations. Seventeen complete linear/pure-square/pure-fourth-power eliminations force the endpoint sector at depths one through three and the full stabilizer at every larger finite depth; the separate depth-at-least-18 argument closes the remaining integers. All terminal spaces satisfy the original equations identically, and 47,500 entries are independently checked by dense F4 arithmetic. The scope is all characteristic-two fields and all nilpotent generators for a monomial multiplier; arbitrary rational multipliers of a higher-nilpotency generator remain open. The new evidence is the focused checker and documentation replay, not a new whole-suite or factoring-score claim.

The filtration checkpoint enlarges constant zero-profile-preserving gauges from exact energy blocks (dimension 69) to the energy-filtration group (dimension 347). The regular joint infinity operator has invertible constant term, so left multiplication preserves minimum entry valuation and gives the exact smallest infinity order `-8-max(H_r-H_c: (U^-1)_rc nonzero)`. Its opposite filtration intersects the zero-preserving one exactly in the energy-block group. The proof is global over every characteristic-two field. Independent controls regenerate the 347-variable tangent kernel and all 834 strict-filtration F4 elementary seed endpoints and exact entry poles. Spending a worse pole budget and source-specific zero-profile preservation outside this filtration remain open; no numerical screen is used as a whole-group cost proof.

## 2026-11-01: Single-bit weight-one coefficient decoder and session handoff

Exact sources `546637a` and `4c2f733` give the ordinary precision reductions and the binary decoder in [the weight-one note](../research/weight-one-divisor-coefficients.md), Section 4.1. For every odd prime level, one bit per cyclotomic coordinate suffices; the full algebra modulo two retains essential nilpotent information. The proof establishes injectivity by cyclotomic root multiplicities, and the implementation uses polynomial-time binary Gaussian elimination. Full integral normalization needs only modulus four. The paired CRT reconstruction credit to Discussion reply 19324705 and the Garg--Schost interpolation attribution are retained.

The new finite controls check 520 binary conversions from 910 supplied coefficient vectors at eight prime levels. Existing exact and division-based conversions and the separate paired support reconstruction remain checked. The finite coefficient provider enumerates known divisors; no faster composite-index coefficient engine or factoring exponent above `1/2` is supplied.

Full `make verify` passed on integration head `0c4777558d85490376ada8b2e737e82820ed0c26`: 530 Markdown files, 59 unit tests, five crossing tests and the default research gates. The tested source fingerprints matched after completion. Subsequent concurrent integration preserved the coefficient decoder, checker and reconstruction dependencies byte-for-byte; final focused and documentation gates passed. The final tree includes additional concurrently integrated default checks that were not rerun as a complete suite by this session. For the class-number citation conflict, retained the concurrent published-reference repair `645fad8` rather than duplicating the local derivation from `62b1aaa`.

The owner requested wrap-up and paused the campaign. The coefficient scope and the unfinished Kronecker/Fay exploration are released in Discussion #51. The latter has no new repository theorem or fast binary-index extraction algorithm. No new Lean build or independent mathematical review is claimed. The complete factoring score remains `1/2`.

## 2026-11-01: PR #29 Ward and crossing checkpoint integration

Exact source head `7154cec09b982b54a0fcc0bb7a303318fc1f90da`, retaining original source attribution and AI disclosures. Restored 434 missing files and indexed 83 research notes. Preserved main's current result/evidence wording. The [source README](https://github.com/1337isnot1337/integer-factoring/blob/7154cec09b982b54a0fcc0bb7a303318fc1f90da/README.md) retains its original detailed inventory. Polynomial entry algorithms, conditional endpoint transfers and construction-specific barriers are not represented as a new numerical factoring bound.

Restored the 53-command `verify-constructions` default dependency, the optional `verify-native-crossings` target, and 25 additional optional post-half commands under `verify-pr29-extended`, preserving the previous post-half checks. Full `make verify` passed on the integration tree: 678 Markdown files, 63 unit tests, five crossing tests, all 53 new default construction commands and all existing default gates.

Not newly run: the larger native and optional extended screens or Lean targets. The source implementations and evidence for those remain available without a fresh-run claim. Ordinary `git diff --check` reports 21 source-inherited blank-at-EOF warnings; those bytes were preserved to retain source fingerprints. The same whitespace check with only blank-at-EOF warnings disabled passed. No other whitespace errors were ignored. Integration is not independent mathematical acceptance or a PR approval review; the complete time exponent remains `1/2`.

<!-- beyond-half-review-20261101:start -->

## November 1: smooth-norm processing and constructive selector checkpoint

The [research guide](../research/beyond-half-research-checkpoint.md) connects the
ordinary arithmetic reductions, constructed power selector, adaptive 144-bit
factorization, exact channel improvements, and eight-letter nonorthogonal
crossing. The complete factoring exponent remains 1/2. No new Lean formalization
or Lean rebuild is claimed.

The [verification record](../evidence/beyond-half-verification-20261101.json)
records the 195-command baseline suite at d7913c7, all 30 verification additions
at d42e464, focused core replays, source hashes, and later integration checks.
The full 144-bit search was repeated, including both its exhausted first phase
and 154,884-query continuation. Independent agents reviewed the scope and replayed
the accepted norms, exact relations, partial extraction, and both crossing cuts.

Uniform generation and represented-norm existence use the pinned Cureton--Thorner
all-number-field finite-order Hecke strip. A fast deterministic norm selector,
a uniform short power-window bound, and a fast dynamic prime-closure bound remain
unproved. Small experimental bases are not presented as satisfying the enormous
uniform analytic cutoff. Partial failure is not a primality certificate, and
passing finite controls is not a proof of a stronger worst-case exponent.
<!-- beyond-half-review-20261101:end -->

## 2026-11-01: PR #30 recovery and Smith integration

Exact source head `28105165edb8e9a0427468dc833f3c09bc0a7fb7`, with original source attribution and AI disclosures retained. Restored the 96 missing files and retained its distinct integral Schur derivation verbatim as `research/integral-schur-crossing-pr30.md`, rather than overwriting the broader PR #23 derivation. The two new dependent notes point to that companion; the existing Schur note cross-links it. The 22-note checkpoint index identifies the source. The k6 status sentence now correctly qualifies its separate `1/30` endpoint as k6; current half-rate evidence wording is unchanged.

Restored parity/recovery/constructive default targets and retained optional native, Schur-rank, exterior and density targets. The shared `verify-polynomial-native` target keeps its previous command and gains the source checker. Full `make verify` passed on the integration tree: 710 Markdown files, 63 unit tests, five crossing tests, the 21 new parity/recovery/constructive commands and all existing default gates. Ordinary `git diff --check` passed.

Not newly run: the larger optional native, Schur-rank, exterior/density or Lean targets. This integration preserves ordinary construction proofs and their finite controls; it is not independent mathematical acceptance, a PR approval review, or an exponent above `1/2`.

## 2026-11-01: PR #31 scalar transfers and construction controls

Exact source head `b1ba452fa4a3d19fe6fea3822c53613f3b582104`, preserving original attribution and AI disclosures. Restored 119 missing files and indexed 32 notes. The fixed-prime scalar transfer and the parallel PR #30 fixed-characteristic recovery derivation share the power-tuple and cyclic-orbit mechanism; both proofs remain available. The certified divisor oracle maximizes the divisor count below an arbitrary bound, rather than computing the divisor count of an arbitrary supplied integer. Its stated Nicolas small-benefit and per-record gap inputs were checked against the primary paper; this is not independent acceptance of the full algorithm or logspace proof.

Preserved all previous default recipes and added 27 construction controls and six frontier controls through separate dependencies. Full `make verify` passed on the integration tree: 752 Markdown files, 63 unit tests, five crossing tests and all default commands. The README subsequently gained a parallel-proof cross-reference and the documentation gate was rerun. The all-field joint-evaluation algebra checker also passed a separate focused replay: 632 polynomial members, 942 closure controls and 3,620 coefficient memberships.

The optional large divisor replay, larger searches and Lean targets were not newly run. Ordinary `git diff --check` reports three trailing-space warnings in `patches/certified-divisor-cli.patch`: these are the required blank context lines of the retained source patch, not implementation whitespace changes. No other whitespace warnings were reported. Integration is not a PR approval review or independent mathematical acceptance. No new fast evaluator or exponent above `1/2` is claimed.

## 2026-11-01: Combined exact-source joint-gauge verification

Full `make verify` and all eight new focused joint-gauge targets passed together in a detached, isolated worktree at `21c7cfb3d507a388ccc788ac301f8dc31d43ebfe`. The full gate checked 764 Markdown files, 63 unit tests, five crossing tests and all default research commands. The focused phase covers the all-rank Riccati reduction, formal polynomial tangent theorem, monomial lattices, arbitrary projector depths, single-evaluation/global rational theorem, all-rank endpoint loops, higher-nilpotency monomial depths and opposite-filtration rigidity. Exact commands, source fingerprints and the local log hash are retained in [the verification record](../evidence/joint-research-verification.json). The verification worktree was clean after completion. Later main integrations are outside this whole-suite pin. This is local author computation, without new Lean, independent mathematical review, leaderboard submission or an exponent above `1/2`.

## 2026-11-01: PR #46 constructor and physical-code integration

Exact source head `8d3a2da11eca5e1a7d064b645c698719ac141098`, retaining source attribution and AI disclosures. Restored 448 missing files and indexed 92 research notes. Current main files and evidence-status wording were preserved. The 48 source verification targets, containing 93 optional commands, remain available without replacing existing targets or making large searches mandatory defaults.

The full default `make verify` passed on the integration tree: 851 Markdown files, 63 unit tests, five crossing tests and all default research commands. Fourteen focused commands also passed: the exact saved-profile 96-letter pressure and rational-type certificates, six prime-constructor checks, two bilinear Cauchy controls and three finite-crossing transfer checks. The 96-letter pressure replay is exact arithmetic on the saved complete profile, not a new construction of the large matrix.

Three additional portable archive commands passed. The root archive retains 6,722 physical candidates and replays nine winners plus seven other scopes. The final archive reconstructs 111 entrywise and 168 finite-frame witnesses, including 182,016 endpoints. The unequal-code replay checks 18 tensors and 1,037,568 endpoints, with 108 numerical cost calculations and 18 exact rational comparisons; none passes the strict threshold. Numerical optimization is not an optimality proof. Ordinary `git diff --check` passed.

The remaining optional targets, the approximately 25-minute full 96-letter profile construction and Lean were not newly run. The complete factoring exponent remains `1/2`; the internal pressure bound `rho < 0.710005874881191` is still above `1/sqrt(2)`. Integration is not independent mathematical acceptance or a PR approval review.

## 2026-11-01: Personal F4 hyperoval crossing certificate

The explicit rank-three F4 hyperoval embedding gives exact typed rho < 710005861/10^9, improving the pinned coordinate q96 upper bound 0.710005874881191. It does not improve the complete factoring exponent 1/2. Both full 9,216-channel profiles were constructed in the recorded approximately 1,551-second run, with all source pins matching and 35,973,996 occupation-block endpoint entries checked. A separate subset-determinant implementation checked 76 complete F4 entries using 107 independent determinants.

Exact outward arithmetic passes at 192 and 224 bits, including an explicit positive type of denominator 2^36 and minimum count 2,813,121. A smaller same-base proposed bound fails as a negative control. The candidate does not dominate the coordinate profile at every scale; exact rational sign controls retain that limitation. Focused entry and pressure gates passed; docs passed on the local package. These quick gates join the defaults, but the large profile construction does not. No new Lean run or independent mathematical acceptance is claimed. The full default suite was not yet rerun after this package. Substantial OpenAI assistant support.

## 2026-11-01: Complete symmetric binary coefficient-envelope exclusion

The literal q3, theta4 fused seed was matched entry by entry to a fresh run of the existing finite fusion constructor, including invariance and the strict endpoint. Its symmetric two-copy scalar-projected coefficient envelope has binary dimension 21. The complete 2^21-vector screen finds 17 nonzero tensors with rank product at most 501, spanning dimension six and annihilated by a functional detecting the crossing. An independent column-elimination pass on all 2^20 vectors of that functional's affine-one hyperplane finds none below the threshold. Python separately checks all 34 retained full-matrix ranks. Exact integer comparisons put sqrt(502) above the nine-letter critical cost.

`make verify-fused-scalar-projection` passed after evidence generation. This excludes every binary linear decomposition in the stated symmetric envelope, rather than selected polynomial coefficients. It does not settle the ordered 33-dimensional space, extension-field decompositions, higher powers or physical type restrictions. The whole repository suite was not newly rerun for this checkpoint. No new Lean, independent review or exponent above `1/2` is claimed. Substantial OpenAI assistant support.

## 2026-11-01: All-length primary-character code obstruction

An explicit simultaneous primary decomposition of the q3 theta4 opposite-cut coefficient operators has block dimensions (1,2,1,1,2,2), with all 54 within-block radical pair products zero. Ten complement-rank certificates prove weighted character-code distance three over every extension of F16. The ordinary tensor-code argument then gives opposite-cut rank at least 3^m for every nonzero vector in the full unprojected coefficient tensor-product span, at every length m. Rank subadditivity on the other cut forces every crossing decomposition in that span to have cost at least 3^(3m/2), hence rho at least 3/4.

The focused primary-code checker passed. Its exact finite algebra supports the stated all-length proof; scalar projection and physical type restrictions are explicitly outside that theorem. The earlier binary envelope checker addresses a separate finite projected space. No new factoring bound, new Lean or independent review is claimed. Substantial OpenAI assistant support.

## 2026-11-01: Global energy-refined Newton ceiling

A constant invertible diagonal over F16, scalar on each joint-energy block, commutes with every rational energy gauge. The full symbolic polynomial det(I+z D D_H A) therefore yields invariant ceilings for all zero Smith prefixes, without assumptions on gauge degree, polynomiality or endpoint behavior. The computed Newton slopes are (8,10,10,12,14,14,16,16,16,16,18,18,18,18,20,20,20,20,20,22,22,22,22,24,28). Constant-field extension makes the proof apply to every characteristic-two field. This improves the ordinary untwisted Newton ceiling but does not establish attainability or control the other cut.

The focused checker passed after generation: complete symbolic coefficient polynomials, 256 independent finite-field determinant evaluations and agreement of the top coefficient with a separate fraction-free determinant. The hypothetical pressure comparison in the note assumes the reference infinity profile and is explicitly not an attained bound. A separate finite search of 7,243 strict square-zero monomial zero-loop candidates found no improvement at inherited heights; it is exploratory and not a uniform exclusion. The full suite at 38bb52b is still running; no new full-suite pass is claimed for this source. No new Lean or exponent above `1/2`.


## 2026-11-01: explicit cyclotomic character decoder

Mathematical and implementation source `e51a9c1aabd8a428daa4612545f1b11610af67b9`; compatibility source `b2ea123`. The [ordinary proof](../research/cyclotomic-character-decoder.md) constructs an inverse of the existing one-bit character map by Galois substitutions, odd-order Fourier inversion and exact nilpotent divisions. It avoids a stored Gaussian system. The previous root-multiplicity proof supplies the exact kernel calculation used for the omitted-character witnesses; the earlier reduction and sparse reconstruction retain their attribution. This checkpoint uses OpenAI Codex assistance.

[Focused evidence](../evidence/cyclotomic-character-decoder.json) compares the new decoder with the older integer cyclotomic arithmetic and Gaussian decoder on 1,036 binary vectors at 17 prime levels, plus 120 general even-order vectors, 504 actual divisor vectors, 114 omission witnesses and 11 rejected invalid inputs. The full nilpotent algebra is retained, including length 128. The sampling necessity result is scoped to arbitrary paired binary vectors, not actual divisor inputs or factoring lower bounds.

The [verification receipt](../evidence/character-decoder-verification.json) records a clean frozen worktree at the exact source above. `make -j2 --output-sync=target verify` passed: 335 Python command blocks, 63 unit tests, five crossing tests and local links in 890 Markdown files. The source hashes and compressed local log are retained. The initial exploratory run on the changing checkout was terminated before completion and is not counted as a pass. Eight targeted overlap-tool tests and the complete live remote overlap scan also passed after replacing unsupported `gh api --slurp` with incremental parsing of paginated JSON arrays.

The final follow-up clarifies a reused symbol and attribution in the note and adds this receipt; all five pinned implementation/test sources are unchanged. Later concurrent main changes are outside the frozen full-run claim. The GitHub Actions run for the source was cancelled; no completed CI pass or new Lean build is claimed. The fast coefficient engine remains open, and the complete factoring exponent stays `1/2`.

The concurrent restart-option change in `search_f4_recentered_crossings.py` changes only its `main` function. An exact AST comparison against 38bb52b, excluding that CLI function, confirms that all imported module code is identical. The affected joint evidence records retain their entire mathematical payload and refresh that source fingerprint. The new graded-Newton focused replay passes on the integrated source. Complete focused recomputations are also running as a redundant integration check; they are not yet recorded as completed passes here.

## 2026-11-01: Common-core natural-pressure floor

A complete 92-label core profile and exact Toeplitz-rank partition monotonicity give a lower channel sum for every stated 96-label boundary retraction over characteristic two. An integer degree-2827 derivative, exactly two Descartes sign variations, exact root-box signs and a direct integer-power minimum-box comparison prove rho >= 17/24 at every natural pressure scale and physical type. This excludes only the stated natural-pressure method and its boundary-supported cone-weight corollary, not general crossing decompositions or factoring algorithms. The core construction completed in about 516 seconds; its full source fingerprints match. The fast all-scale certificate and local docs pass. No new Lean or independent acceptance.

Follow-up verification of the earlier hyperoval package: the full default suite passed at exact published 85b5ee5c7ef0261884fb7875930e0c1750d416f7, including 880 Markdown files, 63 unit tests, five crossing tests and all defaults. Later concurrent main changes are outside that complete-run pin. Substantial OpenAI assistant support.


## 2026-11-01: theta-section engine and Tate-coordinate construction

Theta source `48ecca8` and markup-only repair `4d73ddd` construct signed divisor coefficients by a theta quotient with a completely charged binary section recurrence. Its current cost is exponential. A uniform identity-submatrix proof forces common depth-k linear section spaces to have dimension at least floor(2^k/36), while the q=3 coefficient has a polynomial-time square test. This is a representation restriction, not a factoring lower bound.

The [frozen full-suite receipt](../evidence/theta-section-verification.json) records a clean completed `make -j2 --output-sync=target verify` at `4d73ddd662fa8c8cccd32dc783019a5687ec0a1c`: 337 Python command blocks, 63 unit tests, five crossing tests and 919 Markdown files. An earlier launch before worktree checkout completion failed immediately and is explicitly excluded. Later main integrations and the new Tate compiler are outside that completed whole-suite pin.

The [Tate-coordinate proof](../research/tate-coordinate-coefficients.md) constructs the intrinsic degree-(q-1)/2 nodal algebra and a polynomial-size rational-diagonal compiler for the new variable v. Its coefficients have poly(q,log(n+1)) cost. At odd original indices the full target is X_q(Delta(t)); the Delta composition remains open. Its algebra trace has the separately proved original-index square-test formula and Artin--Schreier identity. The remaining trace-zero part is not supplied. Standard Tate uniformization, division-polynomial, Jacobi and Furstenberg results retain their primary-source attribution. Substantial OpenAI Codex assistance.

[Exact focused evidence](../evidence/tate-coordinate-compiler.json) includes 65 nodal identities, nine nodal algebras and embeddings through prime 31, 221 Hensel/diagonal coefficient controls, finite Delta compositions including 108 odd original indices, 128 small binomial comparisons, 30 internal-index controls at q=3 through 1,024 bits and 54 large original-index trace controls. The large internal controls do not compute the full original-index divisor data. The original-index large controls check only one trace. No new Lean or independent acceptance, and no complete factoring exponent above 1/2.


## 2026-11-01: Tate full verification and nonlinear phase-five checkpoint

The [Tate verification receipt](../evidence/tate-coordinate-verification.json) now records a completed clean frozen default run at `c87e49f52f47f93108ed0d12b72ce95536ef236a`: 338 Python command blocks, 66 unit tests, five crossing tests and 939 Markdown files. All seven pinned source hashes match the checked files retained here. The local compressed log is preserved. Later phase-five files and concurrent main additions are outside that complete-run pin.

Phase-five source `60c0cd9` proves `X=B(c+V)`, `B=1+V+V^2`, with V the sum of the positive square series at scales one and five. Four binary quadratic-form parities give an exact nonlinear original-index evaluator with square-root-time cost, still exponential in bit length. The [focused receipt](../evidence/q5-quadratic-coefficients.json) compares 8,192 literal odd-index coefficients and records 1,414,130 charged quadratic-enumeration steps, 64 Tate prefix controls, 946 distinct-prime semiprimes, 12 large cheap controls and 8 large controls with an explicitly supplied divisor-count residue. The latter do not implement an unrestricted fast tau oracle. Focused replay and the documentation gate pass.

The independent q=5 divisor-class derivation in [Discussion reply 19327344](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19327344) is credited, with its alternative coefficient basis made explicit. Its reported larger replay is not counted as a new local check. No new Lean build, independent acceptance of a factoring improvement, or exponent above 1/2.

## 2026-11-01: Zero-loop reachability and a coupled pencil

The new zero-loop checker proves over every characteristic-two field that no strict gauge I+P/t, with any constant square-zero energy-block P, can keep the first zero order eight while raising the second to ten. Necessary equations reduce 69 to 41 to 36 generator parameters. A general rank-one leading matrix and its tangent add 97 coordinates; exact affine/Frobenius reduction goes 133 to 60 to a constant contradiction. This proof permits a changed leading coefficient and has 13,125 independent dense F4 entry controls. The same checker verifies the identities for an eleven-parameter polynomial pencil coupling poles at zero and one. Its finite 198-case inherited-height search is exploratory and found no improvement. No new Lean, independent mathematical review or factoring exponent above `1/2` is claimed.

## 2026-11-01: Complete all-field obstruction for the coupled pencil

The full eleven-parameter two-place pencil now has an ordinary all-field obstruction: every nonzero q of degree at most eleven with q(1)=0 has typed Smith rho at least 8/11 at the inherited physical heights. The zero polynomial is the original reference and is excluded from that floor. The normalized zero matrix is a constant rank-two perturbation S+qJ. Twenty-five symbolic minimum-minor witnesses show that no polynomial q improves its zero profile at any scale. Twenty-five infinity witnesses, unique-valuation terms and unchanged finite jets give eleven degree-specific partition lower bounds. After upward integer rounding, eleven integer pressure polynomials have positive Bernstein certificates, every scaled coefficient at least 22. The rational AM-GM bound proves the rho floor without numerical optimization; 2*8^2>11^2 puts it above 1/sqrt(2).

The focused zero-prefix and complete-pencil gates pass. The new gate reconstructs all fifty symbolic minors and checks them with 800 independent F16 Gaussian determinants, then recomputes all integer pressure certificates. Other heights, actual coefficient recombinations and other generator directions remain open. The full default suite was not rerun for this follow-up. No new Lean formalization, independent mathematical acceptance or improved factoring score is claimed. Substantial OpenAI assistant support.

## 2026-11-01: Coefficient checkpoint full gate and fingerprint replay

The complete default suite plus both fused coefficient-space targets passed in a clean detached worktree at exact source `38bb52ba47f5195d0c43b2b6406449e543bb2027`: 895 Markdown files, 63 unit tests, five crossing tests and 335 default Python command blocks. The [receipt](../evidence/coefficient-research-verification.json) retains exact commands and the local log hash. This full-suite pin precedes the Newton and zero-prefix additions and later main integrations. All twelve joint-gauge records were also recomputed after the CLI-only fingerprint change; every mathematical payload matched exactly. The [comparison receipt](../evidence/joint-fingerprint-refresh.json) records that separate result.
