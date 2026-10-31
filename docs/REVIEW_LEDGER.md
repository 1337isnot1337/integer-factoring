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
