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

## Maintaining this record

For each later decision, record the exact reviewed head, changed-file coverage,
checks actually run or reused, accepted scope, exclusions, source attribution,
retaining commit or archive tag, and actual GitHub disposition. A snapshot or
successful finite replay alone must not be recorded as mathematical acceptance.
When a PR changes, retain the historical checkpoint and review the delta.
