# Reviewing contributions

Independent reproductions, mathematical audits, and attribution corrections
are welcome. You can review one claim or one checker without reviewing an
entire research branch. Record the scope you actually checked.

## A focused review

1. Read the [README](../README.md), relevant note, source pins, and
   [ledger](REVIEW_LEDGER.md). Check the live PR head and its base, including
   whether the PR is stacked on another branch. Record the exact commit.
2. Identify the claim and inherited assumptions. Compare every changed file
   in your review scope. Note files or dependencies you did not inspect.
3. Run the relevant reproduction commands in an isolated checkout. Record
   environment and outcomes. Prefer exact arithmetic for exact claims; check
   a meaningful negative control where it can expose the asserted failure mode.
4. For a bound, inspect uniform construction, all levels, bit costs, parameter
   choices, reduction, recovery, and finite exceptions. Finite replay alone
   cannot supply these arguments. For Lean changes, use the matching pinned
   toolchain and check the theorem statement, dependencies, and axiom report.
5. Leave a specific finding: what passed, what failed, what remains unknown,
   and the next action. Credit the author and any reused review evidence.

## Evidence and disposition

| Description | What it establishes |
| --- | --- |
| Inventoried | Head and changed files recorded; no correctness review implied |
| Finite replay passed | The stated finite checks passed in the recorded environment |
| Scoped review | The listed argument or implementation was reviewed; exclusions remain explicit |
| Lean kernel checked | The named theorem checked against its stated assumptions and axioms |
| Incorporated | Identified work is present on `main`, with source and credit retained |
| Archived | Work is retained at a permanent reference; it is outside the active queue |

These descriptions can apply together. Incorporation or archival does not
certify a new exponent. The [ledger](REVIEW_LEDGER.md) records specific decisions,
with links to the evidence and the remaining obligations.

## Maintainer closure and branch cleanup

Before closing a PR as incorporated, recheck its head and verify coverage of
all changes. Link the retaining commit and any adaptations; leave unresolved
work open or name exactly what was excluded. A changed head needs a new review.
Do not close by title, age, shared topic, or green CI alone.

An obsolete approach may instead be archived without correctness acceptance.
Explain why it leaves the active queue, preserve its exact head under an archive
tag, and credit its source. Unique useful tests or unresolved proposals can stay
open even when their original headline is superseded.

Delete a remote branch only after checking that its exact head is merged or
archived, its PR is resolved, and no open PR uses it as a head or base. Recheck
the live branch SHA immediately before deletion and use an expected-head guard.
Preserve local worktrees and ongoing work. GitHub issue/PR closure is a separate
step from evidence review and should be verified after the action.

This workflow adapts the contribution, exact-head review, scoped-result, and
credit practices of
[CrocSwap/integer-mult-bounds at `2485805a`](https://github.com/CrocSwap/integer-mult-bounds/blob/2485805a6d47f36ea83ddde6327d453ac82d44dc/CONTRIBUTING.md).
Its larger batch machinery is unnecessary for this repository's current queue.
