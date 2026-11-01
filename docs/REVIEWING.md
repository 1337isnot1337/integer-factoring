# Reviewing contributions

Review helps improve the repository, but independent mathematical review is
not a prerequisite for integrating useful work. Review one claim, one checker,
or one implementation at a time and record exactly what was inspected.

## A focused review

1. Read the [README](../README.md), relevant note, source pins, and
   [ledger](REVIEW_LEDGER.md). Check the live PR head and base, including
   whether the PR is stacked. Record the exact commit.
2. Identify the claim and inherited assumptions. Compare the changed files in
   the scope being reviewed. Note files or dependencies not inspected.
3. Run relevant reproduction commands in an isolated checkout when feasible.
   Record environment and outcomes. Prefer exact arithmetic for exact claims;
   check a meaningful negative control where it can expose the asserted failure
   mode.
4. For a bound, inspect uniform construction, all levels, bit costs, parameter
   choices, reduction, recovery, and finite exceptions. Finite replay alone
   cannot supply these arguments. For Lean changes, use the matching pinned
   toolchain and check the theorem statement, dependencies, and axiom report.
5. Leave a specific finding: what passed, what failed, what remains unknown,
   and the next action. Credit the author and any reused review evidence.

## Evidence and integration

| Description | What it establishes |
| --- | --- |
| Inventoried | Exact head and changed files recorded; no correctness review implied |
| Finite replay passed | The stated finite checks passed in the recorded environment |
| Scoped review | The listed argument or implementation was reviewed; exclusions remain explicit |
| Lean kernel checked | The named theorem checked against its stated assumptions and axioms |
| Incorporated | Identified work is present on `main`, with source and credit retained |
| Independently reviewed | The stated mathematical or implementation scope received independent review |

These descriptions can apply together. Incorporation or merging does not
certify a new exponent. The [ledger](REVIEW_LEDGER.md) records integration and
review scope separately.

## Merge policy

Merge every useful PR by default, including partial or ongoing research,
negative results, finite certificates, tests, documentation, and maintenance.
Do not leave a PR open because it lacks a stronger factoring exponent, Lean
formalization, a complete argument, or independent mathematical review.
Resolve stacking and merge conflicts, preserve attribution, state limitations,
and run relevant checks when feasible. If a PR is still being developed, merge
the useful checkpoint and continue on `main`.

A merge is an integration decision, not a correctness verdict. The README and
research notes must continue to distinguish proposed, finite, conditional,
formally checked, and independently reviewed claims.

## Branch cleanup

Prefer direct commits to `main`; use short-lived branches only when isolation
or external contribution requires them. Delete a remote branch only after
checking its exact head, confirming no open PR uses it as a head or base, and
rechecking the live branch SHA immediately before deletion. Preserve local
worktrees and ongoing work. Do not delete branches as part of routine PR
merging unless they are demonstrably unused.

GitHub's automatic head-branch deletion is enabled for merged PRs. The
[branch cleanup workflow](../.github/workflows/cleanup-branches.yml) also runs
every six hours to remove remote branches whose exact tips are ancestors of
the current default branch and whose latest commits are at least 24 hours old.
It preserves the default branch, protected branches, and all open PR heads and
bases. It rechecks live tips and PR references, then uses an atomic push with
exact SHA leases so concurrent updates abort deletion. Unique unmerged commits
are retained; squash-merged branches are handled by GitHub's merge setting.
Local branches, worktrees, and uncommitted files are outside this automation.

Preview with `python3 tools/cleanup_branches.py`; add `--apply` to delete.
The workflow's manual Run workflow button previews by default; check `apply`
to delete eligible branches. Both commands require authenticated Git and `gh`.

This workflow adapts the contribution, exact-head review, scoped-result, and
credit practices of
[CrocSwap/integer-mult-bounds at `2485805a`](https://github.com/CrocSwap/integer-mult-bounds/blob/2485805a6d47f36ea83ddde6327d453ac82d44dc/CONTRIBUTING.md).
