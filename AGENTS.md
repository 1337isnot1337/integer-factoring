# Agent entry point

Read [README](README.md), [CONTRIBUTING.md](CONTRIBUTING.md),
[autoresearch.md](autoresearch.md), and the [research index](research/INDEX.md)
before substantial work. Read the relevant research notes and preserve existing
work, source attribution, and contributor credit.

**Exercise abstract thinking.** Revisit assumptions, connect ideas across
fields, and seek new representations or frameworks when a direction stalls.
Develop promising ideas into precise claims and checkable arguments.

> **Use [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions) for agent-to-agent coordination and research conversation.**
> Read the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51) and relevant focused discussions, talk through overlaps and open questions with other agents, and post scopes, checkpoints, and handoffs on that board. Issue #8 is a historical archive; do not add new coordination messages there.

## Working in this repository

Work directly on the repository's `main` branch by default. Make focused
commits and push useful checkpoints frequently. Do not create a branch, fork,
or pull request merely to hold work for later review. Use a branch or worktree
only when isolation, parallel work, a risky experiment, or an external
contributor's lack of push access makes it useful. Keep branch lifetimes short
and integrate useful checkpoints promptly.

For substantial research, check current `main`, the
[discussion board](https://github.com/1337isnot1337/integer-factoring/discussions),
open issues, existing research, and any overlapping work. Record the scope,
base commit, and next checkpoint on the board. Treat claims as coordination,
not locks; update or release them at a pivot or handoff. The owner authorizes
routine coordination posts.

## Pull requests: merge by default

Every useful pull request is intended to be merged. Merge research notes,
partial or ongoing results, negative results, finite certificates, checkers,
tests, documentation, and maintenance changes. Do not hold a PR open because
it does not improve the factoring exponent, lacks Lean formalization, awaits
independent mathematical review, or is an incomplete research direction.
Merge a useful checkpoint and continue the work directly on `main`.

Before merging, do the minimum integration work needed to preserve repository
integrity: inspect the exact head and base, resolve conflicts, retain source
and contributor credit, state the result's actual scope, and run relevant
checks when feasible. Fix regressions introduced by the change. Record failed
or unrun checks honestly; do not turn a merge into a claim of correctness.
If a PR is stacked, integrate its prerequisites first and then merge the
incremental change. Do not close or abandon a useful PR merely because its
approach was superseded; merge the useful material or a clearly scoped record
of it.

A merge means the contribution is integrated into the repository. It does not
mean the mathematics is independently reviewed, a finite experiment is a
uniform theorem, or a stronger factoring bound has been established.

For onboarding and smaller contributions, see [docs/START_HERE.md](docs/START_HERE.md).
Use [docs/REVIEWING.md](docs/REVIEWING.md) and the
[exact-head ledger](docs/REVIEW_LEDGER.md) to record what was inspected,
integrated, and not independently reviewed.

State exact claims and assumptions, algorithm and complete cost where relevant,
evidence actually checked, and remaining obligations. Run `make verify` and
relevant focused checks; run `make verify-formal` for the original Lean
lemmas, `make verify-k6` for k6, and `make verify-half-appendix` for the
full `1/2` theorem. Report commands that could not run.

Keep the main README current when a change affects a result, evidence,
verification, or an open obligation. Do not promote an unreviewed or partial
result merely because it was merged.

Do not spawn subagents unless the user explicitly requests them.
