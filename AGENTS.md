# Agent entry point

Read [autoresearch.md](autoresearch.md), [CONTRIBUTING.md](CONTRIBUTING.md),
the [README](README.md), and the [research index](research/INDEX.md) before
substantial work. Read the relevant research notes for your task. Follow the user's assigned task and preserve existing work in progress.

Before starting substantial research or a new direction, check current `main`,
the [coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8),
open Issues and pull requests, and existing research for overlap. Post your
name or agent, scope, base commit, branch or worktree, and next checkpoint on
the board. Update the claim at a pivot or handoff, link the evidence and PR,
and explicitly release it when you stop. Coordinate overlapping work; old
comments are not permanent locks. Use a focused branch or worktree and retain
source attribution. The owner authorizes these routine coordination posts.

For onboarding and smaller contributions, see [docs/START_HERE.md](docs/START_HERE.md).
Small corrections and reproduction reports can go directly to an issue or PR.
For review and queue cleanup, use [docs/REVIEWING.md](docs/REVIEWING.md) and
the [exact-head ledger](docs/REVIEW_LEDGER.md). Record source, evidence, credit,
and unresolved scope separately from GitHub merge or closure state.

All pull requests must be ready for review. Do not create draft PRs or convert
open PRs to drafts. Keep ongoing work reviewable and describe incomplete work
and remaining obligations in the PR body; ready for review does not mean ready
to merge. See [CONTRIBUTING.md](CONTRIBUTING.md#pull-request-status).

State the exact claim and its assumptions, the algorithm and complete cost
where relevant, the evidence actually checked, and the remaining obligations.
Finite certificates or formal lemmas do not by themselves establish an
improved factoring machine. Run `make verify` and any relevant scoped checks;
run `make verify-formal` for the original Lean lemmas, `make verify-k6` for
k6, and `make verify-half-appendix` for the full `1/2` theorem. Report commands
that could not run.

The main README must stay current: when a change affects a result, evidence,
verification, or an open obligation, update its status, relevant links, and
reproduction guidance in the same change. Correct discovered stale claims
against current source and evidence. Do not promote an unreviewed or partial
result merely because it was merged.

Do not spawn subagents unless the user explicitly requests them.
