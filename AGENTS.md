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

## Research priority: cross the barrier, not polish a blocked family

The primary goal is a stronger **complete factoring bound**, not smaller
constants in an intermediate construction. Prioritize abstract, structural
approaches that can change the asymptotic result. Use computation to test
those ideas, rather than set the research agenda.

Before a substantial search or extension, write a short research question:

- What mathematical mechanism can improve the complete bound?
- Which known obstructions apply, if any? Which exact assumptions will change?
  If coverage is uncertain, identify the uncertainty instead of claiming an escape.
- What lemma, construction, or diagnostic experiment would resolve the next
  question? What result justifies further work? What result requires a different approach?

A speculative mechanism is allowed. It need not be proved before exploration.
But it must be more specific than changing a field, increasing a cutoff,
enlarging a search, or trying another parameter without a structural reason.
Read the relevant theorem's hypotheses. A method-specific barrier is not a
universal impossibility result.

**Stop optimizing a family once a proved barrier excludes the target under
its current assumptions.** Do not launch more local sweeps, larger replays,
or finer certificates merely to improve its last decimal places. A bounded
experiment is appropriate if it tests a specific way to escape the barrier,
checks a suspected error, or establishes a reusable structural lemma. Record
that purpose first. Finish necessary verification of an already claimed
result, preserve useful work, then redirect the main effort.

For the balanced two-cut route, the relevant crossing target is
`rho < 1/sqrt(2)`, with a valid uniform transfer and all costs included.
A decrease above that threshold is not a better factoring score. Other
routes can have different targets. State their actual transfer criterion.
The protected-core barriers already exclude boundary-only optimization of
the retained 96-letter family under their stated natural-scaling assumptions.
See [the protected-core results](research/jensen-tail-protected-core-floor.md)
and [the refinement's scope](research/f8-hyperoval-retraction-refinement.md).

Favor new representations, changes to the protected structure, collective
constructions, or reductions that preserve useful arithmetic correlations.
These are examples, not claims that those routes work. Do not replace one
parameter-search loop with an endless sequence of narrower negative results.
Use each obstruction to identify the next assumption to change.

At each research checkpoint, state whether the complete bound changed,
what structural fact was learned, and why the next step can still matter.
Do not use commit count, test volume, or tiny coefficient gains as substitutes
for progress toward the goal. This priority does not change the policy of
integrating useful partial results, corrections, or negative results.

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
