# Contributions and review

New here? Use the [quickstart](docs/START_HERE.md). Contributions include proof
corrections, independent replays, small checkers, documentation, source attribution,
and well-scoped negative results. A stronger factoring exponent is not required.
You may open an issue or share a partial finding before preparing a PR.

For research, read [autoresearch.md](autoresearch.md), the relevant notes in
[research](research/INDEX.md), and [AGENTS.md](AGENTS.md) when using an agent.
Use one focused change and identify its predecessor, source revision,
assumptions, proof or construction, and precise scope. Preserve attribution and
the upstream-compatible Apache-2.0 license; disclose substantial AI assistance.

Before substantial research, check and post your intended scope on the
[coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8).
Link your branch or worktree and base commit, report meaningful changes and
handoffs, and release your claim when finished. Coordinate overlapping lanes.
Small corrections, documentation fixes, and reproduction reports do not need a
prior research claim. Existing claims are advisory; check their latest update.

## Sending a contribution

1. Fork the repository, clone your fork, and create a focused branch, for example
   `git switch -c contribution/my-change`. Record `git rev-parse HEAD` as your base.
2. Make the smallest useful change. Keep pinned upstream sources and existing
   evidence intact; record fresh runs separately. For a research note, use the
   optional [finding outline](docs/START_HERE.md#writing-a-research-finding).
3. Run the checks below that cover the change. Include the actual commands,
   outcomes, environment, and unrun checks in the PR template.
4. Open a PR against `main`; use a draft for ongoing work. Link predecessor work
   and any coordination claim. A reviewer may request a smaller independently
   reviewable checkpoint from a large research branch.

For a bug or mathematical correction, identify the exact file, statement, and
commit, then give the smallest counterexample or reproduction you have. For an
independent replay, report failures as well as successes. The
[issue templates](https://github.com/1337isnot1337/integer-factoring/issues/new/choose)
provide short outlines; leave unknown fields explicitly unknown.

## Claims and checks

For an asymptotic factoring claim, give a uniform algorithm and account for
the full bit-operation cost, including parameter choices, recursion, reduction,
precision, and finite exceptions as applicable. Explain which parts are proved,
implemented, formally checked, or still conditional. Finite tests and checked
lemmas support only their stated scope. A merged pull request does not itself
constitute independent proof review or establish a new bound.

Run `make verify` and relevant focused checks. Use `make verify-formal` for
the original local lemmas. Use `make verify-k6` for the k6 proof and
`make verify-half-appendix` for the full `1/2` theorem. Run
`make verify-docs` for documentation changes. List any unrun commands.
Include reproduction commands, actual outcomes, remaining obligations, and
contributor credit in the pull request.

The default gate uses Python 3 and its standard library, plus Git and Make.
Lean and optional numerical packages are needed only for their respective
checks. A documentation contribution does not require installing Lean.

Keep the README current with changes to results, evidence, verification commands,
and open obligations. Update the relevant links and status in the same pull
request; distinguish proposed improvements from completed and reviewed results.

## Review and credit

The [review guide](docs/REVIEWING.md) describes evidence and integration checks.
The [ledger](docs/REVIEW_LEDGER.md) records exact inspected heads and dispositions;
GitHub is the live queue. Passing CI, merging, formal kernel checking, and
independent mathematical review have separate meanings.

Retain individual credit when adapting, combining, or superseding contributions.
Update [CONTRIBUTORS.md](CONTRIBUTORS.md) when importing another contributor's
work; cite the PR and exact source. A useful result can be credited without
advancing the current exponent or being imported into `main`.
