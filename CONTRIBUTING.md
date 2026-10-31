# Contributions

Contributions include proof corrections, independent replays, small checkers,
documentation, source attribution, maintenance, and well-scoped negative or
partial research results. A stronger factoring exponent is not required.

## Default workflow: work directly on main

Clone the repository and work on `main`:

```sh
git clone https://github.com/1337isnot1337/integer-factoring.git
cd integer-factoring
git switch main
git pull --ff-only
```

Make focused changes, run the relevant checks, and commit and push useful
checkpoints frequently. Do not create a fork, branch, or pull request just to
hold work for review. Use a short-lived branch or worktree only when isolation,
parallel work, a risky experiment, or lack of push access makes it useful.
External contributors without push access can submit a PR.

For substantial research, read [autoresearch.md](autoresearch.md), the relevant
notes in [research](research/INDEX.md), and [AGENTS.md](AGENTS.md) when using an
agent. Check the [coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8)
for overlap, record the base commit and scope, and update or release the claim
when the work changes or stops. Coordination claims are advisory, not locks.

## Pull requests: merge every useful contribution

Pull requests are an integration mechanism, not a waiting room. Every useful
PR should be merged by default, including ongoing research, partial arguments,
negative results, finite experiments, documentation, tests, and maintenance.
Do not defer or reject a PR merely because it does not improve the factoring
exponent, is not fully complete, lacks Lean formalization, or has not received
independent mathematical review. Merge a useful checkpoint and continue
research on `main`.

Before merging, inspect the exact head and base, preserve source attribution
and contributor credit, resolve conflicts, and run relevant checks when
feasible. Fix regressions introduced by the PR. State failed and unrun checks
plainly. If the PR is stacked, merge its prerequisite first and then integrate
the incremental change. Do not leave useful PRs open indefinitely; do not
close useful work as superseded without incorporating its unique material or
a clear record of it.

A merge records integration only. It is not mathematical acceptance,
independent review, Lean verification, or evidence of a stronger factoring
bound. Keep claims, finite evidence, formal verification, and review status
separate in the README, research notes, and review ledger.

## Research and evidence

Use one focused change and identify its predecessor, source revision,
assumptions, proof or construction, and precise scope. Preserve attribution and
the upstream-compatible Apache-2.0 license; disclose substantial AI assistance.

For a bug or mathematical correction, identify the exact file, statement, and
commit, then give the smallest counterexample or reproduction available. For
an independent replay, report failures as well as successes. The
[issue templates](https://github.com/1337isnot1337/integer-factoring/issues/new/choose)
provide short outlines; leave unknown fields explicitly unknown.

For an asymptotic factoring claim, give a uniform algorithm and account for
the full bit-operation cost, including parameter choices, recursion, reduction,
precision, and finite exceptions as applicable. Explain which parts are proved,
implemented, formally checked, or still conditional. Finite tests and checked
lemmas support only their stated scope.

Run `make verify` and relevant focused checks. Use `make verify-formal` for
the original local lemmas, `make verify-k6` for the k6 proof, and
`make verify-half-appendix` for the full `1/2` theorem. Run
`make verify-docs` for documentation changes. List any unrun commands.
The default gate uses Python 3 and its standard library, plus Git and Make.
Lean and optional numerical packages are needed only for their respective
checks.

Keep the README current with changes to results, evidence, verification
commands, and open obligations. Retain individual credit when adapting,
combining, or superseding contributions. Update [CONTRIBUTORS.md](CONTRIBUTORS.md)
when importing another contributor's work; cite the PR and exact source.

## Review and credit

The [review guide](docs/REVIEWING.md) describes lightweight integration checks.
The [ledger](docs/REVIEW_LEDGER.md) records exact inspected heads, integration,
and review scope. Passing CI, merging, formal kernel checking, and independent
mathematical review have separate meanings.
