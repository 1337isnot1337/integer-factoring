# Contributing research effort

The goal is to improve deterministic integer factoring bounds with precise,
reproducible evidence. Start with the [README](README.md), the pinned upstream
revision in [upstream.lock.json](upstream.lock.json), and the
[research index](research/INDEX.md). The [half-rate theorem](formal-appendix/HalfTheorem.lean)
proves the exponent `1/2` in Lean, including the finite-alphabet machine.
The theorem composes pinned source results. Big O Speedrun independently rebuilt
commit `e91442d3` with Lean and nanoda; the [README](README.md) links the report
and verifier run. The k6 variant retains its separate Lean theorem at exponent
`1/30`; its sharper ordinary analysis reaches `log(12/11)/log12`. The older
combined rate in the [integration note](research/integration.md) remains a
separate candidate. Its generalized evaluator remains open and is not needed
by the `1/2` proof.

Before revisiting sparse operators, rectangular layouts, or small crossing
blocks, read [post-half-directions.md](research/post-half-directions.md).
It records proved restrictions and exact counterexamples found while trying
to improve on `1/2`, with reproducible finite rank checks. No stronger
exponent is established by those experiments.

## Workflow

1. Start from current `main`. Check the
   [discussion board](https://github.com/1337isnot1337/integer-factoring/discussions),
   [open issues](https://github.com/1337isnot1337/integer-factoring/issues),
   [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls),
   and relevant research notes for overlap. With GitHub CLI access,
   `python3 tools/find_overlap.py KEYWORD` searches the board, open work, and
   this checkout's research Markdown. Read the stated coverage limitations.
   The tool's board search still covers historical issue #8; check current
   [Discussions](https://github.com/1337isnot1337/integer-factoring/discussions)
   separately for active agent conversations and claims.
2. Work directly on `main` by default. Make focused commits and push useful
   checkpoints frequently. Use a short-lived branch or worktree only for
   isolation, parallel work, a risky experiment, or when direct push access is
   unavailable. Do not create PRs merely to park ongoing work.
3. Pick a precise question and record its assumptions, source revision, and
   success or stopping condition. For substantial work, post the scope, base
   commit, and next checkpoint on the discussion board; update or release
   the claim when the work changes or stops.
4. Explain the construction or obstruction and its scope. For a proposed bound,
   account for all levels, bit costs, parameter restrictions, and the link to
   a uniform factoring machine. Separate conditional or finite evidence from
   a complete asymptotic theorem.
5. Run `make verify` and relevant focused checks. Use `make verify-formal`
   for the original local lemmas, `make verify-k6` for the certified `1/30`
   machine, and `make verify-half-appendix` for the certified `1/2` machine.
   Each Lean target needs its pinned toolchain and upstream dependencies.
   Record actual results and commands that did not run.
6. Write a concise handoff with sources, contributor credit, limitations,
   and the next concrete question. A useful partial result or counterexample
   is welcome; do not describe it as a completed factoring improvement.
   Update the main README in the same change when status, evidence,
   verification, or open obligations change.

## Pull requests

Use a PR only when a contributor cannot push to the repository, a separate
integration review is useful, or branch isolation is necessary. Every useful
PR is expected to merge. Merge partial and ongoing research, finite results,
negative results, checkers, and documentation; do not wait for a stronger
exponent, Lean formalization, or independent mathematical review. Preserve
the exact scope and limitations in the research note and ledger. A merge means
the work is integrated, not that its mathematical claim has been accepted.

For an unassigned contribution, the [quickstart](docs/START_HERE.md) lists
smaller reproduction, correction, and documentation tasks alongside research.
A useful handoff should include an open next question, its relevance, a feasible
first step, and any blocker. Follow an explicit task instead of selecting a new
one.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the contribution policy and
[docs/REVIEWING.md](docs/REVIEWING.md) for integration and review scope.
