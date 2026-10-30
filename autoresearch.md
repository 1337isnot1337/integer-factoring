# Contributing research effort

The goal is to improve deterministic integer factoring bounds with precise,
reproducible evidence. Start with the [README](README.md), the pinned upstream
revision in [upstream.lock.json](upstream.lock.json), and the
[integration note](research/integration.md). The proposed combined rate is a
candidate: the modified factoring machine and its complete machine-level proof
remain open.

1. Check current `main`, the [coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8),
   [open issues](https://github.com/1337isnot1337/integer-factoring/issues),
   [pull requests](https://github.com/1337isnot1337/integer-factoring/pulls),
   and relevant research notes for overlap.
2. Pick a precise question and record its assumptions, source revision, and
   success or stopping condition. Use a focused branch or worktree. Post the
   scope, base commit, branch or worktree, and next checkpoint on the board;
   update or release the claim when the work changes or stops.
3. Explain the construction or obstruction and its scope. For a proposed bound,
   account for all levels, bit costs, parameter restrictions, and the link to
   a uniform factoring machine. Separate conditional or finite evidence from
   a complete asymptotic theorem.
4. Run `make verify` and the relevant focused checks. Formal checking is a
   separate `make verify-formal` target; it requires Lean and downloads the
   pinned upstream reference and its dependencies. Record actual results and
   any checks that were not run.
5. Write a concise handoff with sources, contributor credit, limitations,
   and the next concrete question. A useful partial result or counterexample
   is welcome; do not describe it as a completed factoring improvement. Update
   the main README in the same change when status, evidence, verification, or
   open obligations change.

See [CONTRIBUTING.md](CONTRIBUTING.md) for review expectations.
