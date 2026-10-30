# Contributing research effort

The goal is to improve deterministic integer factoring bounds with precise,
reproducible evidence. Start with the [README](README.md), the pinned upstream
revision in [upstream.lock.json](upstream.lock.json), and the
[research index](research/INDEX.md). The [half-rate theorem](formal-appendix/HalfTheorem.lean)
proves the exponent `1/2` in Lean, including the finite-alphabet machine.
The theorem composes pinned source results and is not externally reviewed. The k6 variant
retains its separate Lean theorem at exponent `1/30`; its sharper ordinary
analysis reaches `log(12/11)/log12`. The older combined rate in the
[integration note](research/integration.md) remains a separate candidate.
Its generalized evaluator remains open and is not needed by the `1/2` proof.

Before revisiting sparse operators, rectangular layouts, or small crossing
blocks, read [post-half-directions.md](research/post-half-directions.md).
It records the proved restrictions and exact counterexamples found while
trying to improve on `1/2`, with reproducible finite rank checks. No stronger
exponent is established by those experiments.

1. Check current `main`, the [coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8),
   [open issues](https://github.com/1337isnot1337/integer-factoring/issues),
   [pull requests](https://github.com/1337isnot1337/integer-factoring/pulls),
   and relevant research notes for overlap. With GitHub CLI access,
   `python3 tools/find_overlap.py KEYWORD` searches the board, open work, and
   this checkout's research Markdown. Read the stated coverage limitations.
2. Pick a precise question and record its assumptions, source revision, and
   success or stopping condition. Use a focused branch or worktree. Post the
   scope, base commit, branch or worktree, and next checkpoint on the board;
   update or release the claim when the work changes or stops.
3. Explain the construction or obstruction and its scope. For a proposed bound,
   account for all levels, bit costs, parameter restrictions, and the link to
   a uniform factoring machine. Separate conditional or finite evidence from
   a complete asymptotic theorem.
4. Run `make verify` and the relevant focused checks. Use `make verify-formal`
   for the original local lemmas, `make verify-k6` for the certified `1/30`
   machine, and `make verify-half-appendix` for the certified `1/2` machine.
   Each Lean target needs its pinned toolchain and upstream dependencies.
   Record actual results and commands that did not run.
5. Write a concise handoff with sources, contributor credit, limitations,
   and the next concrete question. A useful partial result or counterexample
   is welcome; do not describe it as a completed factoring improvement. Update
   the main README in the same change when status, evidence, verification, or
   open obligations change.

See [CONTRIBUTING.md](CONTRIBUTING.md) for review expectations.
