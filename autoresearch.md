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
   It scans all accessible discussion titles/bodies, comments, and replies,
   including closed discussions, with pagination at each level. Historical
   issue #8 comments are listed separately. A failed remote request makes the
   search incomplete (nonzero exit), rather than reporting no overlap.
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

## Retained fixed-characteristic checkpoint from PR #22

The [fixed-characteristic transfer](research/fixed-characteristic-transfer.md)
and [weighted constructor](research/weighted-prefix-constructor.md) establish
a separate route through cyclic polynomial algebras over one fixed prime
field. They transfer a weighted evaluator at cost `F(n)poly(n,q)` to factoring
at cost `F(n)2^{O(n/log n)}`. A faster modular evaluator remains the missing
step; these notes do not improve the current exponent. The
[modular crossing construction](research/past-half-pieces.md) now closes
the same half-rate bound over characteristic two, with scale one, an
invariant quotient, and a direct ring-operation count. Its finite 17-piece
recombination improves rank product `36` to `27` on three labels, without
an asymptotic improvement. Inspect its exact hypotheses before extrapolating.

The [folded Pfaffian formula](research/modular-pfaffian.md) evaluates all
final characteristic-two pieces in polynomial time in their label weight.
For uniformly computable labels at any polynomial weight in `lambda`,
actual rank factorizations now satisfy the required `2^{O(lambda)}`
construction cost. A proved uniform rank-width bound below the half-rate
coefficient remains open. Nine exact high-weight label tests found no
improvement. The [sparse-grid analysis](research/past-half-sparse-grids.md)
also rules out improving the declared central-channel coefficient by
arbitrarily changing regular grid densities; it makes no general rank
or factoring lower-bound claim.

For labels with a shared complete binary block, the
[exact Schur reduction](research/fermi-sea-schur.md) computes the effect
of that background: a common scalar convolution of the original piece
profile. Separated coefficient supports preserve every individual rank.
Read its fixed-cardinality and active-mode restrictions before applying
it to other families. The modular note also explains why type aggregation
of the finite 17-piece certificate retains a pure type too large to beat
the half-rate coefficient.
