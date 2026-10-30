# Contributions and review

Start with [AGENTS.md](AGENTS.md), [autoresearch.md](autoresearch.md), and the
relevant notes in [research](research/). Use one focused change and identify
its predecessor, source revision, assumptions, proof or construction, and
precise scope. Preserve attribution and the upstream-compatible Apache-2.0
license; disclose substantial AI assistance.

Before substantial research, check and post your intended scope on the
[coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8).
Link your branch or worktree and base commit, report meaningful changes and
handoffs, and release your claim when finished. Coordinate overlapping lanes.

For an asymptotic factoring claim, give a uniform algorithm and account for
the full bit-operation cost, including parameter choices, recursion, reduction,
precision, and finite exceptions as applicable. Explain which parts are proved,
implemented, formally checked, or still conditional. Finite tests and checked
lemmas support only their stated scope. A merged pull request does not itself
constitute independent proof review or establish a new bound.

Run `make verify` and relevant focused checks. Use `make verify-formal` for
the original local lemmas and reference envelope. Use `make verify-k6` for
the k6 source patch or its complete rate proof. List any unrun checks.
Include reproduction commands, actual outcomes, remaining obligations, and
contributor credit in the pull request.

Keep the README current with changes to results, evidence, verification commands,
and open obligations. Update the relevant links and status in the same pull
request; distinguish proposed improvements from completed and reviewed results.
