# Agent entry point

Read [autoresearch.md](autoresearch.md), [CONTRIBUTING.md](CONTRIBUTING.md),
the [README](README.md), and the relevant research notes before substantial
work. Follow the user's assigned task and preserve existing work in progress.

Before starting a new direction, check current `main`, open Issues and pull
requests, and the existing research for overlap. Use a focused branch or
worktree when sharing a checkout. Record the source revision and attribution.

State the exact claim and its assumptions, the algorithm and complete cost
where relevant, the evidence actually checked, and the remaining obligations.
Finite certificates or formal lemmas do not by themselves establish an
improved factoring machine. Run `make verify` and any relevant scoped checks;
run `make verify-formal` for changes to formal proofs when the toolchain is
available. Report checks that were not run.

When a change affects a result, evidence, verification, or an open obligation,
update the README's current status, relevant links, and reproduction guidance
in the same change. Check claims against the actual source and evidence; do not
promote an unreviewed or partial result merely because it was merged.

Do not spawn subagents unless the user explicitly requests them.
