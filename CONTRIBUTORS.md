# Sources and contribution credit

Credit follows the particular argument, code, certificate, or review used.
Useful independent reproductions, corrections, partial findings, and superseded
work retain credit even when they do not advance the selected exponent.

## Source foundations

- **Isaac Rudich and Louis-Martin Rousseau:** factoring reduction,
  multiplication networks, and source machine/counting infrastructure in
  [SubexponentialDeterministicFactoring](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring).
- **Claude / Anthropic:** the Appendix B clique evaluator in
  [formal-math](https://github.com/anthropics/formal-math), and the proof/code
  generation credited by the factoring reference.
- The exact upstream revisions are recorded in
  [upstream.lock.json](upstream.lock.json). Retain source-specific authorship
  and license notices when using these works.

## Repository work

- **Cinna Davis / 1337isnot1337:** repository research, integration,
  formalization, verification tools, and maintenance. These files used
  substantial OpenAI assistant support; per-result notes retain detailed scope
  and attribution.
- **Nil / BOS Reviewer, through nmamano:** the independent Big O Speedrun
  verifier rerun of the half-rate endpoint with Lean and nanoda, reported in
  [issue #28](https://github.com/1337isnot1337/integer-factoring/issues/28)
  and [run 41108289491](https://github.com/bigospeedrun/verifier/actions/runs/41108289491).

The [research index](research/INDEX.md), individual PRs, and
[review ledger](docs/REVIEW_LEDGER.md) retain the provenance and status of
specific contributions. Archived PRs remain citable at their exact heads;
archival does not assert that their claims were independently verified.

## Adding or correcting credit

When importing or adapting another person's work, add their preferred name or
handle, contribution scope, PR or commit, and dependencies here and in the
affected note. Keep credit for parallel or superseded work. An attribution
correction can be submitted directly as an issue or PR. Avoid including personal
contact details without the contributor's request. Contributions use the
[Apache-2.0 license](LICENSE), with upstream notices preserved.
