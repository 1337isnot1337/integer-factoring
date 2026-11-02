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
  Sol worker 1's [near-height smooth-descent deduction](research/near-height-smooth-descent.md)
  reuses the credited norm-word decoder and polynomial form arithmetic;
  the [prescribed-target counterfamily](research/prime-power-descent-obstruction.md)
  combines classical Bonferroni/CRT estimates with actual prime-power forms.
  The coordinating worker contributed its two semiprime controls in
  [Discussion 51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359827);
  the checker retains both original receipt hashes and comment-level credit.
  The coordinated [event-boundary recurrence and Mellin contribution](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19356480)
  supplied the finite group-algebra identity, transform multiplier and original
  600-control panel integrated into the [event-word boundary note](research/event-word-boundary.md).
- **Sol worker 2 (OpenAI-assisted repository campaign):** the
  [scoped semiprime-field inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360111)
  and [independent all-component controls](evidence/semiprime-field-independent-review/provenance.json),
  retained byte-for-byte with their original source and record hashes. The
  inspection does not claim a full independent rederivation or a coefficient
  engine; its limits remain in the result note and review ledger.
- **Homelab worker 4 (DeepSeek campaign, under Cinna Davis):** the
  one-hot color-class cube construction, specialized profiles and
  vertex-letter grading corrections, retained at exact source
  [`1772970`](https://github.com/1337isnot1337/integer-factoring/commit/1772970a7dbaab83eeb0607014a1111e46718d86).
  The [integration note](research/color-class-cube-profile.md) preserves
  that source and distinguishes the later generic sector obstruction.
- **Nil / BOS Reviewer, through nmamano:** the independent Big O Speedrun
  verifier rerun of the half-rate endpoint with Lean and nanoda, reported in
  [issue #28](https://github.com/1337isnot1337/integer-factoring/issues/28)
  and [run 41108289491](https://github.com/bigospeedrun/verifier/actions/runs/41108289491).
- **orda-aska:** ordinary multiplier-network treewidth and relation-cycle
  analysis in [PR #73](https://github.com/1337isnot1337/integer-factoring/pull/73),
  including the balanced-separator lower bound and finite checker scripts.
  The contribution discloses AI research-agent assistance; its hand proof
  requires independent human review.
  The same contributor supplied the Fourier window-duality note and saved
  arithmetic experiments in [PR #75](https://github.com/1337isnot1337/integer-factoring/pull/75),
  with Grok Bot assistance disclosed in the source. The integration retains
  that credit and data, with explicit domain and quantifier repairs.
  The [graph/lattice follow-up](research/leading-power-window-graphs.md),
  source [`a427442`](https://github.com/1337isnot1337/integer-factoring/commit/a427442),
  supplies the weighted principal-word norm gap and injective-ball lemma
  used by the later event-separation and word-boundary deductions.

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
