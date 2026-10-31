# Start contributing

You can contribute a small correction, reproduce a check, review one argument,
or develop a construction. You do not need a new factoring bound to help.

## Run the basic checks

Install Python 3, Git, and Make, then run:

```sh
git clone https://github.com/1337isnot1337/integer-factoring.git
cd integer-factoring
git rev-parse HEAD
python3 --version
make verify
```

This runs the default Python tests, exact finite certificates, and local
Markdown-link check using the standard library. It does not rebuild the Lean
proofs. No GitHub account is needed to clone or check the public repository.
For documentation changes, `make verify-docs` is the focused check. Optional
rank experiments can use NumPy; the default checks do not require it.

To send changes, fork the repository on GitHub, clone your fork, and create a
branch with `git switch -c contribution/my-change`. See
[CONTRIBUTING.md](../CONTRIBUTING.md) for the PR workflow and checks.

## Choose a starting point

| Interest | First contribution | Read first |
| --- | --- | --- |
| Reproducibility | Run one checker at a recorded commit; report environment, output, and any mismatch | [Verification commands](../README.md#verification), [tools](../tools/) |
| Documentation | Fix a confusing definition, missing prerequisite, or stale link | [README](../README.md), [research index](../research/INDEX.md) |
| Mathematics | Audit one lemma, supply a counterexample, or close a stated obligation | The relevant note in the [research index](../research/INDEX.md) |
| Lean | Rebuild a pinned proof or review one dependency or axiom report | [Half-rate proof map](../formal-appendix/README.md) |
| Research | Develop an effective family or prove an obstruction in a clearly stated model | [Research workflow](../autoresearch.md), [directions beyond half](../research/post-half-directions.md) |

Browse [starter issues](https://github.com/1337isnot1337/integer-factoring/issues?q=is%3Aissue%20is%3Aopen%20label%3A%22good%20first%20issue%22),
[help wanted](https://github.com/1337isnot1337/integer-factoring/issues?q=is%3Aissue%20is%3Aopen%20label%3A%22help%20wanted%22),
and [open PRs](https://github.com/1337isnot1337/integer-factoring/pulls).
You can also [report a finding](https://github.com/1337isnot1337/integer-factoring/issues/new/choose).
For substantial research, check and post a narrow scope on the
[coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8)
before starting. Small corrections and reproduction reports can go directly
to an issue or PR.

For ongoing multi-agent research direction, the
[Astra supervisor prompt](RESEARCH_SUPERVISOR_PROMPT.md) provides a
copy-ready brief for Dot, including task selection and mathematical guidance.

## Understand the evidence

The [README result table](../README.md#results-and-evidence) describes the
results on `main`. The complete half-rate Lean theorem is
[`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](../formal-appendix/HalfTheorem.lean).
The earlier k6 proof is a separate endpoint. No exponent above `1/2` is
established here.

A finite rank computation checks its declared finite cases. An obstruction
applies to its stated construction class. A proposed asymptotic improvement
also needs a uniform algorithm and complete bit-cost accounting. Open PRs
can contain later, unreviewed work; a merged PR is not mathematical acceptance.
The [review guide](REVIEWING.md) explains these distinctions.

## Writing a research finding

Use a focused note under `research/`; follow existing naming and checker
conventions. A useful outline is:

1. **Question and scope:** exact claim, assumptions, and what is excluded.
2. **Predecessors:** pinned source commits, related notes or PRs, and credit.
3. **Argument:** proof, construction, or counterexample; complete costs for a bound.
4. **Replay:** commands from the repository root, dependencies, input sizes,
   actual outputs, and whether each check is finite or uniform.
5. **Open next questions:** one concrete next step, why it matters, and blockers.

Keep checkers under `tools/`, meaningful tests under `tests/`, and compact
certificates under `evidence/` when appropriate. State randomness, search bounds,
and stopping criteria for experiments. An unsuccessful bounded search is a
partial result; explain exactly what it covered. Preserve frozen evidence and
update the research index and README when their status or commands change.
