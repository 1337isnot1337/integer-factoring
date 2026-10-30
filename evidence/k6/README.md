# Completed k6 verification

The modified factoring program and the closed exponent-`1/30` theorem passed
a fresh-source Lean build on 2026-10-30. Lake completed **2,252 jobs**, including
the modified reference `Solution` and all eight `FactoringK6` modules. A separate
GitHub Actions build of the same proof sources also passed all 2,252 jobs.

## Inputs and result

| Item | Recorded value |
| --- | --- |
| Reference commit | `c91cc3f4947caca8977e81a5d14244f73b5ae2b1` |
| Lean | `4.35.0`, commit `eb6cafb7290e8bd5f2937522d8d804f58f02bafc` |
| Mathlib commit | `56a005855ea17e0b0958f73ec801ae140b0a1679` |
| Patched reference paths | 26 |
| New proof modules | 8 |
| Protected inputs | `Challenge.lean`, `FactoringAux/Definitions.lean`, `lean-toolchain`, `lake-manifest.json` |
| Final declarations | `Factoring.K6.factoring_root_rate`, `Factoring.K6.factoring_real_rate`, `Factoring.K6.certified_improvement` |
| Axioms of each final declaration | `propext`, `Classical.choice`, `Quot.sound` |
| Local Lake result | `Build completed successfully (2252 jobs).` |
| Fresh endpoint command | `lake env lean FactoringK6/Certified.lean`, exit 0 |
| Final source/log checker | `tools/check_k6_receipt.py`, exit 0 |

The [success receipt](verification.json) binds all captured source files,
the patch, the six original logs, the axiom receipt, and the final checker by
SHA-256. The [source manifest](sources.json) records the copied proof files,
patched reference files, and protected reference inputs. The
[separate source and semantic audit](../../research/k6-proof-audit.md) also
checked the other 129 reference source files and all nine pinned dependencies.

## Retained evidence

- [Fresh endpoint output](Certified.log) and [parsed axiom receipt](axioms.json).
- [Complete local Lake build log, gzip compressed](build.log.gz).
- [Patch applicability log](patch-check.log) and [patch application log](patch-apply.log).
- [Toolchain output](toolchain.log) and [dependency cache output](cache.log).
- [Hosted run metadata](hosted-checkpoint.json) and
  [its build/endpoint excerpt](hosted-checkpoint-excerpt.log).

`build.log.gz` decompresses to the exact `build.log` whose hash appears in
`verification.json`; no warning or failure lines were removed. The source
build contains upstream deprecation and lint warnings, but no failed target,
`sorry`, or unexpected final axiom.

## Reproduce

With Lean 4.35.0 and Lake on `PATH`, run from the repository root:

```sh
make verify-k6
```

The driver creates a fresh source tree from the pinned reference, applies the
actual program patch, checks protected definitions and dependency revisions,
builds `Solution` and `FactoringK6`, and runs the endpoint source a second time.
It then rehashes the captured inputs and requires all three endpoint reports
in both Lean passes. Failed runs invalidate previous success receipts.

Mathlib's pinned dependency cache may be reused. Compiled files from an earlier
factoring variant are not reused in the fresh source tree.

The recorded local and hosted builds used the
[checkpoint driver](https://github.com/1337isnot1337/integer-factoring/blob/71c9dc0946e11f117484725eb2e6a154d29d5d45/tools/check_k6.sh).
After the local build finished, the final receipt checker independently
validated its sources and logs and emitted `verification.json`. The current
driver includes that same checker automatically. Its ten regression tests
cover stale sources and receipts, incomplete proofs, unexpected axioms, and
partial receipt-write failure.

`make verify` also passed on the final source tree: 40 discovered Python tests,
five crossing tests, the exact coefficient-basis certificates, and the schedule
parameter checks. Those finite checks support their stated local scope; the
Lean endpoint establishes the complete asymptotic machine bound.

The completed hosted checkpoint is
[run 41000040607](https://github.com/1337isnot1337/integer-factoring/actions/runs/41000040607)
at commit `71c9dc0946e11f117484725eb2e6a154d29d5d45`.
Current verification runs are available in the
[k6 workflow](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify-k6.yml).

The [mathematical note](../../research/certified-k6.md) states the exact theorem
and score implication. These checks establish the claimed asymptotic bound
relative to the pinned proof stack. They do not establish practical performance
or a leaderboard submission/acceptance.
