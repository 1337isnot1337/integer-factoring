# Degree-twelve factoring schedule

This directory contains a modified uniform factoring program and its complete
machine-time proof. The program uses

\[
K_0=48,\qquad G_d=5K_d^6,\qquad K_{d+1}=25K_d^{12}.
\]

The intended closed endpoint is `Factoring.K6.certified_improvement` in
[`FactoringK6/Certified.lean`](FactoringK6/Certified.lean). It states the
strict exponent comparison and the existence of one deterministic,
finite-alphabet factoring machine with

\[
T(n)\le C\exp\!\left(Cn/(\log n)^{1/30}\right)\quad(n\ge2).
\]

The input/output and time definitions remain the pinned reference's
`FactoringMachine`, `FactorsWithin`, binary `encodeNat`, and `Nat.minFac`.

## Verification state

The closed theorem and all eight modules passed an isolated local Lean 4.35.0
build. The three final axiom reports contain only `propext`,
`Classical.choice`, and `Quot.sound`. A fresh-source verification and the
hosted reproducibility run are in progress; the recorded final evidence will
replace this checkpoint status once those gates finish.

## Reproduction

From the repository root, with Lean 4.35.0 and Lake on `PATH`:

```sh
make verify-k6
```

The command archives the pinned public reference, gives the disposable source
tree its own Git root, applies `reference.patch`, rejects skipped patches,
checks reverse applicability, and verifies that the original challenge,
factoring definitions, toolchain, and dependency manifest are unchanged.
It builds the modified `Solution` and all `FactoringK6` modules, then checks
the endpoint source again and inspects its actual axiom reports.

The source snapshot is newly created on every run. Dependency source revisions
are checked against the pinned manifest before and after the build. The
verification does not reuse compiled files from another variant. Logs and
source hashes are written to `build/k6/`.

## Proof modules

| Module | Obligation |
| --- | --- |
| `ChainEnvelope` | Sharper volume contraction and explicit bounds at every recursion level |
| `BottomException` | Preserve the original base count 48 and pay for its one exceptional transition |
| `BudgetBounds` | Bound the top construction, bookkeeping, padding, and chain overhead |
| `RateArithmetic` | Integer logarithms, thirtieth roots, cutoff inequalities, and exact comparison certificates |
| `Quantitative` | Bound the actual uniform engine exponent, including all finite input exceptions |
| `MachineBridge` | Carry the cost through the counting program, compiler, and finite-alphabet simulator |
| `RateInterpretation` | Convert the exact integer bound to the real logarithmic-power statement |
| `Certified` | Discharge all schedule hypotheses and state the closed machine theorem |

The upstream source and its factoring/compiler infrastructure are due to
Isaac Rudich and Louis-Martin Rousseau. The source revision and compatible
Apache-2.0 attribution are retained in the file headers and lock file.
