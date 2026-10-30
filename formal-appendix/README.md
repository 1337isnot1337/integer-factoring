# Lean proof of the half-rate factoring theorem

The `formal-appendix/` files prove the pinned reference's `Factoring.SqrtLogLogConjecture`. The main endpoint is [`HalfTheorem.lean`](HalfTheorem.lean), theorem `IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`.

## Proof map

| Stage | Modules | Role |
| --- | --- | --- |
| Factoring and prefix correctness | [`SelectedFactor.lean`](SelectedFactor.lean), [`SelectedPrefix.lean`](SelectedPrefix.lean), [`PureFactor.lean`](PureFactor.lean), [`ReferenceStatement.lean`](ReferenceStatement.lean) | Connect the program to the reference output specification. |
| Network and table construction | [`CliqueNetwork.lean`](CliqueNetwork.lean), [`CliqueRows.lean`](CliqueRows.lean), [`TopEvaluator.lean`](TopEvaluator.lean), [`TopGrouping.lean`](TopGrouping.lean), [`SelectedRow.lean`](SelectedRow.lean), [`SelectedTable.lean`](SelectedTable.lean) | Define the grouped clique input, row values, and written tables. |
| Search and encoding | [`GroupingPredicate.lean`](GroupingPredicate.lean), [`GroupingSearchPoly.lean`](GroupingSearchPoly.lean), [`SelectedGateCode.lean`](SelectedGateCode.lean), [`WordEncoding.lean`](WordEncoding.lean) | Prove executable selection and exact encodings. |
| Padding and rate | [`ExplicitPad.lean`](ExplicitPad.lean), [`MasterPad.lean`](MasterPad.lean), [`HalfSchedule.lean`](HalfSchedule.lean), [`DiscreteRate.lean`](DiscreteRate.lean) | Charge input and evaluator costs and prove the square-root-log rate. |
| Machine compilation | [`PaddedWordMachine.lean`](PaddedWordMachine.lean), [`PaddedMachineRate.lean`](PaddedMachineRate.lean), [`HalfTheorem.lean`](HalfTheorem.lean) | Build the finite-alphabet machine and prove the reference proposition. |

The remaining modules prove the individual arithmetic, gate, cluster, width, and polynomial-time facts used by these stages. The final result does not rely on the unrelated placeholder theorem in the unchanged reference challenge file.

## Reproduce

From the repository root:

```sh
make verify-half-appendix
```

The script [`tools/check_half_appendix.sh`](../tools/check_half_appendix.sh) fetches pinned versions of two source developments. It builds the Appendix B modules with Lean **4.33.1** and checks the bridge against the factoring source. This build has a different toolchain from the k6 proof (**4.35.0**).

The script audits the axiom reports of supporting declarations and the final theorem. The accepted axioms are `propext`, `Classical.choice`, and `Quot.sound`. The last successful PR verification reported **262 axiom reports**. This is a formal machine theorem, not a practical factoring benchmark.

See the [written composition](../research/half-rate-factoring.md) and the [source/encoding audit](../research/half-evaluator-audit.md) for the mathematical overview.
