# Source audit of the sixth-power grouping variant

## Scope and result

This final audit uses the clean reference checkout at commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. The audited variant has
`G = 5 K^6`, `K_next = 25 K^12`, and the original bottom group count
`K0u = 48`.

The final [`reference.patch`](../variants/k6/reference.patch) has SHA-256:

```text
cc098155c5505cb3e04212c673ab1cc86e3a6639b6b870389c25e20f8ae22e6e
```

The source review found no change to the factoring problem, binary
encoding, machine model, finite-alphabet requirement, source operation
semantics, compiler, or simulator. The exact file comparisons and hashes
are recorded in [`evidence/k6_source_audit.json`](../evidence/k6_source_audit.json).

The final comparison covers all 146 original Lean files: 25 changed and
121 remain byte-for-byte identical. The patch adds
`FactoringAux/GroupingThreshold.lean`, for 26 patched files in total.
All 43 explicitly protected files remain identical to the reference.
The new quantitative proof modules are recorded separately in the audit
manifest, including their source hashes.

## Completed kernel verification

The patched source's `Solution` target completed successfully with 1,948
Lake jobs and exit status zero. The closed `FactoringK6.Certified` endpoint
also compiled with exit status zero against the final variant with initial
group count 48.
Its three endpoint declarations are:

- `Factoring.K6.factoring_root_rate`;
- `Factoring.K6.factoring_real_rate`;
- `Factoring.K6.certified_improvement`.

Each endpoint's axiom report contains exactly `propext`,
`Classical.choice`, and `Quot.sound`. The patched source's
`Factoring.subexponentialRate` has the same axiom set. The audit manifest
records the build completion marker, log hashes, endpoint reports, and
checked source hashes. These checks can be reproduced with
[`tools/check_k6.sh`](../tools/check_k6.sh).
The [local verification record](../evidence/k6-verification.json) separates
these local checks from the fresh `k6-rate` CI build. The reproduction script
verifies all patched and protected source hashes before any Lean command.

## Grouping correctness at both uses

`FactoringAux/EngineGroup.lean` changes `good_exists` to invoke
`Factoring.Improvements.network_grouping_five_sixth`. This theorem gives
the complete original grouping interface for `g = 5 K^6`:

- Each group has at most `2g/K` vertices, expressed without division.
- Each pair width satisfies `K^3 Y ≤ K U + g P`.
- Coarsening and writing a sufficiently wide clique preserve its value.

The recursive row and column networks reach this theorem through
`value_nextT`. Their row-degree bounds and total-volume bound remain the
original generic bounds. The top multiplication network reaches the same
theorem through `EngineTop.top_rep` and `UniTopSpec.top_value`, using
`clusterNet_fits` and `clusterNet_total_le`.

The new existence proof introduces no search oracle. The program still
enumerates the `K^g` coloring codes and selects `firstGood`. Its code-space
bound is computed with the new `g`, and the existing correctness proof
links the first successful code to the complete grouping guarantee.

## Runtime parameter construction

The mathematical schedule and its runtime expressions agree:

| Component | Mathematical value | Runtime expression |
| --- | --- | --- |
| Bottom count | `48` | The original literal `K0u` |
| Grid side | `5 K^6` | `UniParams.gE`, `UniSet.GEx` |
| Next group count | `25 K^12` | `UniParams.kbE` |
| Number of coloring codes | `K^(5 K^6)` | Level and top parameter blocks |
| Parameter budget | `(K+2)^(5 K^6+50)` | `UniSet.lvlV` and its bound proofs |

The level and top specifications, expression representations, scalar
parameter blocks, and value bounds were updated consistently. The syntax
of the uniform program remains fixed. Input length controls runtime loops
and stored parameters.

Some conservative analytical bounds still mention `3 K^8`. These are
supported by the proved inequality `Gs_le_eighth`, which shows
`5 K^6 ≤ 3 K^8` for the applicable group counts. They are upper bounds in
proofs, rather than obsolete runtime parameter constructors.

The successor-count budget changes from `9^40 (K^40)^16` to
`25^40 (K^40)^12`. The assembly accounts for this through the value bound
`2^(12 Z+186)` and the checked estimate `25^40 ≤ 2^186`.

The final source fix updates the monotonicity proof
`UniCostPath.lvlV_mono`: its power comparison uses exponent 6, matching
`lvlV(K) = (K+2)^(5 K^6+50)`. This changes the proof of the parameter
bound. The program and machine-semantics definitions are unchanged by
that fix.

## Counting and machine semantics

The crossing tensors, coefficient construction, grid contraction, choice
multiplicity, and modular zero-test machinery retain their reference
definitions. This certified variant uses the original tensor identity;
the separate rank-56 proposal remains a research candidate.

The exact definitions in `Challenge.lean` and
`FactoringAux/Definitions.lean` are unchanged. They specify `Nat.minFac`,
Mathlib's binary input and output encoding, a `TM2ComputableInTime`
machine, and finiteness of every stack alphabet. `Solution.lean` is also
unchanged.

The expression evaluator and compiler (`Ex.lean`, `ExCmd.lean`), copied
machine and source-language implementations, `SProg` compiler, simulator,
factoring driver, and `FactoringAux/UniFinal.lean` are byte-for-byte identical to
the reference. The revised bounds therefore retain their operation-cost
accounting, including array access, arithmetic, loops, table construction,
and simulation overhead.

The placeholder theorem in `Challenge.lean` is part of the unchanged
benchmark statement. No proof module imports it. The completed final
theorem reports contain only the three standard axioms listed above.

## Bottom transition and the quantitative rate

The variant retains the original 48-group brute-force program.
`FactoringK6.BottomException` proves that every positive group height is
already at least 600. All transitions above the bottom therefore satisfy
the sharper `23/25` contraction. The bottom transition uses the original
nonexpansion estimate and the existing additive padding.

The resulting volume budget uses `(23/25)^(ell-1)` when `ell ≥ 1`.
The formal weighted bound absorbs this one missing contraction into the
fixed coefficient `25 C`; with the denominator estimate's `C = 3`, this
becomes 75. The logarithmic saving exponent remains `1/30`. The finite
range with `ell = 0` is covered by a finite constant in the global bound.

The original `UniPathSpec.lean` was restored byte-for-byte before the final
build with initial group count 48. Experimental proof-resource changes for a 600-group bottom
are absent from the delivered source patch.

## Interpretation of the strict comparison

`RateInterpretation.baseline_log_rate_lt_one_thirtieth` proves

\[
\frac{\log(12/11)}{\log16}<\frac1{30}.
\]

The exponent `1/30` is a stronger proved upper-bound saving exponent than
the limiting exponent from the pinned reference's degree-sixteen schedule
and contraction analysis. A larger saving exponent decreases the
asymptotic upper-bound expression `exp(O(n/(log n)^alpha))`.

`RateInterpretation.exists_real_rate_bound` proves the conversion from
the exact integer time bound to the conventional real expression.
The closed theorem applies it to one finite-alphabet machine computing
`Nat.minFac` correctly on every natural-number input. It proves that some `C > 0`
satisfies

\[
M.\operatorname{time}(n)\le
C\exp\!\left(\frac{Cn}{(\log n)^{1/30}}\right)
\qquad\text{for every integer }n\ge2.
\]

The comparison concerns proved upper bounds. It supplies no lower bound
on the reference machine's actual running time and no measured practical
speedup for particular input sizes. Its scope is the pinned reference
analysis; universal novelty and superiority to every integer-factoring
algorithm are not asserted.
