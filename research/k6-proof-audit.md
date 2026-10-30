# Independent audit of the k6 factoring variant

This is a review by a separate agent within the same working session, not
external peer review. The reviewing agent did not author the k6 Lean modules
or the reference patch. Its work consisted of source-level mathematical and
semantic review, checking the scoring interpretation, and independently
comparing source and log hashes. The implementation agent ran the builds;
this review inspected their recorded outputs. The completion addendum below
records the successful fresh-source verification.

## Scope and result

I reviewed all eight Lean files under `variants/k6/FactoringK6/` and every hunk of
`variants/k6/reference.patch`. I also checked the reference problem definitions
against Mathlib's machine and time definitions. I found no mathematical gap,
undischarged hypothesis, weakened output requirement, or change in the machine
model in this proof package.

The reference revision is
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. The audited patch SHA-256 is
`cc098155c5505cb3e04212c673ab1cc86e3a6639b6b870389c25e20f8ae22e6e`.
The initial local variant build completed 2,238 jobs successfully. Its three
endpoint theorems report only `propext`, `Classical.choice`, and `Quot.sound`.
This report distinguishes that successful build from the separate clean
reproduction gate described below.

## Exact endpoint

`Factoring.K6.factoring_root_rate` proves that there exist one
`FactoringMachine` and one natural constant C such that, for every input length
n,

    M.time n <= C * 2^(C * (n / (nthRoot 30 (log_2 n) + 1) + 1)).

Here the logarithm, root, and division are natural-number operations.
`FactorsWithin` includes finiteness of every stack alphabet.

`Factoring.K6.factoring_real_rate` proves that there exist one such machine M
and a positive real C such that, for every natural n >= 2,

    M.time n <= C * exp(C*n / (ln n)^(1/30)).

The type `FactoringMachine` is unchanged:

    TM2ComputableInTime encodeNat encodeNat Nat.minFac.

Mathlib's `TM2ComputableInTime.outputsFun` counts time at the length of the
encoded input. Consequently n is the input bit length, not the integer being
factored. Correctness is already part of the machine's type. It applies to
every encoded natural input, rather than a selected family of composites.

The machine is uniform: the existential quantifier for M precedes all input
lengths, and the source program computes its recursion depth from the input.
There are finitely many stacks, states, labels, and symbols. The compiler and
simulation bridge remain the original reference machinery.

## What changed in the algorithm

The original schedule uses G = 3 K^8 and K_next = G^2 = 9 K^16. The variant
uses G = 5 K^6 and K_next = 25 K^12. This is an actual change to the uniform
program: `kbE`, `gE`, `GEx`, `lvP`, and the top and level search/table
parameters all use the new constants. It is not a reinterpretation of the old
program's cost.

The new `GroupingThreshold` module derives a sufficient condition

    g > K^2 + 4 D^2 K^4

for grouping at accuracy D, from the reference's second-moment estimates.
Setting D = K yields the threshold K^2 + 4 K^6. For K >= 2,

    K^2 + 4 K^6 < 5 K^6.

The argument unions the bad size events and the bad pair-weight events; the
total number of bad colourings is strictly less than K^g. It retains the
original network value identity. The executable program still enumerates
codes deterministically and uses `firstGood`. The proof of existence does not
introduce random bits or an oracle into the program.

The patch touches 26 reference paths. It does not edit `Challenge.lean`,
`FactoringAux/Definitions.lean`, the encodings, `Nat.minFac`, the compiler, or
the TM simulation. All changed program expressions have corresponding updated
value and cost proofs.

## Quantitative proof chain

| Component | Audited responsibility |
| --- | --- |
| `GroupingThreshold` in the patch | The smaller grouping size still admits a good grouping and preserves value. |
| `ChainEnvelope.lean` | The volume contracts by 23/25 for group counts at least 600; recursive costs fit a single explicit envelope. |
| `BottomException.lean` | The actual initial count 48 is handled, losing only one contraction and a fixed factor 25. |
| `RateArithmetic.lean` | The degree-12 schedule and the 23/25 contraction imply a thirtieth-root logarithmic saving, with floors. |
| `BudgetBounds.lean` | Top construction, path counts, modulus/value sizes, and polynomial terms obey the same weighted budget. |
| `Quantitative.lean` | Combines the budgets into a global bound for `Zu`, including the finite range where depth is zero. |
| `MachineBridge.lean` | Reconstructs the actual source-program specification and invokes the unchanged finite-alphabet TM simulation. |
| `RateInterpretation.lean` | Converts the exact natural-number rate to the stated real-log rate. |
| `Certified.lean` | Discharges the schedule and initial-value hypotheses and gives the closed machine theorem. |

### Contraction and the bottom exception

The reference transition bound is

    12 K N_next <= (11 K + 24) N + 24 K G^4.

For K >= 600, the coefficient satisfies

    11/12 + 2/K <= 23/25.

This gives the new envelope's contraction. The actual machine starts at K = 48,
so applying that contraction indiscriminately would be invalid. The package
does not do that: every positive schedule height has K >= 600, while the sole
bottom transition uses the weaker nonexpansion bound. The resulting main
budget is proportional to

    (23/25)^(max(ell-1, 0)) * N_top.

The lost contraction contributes only a fixed factor. The chain proofs retain
the additive G^4 padding, level charges, and logarithmic bookkeeping.

### Depth and the exponent 1/30

The schedule obeys

    2 K_d <= (2 K_0)^(12^d).

The stopping rule K_ell^40 <= n < K_(ell+1)^40 therefore bounds log_2 n by a
constant times 12^ell. The exact integer certificate

    12 * 23^30 < 25^30

ensures that the thirtieth root grows more slowly than (25/23)^ell. The checked
weighted inequality is

    (nthRoot 30 (log_2 n) + 1) * 23^ell <= 3 * 25^ell.

It combines with the contraction and all overhead bounds to yield

    (nthRoot 30 (log_2 n) + 1) * Zu(n) <= constant * n

in the positive-depth range. A finite sum bounds all remaining inputs. These
constants are independent of n, even though they are very large.

### Actual machine time

The bridge uses the source program `drvC (uniC UP)`. The new value bound is

    2^(12*Zu(n) + 186),

using the exact estimate 25^40 <= 2^186. It proves the source specification,
cost and value requirements, then calls the original TM compiler/simulator.
That simulation introduces only a fixed polynomial in the proved source
budget, absorbed by the constant in the exponential. No arithmetic operation
on unbounded integers is silently charged unit cost in the final theorem.

### Natural and real logarithms

For n >= 2, put D = nthRoot 30 (log_2 n) + 1. The proof establishes both

    (ln n)^(1/30) <= D <= n

and

    floor(n/D) + 1 <= 2*n/(ln n)^(1/30).

It then converts powers of 2 to the real exponential. This handles the
quotient rounding and the added 1 rather than discarding them informally.

## Score and strict comparison

The official BigO source for this task assigns score alpha to a bound of the
form

    2^(b/(log b)^(alpha-o(1))),

meaning the bound with alpha-epsilon eventually holds for every fixed positive
epsilon. Its reference entry has score 0 for the stated subexponential result.

The checked C*exp(C*b/(ln b)^(1/30)) theorem meets the score-1/30 definition:
after taking logarithms, dividing by b/(ln b)^(1/30-epsilon) makes the main
term C/(ln b)^epsilon tend to zero. The logarithm of the prefactor also tends
to zero under that division. Fixed polynomial factors and logarithm-base
conventions are absorbed in the same way.

The stronger exponent extracted analytically from the unchanged reference
recurrence is

    alpha_ref = ln(12/11)/ln(16) = approximately 0.031382720520964705.

The endpoint `certified_improvement` additionally proves

    alpha_ref < 1/30.

It derives this from the exact certificate 12^30 < 16*11^30. The new exponent
is approximately 0.033333333333333333. The comparison theorem itself proves
the numerical inequality and the new machine bound; it does not supply a
separate formal runtime theorem for the unchanged reference at alpha_ref.

Primary scoring source:
https://github.com/bigospeedrun/site/blob/5f428921611f2f0b740e4553a1c2cf92cd6165b2/data/integer-factoring.json

## Precise limitations

1. This is a strict improvement in the proved asymptotic upper bound and the
   stated scoring rule. It is not a lower bound on the reference program's
   actual time, nor proof that the new program runs faster on every input.
2. No practical runtime advantage, usable crossover point, or implementation
   benchmark is established. The constants and the small-input cutoff are
   enormous.
3. The sqrt-log conjecture would require exponent 1/2. Exponent 1/30 does not
   prove it, and does not prove polynomial-time factoring.
4. The official leaderboard's acceptance or review status does not follow
   automatically from a local theorem and must not be claimed without a
   separate submission and response.
5. The proof is checked relative to Lean 4.35.0, the pinned Mathlib dependency,
   and the standard listed axioms. This review audits the changed mathematics
   and the connection to the reference's exact machine definitions; it does
   not independently rederive the complete original compiler and reference
   proof in prose.

## Clean-reproduction issue found during this audit

The first run of the new clean verification script unpacked the reference
inside the enclosing project's Git worktree. In that location,
`git apply --no-index --check` exited successfully but skipped every one of the
26 patches. A direct hash comparison with the successfully compiled variant
detected this. The fresh tree still had G = 3 K^8 and K_next = 9 K^16.

This was reported to the main agent immediately. It affects that attempted
reproduction, not the initial successfully compiled k6 tree. The clean
reproduction must establish its own Git root, or use a directory outside the
enclosing worktree, and verify that the patch was actually applied before
compilation. Final promotion requires matching all patched source files and
all eight endpoint modules to the audited package, a successful clean build,
and the three endpoint axiom reports.

The script was subsequently corrected to initialize a separate Git root,
reject skipped patches, and check reverse application. I independently
compared the replacement clean tree `build.UAZKv1E8` with the successfully
compiled initial variant and the new `sources.json` manifest. All 26 patched
reference files, all four protected reference files, and all eight
`FactoringK6` files match byte for byte. The patch hash also matches. The clean
tree now contains exactly the audited and previously compiled source; its
fresh build was still running at the time of this report update.

## Completion addendum — 2026-10-30

The corrected fresh-source build `build.UAZKv1E8` has now completed
successfully. I inspected `build/k6/build.log`, `Certified.log`, `sources.json`,
and `verification.json` after completion. The build log records all eight
`FactoringK6` modules and the modified reference `Solution` as built, and ends
with `Build completed successfully (2252 jobs).` The recorded toolchain is
Lean 4.35.0.

The endpoint source was then checked again. The
[endpoint axiom log](../evidence/k6/Certified.log) reports exactly
`propext`, `Classical.choice`, and `Quot.sound` for each of:

- `Factoring.K6.factoring_root_rate`
- `Factoring.K6.factoring_real_rate`
- `Factoring.K6.certified_improvement`

The [verification receipt](../evidence/k6/verification.json) records
`status: success` at `2026-10-30T04:04:21.868615+00:00`. Its reference revision
and patch hash are the ones audited above. I independently rehashed the
receipt's 38 source files in the build tree: 26 patched reference files, four
protected reference inputs, and eight new proof modules. I also rehashed the eight
proof files in the repository, all six recorded logs, the source manifest,
the patch, and the receipt checker. Every value matched the receipt.

In particular, the source manifest SHA-256 is
`8523161194d43ff499dee16d2b9e4325630fecfd88d9256f704d6b301949e15a`.
The endpoint log SHA-256 is
`ee7fb20cd65890b3def07e8291beb2c53feae8c89b7b8becc457c5c627eec91f`.
These checks bind the successful build and endpoint reports to the same
source package reviewed in this audit.

The local clean-reproduction gate is therefore complete. The earlier
skipped-patch attempt is retained above as an audit finding, not as evidence
for the theorem. The result supported by the completed gate is one uniform
finite-alphabet factoring machine with the stated exponent-1/30 upper bound
and the exact strict numerical comparison to `alpha_ref`. The score
implication and limitations stated above remain unchanged. This addendum
records local verification; it does not assert an external review or a
leaderboard acceptance.
