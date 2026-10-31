# Global code-coset filters with unrestricted multiplier degree

## Result and scope

The exact reduction in
[`boolean-multiplier-automata.md`](boolean-multiplier-automata.md) permits a
growing code-coset support to use an arbitrary polynomial multiplier, over any
characteristic-two coefficient extension and with unrestricted ramification.
This note records completed exploratory searches using that complete
multiplier envelope.

All 15 recorded scale-and-charge searches chose the baseline. No strict
crossing was found. The result is not a proof of a pressure floor over all
codes, all assignments of physical modes, or all scale-and-charge values.

These filters are global. A single syndrome equation can involve 64, 128, or
256 physical modes. They are not limited to repeating a small independent
tile, and the calculation retains every actual physical mode energy.

## 1. Explicit support families

The source manifest contains 15 full-rank binary parity-check matrices with
eight columns, at ranks 2 through 6. It includes systematic matrices, an
extended Hamming column set, and two graph incidence column sets. The exact
columns are part of the evidence, so the named constructions are not needed
to identify the family.

At physical length \(n\), each source column receives a group of physical
modes. The active modes are the largest-energy 8, 16, \(n/2\), or \(n\) modes.
The assignments are contiguous, interleaved, reflected hot/cold, and seeded
shuffled. Equal nonzero column types are merged only after their actual mode
energy lists have been assigned. All remaining modes are explicit zero-column
spectators.

Every syndrome is included. The support counts at each physical length are:

| Rank | Source matrices | Mode assignments | Supports including every syndrome |
|---|---:|---:|---:|
| 2 | 1 | 16 | 64 |
| 3 | 4 | 64 | 512 |
| 4 | 7 | 112 | 1,792 |
| 5 | 2 | 32 | 1,024 |
| 6 | 1 | 16 | 1,024 |
| Total | 15 | 240 | 4,416 |

The three physical lengths therefore account for 13,248 specified supports.
Every support is evaluated under its complete single-support multiplier
envelope. There is no multiplier degree or coefficient-field-size cap.

## 2. Exact compression and cost accounting

For a support \(g=1_{Hk=b}\), the parity of occupations in each group of equal
nonzero columns determines its homogeneous-layer pattern. Group parities are
independent because the groups contain disjoint physical modes. Their exact
probabilities are

\[
 \Pr(z_h=1)=\frac{1-\prod_{i:h_i=h}\tanh(t\lambda_i/2)}2,
\]

with \(\lambda_i+J\) in place of \(\lambda_i\) for the charged measure. The
implementation evaluates both parity probabilities separately, using
\(\operatorname{expm1}\) to preserve small tails when the product is close to 1.

The charged pattern distribution supplies the zero-place expectation. The
uncharged syndrome probability supplies the entire infinity-degree cost of
\(f_k=1+R g(k)\). Both costs are added to the full Gaussian pressure. Spectator
modes remain in that Gaussian term. No average-bias approximation or deletion
of an infinity cost is used.

At each trial scale and charge, the multiplier is optimized by the finite
automaton: continuous F4 cutoffs use their explicit minima, and each binary
periodic tail uses its two nearest admissible divergence positions. This is
the exact mathematical lower envelope; its numerical evaluation uses floating
arithmetic. Any potential strict improvement would still require a literal
finite polynomial and an exact pressure certificate.

## 3. Completed genuine-ratio searches

The outer search used \(\tau=nt\), with
\(-4\le\log\tau\le\log15\) and \(-6\le J/n\le6\). A fixed grid preceded
the seeded differential-evolution search. The objective was the positive root
\(\rho\) of the full pressure. Absolute pressure at vanishing entropy was not
used as the optimization objective.

Every rank returned the same baseline value at a given physical length:

| Physical length | Best recorded \(\rho\) | Full pressure at \(1/\sqrt2\) | Physical entropy at the root |
|---|---:|---:|---:|
| 64 | 0.7075072433080506 | 0.0195069230211 | 24.3623484370 |
| 128 | 0.7072351467479230 | 0.0109546875847 | 42.6737441957 |
| 256 | 0.7071468886858093 | 0.00605187413939 | 75.4477996850 |

The complete mode assignments, automaton hashes, outer-search parameters,
selected multiplier descriptions, and pressure roots are stored in
[`global-code-multiplier-searches.json`](../evidence/global-code-multiplier-searches.json).
All selected multiplier descriptions are the baseline.

## 4. Reproduction and independent replay

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_boolean_multiplier_automata.py
OPENBLAS_NUM_THREADS=1 python3 tools/search_global_code_multipliers.py \
  --output evidence/global-code-multiplier-searches.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_global_code_multipliers.py
```

The first command verifies the exact algebra and unbounded-tail reduction.
The last command regenerates all support assignments and automaton hashes,
checks that each assignment and its spectators cover the full physical energy
list, and replays all 15 pressure roots. All replayed physical entropies are
positive. The largest absolute root residual is below \(4.3\times10^{-10}\).

The replay record is
[`global-code-multiplier-replay.json`](../evidence/global-code-multiplier-replay.json).

The next larger constructive family assigns an independent polynomial to each
syndrome. Those joint coefficient vectors leave the one-support form
\(1+R g_b\), so this complete multiplier reduction does not close that family.
