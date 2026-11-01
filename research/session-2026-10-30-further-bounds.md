# Further-bound research session: 30 October 2026

## Outcome and scoring target

The strongest complete deterministic factoring exponent in this branch is
still **one half**. This session adds ordinary proofs, exact algorithms,
new crossing constructions, and reproducible certificates. It does not
establish an exponent above one half. No new Lean formalization was added.

The [challenge specification](https://github.com/bigospeedrun/site/blob/main/data/integer-factoring.json)
scores one deterministic algorithm that returns the smallest prime factor
on every input. A positive score `alpha` means that, for every fixed
`epsilon>0`, its time is eventually at most

```math
2^{b/(\log b)^{\alpha-\epsilon}}.
```

The specification checked during this session lists this repository's
one-half entry in the **claimed** tier. The work here makes no claim of
new external acceptance. Polynomial factors and improved finite constants
alone do not change the score.

The branch is `research/further-bounds-20261030`, based on
`e91442d3162cae7f775e008efb2ea998c4478344`, and is collected in
[PR 26](https://github.com/1337isnot1337/integer-factoring/pull/26).
The [research index](INDEX.md) links every theorem and experiment.

## 1. Positive algorithmic improvements

### The optimal divisor-count parameter

The [exact divisor oracle](exact-divisor-count-parameter.md) computes

```math
D(n)=\max_{1\le N<2^n}\tau(N),\qquad
q_{\mathrm{exact}}(n)=\lceil\log_2 D(n)\rceil.
```

Thus `q_exact` is the smallest possible integer parameter in any bound
`tau(N)<=2^q` depending only on this input-length bound. This is a proved
finite optimum, rather than a tested upper bound.

| Bit bound | Maximum divisor count | Earlier split parameter | Optimal parameter |
| ---: | ---: | ---: | ---: |
| 8 | 20 | 5 | 5 |
| 16 | 120 | 9 | 7 |
| 32 | 1,920 | 17 | 11 |
| 64 | 184,320 | 33 | 18 |
| 128 | 318,504,960 | 65 | 29 |

The proof reduces the search to nonincreasing exponent vectors on the
first primes. Elementary prime bounds and a partition generating function
give `2^O(sqrt(n)) poly(n)` bit time. Compressed partitions and recomputed
multiplication bits give a complete construction using
`O(sqrt(n) log(n+2))` working bits. The practical Python oracle uses ordinary
integers and polynomial workspace; the complete compressed low-space
search is specified by the proof, and its multiplication-bit primitive is
implemented and independently checked.

The separate [split bound](optimized-divisor-count-bound.md) improves the
earlier auxiliary leading coefficient nine to one. The exact oracle has
the same asymptotic `q_exact(n)~n/log2(n)` and improves finite values.
Strict count recovery still needs `2^b>D(n)`: when `D(n)` is a power of
two, that exponent is one larger than `q_exact`.

### Field and jet linear algebra with small auxiliary memory

The [field-rank construction](polylog-space-rank-recovery.md) computes rank,
consistency, and requested affine coordinates using
`m^O(log m) poly(log q)` time and
`O(log m (log m+log q))` auxiliary bits, in addition to the stated entry
oracle resources. Polynomial weighting isolates a nonzero minor;
recomputed prefix bases avoid retaining pivot lists.

The [jet-ring construction](polylog-space-jet-factorization.md) works over
`F_q[t]/(t^H)`, including nonunits. If `R_H` and `R_(H-1)` are the two
appropriate Toeplitz field ranks, its minimum inner factor dimension is
exactly `R_H-R_(H-1)`. Selecting generators from the full image modulo its
radical avoids a false prefix-basis argument. Output entries can be
queried without retaining the complete factor arrays.

These are resource improvements to auxiliary operations. In particular,
an exponentially large matrix can still require `2^O(lambda^2)` time for
this small-space query algorithm; the notes do not silently replace that
cost by `2^O(lambda)`.

## 2. A more general route from certificates to a score

The [vanishing-slack transfer](vanishing-slack-score-transfer.md) allows
width `4 rho lambda+o(lambda)` and gives score
`min(log2(1/rho),delta)` under its stated entropy and effective-construction
hypotheses. The [gap-free finite-seed theorem](gap-free-seed-score-transfer.md)
uses bounded search to attain an enumerable certificate infimum without
knowing a convergence modulus or how quickly improving seeds appear.

[Exact finite-field specialization](finite-field-rank-specialization.md)
preserves finite rank profiles, including those originally over
characteristic zero. [Input-wise field selection](universal-field-seed-selection.md)
chooses one field for the entire input and charges all validation costs.
Size filters precede expensive validation, so malformed records cannot
force enormous field searches. Prime-field seeds give the same infimum
after amplification and projection.

For a finite crossing on `d` labels, write

```math
S=\sum_i\sqrt{R_iC_i},\qquad
\rho=\frac{\log_2 S}{2\log_2 d}.
```

The numerical target for beating one half through this method is
`S<d^sqrt(2)`, together with the transfer's effective identity and
factorization obligations. No supplied new family reaches this target.

The separate [verified universal-search argument](universal-verified-factoring-search.md)
explains why a single machine can existentially attain the supremum of
achievable scores. It is a uniformity result, not a numerical improvement.

## 3. Trace kernels: a disproved obstruction and an open synthesis problem

The [binary trace construction](trace-bilinear-conservation-kernels.md)
uses `F_(2^t)` to give individual conservation tensors with both cut ranks

```math
r=c=(2^t-1)2^k.
```

Their ordinary kernel rank is at least `binom(t,k)`. This disproves a
proposed general inequality `rc>=d rank(K)^2`; that inequality therefore
cannot rule out this route to better pieces.

A [complete coefficient certificate](trace-kernel-small-certificate.md)
sharpens the counterexample to `t=19,k=7`:

```math
d=524288,\qquad \operatorname{rank}K=94183,\qquad
r=c=67108736.
```

The exact comparison is

```math
r c=4503582447517696
 <4650663930232832=d\,(\operatorname{rank}K)^2.
```

All 153 smaller parameter pairs in this trace-product family are excluded
by the proved support bound. Thus nineteen is minimal **within that
family**. The certificate does not materialize the enormous four-leg
tensor: it checks the complete multiplicative coefficient support and
independent field identities.

[Restricted-label kernels](trace-kernel-restriction-and-composition.md)
have an exact rank formula through syndrome-image intersections. This
allows precise rank evaluation on arbitrary label subsets. The same note
proves low-degree composition cancellations on the full additive label cube
in characteristic two.

A separate [two-patch construction](two-patch-trace-crossing.md) supplies
an actual complete crossing with polynomial degree exactly `2t-1`.
Its five tested raw coefficient rank sums exceed `d^2`, so it does not
convert the good individual trace pieces into a competitive decomposition.
At widths three and four, [complete binary-span optimization](two-patch-span-recombination.md)
proves the stronger exact minima `d^2`: arbitrary binary recombination
and any number of pieces cannot improve the one-piece crossing in those
two spans. All 320 span tensors and 69,312 potential inequalities are
checked. This is specific to F2 coefficients and the specified field models and trace bases.

A [tagged-coordinate follow-up](tagged-two-patch-span-controls.md) also
completes every single weight-two tag at those widths: fourteen profiles,
4,864 tensor evaluations, and 2,190,592 exact potential inequalities.
Every optimum is again `d^2`. The [bivariate follow-up](bivariate-two-patch-span-control.md) proves
that the containing coefficient span has dimension exactly `4t-1` at
every width, over every characteristic-two coefficient field. Its
complete `t=3`, tag `(0,0)` binary span also has exact optimum sixty-four:
all 2,048 tensors, 4,096 dense cut comparisons, and 4,192,256 potential
inequalities are checked. Other bivariate tags and widths remain open.

The remaining constructive task is to synthesize the full crossing while
controlling total rank cost, possibly using restrictions or different
trace bases. A large ordinary kernel rank alone does not accomplish this.

## 4. Observable tables: a wider interface and exact exclusions

The [observable realization](observable-transition-realization.md) lifts
a finitely certified matrix table to diagonal rank-one auxiliary pieces.
It includes examples without a common left or right crossing eigenvector.
The [word transfer](diagonal-rank-one-word-transfer.md) gives exact ranks
for individual auxiliary words without cross-half commutation.

A [linear detector](linear-rank-detectors-for-observable-tables.md) then
gives lower bounds on the resulting word-rank spectral base. A [universal choice](local-alphabet-and-pairing-spectral-barriers.md) yields

```math
\operatorname{spr}(G_{\mathrm{live}})
\ge \frac{2d^3}{d^2+1}.
```

This excludes improving constructions on two, three, or four external
labels **within this spectral method**. In characteristic two, the scalar
entry space spanned by the three standard pairings is excluded at every
alphabet and every auxiliary dimension. These results do not lower-bound
the true ranks of grouped pieces, where cancellation can help.

[Finite functional certificates](finite-functional-rank-detector-certificates.md)
extend the detector method to specified larger entry spaces. Each
certificate is checked over its stated coefficient field; a supplied
counterexample explains why F2 validity does not imply extension-field
validity. On five labels, [all six overlap orbits](overlapping-subalphabet-spectral-certificates.md)
of two or three restricted three-label pairing corrections have spectral
lower bounds at least thirteen over F2. This covers all 165 labelled
families and any finite auxiliary dimension.

A further [twisted-pairing catalogue](twisted-pairing-observable-preflight.md)
contains 7,500 distinct index choices, with possible repeated scalar
spaces. All 112,500 tensor evaluations give a crossing-coefficient
detector of weight eleven. This excludes every auxiliary size for this
F2 catalogue before auxiliary matrices are searched. It does not exclude
arbitrary twisted pairings.

Consequently, further auxiliary-state searches should first test the
scalar entry space. A detector can exclude an entire space before any
large matrix enumeration. Escaping a detector is only a necessary step;
it does not itself produce a crossing or an upper bound.

## 5. Noncommutative curves and exact channel growth

[Mismatched mode weights](mismatched-mode-weyl-crossings.md) give coefficient
algebras outside the earlier additive family. A [short commutator](mismatched-commutator-projector.md)
constructs an all-energy projector and, at energy 63, a paired-sector
tensor with cuts `(361,1444)`. It is proved outside the polynomial algebra
of the constant operator at that energy.

[Partial-mode curves](partial-mode-weyl-sweeps.md) isolate adjacent exchanges
while retaining a complete endpoint. Using every adjacent edge once gives
degree `O(h log h)` and a unique intermediate occupation path for each
entry. This is an effective noncommutative interface; the current raw
profiles and constant powers still fail the finite target.

For channel graphs, [integer base walks and stationary flows](exact-channel-type-growth.md)
give the exact limsup rate of matching-type counts, including connectivity
and congruence restrictions. A [finite carry graph](finite-carry-channel-feasibility.md)
decides exact positive-support feasibility with finite certificates for
both outcomes. These results make the count bound computable to arbitrary
rational precision. They do not identify actual tensor ranks with counts.

## 6. Verification, finite screens, and preserved limitations

The main reproduction commands are:

```sh
make verify
make verify-search-native
python3 tools/check_exact_divisor_maximum.py --large
python3 tools/check_trace_weight_restriction.py --full --output /tmp/trace-restriction.json
```

The integrated-run receipt is recorded separately in
[the session verification evidence](../evidence/session-2026-10-30-verification.json).
The receipt identifies the exact tested source tree. Documentation added
after that tree receives its own link and whitespace checks.

Earlier focused verification includes the full energy-32 Weyl orbit and
protected minor, all 16,824 recorded site-weight rank pairs, all 33 selected
eight-feature profiles, the exact trace coefficient support, and the
new divisor, restriction, observable, overlap, and partial-mode packages.
The complete finite screens and their precise scopes are linked in the
[index](INDEX.md); none is extrapolated to an unproved all-size statement.

The width-18 run-count experiment is explicitly **incomplete**: 3,367 of
4,583 types were computed. Its restartable checkpoint is preserved in the
[polynomial-shear note](polynomial-shear-diagnostics.md). There is no
claimed width-18 final profile. The experiment was stopped after a pinned
theorem from [PR 27](https://github.com/1337isnot1337/integer-factoring/pull/27)
excluded its whole construction family; the note cites the precise source.

No local Lean build was performed. The new results rely on the ordinary
proofs and exact checkers stated in their notes. The strongest complete
factoring exponent and the challenge entry remain one half.
