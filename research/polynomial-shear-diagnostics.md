# Polynomial shear controls and a partial long-profile checkpoint

Work over a coefficient field of characteristic two.

## 1. A targeted family beyond constant conjugation

These controls test whether polynomial conjugation can improve the
coefficient ranks of the elementary binary pairing pencil. The endpoint
identity is exact, but all 27 tested coefficient decompositions fail the
finite rank-sum target.

Use normalized row order `(x,y)|(v,u)`, so the crossing is the identity.
Let `C_i` swap the two physical bits at site `i`, and put

```math
F_t(z)=\prod_{i=1}^{t}\bigl(I+(1+z)C_i\bigr),
\qquad F_t(1)=I.
```

On two consecutive sites, use the nilpotent shear

```math
K=|00,11\rangle\langle01,10|
   +|11,00\rangle\langle10,01|,
\qquad K^2=0.
```

For a fixed binary polynomial `p`, the gate `H_i(z)=I+p(z)K_i` satisfies
`H_i(z)^2=I`. Thus, for an ordered list of gates,

```math
\widetilde F_t(z)=G(z)F_t(z)G(z)^{-1},
\qquad G(z)=H_s(z)\cdots H_1(z),
```

is polynomial and has `tilde F_t(1)=I` as its exact endpoint. Unlike constant
conjugation, this operation can change individual coefficient row ranks.
With each adjacent edge used at most once, the degree is at most
`t+2(t-1)deg(p)`. There are therefore `O(t)` coefficients for fixed `p`.

The complete polynomial tensor also has bounded matrix-product auxiliary
dimension. At each cut between consecutive sites, at most one gate from
`G` and one from `G^(-1)` cross the cut. Each two-site gate on local
physical dimension four has an operator factorization with at most 16
terms. Multiplying these two factorizations gives auxiliary dimension at
most 256, independent of `t`. This is an endpoint and description bound;
it is not a useful rank bound for the resulting coefficient pieces.

## 2. Complete finite outcome

The screen uses `p=z,1+z,z+z^2`, lengths `t=2,4,6`, and three gate orders:

* `chain`: every adjacent edge in increasing order;
* `pairs`: disjoint edges starting at positions `0,2,4,...`;
* `layers`: those even edges followed by the odd edges.

These are exactly 27 profiles. Every profile has complete integer row
and column rank lists in `polynomial-shear-verification.json`. The
endpoint is checked directly. An independent implementation constructs
the full polynomial row matrices at lengths two and three, applies the
literal polynomial shears, and ranks both full four-leg flattenings.
All 18 comparisons with the parity-block engine agree.

For each complete profile the checker computes the integer lower bound

```math
L=\sum_j\lfloor\sqrt{r_jc_j}\rfloor.
```

Every tested profile has `L>4^t=d^2`. Thus every one fails the stricter
target `S<d^sqrt(2)` by exact integer arithmetic. No floating-point
threshold comparison is needed.

For orientation, the smallest observed rank sum at each tested length is:

| Length | Polynomial and order | Integer enclosure for `S` | Decimal diagnostic |
|---:|---|---:|---:|
| 2 | `1+z`, any listed order | `[24,26]` | 24.92747 |
| 4 | `1+z`, disjoint pairs | `[497,502]` | 498.99299 |
| 6 | `1+z`, disjoint pairs | `[9249,9257]` | 9254.43752 |

The test rules out only these literal finite profiles. It does not
exclude recombinations, different gates, or all polynomial conjugations.

Run from the repository root:

```sh
python3 tools/check_polynomial_shear_probe.py --output /tmp/polynomial-shear-verification.json
```

The checker regenerates all 27 profiles, validates the complete parameter
set, records a hash of the exact rank lists, and rereads the final JSON
after its atomic write. `polynomial_shear_probe.py` is the standalone
profile engine; its optional timings are not part of the exact-rank hash.

## 3. A deliberately incomplete width-18 computation

The run-compressed count/run profile at width 18 was stopped after
**3367 of 4583 run types**. The remaining 1216 types were not computed.
There is no completed width-18 profile and no reported width-18 rank
exponent from that run.

The preserved file
`two-state-run-count-bigrading-t18.json.checkpoint` is an atomic partial
accumulation for parameters

```text
a_S=0, a_C=10, b_SS=0, b_SC=1, b_CS=0, b_CC=0.
```

Its SHA-256 is
`7ffeef031a0371e9e30015b09f06f3a060489b346c0c1b835e56331ab819d532`.
The restartable wrapper is `two_state_run_profile_checkpoint.cpp`, with
SHA-256
`d8beb24d0316ffcb8102da937fce44f9f55062a256291cb3ff8de47a2506e779`.
It validates the saved width and costs and replaces the checkpoint
atomically after each completed type. The rank cache is not serialized,
so recomputation counters can differ after a restart.

An optional explicit restart is:

```sh
g++ -O3 -std=c++17 tools/two_state_run_profile_checkpoint.cpp -o /tmp/run-profile-checkpoint
cp evidence/two-state-run-count-bigrading-t18.json.checkpoint /tmp/run-count18.json.checkpoint
/tmp/run-profile-checkpoint profile 18 /tmp/run-count18.json 0 10 0 1 0 0
```

The partial accumulation is retained for reproducibility, not as a proposed
next score improvement. The full two-state word family is already excluded
by
[the pinned diagonal-word theorem](https://github.com/1337isnot1337/integer-factoring/blob/44832e5316a4bd8bd6caebe962c31aa6897fb77d/research/diagonal-word-rank-sum-barrier.md).
In the original row order, its two symbols are `C=I` and
`S=E+C=I+E=N`; every count/run coefficient belongs to
`tensor(span{I,N})`. That theorem permits arbitrary ordered coefficients,
groupings, and numbers of pieces, and proves `S>=d^(3/2)`. It therefore
applies independently of the unfinished finite computation. This inclusion
does not by itself cover the larger three-symbol or eight-feature families.
