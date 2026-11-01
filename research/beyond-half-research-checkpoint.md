# Research checkpoint beyond the half-rate theorem

This checkpoint connects the ordinary proofs and exact programs developed
on October 31--November 1, 2026. It introduces no new Lean formalization.
The strongest complete factoring result remains saving exponent **1/2**.
The [existing theorem](half-rate-factoring.md) and its machine statement
are unchanged.

The [scoring rule](https://bigospeedrun.com/integer-factoring) requires one
deterministic algorithm for every input, with all bit operations charged.
A successful finite factorization, a faster certificate, or an existence
proof for short relations does not establish a larger score.

## 1. Arithmetic: the remaining task is selection

The arithmetic route now has effective constructions for every stage after
a suitable represented smooth integer has been selected. Its uniform
generation and existence results use the
[pinned Cureton--Thorner finite-order Hecke theorem](polylog-prime-form-generators.md).
The exact algebraic decoders, relation checks, and finite factorization
certificates do not need that analytic input.

A further [consecutive-power selector](smooth-leading-power-selection.md)
now constructs its own candidate norms and coordinates. Its
[adaptive continuation](adaptive-power-window-selection.md) factors a fresh
144-bit semiprime after 8,232,730 norm tests, including the exhausted first
phase. This finite search uses no intermediate norm-factoring oracle.
A uniform bound on the required exponent window remains open.

| Stage | Established result | Main reference |
| --- | --- | --- |
| Generate prime forms | An explicit polynomial cutoff generates the full proper class group, including nonfundamental orders. Modular roots and primes are constructed deterministically. | [Cutoffs](explicit-prime-form-cutoffs.md), [factor base](prime-form-root-construction.md) |
| Prove witnesses exist | Powered small-bias prime samples give a short primitive smooth norm in every class. Their product support is used only in the proof. | [Powered prime norms](powered-prime-smooth-norms.md), [density](smooth-band-density.md) |
| Recover a representation | Deterministic orientation matching reconstructs coordinates from an unsigned factored norm. The relaxed decoder also handles imprimitive norms and inert square content. | [Primitive decoder](factored-norm-orientation-search.md), [arbitrary decoder](arbitrary-factored-norm-representation.md) |
| Build independent relations | A height bound alone suffices. Different target powers and logarithmic column weights give a nonsingular designated minor. | [Weighted relations](height-only-smooth-norm-relations.md) |
| Extract a factor | Repeated binary relation dependencies and exact halving find a useful ambiguous class. No relation determinant, class order, or annihilator is computed. | [Relation halving](determinant-free-ambiguous-extraction.md) |
| Select successful integers | A uniform deterministic bound at the desired score is still missing. Every failed candidate and retry must be counted. | [Integer-only interface](selected-norm-factoring-interface.md), [adaptive finite search](adaptive-height-norm-selection.md) |

### A precise contract that would suffice for score one

Let `Delta=-4N`, and let `b=bits(abs(Delta))`. Fix an integer `d>=2`.
The proved parameter family has

$$Y=b^{32d+o(1)},\qquad
  \log_2 H=(C_d+o(1))b,\qquad
  C_d=1+\frac{1}{2(d-1)}.$$

There are polynomially many prime-form generators. For a split base
prime `p`, let `g_p` be its chosen oriented class, and let `K_p` be the
least positive integer such that `p^{K_p}>H`.

A selector must cover the designated targets by returning an odd
`Y`-smooth integer `n<=H`, coprime to `N`, represented by `g_p^{K_p}`.
It need not supply a factorization, orientations, coordinates, or a
primitive representation. Uniform existence follows from the stronger
primitive witness theorem.

The arbitrary-norm decoder recovers a primitive represented divisor
`m=n/c^2` and its signed prime word `v`. Its logarithmic weight satisfies

$$\sum_q |v_q|\log q=\log m\le\log H<K_p\log p.$$

This gives the independent relation `K_p e_p-v`. Once every target is
covered, the full relation matrix and proved full generation allow the
halving extractor to return a proper factor of an odd integer with at
least two distinct prime factors. Primality, perfect powers, and recursive
splitting are handled as in the [arithmetic framework](arithmetic-framework.md).

For a norm with `B` bits, mixed-radix matching uses

$$\sqrt{(B+1)\tau(n)}\,\mathrm{poly}(b+B)
  \le 2^{(1/2+o(1))B/\log_2 B}\,\mathrm{poly}(b).$$

Consequently, if a selector uses time `T_sel` and emits `s` written
candidates, a valid total processing bound is

$$T_{\rm sel}
 +s\,2^{(C_d/2+o(1))b/\log_2 b}\,\mathrm{poly}(b)
 +\mathrm{poly}(b).$$

The count `s` includes candidates that are rejected. A deterministic
selector with total cost `2^O(b/log b)` would therefore establish score
one through this interface. Such a selector has not been supplied.
Coordinates would make post-selection processing polynomial for this
polynomial base and preserve any fixed proposed score.

The smaller-support parameter choice gives
`log_2 Y=(1/sqrt(2)+o(1))sqrt(b)`,
`Omega(n)<=(sqrt(2)+o(1))sqrt(b)`, and
`log_2 n<=b+O(sqrt(b)log b)`.
The resulting per-candidate construction and primitive decoding cost is
`2^((1/sqrt(2)+o(1))sqrt(b))`. Processing the complete generator list
and relation matrix must still be charged; this is not a selection bound.

A [uniform quasipolynomial-base specialization](quasipolynomial-base-norm-processing.md)
sharpens the leading coefficient without choosing an arbitrarily large
fixed polynomial degree. Set `q=ceil(log_2 b)` and `log_2 Y=32q^2`.
Then `log_2 H=(1+o(1))b`, and a full batch of supplied unsigned norms
can be processed in

$$2^{(1/2+o(1))b/\log_2 b}.$$

The entire powered-prime word support can also be enumerated within
that scale. Recovering the corresponding smooth pair remains open;
the support enumeration does not by itself select represented norms.

### Partial relations and adaptive targets

The extractor can use independent rectangular relation matrices before
the full list is available. If full generation is known, failure with
row deficiency `f` implies `omega(N)<=f+1`. It does not certify that an
arbitrary failed partial search has exhausted all useful two-torsion.

The adaptive selector may use any target power `k<=K_p` when its
recovered primitive norm satisfies `m<=H` and `m<p^k`. The bounded arbitrary-norm
decoder checks this condition across all possible contents. Filtering
only its first representation would be incomplete.

### What the finite factorizations establish

The searched fixtures have 60, 92, and 124 bits. The 124-bit fixture
contains 942 relation rows. Its integer-only replay receives target
primes and unsigned integers, independently reconstructs each relation,
and obtains a proper factor.

The [partial-search record](partial-norm-search.md) charges the actual
coordinate searches, including failed finite targets. The
[height-only consumer](height-only-selected-integers.md) and
[arbitrary-content consumer](arbitrary-factored-norm-representation.md)
exercise the relaxed contracts. These experiments validate particular
inputs and implementation behavior. Their practical cutoffs are not the
enormous cutoffs in the uniform analytic theorems.

### A constructed selector, including failed windows

The newer power selector reduces consecutive powers of each prime form,
tests the leading coefficient by an exact modular primorial power, and
divides successful values only by the known base primes. CRT composition
with an HNF fallback accelerates the consecutive powers. The coordinate
witness is automatically `(1,0)` in the reduced target form.

| Control | Base generators | Accepted rows | Norm candidates tested |
| --- | ---: | ---: | ---: |
| 60-bit input | 18 | 18 | 2,610 |
| 92-bit input | 85 | 82 | 218,738 |
| Original 124-bit input | 283 | 274 | 3,650,183 |
| [Fresh fixed-seed 124-bit input](fresh-leading-power-control.md) | 262 | 248 | 2,301,916 |
| [Fresh fixed-seed 144-bit input](adaptive-power-window-selection.md) | 513 | 499 | 8,232,730 |

These counts include all rejected candidates in the stated executions.
The 144-bit search first exhausted a 50,000-step window with 494 rows.
Continuing only unfinished targets added five rows and a factor after
154,884 more candidates. Doubling windows needs fewer than `2mW` norm
tests whenever all `m` targets have a successful window of length `W`;
base reconstruction and retained-state validation are also charged.

The existence theorem for represented smooth values does not imply that
a reduced leading coefficient is smooth. It therefore does not prove
the short-window or sparse-exponent-list hypothesis needed by this selector.
Finite success is not used to claim a worst-case bound.

A separate [dynamic leading-norm constructor](dynamic-leading-norm-closure.md)
factors smaller norms recursively and grows its prime base. Fresh unit
elimination keeps every newly retained prime below the square root of the
leading-norm height. [Odd-pivot elimination](odd-pivot-two-primary-elimination.md)
extends this compression while retaining genuine full-group identities.
The number of introduced and pending primes still needs a uniform bound.

## 2. Tensor certificates: constructive improvements

The new optimization changes channel bases at each declared lattice
state. It preserves the physical factorization required by the existing
all-length transfer theorem.

| Certificate | Exact improvement on the stated companion example |
| --- | --- |
| [State-dependent channels](state-adaptive-smith-channels.md) | At `s=2`, `B<=41/3<96/7`, improving the fixed-channel certificate. |
| [Constant basis optimization](matroid-adaptive-smith-channels.md) | At `s=5`, `B<=24073/1875`; the old fixed representation stays above `12.9` at every scale. |
| [Complete regular bases](dvr-adaptive-smith-channels.md) | `B<=32053/2500=12.8212`; a strict improvement over every constant-basis policy is proved on one side. |
| [Subspace oracle](subspace-smith-channel-oracle.md) | Linear kernels and weighted matroid bases replace exponential primitive-jet enumeration. The local oracle is polynomial in the explicitly written state family and precision. |
| [Blocking hierarchy](block-subspace-channel-hierarchy.md) | Two-site and three-site certificates give per-site bounds below `12.542` and `12.386`, respectively. Divisible blocking and precision enlargement have proved monotonicity. |

An [exact separation theorem](channel-pressure-separation.md) brackets
the one-site all-state optimum at `s=5` between approximately
`12.8203599666` and `12.8212`. Its lower certificate concerns the stated
global policy problem, including reducible cases. It need not control
every chosen boundary start.

These are improvements to the certificate machinery. They do not yet
beat the best factoring construction or give a larger saving exponent.

### An actual eight-letter crossing

The [nonorthogonal reflection retraction](eight-letter-reflection-retraction.md)
is an explicit characteristic-two crossing of actual degree 32. Its
rational common type and exact scale `s=(73/71)^18` certify `rho<18/25`.
It has strictly lower two-cut pressure than its matched coordinate
restriction for every finite `s>=3/2`. Both complete cut profiles,
the swap endpoint, and the scalar restorations are checked exactly.

This is a constructive gain within that family. The certified trial rate
`0.72` remains above `1/sqrt(2)`, and the unrestricted ambient reflection
seed remains numerically stronger. No global optimum or improved
factoring exponent is asserted.

## 3. Scoped obstructions guide the next construction

The strongest new alphabet result is the
[all-degree seven-letter adjugate floor](adjugate-seven-letter-pressure-floor.md):
bare ungauged local Smith pressure is greater than `709/1000`, hence
greater than `1/sqrt(2)`, over every field. The stated constant auxiliary
contraction extension has the same pressure conclusion. Independent
endpoint gauges, adaptive cancellation, and direct coefficient-rank
methods are outside its hypotheses.

An explicit eight-letter formal profile satisfies the current algebraic
constraints and has pressure below `707/1000`. No polynomial crossing
realizing that profile is known. A feasible relaxation is a construction
target, not a new algorithm.

Further exact results delimit particular constructions:

- [Smith orders bounded by the degree](bounded-order-smith-pressure-floor.md)
  force `rho>442/625=0.7072` through ten letters at every degree.
  This includes the stated eight-letter Bell polynomial paths. The
  general eight-letter formal profile has larger orders and remains
  outside this hypothesis.
- [Arbitrary binary coordinate codes](arbitrary-binary-code-pressure-floor.md)
  retain the elementary `rho>=3/4` floor, even with the stated independent
  diagonal gauges and paired embeddings.
- [Protected scalar subtensors](protected-subtensor-pressure.md) and
  [bounded boundary extraction](boundary-extracted-scalar-pressure.md)
  explain why the specified auxiliary extensions retain their scalar
  full-word pressure obstruction.
- [Clifford compression](clifford-compression-pressure-floor.md) proves
  low-order layers, near-full pressure floors, and 57 exact exclusions.
  The [threshold examples](clifford-threshold-sharpness.md) show where the
  generic layer argument stops.
- [Homogeneous nonorthogonal restrictions](reflection-homogeneous-restrictions.md)
  strictly improve a displayed coordinate deletion, while the
  [finite-family certificate](homogeneous-restriction-pressure-floor.md)
  proves `rho>=5/7` for all 120 restrictions at every scale.
- [Every coordinate code of two reflection sites](reflection-coordinate-code-pressure-floor.md)
  has the same `5/7` floor at the inherited gauges. The exhaustive
  classification covers 65,519 codes in 33,269 site-exchange orbits.

None of these is a lower bound on integer factoring. Each note states
the interface, fields, gauges, finite families, and optimization variables
covered by its conclusion.

## 4. Next work that could change the score

The direct arithmetic target is a deterministic selector for the
designated prime-form powers, with every unsuccessful candidate and all
setup work included in its complexity. The contract above permits much
more freedom than primitive oriented witnesses. The
[universal-list coverage bound](norm-selector-coverage.md) prevents one
short fixed list from covering every class at these heights. It does
not exclude a short list tailored to the polynomially many required
generator-power targets, even if those targets share the same list.

The constructive tensor target is an actual crossing or an exact
coefficient method outside the proved pressure restrictions. The
eight-letter formal profile, larger nonorthogonal restrictions, and
adaptive cancellation provide concrete tests. A new construction must
include its endpoint identity, both cut bounds, and the complete
uniform transfer cost before it can support a factoring claim.

The source-scope checks also matter. The
[recent ideal sampler](rigorous-ideal-sampling-scope.md) is randomized
under ERH. Its near-uniform support cannot be exhaustively enumerated
within the desired deterministic bound. A stronger numerical zero-free
result for Dirichlet characters alone does not replace the varying-field
Hecke input used here.

## Verification and status

Every linked packet includes an ordinary proof, an executable checker,
or a precisely scoped source audit. Evidence records include exact
arithmetic controls and source hashes where applicable. The research
index links the full packet set and replay commands.

The [verification record](../evidence/beyond-half-verification-20261101.json)
pins the executed suites, source hashes, finite replays, and analytic
dependency. `make verify-beyond-half` runs the core checkpoint checks.

Independent agent review and passing finite checks do not turn the
uniform selector obligation into a theorem. This checkpoint records
what is proved, what the programs consume, and the remaining cost that
must be bounded to improve the score.
