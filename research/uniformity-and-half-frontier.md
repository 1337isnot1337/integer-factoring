# What a universal ceiling at one half would imply

The project has established `1/2` and proved several barriers for specified
crossing families. Those statements do not prove a ceiling for all
deterministic factoring algorithms. Standard derandomization consequences
help make the distinction precise:

- Assuming `P=BPP`, there is a uniform deterministic factoring algorithm
  in `exp(O(sqrt(n log n)))` bit operations.
- Unconditionally, advice depending on the bit length can select a family
  of successful random tapes giving that time bound for every input of
  the length. Constructing this advice efficiently is unresolved.

The first assertion is conditional; the second uses a different computation
model from the challenge. Neither is a new qualifying score, an explicit
derandomization, or a new Lean result. The established score remains `1/2`.
These are standard padding and probabilistic-method arguments, included to
answer the owner's hard-limit question and identify a uniform construction
target.

## 1. Fix the rigorous randomized baseline and its verification

Use the Lenstra--Pomerance expected-time bound already source-audited in
[arithmetic-framework.md](arithmetic-framework.md), from
[*A rigorous time bound for factoring integers*](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf).
On `n`-bit inputs it is `exp(O(sqrt(n log n)))`. Polynomial factors for
complete factorization, deterministic primality verification, bit-model
simulation, and taking the smallest prime factor fit within that bound.

Fix a sufficiently large integer constant and a time-constructible integer
bound

```math
T(n)=2^{O(\sqrt{n\log(n+2)})}
```

that dominates the expected time of the verified complete-factorization
routine by a factor of at least four. Small inputs can be absorbed into
the same fixed machine and bound. Cap a trial at `T(n)` steps and accept
only a complete list of deterministically verified primes whose product
is the input. A returned answer is correct. Markov's inequality, with
the chosen slack, gives success probability at least `3/4` for every
input. An unsuccessful trial returns failure.

Each capped trial reads at most `T(n)` independent random bits. Fixing a
uniform tape of that length therefore models the whole trial, including
all adaptive random choices. No independence assumption about its internal
arithmetic events is needed for the arguments below.

## 2. Padding proves the conditional uniform bound

Define one decision language on padded inputs

```math
L=\{\langle N,i,1^{T(n)}\rangle:
       \text{bit }i\text{ of }\operatorname{spf}(N)\text{ is }1\},
\qquad n=\operatorname{bits}(N).
```

Reject malformed padding and out-of-range indices deterministically.
Checking the prescribed padding length is polynomial in the presented
input length; the algorithm can reject a length mismatch without writing
the missing padding. For a valid padded input, run the capped randomized
routine and output the requested bit on success, or zero on failure.
This has bounded error and runs in time polynomial in the padded length
`M=Theta(T(n))`. Hence `L` belongs to `BPP`.

If `P=BPP`, this **one fixed language** has a deterministic algorithm
taking `M^c` steps for one constant `c`. On an ordinary input `N`, construct
the padding and ask for its at most `n` output bits. The total bit cost is

```math
O(T(n))+n\,T(n)^{O(1)}
=\exp(O(\sqrt{n\log n})).                          \tag{1}
```

The deterministic machine is uniform across lengths; its constant does
not depend on `N`. For every fixed real `alpha>0`,
`sqrt(n log n)=o(n/(log n)^alpha)`. Thus (1) beats every fixed logarithmic
saving exponent in the challenge's scale, conditional on `P=BPP`.

Consequently, a universal claim that no deterministic algorithm can beat
score `1/2` would imply `P!=BPP`. A barrier for a particular crossing
ansatz has no such implication. This is a conditional logical consequence,
not a proof that `P=BPP` or a claim that the conditional machine is supplied.

For primary context, [Impagliazzo--Wigderson (1997)](https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/IW97/proc.pdf)
prove `P=BPP` under the hypothesis that some language in
`E=DTIME(2^(O(n)))` requires circuits of size `2^(Omega(n))` at all
sufficiently large lengths. Combining that stated hypothesis with the
padding argument gives the same conditional factoring bound. No such
circuit lower bound is proved in this repository.

## 3. A small universal tape family exists without an assumption

For a fixed `n`, choose `n+1` independent uniform tapes of length `T(n)`.
For each fixed input the probability that all trials fail is at most
`(1/4)^(n+1)`. There are at most `2^n` inputs, so the union bound gives

```math
\Pr[\text{some }n\text{-bit input is missed}]
\le2^n4^{-(n+1)}<1.                               \tag{2}
```

Therefore a single family of at most `n+1` tapes succeeds for **every**
input of that bit length. It has `O(n*T(n))` advice bits, and trying all
its members costs `O(n*T(n))` steps, up to the already charged verification
overhead. This is an unconditional nonuniform existence statement.

A different family may be needed for each length. Its existence does
not put the infinitely many families into one finite-alphabet machine.
Exhaustively validating a candidate against all inputs of the length
also incurs a `2^n` factor; that procedure gives no desired uniform bound.

## 4. A quantitative uniform target

A sufficient specialized construction is one deterministic procedure that,
given `n`, emits a tape family with the covering property proved by (2). Charge
its generation time `S(n)`, its tape count `r(n)`, and every trial. The
resulting uniform factoring cost is

```math
S(n)+r(n)T(n)+\operatorname{poly}(n).                \tag{3}
```

For example, producing at most `n+1` such tapes in
`S(n)=2^(O(n/(log n)^beta))` for a fixed `beta>1/2` would improve the
score. The tape-output length `O(n*T(n))` fits within that target.
The missing theorem is a uniform construction with the all-input covering
guarantee. The probabilistic existence proof does not provide it.

This target keeps the actual randomized factoring computation, rather than
requiring a generator that fools every circuit of comparable size. The
arithmetic class-group lane studies concrete sampling and witness-extraction
steps inside that computation. The global
[pole trace](interval-pole-trace.md) and
[Cauchy scalar](pole-cauchy-content.md) are separate arithmetic-specific
targets. None of these open interfaces currently establishes an improvement.
