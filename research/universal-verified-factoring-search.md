# Universal verified search already attains the supremum factoring score

## Statement and scope

There is one deterministic finite-alphabet Turing machine whose factoring
score is at least the score of every deterministic factoring machine.
Consequently, if the supremum of achievable scores is finite, this one
machine attains that supremum under the challenge's `o(1)` convention.

This is an application of universal verified search, not a new
optimality principle or a new numerical factoring bound. It explains
why an existential uniform endpoint does not by itself require a
tensor-specific seed selector. The direct selector additionally gives
a uniform construction inside the tensor method and exposes its local
widths, field choice, and resource costs.

The only external algorithmic input below is deterministic polynomial
time primality testing: Agrawal, Kayal, and Saxena,
[*PRIMES is in P*, Annals of Mathematics 160 (2004), 781--793](https://doi.org/10.4007/annals.2004.160.781).
The argument uses only that theorem, not a particular polynomial degree.

## 1. A polynomially checkable output

On input the binary integer `N>=2`, a valid output record is a finite
list of primes, with repetition permitted, whose product is `N`.
Sorting is optional. If `n` is the input bit length, the list has at
most `n` entries and its total bit length is `O(n)`: the sum of the
binary logarithms of its entries is `log2 N`, and there are at most
`log2 N` entries. An ordinary delimiter encoding still has length
`O(n)`.

One deterministic verifier checks the syntax and the linear output
length bound, checks primality of every entry, and multiplies the
entries to compare their product with `N`. It runs in `n^O(1)`
finite-alphabet steps. Intermediate products can be rejected as soon
as they exceed `N`. A valid record therefore gives the complete prime
factorization, and its smallest entry gives the required least prime
factor. In particular, a record consisting only of `N` is accepted
exactly when `N` is prime.

There is no assumption that candidate programs are correct, total,
or provably correct. Their outputs are accepted only by this verifier.

## 2. A uniform geometric simulation

Fix an effective binary encoding of deterministic finite-alphabet
machines and a universal simulator. At stage `j=1,2,...`, enumerate
every code of length at most `j`. For a code of length `l`, simulate
its program on `N` for `2^(j-l)` simulated steps, starting again from
its initial configuration. If it halts with an output record within
the fixed length bound, run the verifier. Return the smallest entry
of the first accepted record.

The total simulated steps in stage `j` are at most

```math
\sum_{l=0}^{j}2^l2^{j-l}=(j+1)2^j.
```

Universal simulation, copying the input, enumerating programs,
handling their bounded output strings, and verifying their answers
can all be implemented with one fixed polynomial overhead. For
example, a simple tape-configuration encoding and scan-based
simulation suffice; no unit-cost unbounded-word operations are used.
The total work through stage `j` is therefore at most

```math
\operatorname{poly}(n+j)\,2^{O(j)},                    \tag{1}
```

with absolute constants determined by the universal machine and the
verifier. A basic exhaustive factoring program guarantees termination
on every input, independently of any faster program.

Let `P` be a fixed program that outputs a correct complete prime
factorization on every input. Write `l_P` for its fixed description
length and `T_P(n)` for its worst-case running time on `n`-bit inputs.
At the latest, it succeeds at stage

```math
j=l_P+\lceil\log_2(T_P(n)+1)\rceil.
```

Equation (1) gives a fixed integer `c` and a constant `C_P`, depending
on the program but not the input, such that the universal machine has

```math
T_U(n)\le C_P\bigl(n+T_P(n)+1\bigr)^c.              \tag{2}
```

The scheduler does not know `P`, its correctness, its running-time
bound, or `C_P`. These quantities occur only in the proof of (2).
The machine itself is finite, deterministic, and independent of the
score tolerance.

## 3. Least-factor programs and the score endpoint

A correct least-prime-factor program produces a complete prime
factorization by repeated calls on the remaining cofactor and exact
division. There are at most `n` calls. If `T_P` is not monotone, this
wrapper is bounded by a polynomial in `n` times

```math
\max_{1\le m\le n}T_P(m).
```

Every fixed positive-log-saving asymptotic bound in question has an
eventually increasing majorant, so this envelope has the same score.
The finitely many smaller inputs contribute only a fixed constant.
Thus using complete factorization as the verified output does not
exclude any score achievable by a least-factor machine.

For the challenge, score at least `s` means that for every fixed
`epsilon>0`, the same machine eventually satisfies

```math
T(n)\le 2^{\,n/(\log n)^{s-\epsilon}}.             \tag{3}
```

A fixed polynomial overhead in `T(n)+n` preserves this score: a
constant factor multiplying the runtime exponent is absorbed by any
positive extra power of `log n`. Therefore (2) makes `U` at least
as good in score as every correct program.

More explicitly, let `s_*` be the finite supremum of all achievable
scores. Given `epsilon>0`, choose a correct program of score greater
than `s_*-epsilon/4`. Apply its score bound with another
`epsilon/4` of slack, and then use (2). Eventually the constant
polynomial overhead is absorbed by the remaining logarithmic slack,
giving (3) with `s=s_*`. The chosen program and eventual threshold
may depend on `epsilon`; the universal machine does not.

No effective enumeration of correct programs, computable rate of
improvement, or attained optimal finite certificate is required.
If the supremum is infinite, the same argument gives every fixed
finite score to `U`; that conclusion alone does not imply polynomial
time factoring.

## 4. What this observation does not preserve

The argument does not prove a stronger numerical score, identify a
fast practical program, or provide a useful threshold for its
asymptotic bound. Its constants can be enormous. Nor does it preserve
the sharp leading constant in a runtime exponent.

It also does not automatically preserve a simultaneous sublinear
working-space bound. The universal simulator and polynomial-time
primality verifier above were charged by time, and may use polynomial
space in their time or input bounds. A direct construction with a
specified time-space tradeoff remains additional mathematical content.

Pratt certificates offer an elementary polynomial-time verification
alternative when the candidate already outputs suitable certificates.
However, replacing a fast factoring program by a certificate-producing
program requires a bound on certificate generation, including its
primality witnesses. Merely knowing that short certificates exist does
not prove that this transformation preserves the factoring score.
Using deterministic polynomial-time primality testing avoids that
extra obligation.

For the finite-seed application, the universal-search conclusion is
especially simple. Each individual valid finite seed gives an ordinary
factoring program with its own field and constants. Universal verified
search combines all those programs without mixing coefficient fields
inside any contraction. Hence the supremum of their numerical scores
already belongs to one uniform machine. The explicit bounded seed
selector is a construction-level realization of that conclusion.
