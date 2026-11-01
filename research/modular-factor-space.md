# Polynomial-space evaluation of the ambient modular factors

## Result

The original scale-one characteristic-two partition factors admit a
uniform evaluator using polynomial space and `2^(O(lambda))` time.
It stores neither the expanded dual polynomial nor an exponentially large
table of normal-order coefficients.

This is a local factor-evaluation result. A polynomial-space factoring
theorem additionally needs a depth-first evaluator and a space-conscious
nonzeroness test; those are separate algorithmic claims.

Use the normalized quotient construction in
[`modular-crossing.md`](modular-crossing.md), or the exact-half-rate
selection established in
[PR22's Section 7](https://github.com/1337isnot1337/integer-factoring/blob/c12da897fcc6081685e8c663217ce9083db87b4a/research/past-half-pieces.md).
The underlying source operators are from Anthropic's Appendix B at
[`16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
All new assertions here are ordinary proofs and exact finite checks.

## 1. Normal ordering on distinct-odd states

The label `mu` is a distinct-odd partition of weight `h`.
The input `sigma` and output `pi` are distinct-odd partitions of weights
at most `4h`. Put `ell=length(mu)`.

In a nonzero normal-order term, let `eta` be the parts that remain
unchanged from input to output. It is a subset of `pi intersection sigma`
and has size

```math
|\eta|=\frac{|\pi|+|\sigma|-\ell}{2},
\tag{1}
```

where vertical bars in (1) denote particle counts. If this is not a
nonnegative integer, the coefficient is zero.

For each possible `eta`, form a list of signed columns: the created modes
`pi minus eta` and the annihilated modes `sigma minus eta`.
There are exactly `ell` columns. Assigning the operator parts to them
gives a permanent. In characteristic two this equals a determinant,
but the permanent formula is convenient for streaming.
An input mode occurs once and is odd, so its differentiation multiplier
is one in `F_2`. No factorial or repeated-mode multiplier is missing.
The same physical mode may occur once as a creator and once as an
annihilator; these are distinct signed columns.

For the bottom factor, the two scalar entries for an operator part `n`
and mode `k` are

```math
C(n,k)=[(n-1)\subseteq(k-1)],\qquad
D(n,k)=[(n-1)\mathbin{\mathrm{AND}}k=0].
\tag{2}
```

Sum these permanents over all `eta`. This is exactly the source
normal-order coefficient. The normalization `z_mu^-1` is one because
every retained mode is odd and every multiplicity is one.

## 2. Fold the original dual and the angle index

The top factor uses
`a^mu=prod_(n in mu) sum_j binom(n,j)t_j` in characteristic two.
To include angle index `i`, attach `(1+S)^j` to an operator part `j`
and multiply the complete result by `(1+S)^h`; the coefficient of `S^i`
is the required `binom(h+weight,i)` weighting.

For a single original part `n`, the folded creator entry is

```math
\sum_{j\text{ odd}}\binom nj\binom{k-1}{j-1}(1+S)^j
=(1+S)S^{(n\mathbin{\mathrm{AND}}k)-1}.
\tag{3}
```

The folded annihilator entry is

```math
\sum_{j\text{ odd}}\binom nj\binom{k+j-1}{j-1}(1+S)^j
=(1+S)S^{n\mathbin{\mathrm{AND}}(\mathrm{NOT}\ k)}.
\tag{4}
```

All coefficients are reduced modulo two. These identities follow by
summing the allowed odd binary subsets of `n`: for an odd integer `b`,
`sum_(j odd, j subset b) T^j=T(1+T)^(b-1)`.
Substitute `T=1+S`. The creator's allowed subsets are those of
`n AND k`; the annihilator's non-low bits are those of `n AND NOT k`.

Even operator indices act as zero on the quotient. Repeated odd indices
cancel in the multilinear matching sum; equivalently their normal-ordered
fields act as zero. This is the same valid operator-label reduction
proved in the modular construction.

Extract one common factor `(1+S)` from each row. The entire top-factor
profile is therefore

```math
(1+S)^{h+\ell}
\sum_\eta\operatorname{per}M_\eta(S),
\tag{5}
```

where the creator entry of `M_eta` is
`S^((n AND k)-1)` and the annihilator entry is
`S^(n AND NOT k)`.
The coefficient of `S^i` in (5) is the angle-`i` top factor.
The source's alternating sign is one in this characteristic.

Every row exponent is at most `n-1`. Every permanent in (5) has degree
at most `h-ell`, and the complete profile has degree at most `2h`.
There is no expansion into exponentially many dual monomials.

## 3. Stream the permanent and the unaffected subsets

For an `ell` by `ell` matrix over `F_2[S]`, Ryser's identity is

```math
\operatorname{per}M
=\sum_{J\subseteq\{1,\ldots,\ell\}}
       \prod_{r=1}^{\ell}\left(\sum_{j\in J}M_{rj}\right).
\tag{6}
```

The signs in the usual formula disappear over `F_2`.
Enumerate `J` with an `ell`-bit counter. For each subset, compute the row
sums and multiply them sequentially, using coefficient arrays of degree
at most `h-ell`. Add the product into one accumulator. Then move to the
next subset. The scalar bottom factor is the degree-zero case.

The outer loop enumerating `eta` is also sequential. A conventional
combination iterator stores only its current subset. The small matrix
may be stored explicitly; it has `ell^2` polynomial entries, each of
degree at most `h`. Keeping this matrix is still polynomial space.
No results indexed by all subsets are memoized.

Distinct positive odd parts have sum at least the square of their
number. Consequently

```math
\ell\le\sqrt h,\qquad
|\pi|,|\sigma|\le2\sqrt h.
```

There are at most `2^(2sqrt h)` unaffected subsets and `2^sqrt h`
Ryser subsets for each. Ordinary polynomial arithmetic therefore gives
total time `2^(O(sqrt h)) poly(h)` and space `poly(h)`.
For `h=O(lambda^2)`, these are `2^(O(lambda))` time and polynomial space
in `lambda`.

The regular partition counts and canonical label unranking also use
polynomial space. A bounded-multiplicity partition dynamic program
through weight `4h` has polynomially many cells, each with `O(sqrt h)`
bits. It supports count queries and lexicographic unranking without
storing the exponentially many partitions. The factor algorithm only
retains the chosen input, output, and operator label at a time.

All these computations occur over the prime field. When factors are
used over a commutative `F_2`-algebra, the returned coefficients are still
bits; no large coefficient-algebra element is stored by this local routine.

## 4. Reproduction

Run:

```sh
python3 tools/check_modular_factor_space.py \
  --output evidence/modular-factor-space.json
```

The production factor routines in the checker have no memoization or
exponential tables. The verification section compares them with the
independent source recursion, including every angle index for its finite
cases. It separately checks the folded entries (3)--(4) against direct
binomial sums. Omitting the common angle factor supplies a negative
control. The finite checks verify the formula's transcription; the
space bound is the streaming argument above.

The source operators are prior work by Anthropic. This streamed-factor
analysis and checker used substantial OpenAI assistant support.
