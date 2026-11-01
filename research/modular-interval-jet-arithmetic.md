# Precision-bounded Newton arithmetic for long interval jets

The [unit-progression jet formula](normalized-interval-jets.md) computes
products over a binary-encoded interval without listing its elements.
Its original implementation forms exact elementary symmetric functions
before reduction. The new
[modular_interval_jets.py](../tools/modular_interval_jets.py) retains
only a rigorously sufficient number of bits during the same computation.

For a target precision `p` and degree `d`, all Newton coefficients use
at most

\[
 P=p+v_2(d!)=p+d-s_2(d)
 \tag{1}
\]

bits. Here `s_2(d)` is the number of ones in the binary expansion of
`d`. In the interval-product application `d<p`, so `P<2p`. The endpoint
length still enters through its binary representation; no iteration
depends on its numerical value.

## 1. Modular power sums from binomial coefficients

For `0<=r<=d`, use the integer identity

\[
 \sum_{j=0}^{L-1}j^r
 =\sum_{i=0}^r
     \left\{\begin{matrix}r\\i\end{matrix}\right\}
      i!\binom L{i+1}.
 \tag{2}
\]

The bracket is a Stirling number of the second kind. Equation (2)
follows by expressing `j^r` in the falling-factorial basis and summing
`binom(j,i)` over the interval.

Compute the small Stirling table, factorials, and binomial values
modulo `2^P`. To update `binom(L,i)` from `binom(L,i-1)`, remove powers
of two separately from the numerator `L-i+1` and denominator `i`.
Retain the binary exponent and invert only the odd denominator. Keep
the unit to `P` bits even when the whole binomial is currently zero
modulo `2^P`; its valuation can decrease at a later update.

This yields all power sums with `O(d^2)` operations on `O(P)`-bit
coefficients, plus arithmetic on `L` and the small indices. A general
interval is the difference of two prefix sums. No large exact power
`L^r` is required.

## 2. Newton division has an explicit precision invariant

Let `S_j` be the interval power sums, and `e_j` its elementary
symmetric functions. Newton's identity is

\[
 k e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}S_j.
 \tag{3}
\]

At step `k`, retain `e_k` modulo

\[
 2^{P-v_2(k!)}.
 \tag{4}
\]

Every earlier coefficient on the right of (3) is known to at least
`P-v_2((k-1)!)` bits. Thus the numerator is known to that precision.
Its residue is exactly divisible by `2^v_2(k)`. Divide that power of
two in the ordinary integers, then invert the odd part of `k` at the
remaining precision. This produces (4).

By (1), every coefficient through degree `d` retains at least `p`
bits. This is a finite-precision induction, not division by a nonunit
inside `Z/2^p Z`.

## 3. Use in unit progression products

For a unit jet `u` and an even jet `v`, write

\[
 \prod_{j=\ell}^{r-1}(u+jv)
 =u^{r-\ell}\sum_{q\ge0}e_q(\ell,\ldots,r-1)(v/u)^q.
 \tag{5}
\]

If every coefficient of `v` has valuation at least `a>=1`, only
`q<=floor((p-1)/a)` is needed modulo `2^p`. If `v` also has formal
order `b>=1` and the jet degree is `D`, only `q<=floor(D/b)` survives.
The implementation uses both bounds, the interval length, and the
modular elementary coefficients above.

The interface agrees with `unit_progression_jet` in the original jet
module. The polynomial inverse-moment engine uses this implementation.
Its complete degree-eight-by-degree-eight bilinear control at matrix
dimension `2^200+13` now runs at 627 bits of normalized precision,
enough for 12 bits of absolute precision after the known `2^(-615)`
scale. The numerical dimension is never traversed.

This is an arithmetic and implementation improvement. It does not
evaluate the unresolved Hadamard-square trace or change the factoring
exponent.

## 4. Independent checks

[check_modular_interval_jets.py](../tools/check_modular_interval_jets.py)
compares 1,000 modular elementary computations with exact integer
Newton sums. It compares 1,000 random jet products both with the
original implementation and with literal multiplication of every
small-interval factor. A further control uses 302-bit endpoints,
degree 80, and precision 83, checked against exact integer arithmetic.
Its working precision is 161 bits.

The results are in
[modular-interval-jets.json](../evidence/modular-interval-jets.json).
