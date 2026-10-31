# Upper bounds for canonical Cauchy-packet visibility

## Results and remaining question

Two upper bounds distinguish canonical inputs from arbitrary rational
realizations of the same polynomial class.

* At packet length four, the largest possible minimum odd-center content
  among **all binary tensor pairs with odd original contraction** is 11.
  At length two it is 2. These are exact finite exhaustive results.
* At every packet length `L`, a zero-state indicator row tensor and
  nonnegative canonical `p`-bit column entries have an endpoint value
  visible using only `p+O(L)` bits, whenever their original value is odd.
  The endpoint value itself need not be odd; the claim is that it is
  nonzero at the stated precision.
* For arbitrary nonzero nonnegative canonical tensors, a suitable
  simultaneous Boolean translation makes one endpoint value positive and
  visible with `p_A+p_B+O(L)` bits. However, even all such translated
  endpoint values do not determine the original parity.

The first result matches the
[alternating family](canonical-alternating-packet.md) at `L=4`. It supports
the question whether a comparable `O(L log L)` visibility bound holds for
all canonical binary packet tensors. Neither result proves that question,
and neither supplies a local-data reconstruction algorithm for general
tensors. The factoring score remains `1/2`.

The deformation and canonical-creation convention are from
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`. We use the packet and
local-content notation in the [realization note](cauchy-packet-realization.md).

## 1. A precise finite upper statement

For `L=4`, let `A,B` range independently over all functions from the
sixteen Boolean states to `{0,1}`. If

```math
F(0)=\sum_X A(X)B(X)\equiv1\pmod2,
```

then at least one coefficient of at least one of the four local
polynomials `F(t+8z)`, `t=-3,-1,1,3`, is nonzero modulo `2^12`.
The alternating point tensors make every such coefficient divisible by
`2^11`, so 11 is the exact maximum of the minimum contents.

The finite input space has `2^32` tensor pairs. Exactly

```math
(2^{16}-1)2^{15}=2147450880
```

have odd original contraction. The optional C++ checker covers all of
them with an exact meet-in-the-middle calculation.

### Why the projection filter is exhaustive

The Python wrapper constructs all 116 coefficient matrices
`K_r[Y,X]` modulo `2^32`: four centers times 29 possible degrees.
Two fixed integer linear combinations of these matrices give projections
`H_1,H_2`. A pair whose every local coefficient vanishes modulo `2^12`
must satisfy both projected congruences.

For each of the 65,535 nonzero column tensors `B`, compute `H_1 B` and
`H_2 B`. Split the sixteen bits of the row tensor into two blocks of
eight. Store the 256 possible triples

```math
\left(A_{\rm left}^{\mathsf T}H_1B,
      A_{\rm left}^{\mathsf T}H_2B,
      \langle A_{\rm left},B\rangle\bmod2\right).
```

For each of the 256 right blocks, look up the unique required projected
values and the opposite parity. Every matching left block is retained,
including duplicate keys. This enumerates every pair with odd original
value satisfying the two necessary congruences. Every survivor is then
checked against all 116 original coefficient matrices.

No survivor passes at precision 12. The calculation is repeated with two
different fixed projection seeds, and a precision-11 positive control
is found and independently replayed in Python. The seeds affect only
the amount of work: the projections cannot exclude a genuine zero of
all the local coefficients. There is no probabilistic correctness claim.

At `L=2`, the Python checker directly examines all 120 odd-original
pairs among the 256 possible binary tensor pairs. Their largest minimum
odd-center content is exactly 2.

## 2. An all-size upper bound for a zero-state row

Let `L>=2` be a power of two, set `A` to the indicator of the zero
Boolean state, and let

```math
0\le B(X)<2^p.
```

Suppose `B(0)` is odd, so the original contraction is odd. Put

```math
b_L=\binom{2L-1}{L},\qquad t=L-1.
```

The only nonzero Cauchy entry at this endpoint is
`W_(0,L-1)(t)=-b_L`. If `N_0,N_1` are the sums of `B(X)` on the states
whose last bit is respectively zero or one, then

```math
F(t)=N_0+(1+b_L)N_1.
\tag{1}
```

Here `N_0>=B(0)>=1` and `N_0,N_1<2^(p+L-1)`. In particular (1) is a
strictly positive integer, and

```math
F(t)\le (2^p-1)2^{L-1}(b_L+2)<2^{p+3L-2}.
\tag{2}
```

The last inequality uses `b_L+2<2^(2L-1)`, which follows by including
the other positive binomial coefficients in the binomial sum. Therefore
the constant local coefficient at this one odd center is nonzero
modulo `2^(p+3L-2)`. A sharper finite bound follows by taking one more
than the floor of the binary logarithm of the first upper bound in (2).

This argument explains why
[the eight-wire rational lattice theorem](eight-wire-cauchy-lattice.md)
does not automatically provide a hard low-precision canonical input with
the same fixed row. Favorable odd-denominator rational entries need not
be nonnegative integers of small height.

For a known point-row tensor `A=delta_Y`, simultaneous Boolean translations
of both tensors send that row state to zero and preserve the original
contraction and canonical entries. The translated packet then satisfies
(2). Thus the alternating lower witness is a statement in the fixed
Cauchy coordinates. It does not rule out a useful input-dependent choice
of Boolean coordinates. Extending that idea to general row tensors without
enumerating all `2^L` point states remains open.

## 3. Adaptive visibility for arbitrary nonnegative tensors

Let the two tensor entries lie in `[0,2^(p_A))` and `[0,2^(p_B))`, and
suppose neither tensor is zero. Choose a state `Y` such that

```math
A(Y)>0,\qquad A(Y)\ge A(Y\mathbin\oplus e_0).
```

A global maximum suffices, but only this one neighbor comparison is
required. For binary `A`, every positive entry meets the condition.
Simultaneously translate both tensors by `Y`. This preserves their
canonical entries and their original inner product. In the new coordinates,
`A(0)>0` and `A(0)>=A(e_0)`.

At `t=L-1`, the exact endpoint value is

```math
F(t)=A(0)\sum_XB(X)
 +b_L\bigl(A(0)-A(e_0)\bigr)
       \sum_{X_{L-1}=1}B(X).
\tag{3}
```

Both terms are nonnegative and the first is positive. Moreover,

```math
0<F(t)<2^{p_A+p_B+3L-1}.
\tag{4}
```

Thus the translated packet has a universal linear-height visibility
bound at one endpoint. Finding the required state can still take
exponentially many entry queries for an arbitrary tensor oracle. More
fundamentally, (4) holds even when the original supports are disjoint:
positive endpoint data alone do not determine the original contraction.

## 4. Visibility is not yet reconstruction

An upper visibility bound says that an odd original value cannot produce
all-zero local data at the stated precision. It does not by itself show
how to recover the original parity from arbitrary local data.

For a linear polynomial class, the usual module argument applies the
zero-data statement to a difference of two inputs. Canonical nonnegative
tensor pairs are not closed under that subtraction. A constructive
evaluator improvement needs a decoder, a suitable closure theorem, or
another exact reconstruction identity as well as the visibility estimate.

There is an exact obstruction to decoding from the outer endpoints alone.
Take `A` to be the indicator of zero, and use two designated input bits.
Let `B_even` indicate the states `00,11`, and let `B_odd` indicate
`01,10`, with all other bits fixed at zero. The two column tensors have
identical totals and identical one-bit marginals. Their original
contractions against `A` are respectively one and zero. Formula (3)
therefore gives identical endpoint values for every simultaneous Boolean
translation and every coordinate permutation. Both rank-one endpoints,
`t=L-1` and `t=-L+1`, obey the same obstruction.

The obstruction extends to higher jets. At a rank-one endpoint, only one
row coordinate of `W(t)X` is nonzero. Expand the multilinear polynomial
of the row tensor into monomials. A monomial on `k` row coordinates needs
at least `k-1` positive powers of the local variable. Consequently the
coefficient of local order `j`, as a function of the column state `X`,
has multilinear degree at most `j+1`.

For any `0<=r<=L-2`, let the column tensors instead indicate the even
and odd states on `r+2` designated bits, with the others zero. Their
moments of degree at most `r+1` agree, while their original contractions
against the zero-state row differ. Hence all endpoint jets of orders
`0,...,r` agree exactly, even after arbitrary simultaneous bit flips and
coordinate permutations. Increasing precision cannot distinguish these
data. A general endpoint decoder of this kind must use order at least
`L-1` somewhere.

The current naive height argument for general canonical inputs also
remains quadratic in `L`: products of harmonic sums can require a common
odd denominator with `O(L^2)` bits. Positivity alone does not reduce that
bound to `O(L log L)`. The finite maximum is evidence for a sharper
structural statement, not a proof of one.

## 5. Reproduction

```sh
python3 tools/check_canonical_packet_exhaustive.py \
  --output /tmp/canonical-packet-exhaustive.json
```

The wrapper uses the standard library and an optional C++17 compiler.
The stored
[canonical-packet-exhaustive.json](../evidence/canonical-packet-exhaustive.json)
records the complete finite coverage, both projection passes, the positive
control, adaptive positive-endpoint checks, the finite moment-obstruction
replays, and source hashes.
