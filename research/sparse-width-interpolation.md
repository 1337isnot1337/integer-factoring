# Crossing families from a sparse sequence of widths

## Result and hypotheses

A uniform rank theorem need not hold at every binary label width. It is
enough to construct good families at an unbounded, uniformly generated
sequence of widths with at most polynomial gaps. Tensor products fill the
gaps with subexponentially many pieces, within the effective construction
budget.

Fix a finite field, `1/2<=rho<1`, and positive integer widths
`1=a_0<a_1<...`. Assume that for fixed constants `K>=1,s>=1`,

```math
a_{j+1}<=K a_j^s.
\tag{1}
```

At width `a_j`, assume an effective scale-one crossing family with
`M_j<=(2a_j)^C` pieces, separate binary channel widths at most `2a_j`,
and each total width at most

```math
4rho a_j+C log2(2a_j).
\tag{2}
```

A single uniform algorithm must list the available widths up to any `L`
and construct the families and all factor values for such widths in
`2^O(L)` steps. The sequence and its families cannot be uncomputable advice.
An explicitly computed sequence such as `1,2,4,16,256,...` meets the
sequence-generation requirement; merely postulating some good widths does
not.

**Interpolation theorem.** For every positive width `L`, there is an
effective scale-one family with separate widths at most `2L` and:

- If `s=1`, piece entropy `log2 M=O(log^2(2L))` and total width
  `4rho L+O(log^2(2L))`.
- If `s>1`, put `beta=1-1/s<1`. Piece entropy is
  `O(L^beta log(2L))` and total width is
  `4rho L+O(L^beta log(2L))`.

Construction, index enumeration, and factor evaluation all cost `2^O(L)`.
Combining this with the
[subexponential-piece transfer at PR #23 / 6dd2901](https://github.com/1337isnot1337/integer-factoring/blob/6dd290118668914d1f58a1d23d302963955d1dcb/research/subexponential-piece-transfer.md)
gives the **same** factoring exponent `alpha=log2(1/rho)` as a family at
every width. This conclusion uses that separate ordinary transfer theorem;
it does not follow from the polynomial-piece interface alone.

This is a conditional bridge. It provides neither the rank bound (2) at
an unbounded sequence nor a factoring exponent above one half. Finite
certificates at a bounded list of widths do not meet its hypotheses.
The proof and checker used substantial OpenAI assistant support. No new
Lean formalization is included.

## 1. Greedy decomposition of a label width

Start with remainder `m=L`. Choose the largest available `a_j<=m`, append
it to the block list, and replace `m` by `m-a_j`. The process terminates
because `a_0=1`. The chosen block sizes `b_1,...,b_k` sum to `L`.

By maximality, `a_{j+1}>m`, so (1) implies

```math
b_t>(m/K)^{1/s}.
\tag{3}
```

For `s=1`, if `K>1` the next remainder is less than
`(1-1/K)m`, proving `k=O_K(log(2L))`. An unbounded strictly increasing
sequence with `a_0=1` cannot satisfy (1) with `K=s=1`, so that degenerate
case supplies no counterexample.

For `s>1`, let `gamma=1/s` and `beta=1-gamma`. Concavity gives, for
`0<b<=m`,

```math
m^beta-(m-b)^beta>=beta m^{-gamma} b.
```

Using (3), every step decreases this potential by more than
`beta K^{-gamma}`. Telescoping down to zero proves

```math
k<=K^gamma L^beta/beta.
\tag{4}
```

The argument allows repeated use of an available width. It does not
require a unique digit representation or a bound on the number of
repetitions of each width.

## 2. Tensor the pieces and their explicit factors

Split each `L`-bit label into the fixed blocks of sizes `b_t`. A new piece
index is a tuple `(i_1,...,i_k)` of source piece indices. Define

```math
A_i(x,y,u,v)=\prod_{t=1}^k A_{i_t}^{(b_t)}(x_t,y_t,u_t,v_t).
\tag{5}
```

Distributivity and the source crossing sums give

```math
\sum_i A_i=\prod_t [x_t=v_t][y_t=u_t]=[x=v][y=u].
```

Tensor the row factors and concatenate their binary channel strings;
tensor the column factors in the other arrangement in the same way.
Finite distributivity proves both factorizations. Their widths add, so
they are separately at most `sum_t 2b_t=2L`, and their sum is at most

```math
4rho L+C\sum_t log2(2b_t)
<=4rho L+C k log2(2L).
\tag{6}
```

The piece count is exactly the product of the source piece counts, hence

```math
log2 M=\sum_t log2 M_{b_t}<=C k log2(2L).
\tag{7}
```

Equations (4), (6), and (7) prove the claimed quantitative bounds. No
rank subadditivity estimate or unproved cancellation is used.

## 3. Uniform construction and the factoring corollary

Generate the available widths up to `L`, greedily choose the blocks, and
run the given source constructors. Repeated block sizes may share their
already generated source tables. At most `L` block occurrences are used,
and their total width is `L`; their complete construction fits `2^O(L)`.

The index tuple has `O(k log(2L))=o(L)` bits. Enumerating all its values
therefore costs `2^o(L)` times a polynomial overhead. A factor value is a
product of at most `L` source values. Even materializing every factor
value costs `2^O(L)`: there are at most `2^(4L)` input-and-channel choices
for any one factor and `2^o(L)` pieces. Arithmetic is over a fixed finite
field. Loop counters and binary channel descriptions have `O(L)` bits.
This is a single finite-alphabet uniform algorithm.

For `s=1`, both overheads are `O(log^2(2L))`. For `s>1`, they are
`O(L^beta log(2L))`, where `beta<1`. In either case, for every fixed
`delta>0` and all sufficiently large `L`, the piece entropy is bounded by

```math
C_delta L/(log2(2L))^delta,
```

and the width slack is bounded by `C L/log2 log2(L+4)`. These are exactly
the larger-piece transfer's hypotheses. Choose `delta=1`; because
`alpha=log2(1/rho)<=1`, its exponent `min(alpha,delta,1)` equals `alpha`.
Finite exceptions are handled by enlarging fixed constants and the
explicit source algorithms, as in that theorem.

Thus a future theorem at widths `2^j`, or even `2^(2^j)` with a finite
initial width-one block, can suffice without an all-width rank theorem.
For example, the latter sequence has `s=2,K=1` after its initial step,
which is covered by increasing `K`. This does not amplify a single
finite seed into a stronger rate: that would give bounded block sizes,
linear piece entropy, and would fail the unbounded-sequence hypothesis.

## 4. Finite replay

The all-size result follows from the proof above. The finite checker
exercises greedy decompositions on geometric and polynomial-gap sequences
and verifies tensor crossing sums and multiplicativity of both matrix
ranks on the stored one-bit and two-bit Pfaffian witnesses. It materializes
all 18 product pieces at width three and uses an independently written
Kronecker assembly to compare both arrangements.

```sh
python3 tools/check_sparse_width_interpolation.py
make verify
```
