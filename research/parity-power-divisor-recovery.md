# Factoring by characteristic-two parity and power-divisor queries

## Result and scope

There is a deterministic reduction from integer factorization to parity
counts of multiplication systems with total operand length `O(n)`, where
`n` is the input bit length. It uses `O(n^2)` counting calls for complete
factorization. Every call has only zero-or-one weights and can be evaluated
over the fixed two-element field, or any extension field of characteristic
two. It does not need sparse-polynomial identity tests, a growing coefficient
ring, integer-count reconstruction, or a primality oracle.

The new part is the recovery reduction and its composition of multiplication
networks. The individual multiplication network and its unique-extension
property are due to Isaac Rudich and Louis-Martin Rousseau. We use their
revision
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1),
particularly `FactoringAux/MulNet.lean` and `FactoringAux/MulCount.lean`.

This reduction supplies no new crossing tensor and does not improve the
unconditional factoring exponent `1/2`. It removes a recovery overhead for
characteristic-two evaluators. The earlier
[rough-input sparse recovery](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/rough-input-transfer.md) has
time `exp(O(log A(n)+sqrt(n)))` after rough-input preprocessing. Here a
parity evaluator for the compiled systems with time `A(n)` yields factoring
in `poly(n,A(Cn))` time, for an absolute constant `C`. The evaluator must
cover the compiled systems, not only the original fixed-product query.

The direct recovery stores ordinary binary factors and prefixes. No
sublinear-workspace theorem is asserted by this note.

## 1. Binary digits of prime exponents appear in power-divisor parities

Write

```math
N=\prod_p p^{a_p},\qquad n=\operatorname{bitlength}(N),\qquad N\ge2.
```

For an integer `j>=1`, put

```math
R_j=\prod_p p^{\lfloor a_p/2^j\rfloor},\qquad
\mathcal D_j(N)=\{d\ge2:d^{2^j}\mid N\}.
```

Prime factorization occurs only in this proof. The algorithm does not
compute `R_j` or know the exponents `a_p`.

**Power-divisor parity lemma.** The integers `d` satisfying
`d^(2^j)|N` are exactly the divisors of `R_j`. Consequently,

```math
|\mathcal D_j(N)|
=\prod_p(\lfloor a_p/2^j\rfloor+1)-1.                 \tag{1}
```

This number is odd if and only if the `j`th binary digit of at least
one exponent `a_p` is one.

**Proof.** The divisibility condition is equivalent, separately at every
prime, to `v_p(d)<=floor(a_p/2^j)`. This proves the set identity and (1).
The product in (1) is even if and only if some
`floor(a_p/2^j)+1` is even, which is equivalent to
`floor(a_p/2^j)` being odd. The latter condition is exactly the stated
binary-digit condition. Subtracting one reverses the parity.

Every exponent is at most `n-1`, since
`2^a_p<=p^a_p<=N<2^n`. It therefore suffices to query

```math
1\le j\le J:=\lfloor\log_2(n-1)\rfloor.              \tag{2}
```

When `n=2`, this range is empty. All the parities in (2) vanish if
and only if every prime exponent is zero or one, which is exactly
squarefreeness. If any query is odd, every member of its nonempty
set is a proper divisor of `N`: `j>=1` gives `d^2<=N`, while
`d>=2` gives `d<N`.

## 2. An odd set supports exact prefix recovery

Let a finite set of positive integers below `2^b` have odd cardinality.
Suppose an oracle returns the parity of the members in any binary prefix.
An initially selected prefix with odd count can be refined one bit at a
time. Query its zero child. If that child is odd, retain it. Otherwise
retain the one child, which must be odd because the two child parities
sum to the odd parent parity. After at most `b` calls there is a single
integer in the retained prefix, and it belongs to the set.

The procedure returns a member, not necessarily its least member. For
example, ordinary parity search in the divisor set `{2,3,5}` skips the
even-count prefix containing `2,3` and can return `5`. Returning a proper
factor and recursively factoring both resulting factors is sufficient.

Apply this procedure to an odd set `D_j(N)`. Its queries are

```math
\#\{d\ge2:d^{2^j}\mid N,\ \lfloor d/2^k\rfloor=x\}
\pmod2.                                                   \tag{3}
```

The complete power equation has the unique complementary variable
`e=N/d^(2^j)`. A count that permits `d=1` is changed to (3) by subtracting
`[floor(1/2^k)=x]`. Subtraction is addition in characteristic two.
No actual witness is deleted except the known value one.

## 3. The squarefree case requires only ordinary divisor intervals

Suppose every query in (2) is even. Section 1 has then proved that `N`
is squarefree. Let `r` denote its number of prime factors and set
`B=floor(sqrt(N))`. Because `N>=2` is squarefree, it is not a square.
Complementation `d -> N/d` pairs the divisors with exactly one member
below `sqrt(N)`. Therefore

```math
\#\{d:1<d\le B,\ d\mid N\}
=2^{r-1}-1.                                              \tag{4}
```

For a prime, `r=1` and the set is empty. For a squarefree composite,
`r>=2` and (4) is odd.

One can query this interval using the original divisor-prefix interface.
Write `[0,B+1)` as the standard disjoint union of at most `n` binary
prefix blocks: for each one bit of `B+1`, keep the preceding high bits,
replace that bit by zero, and leave lower bits free. For a block that
contains one, subtract its known divisor-one contribution. Zero is never
a divisor of the positive input.

Scan those blocks for an odd parity. If one is found, Section 2 recovers
a member of that block using only further ordinary prefix counts. It is
at most `B` and greater than one, and hence is a proper factor. If every
block has even parity, (4) proves that the squarefree input is prime.
No outside primality algorithm is used.

Computing `B` and the prefix descriptors takes polynomial bit time.
For example, the usual integer-square-root algorithm uses only integer
comparisons and multiplication. Prefix recovery uses at most `n`
additional parity calls after the scan.

## 4. The complete deterministic factorization algorithm

On a positive input `N>=2`:

1. Compute `n` and the integer bound `J` in (2).
2. Query (1) modulo two for `j=1,...,J`. At the first odd result,
   use (3) and prefix recovery to obtain a proper divisor `d`.
3. If no such result occurs, perform the interval-prefix scan from
   Section 3. An odd block yields a proper divisor `d`. If all its
   blocks are even, declare `N` prime.
4. Whenever a proper divisor is found, factor `d` and `N/d` by the
   same procedure. Collect or stream the resulting prime factors;
   their least value is the required smallest prime factor.

Every recursive call has a smaller positive input. Its termination and
prime declarations follow from Sections 1--3. A complete recursion tree
has at most `n` prime leaves, counted with multiplicity, since every
prime is at least two. It therefore has at most `2n-1` nodes. A node of
bit length at most `n` makes `O(n)` parity calls. The whole procedure
uses `O(n^2)` calls, plus polynomial-time arithmetic and bookkeeping.

It is also possible to keep only the smallest prime found during this
traversal. The algorithm never assumes that the first recovered factor
is prime or is the smallest factor.

## 5. Compile each power-divisor query into short multiplications

Fix `j` from (2). Introduce

```math
x_0=d,\qquad x_i=x_{i-1}^2\ (1\le i\le j),\qquad x_j e=N.
                                                               \tag{5}
```

For every permitted `d`, the intermediate values and `e` are unique.
Thus this system has exactly one solution per power divisor. In
particular, counting system solutions in characteristic two gives the
desired parity, without an extra multiplicity from the squaring chain.

Use the bounds

```math
b_i=\left\lceil\frac n{2^{j-i}}\right\rceil
\quad(0\le i\le j).                                      \tag{6}
```

Every valid solution has `0<x_i<2^b_i` and `0<e<2^n`.
The `i`th squaring network uses a scale sufficient for `b_i` bits,
with both of its inputs restricted to `b_(i-1)` bits. Its output is
restricted to `b_i` bits. The last multiplication uses scale `n`.
The sum of these requested scales is

```math
n+\sum_{i=1}^j b_i<3n+j\le4n.                            \tag{7}
```

This is the relevant length budget. Using `n` bits for every squaring
would instead introduce an unnecessary `log n` factor.

The first squaring has two input copies of `d`. Equality between them
uses exactly `b_0` bit wires. For each intermediate output `x_i` with
`1<=i<j`, copy its `b_i` bits to both inputs of the next square. Copy
the last output's `b_j=n` bits to the first input of the final product.
The number of added equality-wire bits is

```math
b_0+2\sum_{i=1}^{j-1}b_i+b_j=O(n).                       \tag{8}
```

Both source and target layouts may be different. Match their global
bit indices; concatenate the bits joining the same pair of gates in
one fixed order at both ends. This bit refinement works when the two
networks have different digit widths. It introduces no new freely
chosen settings and creates at most two outgoing copies of any
intermediate bit.

High input bits outside the bound in (6) are locally fixed to zero.
The first input's allowed prefix is imposed at its existing input
gates. The known `d=1` setting is subtracted after evaluation as in
Section 2, so no nonlocal test for `d>1` is required in the network.

## 6. Expose an output port in the pinned multiplication network

The source network has input gates `faP` and `faQ`, transform gates
`mid` and `ib`, and a carry chain. At scale `s`, write

```math
S=2^s,\quad L=2S^2,\quad A=2S,\quad h=S,\quad
f=4L+1,\quad M=2^{4L}+1.
```

Each input gate stores `h` base-`2^L` digits. A carry gate has a
coefficient `c` and incoming carry `u`; its product-output digit is

```math
(u+c)\bmod2^L,                                          \tag{9}
```

and its outgoing carry is `floor((u+c)/2^L)`.
In the original `mW`, the digit (9) is compared with the corresponding
fixed digit of `N`. For an intermediate variable output, remove that
comparison and expose the bits of (9) as deterministic port messages.
Keep the first carry equal to zero and the final carry equal to zero.
Also set output bits above the chosen bound `b_i` to zero by local
carry-gate predicates.

The source `good_pinned` proof already forces every `mid`, `ib`, and
carry setting from the two input vectors. Its carry induction uses
only the first carry being zero and the equality of adjacent carry
messages; it does not use the fixed output-digit comparison. That
unique-extension proof therefore remains valid after the port change.
The source `fC_eq_cv` and `val_fC` identify the received coefficients
with ordinary digit-convolution coefficients. The carry recurrence
then gives exactly the binary product at the exposed ports. Conversely,
the unique settings built by `asm` satisfy these port identities for
every admitted input pair.

Thus replacing a fixed product by connected output ports is an exact
relation construction. It does not reinterpret an equality of total
counts as a statement about individual settings. The bit equalities
in (8) impose exactly (5), with the unique witness count already proved.

### Width and setting budgets

Choose the least scale with `L^2>=m` for a requested length `m`.
Then `L^2=O(m+1)`. A source multiplication gadget has

- total binary setting-description length `O(m+1)`;
- total ordered wire width `O(m+1)`;
- maximum single-gate setting length plus incident wire width
  `O((m+1)^(3/4))`;
- polynomial-time local message and weight evaluation.

For the setting bound, the two input layers use `2AhL=2L^2` bits.
The `mid` layer uses at most `2A^2 f`, the `ib` layer at most `A^2 f`,
and the carry layer at most `2A^2 f`. Since `A^2=2L` and
`f=4L+1`, their sum is `O(L^2)`. A noncarry gate has at most
`O(Af)=O(S^3)` incident original bits, and a carry gate has
`O(f)=O(S^2)`. The same estimates bound individual setting lengths.
There are `O(A^2)` original edges of width `f`.

Adding (8) preserves these estimates after summation. An input gate
receives at most its `hL=O(S^3)` stored bits per operand occurrence.
A carry gate exports at most two copies of its `L=O(S^2)` output
bits. No global variable node receives all copies. Equations (7)--(8)
therefore give for the whole compiled query

```math
\sum_v\left(\sigma_v+\sum_w c_{vw}\right)=O(n),\qquad
\max_v\left(\sigma_v+\sum_w c_{vw}\right)=O(n^{3/4}),     \tag{10}
```

where `sigma_v` bounds the bits describing a local setting, and
`c_vw` is the width of its equality message to the neighboring gate.
All descriptions, restrictions, and port functions are uniform and
computable in polynomial time from the query.

The companion [general weighted-network constructor](weighted-network-top-compiler.md)
can now be applied to (10). In particular, its greedy preliminary clustering and the
existing grouping lemma give `P^2 lambda=O(n)`, local setting budgets
`O(n/P)`, and written-table construction `2^O(n/P)` for `P^40<=n`.
This is the same scale needed by the existing clique evaluators.

## 7. Consequence for a characteristic-two evaluator

Suppose a uniform algorithm exactly evaluates the compiled zero-or-one
networks, over a fixed characteristic-two field, in `A(n)` bit steps
when their scale budget is at most `n`. Include network construction,
indexing, factor generation, and field arithmetic in this hypothesis.
Sections 1--6 give

```math
T_{\rm factor}(n)\le\operatorname{poly}(n)\,(1+A(Cn))    \tag{11}
```

for an absolute constant `C`, after replacing `A` by a monotone upper
envelope if necessary. If instead `A` bounds arithmetic circuit size
and generation time, ordinary sequential-tape evaluation gives the
slightly more general statement `T_factor(n)=poly(n,A(Cn))`.

For example, an evaluator with

```math
A(n)=\exp(O(n/(\log(n+2))^\alpha))
```

preserves every fixed `alpha>0`. An evaluator with
`A(n)=exp(O(n^beta))` preserves every fixed `beta>0`; there is no
additional square-root floor from the recovery reduction. These are
conditional statements about an evaluator meeting the stated interface.
They supply no such improved evaluator and establish no new unconditional
factoring score.

The parity test is specific to characteristic two. A field of odd
characteristic does not satisfy the binary-exponent test in Section 1.
No characteristic-independent extension is claimed here.

## 8. Reproduction and limitations

Run the standalone checker:

```sh
python3 tools/check_parity_power_divisor_recovery.py
```

It checks the binary prime-exponent identity against direct integer
divisibility, squarefree certification, exact interval-prefix coverage,
parity witness recovery, and complete factorization against independent
trial division. It also checks the geometric scale budgets and the
variable-output carry/bit-port identities, including unequal layouts and
the two-copy square constraint. Negative controls omit the divisor-one
subtraction, operand-copy equalities, or an output link; each omitted
condition has an explicit failing example.

These finite computations do not implement the large-input clique
evaluator or establish its asymptotic bound. The ordinary proof above and
the separately supplied weighted-network constructor establish the
reduction and its costs. No Lean formalization is claimed. This work used
substantial OpenAI assistant support and retains the source attribution.
