# Rational crossings need no separate polynomial coefficient-height bound

## Result

The grid evaluator can use rational crossing factors whose numerators and
denominators have `2^O(lambda)` bits. A uniform `2^O(lambda)` **bit-time**
algorithm for their values already supplies enough height control.
The polynomial coefficient-height condition of the original integer
interface is stronger than this transfer needs.

This conclusion uses a property of this evaluator: its arithmetic circuit
has small syntactic degree in all scalar inputs. Circuit size alone would
not suffice, because repeated squaring can produce doubly exponential
integers.

**Rational transfer theorem.** Fix `1/2<=rho<1`, `delta>0`, and a constant
`C`. For every positive `lambda`, suppose there is a rational crossing
family with a nonzero rational scale, both exact factorizations, separate
channel widths at most `2lambda`, and

```math
\begin{aligned}
r_i+c_i&\le4\rho\lambda+
 \frac{C\lambda}{\log_2\log_2(\lambda+4)},\\
\log_2\max(2,|I_\lambda|)&\le
 \frac{C\lambda}{(\log_2(2\lambda))^\delta}.
\end{aligned}
\tag{1}
```

Assume one uniform algorithm enumerates the piece set and computes its
metadata, scale, and every piece or factor value in at most `2^(C lambda)` bit
operations, after increasing `C` if necessary. Rational values are
represented by binary numerators and positive binary denominators.
Then one deterministic finite-alphabet machine computes the least prime
factor and the complete prime factorization of every `n`-bit integer in

```math
\exp\!\left(O\!\left(
 \frac{n}{(\log n)^\alpha}+
 \frac{n}{(\log n)^\delta}
\right)\right),\qquad
\alpha=\log_2(1/\rho).
\tag{2}
```

In particular, the exponent is `min(alpha,delta)`.
No additional polynomial bound on factor numerators or denominators is
assumed. Their bit lengths may be exponential in `lambda`.

The arithmetic identities and the information-volume analysis are those
of [subexponential-piece-transfer.md](subexponential-piece-transfer.md).
That note works in one fixed characteristic. The argument here supplies
an exact rational bit-complexity implementation and uses the unweighted
integer divisor counts directly. No sparse-polynomial test is required.
This is a positive extension of the interface; it does not construct a
family with `rho<1/sqrt(2)` or claim a stronger factoring exponent.

## 1. Normalize the scalar before building the circuit

If the crossing sum is `d_lambda[x=v][y=u]`, divide one factor of each
piece by the nonzero rational `d_lambda`. The new sum has scale one.
The ranks and widths do not change. Rational multiplication and inversion
increase the bit length by at most a constant multiple of the original
input lengths, so the stated computation bound remains `2^O(lambda)`.

If the scale is not supplied separately, it can be computed by summing
`A_i(0,0,0,0)` over the enumerated set. The number of terms and each
output length are at most exponential in `lambda`; exact reduced rational
addition consequently still takes `2^O(lambda)` bit operations. The scale
must be proved nonzero; this normalization does not repair a zero scale.

After normalization the evaluator uses only additions and multiplications
of input-table entries and rational scalar factor values. It does not
divide by a table entry. The per-piece normalization is included among
the rational scalar inputs.

## 2. A bounded-degree circuit height lemma

Consider a circuit with at most `S>=2` nodes, counting scalar-input nodes.
Its operations are addition and multiplication. A subtraction, if used,
may be implemented with a scalar input `-1` and multiplication. Suppose
every scalar input has a numerator and positive denominator of at most
`B>=1` bits.

Define syntactic degree by assigning **degree one to every scalar input**,
including the literal constants zero and one. At an addition take the
maximum of the two degrees; at a multiplication take their sum. Let `D`
be an upper bound on the degree of every gate.

Then every gate, evaluated as a reduced rational number, has numerator
and denominator bit lengths at most

```math
1+D\bigl((S+1)B+S\bigr)
=O(D S(B+1)).
\tag{3}
```

### Proof

Temporarily replace each distinct scalar-input node by its own formal
variable. For a polynomial `P`, let `||P||_1` be the sum of absolute
values of its coefficients. Order the nodes topologically. Induction
on the index `k` of a gate of syntactic degree `d` gives

```math
\log_2\max(1,\|P\|_1)\le d k.
\tag{4}
```

For a scalar variable the left side is zero. At an addition, the triangle
inequality increases the maximum preceding logarithm by at most one;
`d>=1` absorbs that one. At a multiplication, coefficient norms multiply
and the two degrees add. This proves (4), and hence `||P||_1<=2^(DS)`.
Every monomial of `P` has total degree at most `D`.

Write the at most `S` scalar inputs as `a_i/q_i`, with `q_i>0`, and put
`Q=product_i q_i`. We have `Q<=2^(SB)`. The denominator of each evaluated
monomial divides `Q^D`, because the exponent of every scalar variable
is at most `D`. Therefore the reduced denominator of `P(a_i/q_i)`
divides `Q^D`.

Also `|a_i/q_i|<=2^B`, so

```math
|P(a_i/q_i)|\le 2^{D(S+B)}.
```

Writing the value over `Q^D` bounds its numerator by
`2^(D(SB+S+B))`. Reduction can only decrease the numerator and denominator
in absolute value. This proves (3).

### Exact implementation

After each rational addition or multiplication, reduce by the greatest
common divisor and make the denominator positive. The unreduced operands
at that one gate have only a constant multiple of the bit length in (3).
Schoolbook integer arithmetic, long division, and the Euclidean algorithm
therefore evaluate a gate in a fixed polynomial of that bound.
Scanning a sequential store to find gate operands adds a fixed polynomial
in `S`; constant-time random access is unnecessary.

Consequently the complete evaluation has cost polynomial in
`S,D,B`. The common denominator `Q^D` is only a proof device; the machine
does not construct it.

The degree convention is essential. If literal one had degree zero,
repeatedly squaring `1+1` would keep formal degree zero while its bit
length doubled at every step. Under the convention above, its syntactic
degree doubles as well.

## 3. The grid circuit has syntactic degree at most `P^6`

Use the schedule

```math
K_0=2^{x_0},\quad
g_j=2^{3j+x_0}K_j,\quad K_{j+1}=g_j^2,\qquad x_0\ge18,
```

and let `P=K_ell`. This section concerns the actual sum-product circuit,
including every rational piece-factor value as a scalar input.
It permits the early volume cutoff from the companion transfer.

At a step from `K_j=g^2` vertices to `K=K_(j-1)` vertices, the row or
column network has the following explicit structure:

- A node value is a product of at most `g` parent-table values and
  `g^2` indicator values.
- An edge value or one of its half-edge factors is a product of at most
  `g+g^2` indicators and scalar piece-factor values.
- Each term in a grouped-table entry is a product of at most `g^2`
  such network values.

These are the products in the source's Proposition 3, equation (35),
and Lemma 8. With scalar leaves assigned degree one, a grouped-table
entry is therefore a polynomial of degree at most

```math
b_j=4g^4=4K_j^2
\tag{5}
```

in parent-table entries and the scalar factors introduced at that step.
Indeed `g+g^2+1<=2g^2` for `g>=2`, and multiplying at most `g^2`
such values with one additional accumulator gives degree at most
`2g^4+1<=4g^4`. This includes product accumulators and padding
indicators. The potentially large **number** of local settings does not
increase the degree: the entry is a sum of those products.

The primary source is
[the pinned paper, Proposition 3 and Appendix A, Lemma 8](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
Its bound `B1=gB+g^2(2lambda)^c` records the same product structure;
we use the explicit products rather than that coefficient-height bound.

At a terminal call at level `k`, brute force sums products of `K_k`
table entries. Compile a nonempty product starting with its first
entry, so no extra multiplication by an initial one is needed.
Its syntactic degree in those entries is `K_k`.

Substituting a polynomial of degree at most `b_j` for a table-input
leaf multiplies the current degree bound by at most `b_j`.
New factor leaves from deeper levels have degree one and satisfy the
same bound. Summing the circuits for distinct piece assignments takes
the maximum degree, not their sum. Thus every branch, and the full
circuit, has degree at most

```math
K_k\prod_{j=k+1}^{\ell}(4K_j^2)
\le P\,4^\ell\prod_{j=1}^{\ell}K_j^2.
\tag{6}
```

Since `K_j>=K_(j-1)^2`,

```math
\sum_{j=1}^{\ell}\log_2K_j\le2\log_2P,
\qquad
\prod_{j=1}^{\ell}K_j^2\le P^4.
```

Also `log2 P>=18*2^ell>=2ell`, hence `4^ell<=P`.
Equation (6) proves

```math
D\le P^6.
\tag{7}
```

This bounds all intermediate gates as well. A full piece entry needed
on an internal edge is computed by its uniform scalar algorithm and
counted as a scalar input. Its generation cost and output bit length
are charged in Section 4. No cancellation of rational polynomials is
used to obtain the degree bound.

## 4. Exponentially long scalar coefficients still fit the runtime

Choose the largest level with `P^40<=n` and use the exact cutoff
`K_j^2 lambda_j<=ceil(sqrt(n))` from the companion transfer.
Its proof gives a uniform envelope, for every visited path,

```math
V_j=K_j^2\lambda_j
\le C_0 n\rho^{\ell-j}.
\tag{8}
```

Since `K_j>=2^(18(j+1))` and `rho^(-j)<=2^j`,

```math
\lambda_j
\le C_0n\rho^\ell\frac{\rho^{-j}}{K_j^2}
=O(n\rho^\ell)
\tag{9}
```

uniformly over all calls, including those that immediately stop.
A value generated in `2^(C lambda_j)` bit operations has at most that
many output bits. Therefore every rational scalar leaf has at most

```math
B\le\exp(O(n\rho^\ell))
\tag{10}
```

bits. The written top-table entries have only `O(n/P)` bits and also
satisfy this bound after increasing a fixed constant.

The circuit-size and circuit-generation analysis uses no field property.
It counts the same finite sums and products over `Q`. The factor
algorithms' complete bit costs `2^(C lambda_j)` are charged in the
local-query bound, just as their finite-field computation costs were.
Hence the companion analysis supplies

```math
S\le\exp\!\left(O\!\left(
 n\rho^\ell+\frac n{(\log n)^\delta}
 +\sqrt n+\frac nP
\right)\right).
\tag{11}
```

All leaves are explicitly generated within this bound. Equivalently,
their maximum individual encoding length is already at most the total
circuit-construction budget; (9)--(10) explain directly why allowing
exponential height in a local label length is harmless here.

Insert (7), (10), and (11) into the height lemma. Since `P<=n^(1/40)`,
`D` is polynomial in `n`. The logarithm of the intermediate bit length
and of the complete exact rational evaluation cost is consequently

```math
O\!\left(n\rho^\ell+\frac n{(\log n)^\delta}
 +\sqrt n+\frac nP\right).
```

Finally `rho^ell=Theta((log n)^(-alpha))` and `n/P<n^(399/400)`.
This proves the runtime in (2).

The top clique computes the integer number of divisors with a prescribed
prefix. Subtract the known contribution of divisor one when appropriate.
Exact rational evaluation preserves this integer, so testing its numerator
for zero gives the prefix decision. The usual greedy prefix search
returns the least divisor greater than one, including `N` when it is
prime. Repeating after division gives the complete factorization with
only polynomially many queries. There is no reconstruction from numerical
approximations or from modular counts.

## 5. Corollary: exact rational rank bounds can be made effective

Suppose instead that the pieces themselves are given by uniform rational
entry algorithms of cost `2^O(lambda)`, with the same piece-entropy bound
and crossing identity. Let `R_i,C_i` be the exact ranks of their two
flattenings. If

```math
\log_2\max(1,R_i)+\log_2\max(1,C_i)
\le4\rho\lambda+
 \frac{C\lambda}{\log_2\log_2(\lambda+4)},
\tag{12}
```

then explicit factor algorithms meeting the theorem can be constructed
by exact Gaussian elimination. No separate closed-form factorization or
polynomial coefficient-height certificate is necessary.

Each flattening is an `M` by `M` rational matrix, with `M=2^(2lambda)`.
Writing every entry takes `2^O(lambda)` total bit operations. If their
numerators and denominators have at most `B=2^O(lambda)` bits, clear all
entry denominators by their product. The resulting integer matrix has
entry bit length at most `O(M^2 B)=2^O(lambda)`.

Exact elimination finds independent columns and pivot rows in time
polynomial in the matrix dimension and this entry bit length. To see
why coefficient growth is controlled, each Schur-complement entry is a
ratio of minors. The elementary determinant bound
`|det H|<=r! max|H_ab|^r` gives polynomial bit length in those same two
parameters. Reduced rational arithmetic therefore has polynomial bit
cost throughout elimination.

Take the pivot columns as one factor and express every column in them
using the invertible pivot-row minor. Cramer's formula gives the same
minor-based height bound for the second factor. The two flattenings can
be handled independently. Zero matrices use a single zero channel.
Pad a rank-`r` factorization to `2^ceil(log2 r)` channels; the two ceiling
costs add at most two bits and are absorbed into the slack in (12).

All dimensions, entries, and exact arithmetic bounds above are
`2^O(lambda)`, so the complete factor generation also has that cost.
The coefficients may have exponentially many bits, which Sections 2--4
now permit. This turns a proved scalable rank profile into a complete
ordinary algorithm under the stated uniform entry and piece-count
hypotheses. A finite rank experiment still does not establish (12).

## 6. Reproduction and limits

Run

```sh
python3 tools/check_rational_height_transfer.py
```

The checker verifies the height lemma on exact rational circuits,
including literal-one leaves, cancellation, shared gates, and nontrivial
denominators. It checks the degree product for all stopping levels of
finite schedules and exercises exact rank factorizations whose entries
are rational. Repeated squaring is a negative control for a size-only
height bound. These finite cases do not replace the all-input arguments.

The recorded output is
[rational-height-transfer.json](../evidence/rational-height-transfer.json).
The grid and top constructor are prior work from the pinned Anthropic and
factoring developments. This degree/height analysis and checker used
substantial OpenAI assistant support. No new Lean formalization was run.
