# Factoring from scale-one crossings over a fixed finite field

## Result and scope

An effective crossing family over one fixed finite field can suffice for
deterministic integer factoring. It need not reconstruct an integer count
by combining counts in different characteristics.

The transfer below uses two separate facts. First, scale-one crossing
identities give an arithmetic circuit over every commutative algebra over
that field. Second, a divisor polynomial has few terms. Its reductions in
a controlled family of cyclic polynomial rings detect whether a permitted
divisor exists.

**Modular transfer theorem.** Fix a prime `p` and an effective family
with parameter `1/2<=rho<1`, as defined in Section 1. Then one deterministic
finite-alphabet machine computes the smallest prime factor of every
`n`-bit integer `N>=2` in time

```math
T(n)\le C\exp\!\left(
\frac{Cn}{(\log n)^\alpha}\right),
\qquad \alpha=\log_2(1/\rho),\qquad n\ge2.
\tag{1}
```

The field and family are fixed parts of the algorithm. The constant is
independent of the input and its length. The endpoint `rho=1/2` gives
`alpha=1`. This is a conditional transfer theorem: it does not construct
a family with `rho<1/sqrt(2)`.

The argument is an ordinary proof. It builds on the grid identities and
operation count in Claude's Appendix A, and the multiplication network
of Isaac Rudich and Louis-Martin Rousseau. It does not reinterpret the
source's integer program as a program on polynomial rings. The ring
evaluator described here is a separate, explicit sum-product algorithm.

## 1. The finite-field hypothesis

Fix an integer `c>=1`. For every positive label length `lambda`, require:

1. A set `I_lambda` with at most `(2lambda)^c` pieces. Each piece
   `A_i` is an `F_p`-valued function of four `lambda`-bit labels.
2. The exact identity

   ```math
   \sum_{i\in I_\lambda}A_i(x,y,u,v)=[x=v][y=u].
   \tag{2}
   ```

3. Two explicit factorizations over `F_p`, in the source's row and
   column arrangements, with respective bit widths `r_i,c_i` satisfying

   ```math
   r_i,c_i\le2\lambda,\qquad
   r_i+c_i\le4\rho\lambda+c\log_2(2\lambda).
   \tag{3}
   ```

4. Uniform algorithms for the index set, widths, piece values, and
   factor values, with cost at most `2^(c lambda)`. Increasing the
   fixed `c` to cover constant overhead is permitted.

The factorization functions take binary channel indices. Padding a finite
channel set by zero entries supplies that convention. A longer width
can be replaced by the direct factorization through the two input labels,
so the separate `2lambda` caps impose no additional rank restriction.

These are hypotheses about explicit algorithms, not only about matrix
ranks. A family over a field that changes with the input is outside this
statement unless a separate uniform cost argument is supplied.

Equation (2) can also start with a nonzero scalar `d_lambda` in `F_p`.
Its inverse in the fixed field normalizes one factor of each piece.
The scale must be a unit. A scalar that becomes zero modulo `p` cannot
be divided out and does not satisfy this hypothesis.

## 2. Base change preserves every crossing identity

Let `R` be any commutative unital `F_p`-algebra. Embed the coefficients
of all pieces and factors into `R`. Finite sums and products commute
with this embedding. Thus (2) and both factorizations hold in `R`.

Now allow the input clique tables to have arbitrary values in `R`.
The grid expansion multiplies (2) at its crossing sites and distributes
the products. For each piece-index assignment, the row and column
factorizations give two presentations of the same term. All these
statements use only finite distributivity and commutativity.

The grouping construction also uses only sums and products of table
entries. Its combinatorial choice depends on channel widths and group
sizes. It does not inspect, compare, or divide table values.

Consequently, the recursion works over `R` even when `R` has zero
divisors. There is no division by an input table value. The integer
algorithm's final division by `d_lambda` to a power is absent because
the scale in (2) is one.

This distinction is essential for the rings used below. For example,
in `F_2[X]/(X^3-1)`, both `1+X` and `1+X+X^2` are nonzero, while
their product is zero.

## 3. Count circuit operations independently of ring representation

Choose an integer `x0>=18` with `2^x0>=c+5`. Use the source schedule

```math
K_0=2^{x_0},\qquad
g_j=2^{3j+x_0}K_j,\qquad K_{j+1}=g_j^2,
\tag{4}
```

and, for the analysis only,

```math
e_0=\tfrac12,\qquad
\zeta_j=2^{-j-5},\qquad
e_{j+1}=(\rho+3\zeta_j)e_j+3/K_j.
\tag{5}
```

**Circuit bound.** There is a fixed constant `C0` such that the algorithm
for a `K_ell`-vertex, `lambda`-bit clique uses at most

```math
S_\ell(\lambda)
=2^{e_\ell K_\ell^2\lambda}
  \bigl(K_\ell(\lambda+1)\bigr)^{C_0K_\ell^2}
\tag{6}
```

ring operations and input-table queries, with structural computation
bounded by a fixed polynomial in `S_ell(lambda)`. This computation includes
piece evaluation over the fixed field. The bound is uniform over `R`.
One deterministic algorithm generates the circuit within that same fixed
polynomial bound.

Here a circuit is a finite sequence of additions and multiplications,
with input-table entries and embedded field elements as its inputs.
No bound on the binary length of an `R` element is asserted by (6).

### The operation-count proof

Use the explicit construction in Appendix A, Theorem A(b).
At level zero, enumerate the `2^(lambda K0(K0-1)/2)` edge-label
assignments. Each term takes polynomially many table queries and ring
operations. This gives (6) at `ell=0`.

At the next level, write `K=K_j`, `g=g_j`, and `N=g^4 lambda`.
The parent has `g^2` vertices. Enumerate the piece-index assignments
to its at most `g^4/4` crossing sites. Their number is at most
`(2lambda)^(c g^4/4)`.

For each assignment, compute the row and column totals from their
integer channel widths. Choose the smaller total. Its bound is the
source's display (41):

```math
U\le(\rho+1/g)N+\frac c4g^4\log_2(2\lambda).
\tag{7}
```

Find a grouping by enumerating all `K^g` colorings and testing its
size and width inequalities. This uses only integer data. The source's
case (ii) grouping proof applies unchanged and gives a child width `L`
with

```math
K^2L\le(\rho+3\zeta_j)N
       +\frac c4g^4\log_2(2\lambda)+K^2.
\tag{8}
```

To answer a query to a grouped table, enumerate at most `2^(2N/K)`
local settings. Multiply and add their input-table and factor values.
Factor generation costs `2^(c lambda)` and polynomial structural
overhead. None of these loops depends on the ring values.

Put `Xi=g^2(lambda+1)`. The index calculations for a term and the tests
for a grouping have fixed polynomial cost in `Xi`. The crude bound
`L<=N+1` gives `K(L+1)<=Xi^(c+5)`. Since `e_j<=1`, (8) and the
induction hypothesis bound one child circuit by

```math
T\le2^{e_j(\rho+3\zeta_j)N}
\Xi^{c g^4/4+K^2+(c+5)C_0K^2}.
```

Inlining a grouped-table query costs at most `2^(2N/K+c lambda)`
times a fixed power of `Xi`. The remaining grouping search contributes
`K^g<=Xi^g`. Sum over the at most `Xi^(c g^4/4)` patterns.
The inequality `c lambda<=N/K` then gives the bound

```math
2^{e_{j+1}N}
\Xi^{c g^4/2+g+K^2+C_4+(c+5)C_0K^2}.
\tag{9}
```

The constants `C4` and the base-case constant depend only on `c`.
For the source's fixed choice of `C0`,
`(c+5)K^2<=g^4/2` absorbs the last term, and the exponent in (9)
is at most `C0 g^4`. This proves (6).

The integer source also tracks a height parameter `B`. The derivation
above counts circuit nodes directly and needs no scalar-height inequality.
Group selection and orientation selection use only the integer widths.
Loop counters, circuit indices, and piece queries have length bounded
by a fixed polynomial in `Xi` and the circuit size. Their construction
has fixed polynomial overhead. No ring value controls an indexing choice.
Bounds on child integer heights are unnecessary for this circuit.

The relevant primary source is
[`eth/paper.pdf`, Appendix A, displays (41)--(44)](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
The numerical induction is also stated explicitly in
[`GeneralPieces/Rec/Count.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/Count.lean).
The expansion and grouping structure appears in
[`GeneralPieces/Rec/Level.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/Level.lean).
These sources prove their stated integer results. The base-change and
circuit interpretation in this section are the new ordinary argument.

## 4. A weighted multiplication network gives the divisor polynomial

For a divisor-prefix query, let

```math
F_{N,k,x}(X)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}X^d.
\tag{10}
```

The reference's multiplication network produces (10) with monomial
weights. Use the source's notation `L=pL(s)`, `A=pA(s)`, `h=ph(s)`.
An input gate `faP(r)` stores digits `u(i)` in positions `r+A i`.
Multiply its original zero-or-one weight by

```math
X^{\,\sum_{i=0}^{h-1}u(i)2^{L(r+Ai)}}.
\tag{11}
```

The product of these additional weights is exactly `X^d`, where `d`
is the integer represented by the first factor's input digits.
All messages, equality tests, settings, and grouping choices are retained.

This assertion does not rely only on an equality of unweighted counts.
The pinned source's `good_card` proves an explicit bijection between
valid complete settings and permitted divisors. Its forward map reads
the `faP` digits; its inverse is `asm(d,N/d)`. Therefore each permitted
divisor contributes its monomial exactly once.

The source is already parameterized by an arbitrary commutative ring:
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
defines the gate weights and digit positions, and
[`MulCount.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean)
contains `good_pinned`, `good_card`, and `mulNet_value`.
Weighting the stated bijection gives (10) over every commutative ring.

To generate a circuit independently of the target ring, give bit position
`j<n` a formal input `z_j`. A local setting receives the product of
the `z_j` whose factor-input bits are one. Positions `j>=n` may receive
one, because every valid divisor is less than `2^n` and has zero bits
there. The resulting circuit computes the multilinear divisor polynomial
in the companion note. Each local monomial takes at most `n` ring
multiplications. Substitution into a cyclic ring happens after this
circuit has been generated.

In a cyclic ring `R_r=F_p[X]/(X^r-1)`, compute every exponent in (11)
modulo `r`. Repeated squaring of two modulo `r` and integer arithmetic
on the original digits suffice. The algorithm does not create a
polynomial of degree `2^n`.

For `P^40<=n`, the existing
[top constructor](half-top-constructor.md) gives the same geometry and
the same bound

```math
P^2\lambda\le4717n.
\tag{12}
```

Its local-setting enumeration takes `2^O(n/P)` operations. Instead of
adding one for each valid setting, add its monomial weight. The weighted
tables require at most this many ring and structural operations, times
a polynomial in `n` and `log r`.
Writing all table entries adds their representation cost. No parameter
in (12) grows with `r`.

## 5. Choose the recursion depth and charge ring operations

Choose the largest `ell` with `P^40<=n`, where `P=K_ell`.
For the finite initial range without such a level, use trial division.
The source's elementary schedule estimates give

```math
P^{40}\le n<P^{400},\qquad
e_\ell\le0.85\rho^\ell,
\qquad
\rho^\ell=O((\log n)^{-\alpha}).
\tag{13}
```

Equations (6), (12), and (13) bound the circuit and construction work by

```math
A(n)=\exp\!\left(O\!\left(
\frac{n}{(\log n)^\alpha}\right)\right).
\tag{14}
```

In detail, the extra exponent is `O(P^2 log(n+2))`, and construction
cost contributes `O(n/P)`. Both are smaller than the exponent in (14),
since `n^(1/400)<P<=n^(1/40)` and `alpha<=1`.

Represent an element of `R_r` by its `r` coefficients in the fixed field.
Addition costs `O(r)` field operations. Schoolbook cyclic convolution
costs `O(r^2)`. The field-operation bit cost is constant because `p`
is fixed. Circuit storage, addressing, and sequential-tape simulation
cost a fixed polynomial in `A(n)`, `r`, and `n`.
Consequently, one weighted query takes

```math
\exp\!\left(O\!\left(
\frac{n}{(\log n)^\alpha}+\log(r+2)\right)\right)
\tag{15}
```

finite-alphabet machine steps. This includes all written tables and
ring elements. No random-access or unit-cost large-integer assumption
is required: a polynomial tape simulation changes only fixed constants
in this exponent.

The algorithm never substitutes `r` for the integer height parameter
in the old padded theorem. Such a substitution would introduce an
unjustified large overhead. It executes the circuit whose size was
bounded before choosing the representation of `R`.

## 6. Sparse nonzeroness and the final factoring algorithm

The sparse divisor-polynomial argument supplies a deterministic family
of cyclic rings with the following properties:

- There are `exp(O(n/log n))` rings.
- Every ring dimension is `exp(O(n/log n))`.
- The polynomial for nontrivial divisors with a given prefix is nonzero
  in at least one ring exactly when that prefix contains such a divisor.
- All ring parameters can be generated within `exp(O(n/log n))` steps.

These facts use the bound `tau(N)=exp(O(n/log n))`, an elementary
prime-supply argument, and collisions of exponent residues. The proof
and finite checks are given in
[sparse-divisor-nonzeroness.md](sparse-divisor-nonzeroness.md).

For each prefix, evaluate (10) in those rings. Remove the known divisor
one by subtracting `X` when it has the selected prefix. Test whether
any resulting coefficient vector is nonzero.

Choose the next bit to be zero if that prefix contains a divisor
greater than one, and one otherwise. The invariant is nonempty because
`N` itself is always an eligible divisor. After `n` bits, this returns
the least divisor greater than one. It is prime, since any composite
divisor has a smaller divisor greater than one.

Combining the ring family with (15) gives

```math
T(n)\le\exp\!\left(O\!\left(
\frac{n}{(\log n)^\alpha}+\frac n{\log n}\right)\right).
\tag{16}
```

Because `0<alpha<=1`, the second term is absorbed by the first.
The polynomial number of prefix steps is absorbed as well. This proves
(1). Inputs zero and one can retain the repository's values `2` and `1`.

## 7. What this changes in the research target

A construction confined to one characteristic is not automatically
blocked by the need to recover integer divisor counts. The sparse
nonzeroness route recovers the required factoring decisions instead.
It uses no randomness, integer-count reconstruction, or field extensions.

The remaining algebraic target is still substantial: an effective family
with `rho<1/sqrt(2)`. The crossing identity, both factorizations, uniform
generation, and width bound must all be proved. A finite rank drop or a
zero crossing scalar supplies none of these missing conditions.

This theorem adds a route from such a modular family to a complete
factoring machine. It does not change the strongest proved exponent
in this repository. The new bridge used substantial OpenAI assistant
support and retains the attribution of both pinned source developments.

The [normalized modular construction](modular-crossing.md) supplies
the exact scale-one interface for each fixed `rho>1/sqrt(2)` with
`rho<1`. That is a concrete use of the theorem below the half exponent.
Its smaller state space does not establish a stronger exponent.

Run `python3 tools/check_modular_crossing.py` and
`python3 tools/check_sparse_divisor_nonzeroness.py` for the finite
identity and sparse-recovery checks. Their recorded outputs are
[modular-crossing.json](../evidence/modular-crossing.json) and
[sparse-divisor-nonzeroness.json](../evidence/sparse-divisor-nonzeroness.json).
These tests do not materialize the asymptotic ring-circuit evaluator;
Sections 1--6 and the cited companion proof establish that conditional
algorithm and its cost.
