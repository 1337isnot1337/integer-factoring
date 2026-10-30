# From a fixed-characteristic evaluator to deterministic factoring

## Result and scope

A faster crossing construction over one fixed finite field can be useful
for integer factoring. It need not work over the rationals or over enough
different prime characteristics to reconstruct an integer count.

This note proves the following transfer statement. Fix a prime `p`.
Suppose one uniform algorithm computes, for every divisor-prefix query and
every prime `q`, the polynomial

```math
F_{N,k,x}(Z)
=\sum_{\substack{d\mid N\\ \lfloor d/2^k\rfloor=x}}Z^d
\quad\text{in }\ R_q=\mathbb F_p[Z]/(Z^q-1)
\tag{1}
```

within `F(n)*(nq)^a` bit operations, where `n` is the bit length of `N`,
`a` is a fixed constant, and `F(n)>=1`. Then one deterministic
finite-alphabet machine computes the smallest prime factor in

```math
F(n)\,2^{O(n/\log n)}.
\tag{2}
```

All finite input exceptions are paid by a fixed fallback. In particular,
if `F(n)=2^{O(n/(log n)^alpha)}`, the resulting saving exponent is
`min(alpha,1)`. Thus any such fixed-characteristic improvement with
`1/2<alpha<=1` would improve the current factoring bound.

The transfer proof is complete. The faster evaluator in its hypothesis
is not supplied by this note. The new modular crossing construction
developed alongside it does not yet establish `alpha>1/2`.
No improved factoring exponent is claimed here.

The weighted network that represents (1), including its construction
cost, is described in [the weighted-constructor note](weighted-prefix-constructor.md).
The sparse polynomial test below is a standard method. In particular,
Bläser, Hardt, Lipton and Vishnoi's
[*Deterministically Testing Sparse Polynomial Identities of Unbounded Degree*](https://cs.yale.edu/homes/vishnoi/Publications_files/BHLV08-ipl.pdf)
uses reductions modulo `Z^q-1`. The contribution here is its explicit
connection to weighted divisor queries and the cost required for factoring.
An elementary proof is included so the transfer needs no black-box
identity-testing theorem.

## 1. Why an ordinary modular count is insufficient

A nonempty set can have cardinality zero modulo `p`. Merely replacing
the integer prefix count by its residue therefore invalidates the prefix
search. In characteristic two, two divisors already give a counterexample.

In contrast, each divisor occurs once in (1). Distinct divisors give
distinct monomials with coefficient one. Consequently the polynomial
over `F_p[Z]` is nonzero exactly when the divisor set is nonempty.
The task is to detect that fact without writing its potentially
exponential degree.

The quotient in (1) can have zero divisors. That causes no problem:
`Z^q-1` is monic, so every residue has a unique representative of degree
less than `q`. Its coefficients can be compared with zero. No inversion
of a nonconstant element of `R_q` is used.

## 2. An exact sparse-polynomial test

Let `f` be a nonzero polynomial over any field. Suppose it has at most
`S` nonzero terms and degree less than `2^n`, with `n>=1`.
For any collection of

```math
M=n(S-1)+1
\tag{3}
```

distinct prime numbers `q`, at least one residue `f mod (Z^q-1)` is nonzero.

**Proof.** Fix one exponent `e_0` in the support of `f`. This is a proof
choice; the algorithm does not need to know `e_0`. For a different
support exponent `e`, a collision with `e_0` modulo `q` requires
`q | (e-e_0)`. The nonzero integer `|e-e_0|<2^n` has at most `n`
distinct prime divisors: the product of `n+1` distinct primes is at
least `2^(n+1)`.

There are at most `S-1` other support exponents. Thus at most `n(S-1)`
primes can cause any collision with `e_0`. One of the primes in (3)
leaves that exponent alone in its residue class. Its coefficient remains
nonzero after reduction, proving the claim.

This only isolates one monomial. Separating every pair would give an
unnecessary quadratic dependence on `S` in the number of primes.
The argument also permits `q=p`. Separability of `Z^q-1` is not required.

An algorithm therefore tests all the residues for the first `M` primes.
It returns nonzero as soon as one residue is nonzero, and returns zero
only after all `M` tests vanish. This is deterministic and exact.

### A completely elementary bound on the primes

The first `M` primes are at most

```math
Q=16(M+1)^2.
\tag{4}
```

Here is a sufficient elementary bound. For any positive integer `m`,

```math
\binom{2m}{m}\ge2^m.
```

For every prime `r`, Legendre's formula gives

```math
v_r\binom{2m}{m}
=\sum_{j\ge1}
\left(\left\lfloor\frac{2m}{r^j}\right\rfloor
      -2\left\lfloor\frac{m}{r^j}\right\rfloor\right)
\le\lfloor\log_r(2m)\rfloor.
```

Each summand is zero or one. Hence
`binom(2m,m) <= (2m)^pi(2m)` and
`pi(2m) >= m/log_2(2m)`.
Apply this with `2m=Q`. Since
`log_2 Q=4+2log_2(M+1)<=4+2M`,

```math
\pi(Q)\ge\frac{8(M+1)^2}{4+2M}\ge M.
```

A sieve or trial division up to (4) finds the required primes in
`poly(M)` bit operations. The proof does not need an estimate for primes
in arithmetic progressions, an irreducible-polynomial construction,
or the prime number theorem.

## 3. A computable sparsity bound from the input length

The required bound on the number of divisors can be obtained without
factoring the input. Put

```math
k=\lfloor\log_2 n\rfloor,\qquad
S(n)=2^{\lceil3n/k\rceil}
\quad(k\ge18).
\tag{5}
```

Then every positive integer `N<2^n` satisfies `tau(N)<=S(n)`.

To prove this, write `N=prod r^{a_r}` and split its prime divisors at
`R=2^ceil(k/2)`. For `r>R`, use `a_r+1<=2^{a_r}` and
`sum a_r log_2 r<n` to obtain

```math
\log_2\prod_{r>R}(a_r+1)
\le\sum_{r>R}a_r
\le\frac{2n}{k}.
```

There are at most `R` smaller primes, and every exponent is less than
`n`. Since `n<2^(k+1)`,

```math
\log_2\prod_{r\le R}(a_r+1)
\le R\log_2(n+1)
\le2^{\lceil k/2\rceil}(k+2)
\le\frac nk.
```

The last inequality follows from
`k(k+2)<=2^floor(k/2)` for `k>=18`.
It holds at `k=18,19`. Increasing `k` by two multiplies the left side
by `(k+4)/k<=2` and the right side by two.
Also `n>=2^k`. This proves the displayed estimate and then (5).

Thus `log S(n)=O(n/log n)`. A trivial bound or trial division handles
the finitely many smaller bit lengths. This argument is a quantitative
version of the classical divisor bound; see also
[Tao's elementary exposition](https://terrytao.wordpress.com/2008/09/23/the-divisor-bound/).
The algorithm computes (5) from `n`, without knowing the primes `r`.

## 4. Exact prefix search in a fixed characteristic

For the prefix query, remove divisor one by subtracting

```math
[\lfloor1/2^k\rfloor=x]\,Z
```

from (1), in each quotient ring. The polynomial before reduction is now

```math
G_{N,k,x}(Z)
=\sum_{\substack{d\mid N,\ d>1\\ \lfloor d/2^k\rfloor=x}}Z^d.
\tag{6}
```

It still has at most `S(n)` terms and degree less than `2^n`.
The sparse test therefore decides whether that prefix contains a divisor
greater than one. It is not a test of the residue of its cardinality.

Start at the empty prefix `b=0`. For `k=n-1,n-2,...,0`, test
`G_(N,k,2b)`. Replace `b` by `2b` when it is nonzero, and by `2b+1`
otherwise. There is always a divisor greater than one under the current
prefix, initially `N` itself. Choosing the first nonempty child preserves
the least such prefix. After exactly `n` steps, `b` is the least divisor
greater than one and is therefore prime.

This proves correctness for every `N>=2`. The conventions for `N=0,1`
can be retained with two explicit branches.

## 5. Complete bit cost and the precise evaluator hypothesis

There are at most `nM` weighted evaluations. Equations (3)--(5) give

```math
M\le nS(n)+1,\qquad
q\le Q=16(M+1)^2,\qquad
\log(nMQ)=O(n/\log n).
```

Each representative in `R_q` uses `q` coefficients of the fixed field.
Addition, multiplication by cyclic convolution, subtraction of the
divisor-one monomial, and a zero test cost `poly(q)` bit operations.
Elementary multiplication of two representatives takes `O(q^2)` field
operations, which is sufficient here. Computing powers of `Z` reduces
the exponent modulo `q`.

Under the stated evaluation bound `F(n)*(nq)^a`, the entire prefix
search, prime generation, ring arithmetic and representation cost are

```math
F(n)\operatorname{poly}(n,M,Q)
=F(n)2^{O(n/\log n)}.
```

All constants and polynomial degrees are fixed. The construction and
search are uniform. If an arithmetic or random-access program is used
as the starting point, its fixed polynomial finite-tape simulation
preserves the displayed exponential rate up to constants. When the
hypothesis is already in bit-machine steps, the bound (2) applies as
written without another simulation exponent.

### A necessary distinction when importing an evaluator

It is sufficient to prove that the evaluator uses `A(P,lambda)` ring
operations over every commutative `F_p`-algebra, with all needed inverses
belonging to the fixed base field. Simulating those operations in `R_q`
then multiplies the arithmetic cost by `poly(q)`. A fixed polynomial
machine simulation must also be charged if the arithmetic model needs it.

It is not sufficient to insert the ring dimension `q` into the entry
bit-length parameter `B` of the existing integer `PartB` theorem. That
theorem has a term of the form `P^2 log(P(lambda+B))`; such a substitution
could put `P^2 log q` into the exponent instead of paying `poly(q)`
outside the operation count. The theorem in this note explicitly
requires the latter cost.

The [weighted constructor](weighted-prefix-constructor.md) has cost
`2^{O(n/P)}poly(n,q)` and preserves the original channel widths.
A modular crossing identity with a nonzero base-field scale is useful
only after the complete recursive operation count is established.
No integer-to-modular evaluator conversion is assumed without that step.

## 6. Verification and remaining target

Run:

```sh
python3 tools/check_fixed_characteristic_transfer.py
```

The checker tests sparse quotient detection over several characteristics,
prefix recovery using a finite exact divisor oracle, ring reductions,
the integer bounds above, and cancellation negative controls.
The divisor oracle in these finite tests is deliberately simple.
It is not an implementation or benchmark of a faster factoring evaluator.

The remaining constructive target is an effective fixed-characteristic
crossing family with a better contraction, together with its recursive
ring-operation bound. This transfer removes the need to lift that family
to characteristic zero or to reconstruct exact integer counts.
It establishes a new route to an improvement, not an improvement itself.
