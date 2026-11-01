# Exact theta sections for the signed divisor coefficients

This note constructs a coefficient engine that does not enumerate divisors.
It turns the signed cyclic divisor data into one coefficient of a theta
quotient and gives an exact binary section algorithm. All intermediate
polynomial supports and ring products are charged. The current implementation
has exponential cost in the binary input length and does **not** improve the
complete factoring exponent, which remains 1/2.

The construction also isolates a restriction on a proposed compression:
already at additive phase order three, the common linear space of depth-k
section states has dimension at least floor(2^k/36). This is a restriction
on common linear representations. The same level-three answer has a direct
polynomial-time square test, so the restriction is not a lower bound for
coefficient algorithms or factoring.

The [earlier coefficient reduction](weight-one-divisor-coefficients.md) and
[sparse signed reconstruction](sparse-divisor-residues.md) supply the recovery
interface. The released Kronecker/Fay lead in Discussion #51 motivated the
theta examination. Jacobi's triple product and coefficient section methods
are classical; no novelty claim is made for them. The implementation, charged
recurrence and common-space witness below were developed with OpenAI Codex
assistance.

## 1. An additive phase algebra

Let q be an odd prime and use

```math
R_q=\mathbb F_2[Z]/\Phi_q(Z),\qquad
\Phi_q(Z)=1+Z+\cdots+Z^{q-1}.
```

This is a different algebra from the even-order multiplicative-character
algebra in the earlier decoder. Here the generator has the **odd additive
phase order q**. We are constructing the signed residue output directly,
not reducing or discarding coordinates of the old supplied character input.

The signed divisor polynomial is

```math
S_N(Z)=\sum_{d\mid N}(Z^d+Z^{-d})\pmod{Z^q-1}.          \tag{1}
```

Its coefficient at zero is zero, opposite coefficients agree, and its
augmentation is zero. Passing to R_q preserves the entire signed vector:
the two degree-below-q lifts differ by Phi_q, whose augmentation is the
odd number q, hence one. Exactly one lift has augmentation zero.

The element 1+Z is a unit in R_q, with explicit inverse

```math
(1+Z)^{-1}=Z+Z^3+\cdots+Z^{q-2}.                         \tag{2}
```

The analogous formula after Z -> Z^a applies whenever q does not divide a.
No field factorization is required, even when Phi_q is reducible.

## 2. The theta quotient

Start in the integral Laurent polynomial coefficient ring, before imposing
any phase relation. Jacobi's triple product gives

```math
\Theta(Z,t)=\sum_{m\in\mathbb Z}(-1)^m Z^m t^{m(m-1)/2}
 =(1-Z)\prod_{n\ge1}(1-t^n)(1-Zt^n)(1-Z^{-1}t^n).
                                                               \tag{3}
```

This is the standard identity in [DLMF 17.8.1](https://dlmf.nist.gov/17.8.E1).
Apply D=Z d/dZ before specialization. Logarithmic differentiation gives

```math
\frac{D\Theta}{\Theta}
 =-\frac{Z}{1-Z}
  -\sum_{n,d\ge1}(Z^d-Z^{-d})t^{nd}.                     \tag{4}
```

Reduce coefficients modulo two and then specialize to R_q. The constant
Theta coefficient is 1+Z, a unit by (2). Thus the quotient is a well-defined
formal series, and for every N>=1,

```math
[t^N]\frac{D\Theta}{\Theta}=S_N(Z)\quad\hbox{in }R_q.    \tag{5}
```

The coefficient at N=0 is Z/(1+Z), not S_0. The engine only accepts positive
indices. Divisors divisible by q contribute zero to (1) automatically.

The order of operations matters. D is not a derivation of the specialized
quotient: D(Z^q)=q Z^q specializes to one, whereas Z^q specializes to the
constant one, whose derivative is zero. The code constructs DTheta from
the integer exponent m before evaluating its phase.

Pair the terms m=s+1 and m=-s in (3). At triangular exponent s(s+1)/2,
the reduced Theta coefficient is Z^(s+1)+Z^(-s). The DTheta coefficient is
Z^(s+1) if s is even, and Z^(-s) if s is odd. There are O(sqrt(N)) terms
through degree N. Formula (3), not an assumed finite q-expansion oracle,
constructs them.

## 3. A complete binary-index section recurrence

For a series A define Lambda_b(A)=sum_j [t^(2j+b)]A t^j, b in {0,1}.
Let sigma act on R_q by Z -> Z^2. This is an automorphism because q is odd.
Write Q_a(t)=Theta(Z^a,t). Characteristic two gives

```math
Q_a(t)^2=Q_{2a}(t^2).                                   \tag{6}
```

If the current target is n=2r+b and its quotient is P/Q_a, then

```math
[t^n]\frac{P}{Q_a}
 =[t^r]\frac{\Lambda_b(PQ_a)}{Q_{2a}}.                   \tag{7}
```

Indeed multiply numerator and denominator by Q_a and use (6); the
denominator only contains even powers of t. Start with P=DTheta and a=1,
and repeat (7) along the binary digits of N from least to most significant.
At target zero the answer is P(0)/(1+Z^a), using (2).

Every step is exact with P and Q_a truncated through the current n: no
coefficient above n contributes to the child through degree floor(n/2).
The denominator phase evolves as a -> 2a mod q and depends only on depth.
The implementation records the full support sizes, candidate term pairs
and ring multiplication count at every step. It never obtains a coefficient
by listing divisors of N.

## 4. Cost and the actual compression gap

The current sparse polynomial implementation allows O(n) numerator terms
and O(sqrt(n)) denominator terms at target n. Thus the sum of candidate
products is bounded by O(N^(3/2)); ring arithmetic and exponent handling
have polynomial cost in q and log(N+1). A conservative binary cost is

```math
(q+\log(N+1))^{O(1)}N^{3/2},
```

with O((q+log(N+1))N+q^2) bits of workspace for polynomial dictionaries,
their indices and the phase table. This is exponential in the input bit
length b=ceil(log_2(N+1)), and is much worse than the accepted factoring
bound. A logarithmic number of recursion levels alone does not make the
algorithm fast.

For the retained target N=65535, the numerator initially has 362 terms.
At q=7 it reaches 11,116 terms and uses 1,576,938 ring products; at q=11
it reaches 14,921 terms and uses 2,236,072 ring products. These are exact
counts for the implemented representation, not asymptotic lower bounds.

One proposed repair is to place every numerator at a given depth in a common
small linear basis. Section 6 proves that a polynomial-size common linear
basis cannot work uniformly in the input length, already at q=3. A nonlinear,
adaptive or parameterized arithmetic representation is outside that result.

## 5. A genuinely cheap level-three coefficient

At q=3 each divisor not divisible by three contributes Z+Z^2=1 in R_3;
each divisor divisible by three contributes zero. Therefore

```math
[t^N]\frac{D\Theta}{\Theta}
 =[N/3^{v_3(N)}\text{ is a square}]\quad\hbox{in }R_3.   \tag{8}
```

The number of divisors of an integer is odd exactly when it is a square,
by pairing d with the complementary divisor. Repeated exact division by
three followed by an integer square root computes (8) in polynomial time
in b. No factorization of N is needed. In particular, the raw polynomial
state growth in Section 4 is not an inherent cost of the coefficient.

This is a single small phase order. It does not supply the collection of
prime-level residue oracles needed by the sparse factoring reduction.

## 6. Exponential dimension for a common linear section space

Let c_N be the binary value in (8), defining c_0=0, and set M=2^k with
k>=6. Consider the finite section matrix

```math
H_k(u,v)=c_{Mu+v},\qquad 0\le u,v<M.
```

We prove, over every field of characteristic two,

```math
\operatorname{rank}H_k\ge\lfloor M/36\rfloor.             \tag{9}
```

Take the odd integers x in (M/2,3M/4) that are not divisible by three.
There are at least floor(M/12) of them. Partition them by
u_x=floor(x^2/M) modulo three; one class I has at least floor(M/36)
members. For each selected root put v_x=x^2 mod M.

The rows u_x are distinct: consecutive eligible odd roots have square
differences greater than M. The columns v_x are distinct too. For any
odd y, its four square roots modulo M are the residues represented by
plus or minus y, optionally shifted by M/2. Exactly one lies in
(M/2,3M/4). This root description follows by factoring z^2-y^2 and
observing that one of z-y,z+y has valuation one.

All selected rows have the same residue A modulo three. Since x is not
divisible by three, v_x=1-MA modulo three. Consequently every cross entry
N_(x,y)=M u_x+v_y is one modulo three. Equation (8) is therefore the ordinary
square predicate on every selected entry.

If N_(x,y)=z^2, then 0<=z<M and z^2=v_y modulo M. Of its four possible
roots, the two below M/2 have row floor(z^2/M) below M/4; the one above
3M/4 has row at least 9M/16. The selected rows lie between these bounds.
Thus z=y, and the distinctness of the selected rows forces x=y. Conversely
N_(x,x)=x^2. The selected submatrix is exactly the identity, proving (9).

Equivalently, the truncated depth-k states

```math
F_{k,v}(t)=\sum_{0\le u<M}c_{Mu+v}t^u,\qquad 0\le v<M,
```

span dimension at least floor(2^k/36). Any common linear representation of
all these states must have that many dimensions, even if its basis depends
on k or its coefficient field is extended.

In (7), the denominator at depth k is independent of the suffix v and is
a unit series. Multiplication by that common denominator is invertible on
truncated series, so the same lower bound applies to a common linear basis
for the corresponding theta numerators. For indices with at most 2k bits,
the required dimension is at least 2^k/36-1.

This excludes a declared linear compression architecture. It does not
exclude sparse formulas depending nonlinearly on the target, bounded-time
quadratic predicates, or other coefficient engines. Formula (8) explicitly
provides such a cheap nonlinear algorithm in this very example. In
particular, neither (9) nor the observed supports imply a factoring lower
bound. Fast coefficient algorithms for **algebraic** series, such as
[Bostan--Christol--Dumas](https://arxiv.org/abs/1602.00545), require their
algebraicity/finite-state hypotheses; those hypotheses cannot be supplied
by assuming this common linear closure.

## 7. Reproduction and next target

Run `make verify-lambert-theta`, or

```sh
python3 tools/check_lambert_theta_sections.py
python3 tools/check_lambert_theta_sections.py --profile-bits 8 10 12 14 16
```

The [engine](../tools/lambert_theta_sections.py) is separate from the
[checker](../tools/check_lambert_theta_sections.py). The
[receipt](../evidence/lambert-theta-sections.json) retains source hashes and:

- 2,056 exact coefficient comparisons for N=1,...,257 at eight prime levels;
- 325 theta-quotient coefficients checked by independent cyclic convolution;
- 293,717 entries of explicit identity submatrices, with largest dimension 468;
- additional literal-divisor and theta-engine checks of small submatrices;
- six closed-form level-three controls through 701 input bits;
- charged numerator profiles through 16-bit targets.

The larger rank controls evaluate the independently proved square predicate;
they do not run the slow general theta engine on those large indices.
All coefficient comparisons retain inputs divisible by the phase prime.
These are ordinary proofs and finite implementation controls, not new Lean
verification or independent mathematical acceptance.

The next construction target is a nonlinear representation that uses the
quadratic theta structure while charging its entire update cost. A fixed or
polynomial-size common linear numerator space is not an available shortcut.
