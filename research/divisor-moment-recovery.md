# A uniform Eisenstein coefficient engine would recover factors deterministically

This is a counting alternative to smooth-norm witness search. After trial
division at the arithmetic subexponential scale, the divisors of N form a
positive measure with at most exp(O(sqrt(b))) support points. A short list
of its odd power sums determines every support point. Exact Hankel
recovery and integer root isolation then give the smallest prime factor.

The power sums are coefficients of level-one Eisenstein series. Therefore
one sufficiently uniform coefficient algorithm would yield a complete
deterministic factoring bound, with its exponent preserved. The coefficient
algorithm is **unproved**. This note supplies the exact reduction and
executable postprocessing, not a faster way to compute those coefficients.
The complete factoring score remains `1/2`.

The recovery is classical Prony/Hankel reconstruction of a finite positive
measure. No novelty claim is made for that method or the divisor-sum
formula for Eisenstein coefficients. The contribution is the uniform
parameter choice, deterministic smallest-factor endpoint, full bit-cost
accounting and independently supplied-moment controls.

## 1. The divisor measure has small support after cheap trial division

Let N>=2 have b binary digits. Put

```math
t=\lceil\sqrt b\rceil,\qquad
D=2^{\lceil b/t\rceil}.
```

Test the integers 2 through min(N,2^t) in order. A first divisor is the
smallest prime factor, so return it. This costs exp(O(sqrt(b))) bit
operations. If no divisor was found, write N=product p_i^(e_i) in the
proof; this factorization is not supplied to the algorithm. Every p_i
is greater than 2^t. Hence

```math
\sum_i e_i<b/t,\qquad
\tau(N)=\prod_i(e_i+1)\le2^{\sum_i e_i}\le D
       =\exp(O(\sqrt b)).                                \tag{1}
```

No squarefree or semiprime assumption is needed. Prime inputs are also
covered: their smallest prime factor is N.

For j>=0 define

```math
a_j=\sigma_{2j+3}(N)
   =\sum_{d\mid N}d^{2j+3}
   =\sum_{d\mid N} w_d\lambda_d^j,
\qquad w_d=d^3>0,\quad\lambda_d=d^2.                    \tag{2}
```

The nodes lambda_d are distinct positive integers, and all weights are
positive. The algorithm requests a_0,...,a_(2D-1) from an exact provider.
It does not enumerate divisors to obtain them.

## 2. Recover the support with rational arithmetic

Let r=tau(N). For k<=D form the Hankel matrix

```math
H_k=(a_{i+j})_{0\le i,j<k}.
```

Its quadratic form is

```math
v^T H_k v
 =\sum_{d\mid N}d^3\left(\sum_{i=0}^{k-1}v_i(d^2)^i\right)^2.
```

For k<=r it is positive definite: a nonzero polynomial of degree below
k cannot vanish at all r distinct nodes. For k>r its rank is r by the
Vandermonde factorization. Exact elimination of H_D thus recovers r.
In particular the first r leading pivots are positive and every
remaining Schur-complement entry is zero. There is no numerical rank
threshold or probabilistic choice of a nonsingular minor.

Solve the nonsingular r-by-r system

```math
\sum_{j=0}^{r-1}a_{i+j}c_j=-a_{i+r},\qquad0\le i<r.     \tag{3}
```

Its unique solution is the coefficient vector of

```math
P(X)=X^r+\sum_{j=0}^{r-1}c_jX^j
    =\prod_{d\mid N}(X-d^2).                            \tag{4}
```

Indeed the right-hand polynomial gives a recurrence for (2), hence
satisfies (3), and H_r is nonsingular. This proves integrality of the
recovered coefficients without assuming it computationally.

Remove the known factor X-1. All remaining roots are squared divisors
between 4 and N^2. A Sturm chain and binary search with half-integer
endpoints isolate the least remaining root in O(b) root-count queries.
The endpoints never equal an integer root. Its integer square root is
the least divisor of N greater than one, and therefore is prime. This
computes the required smallest prime factor without integer-polynomial
factorization, a group-generation hypothesis or a random splitting base.

## 3. The full bit cost of recovery

All requested moments satisfy

```math
\operatorname{bits}(a_j)
 \le1+\log_2D+(2j+3)b=O(Db),\qquad j<2D.
```

Every coefficient of (4) has O(rb+r) bits, since its roots are at most
N^2. Exact elimination and rational solution of (3) have polynomial bit
cost in D and b: their intermediate numerators and denominators are
ratios of minors of matrices of dimension at most D and entries of
O(Db) bits. Hadamard bounds those minor bit lengths polynomially.

For the Sturm step, primitive integer remainders retain coefficient
heights bounded by the relevant Sylvester minors. Exact remainder
computation, sign evaluation at half-integers and the O(b) bisection
steps therefore also have polynomial bit cost in D and b. No numeric
precision promise is involved. A standard deterministic root-isolation
algorithm could be used instead with the same polynomial conclusion.

Thus trial division, recovery and output together cost exp(O(sqrt(b))).
The additional input is exactly 2D coefficient queries, with maximum
power-sum exponent 4D+1. These quantities and their output lengths are
charged; a constant-weight coefficient oracle is not being assumed to
cover this growing family.

## 4. One precise missing coefficient algorithm

Use the q-coefficient-one normalization of Eisenstein series in
William Stein's [*Modular Forms, A Computational Approach*, Chapter 2,
Proposition 2.7 and Definition 2.9](https://wstein.org/books/modform/modform/level_one.html):

```math
\mathcal E_w(q)=-\frac{B_w}{2w}
                 +\sum_{n\ge1}\sigma_{w-1}(n)q^n,
\qquad w\ge4\text{ even}.                              \tag{5}
```

It is a level-one modular form. Therefore

```math
a_j=[q^N]\mathcal E_{2j+4}(q).                          \tag{6}
```

Only N's b-bit binary encoding and the even weight w are supplied to
the coefficient engine. The prime factorization of N is not input.
Every query uses the complete coefficient at the same index N; changing
the weight preserves the full modular form. The moment recovery uses
that symmetry without introducing interval restrictions on its divisors.
The usual divisor-sum formula is an identity, not an algorithm to
evaluate (6) without divisors. Computing all coefficients through N
would also be exponential in b and does not meet the target.

**Conditional transfer.** Fix beta>0. Suppose a uniform deterministic
algorithm computes (6), for all weights queried above, within

```math
2^{O(b/(\log(b+2))^\beta)}(b+w)^{O(1)}                 \tag{7}
```

bit operations, with constants independent of N and w. Then a uniform
deterministic machine computes the smallest prime factor of every N in
time 2^O(b/(log(b+2))^beta).

To prove this, multiply (7) by the number 2D of queries. Both D and
w<=4D+2 are exp(O(sqrt(b))). The recovery cost has the same smaller
scale. For each fixed beta, sqrt(b)=o(b/log^beta(b+2)), so all overhead
is absorbed into the claimed bound. Finite inputs can be handled by
the explicit trial procedure. All arithmetic is exact and its bit cost
has been counted, giving a deterministic finite-alphabet implementation.

In particular a beta>1/2 engine would improve the score. There is no
such engine in this note. Algorithms stated for fixed weight with an
unspecified weight-dependent constant do not establish (7). Neither
does an algorithm requiring the factorization of the coefficient index.

This uses a family of exact divisor moments. Bach, Miller and Shallit's
[*Sums of divisors, perfect numbers and factoring*](https://www.cs.cmu.edu/~glmiller/Publications/Papers/BMS86.pdf),
Theorem 2, gives a randomized polynomial reduction from a single ordinary
divisor sum to factoring. That result is not silently made deterministic
here; (2)--(4) explicitly recover the support from the larger moment
family. The [Cauchy pole-count lane](pole-cauchy-content.md) targets a
different global counting primitive, determinant valuations on factor
intervals. Both remain unresolved computational interfaces.

## 5. Exact implementation and evidence

Run `make verify-divisor-moments`, or

```sh
python3 tools/check_divisor_moment_recovery.py > evidence/divisor-moment-recovery.json
```

The [core](../tools/divisor_moment_recovery.py) accepts an explicit exact
moment callback. It implements positive Hankel elimination, rational
recovery, recurrence checks, primitive Sturm remainders and root
isolation. The uniform wrapper implements the trial/rank bound (1).

The [checker](../tools/check_divisor_moment_recovery.py) enumerates divisor
lists for small inputs and charges that enumeration to its finite
oracle. It compares the recovered polynomial to the independent product
over every squared divisor, with exact and oversupplied rank bounds.
For 60/80-bit controls it takes the existing factor certificates and
independently proves both factors prime by trial division, supplying
four-node moment lists. That source knowledge and expensive trial check
belong to the finite provider, not an algorithmic advance on those inputs.

The [receipt](../evidence/divisor-moment-recovery.json) records 516 core
recoveries and 255 uniform-wrapper controls, including prime powers and
nonsquarefree inputs. The two large controls recover their smallest
prime factors from 16 supplied moments each, through weight 34. None of
these tests implements the missing fast coefficient engine or proves
(7). The ordinary argument has not received independent mathematical
review. This work used substantial OpenAI Codex assistance and no Lean.
