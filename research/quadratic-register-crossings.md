# Complete six-bit nonlinear register catalogue

## Result and scope

Every quadratic Boolean function on five lower register bits, with an
added zero-tail indicator, has been tested as a six-bit feedback rule.
Of the 65,536 functions, exactly 1,952 yield one full 64-cycle.
All ten balanced physical bit partitions of every qualifying cycle,
19,520 systems in total, have exact certificates
`rho_bare >= 17/24 > 1/sqrt(2)`.

This is a complete finite catalogue in its stated function class. It
does not exclude higher Boolean degree, more register bits, arbitrary
physical bases, word-type gauges, or other crossing criteria. It does
not improve the factoring score, which stays `1/2`.

The function class goes beyond maximal linear feedback, so the
[copied-wire obstruction](shift-register-copy-wire-barrier.md) covers
arbitrary feedback only at its stated large register lengths. This
six-bit catalogue lies below that cutoff. Qualification here is by a
literal full-cycle test.

## 1. Explicit nonlinear permutations

Write a state as six bits, with five-bit lower tail `b`. Choose

```math
Q(b)=c+\sum_{i=0}^4c_i b_i+
               \sum_{0\le i<j\le4}c_{ij}b_i b_j\quad\hbox{over }\mathbb F_2.
```

There are `1+5+10=16` independent coefficients. Define

```math
U_Q(s)=((2s)\bmod64)+
 \left(s_5+Q(b)+[b=0]\pmod2\right).                         \tag{1}
```

The output's upper five bits give the old lower tail. Its remaining
bit determines the old top bit uniquely. Thus every map (1) is a
permutation, before any cycle qualification. The indicator is the
degree-five Boolean polynomial `product_i(1+b_i)`, so it is independent
of the quadratic coefficient choices.

For each qualifying full cycle, `I+U_Q` is regular nilpotent of order
64. The [cyclic-algebra constructor](cyclic-algebra-crossings.md) gives
the complete polynomial crossing endpoint, degree 63, and column zero
orders `(0,64,...,64)`. The same complement symmetry as in the
[linear-register catalogue](shift-register-crossings.md) reduces the
balanced bit partitions from twenty to ten.

## 2. Exact certificates for the complete catalogue

For the normalized reversed row polynomial, put
`A0=Phi(U_Q)`, `A1=Phi(U_Q+U_Q^2)`. Literal binary matrix ranks give
the zero-order count `r` and order-one count `k`. The established
per-channel inequalities give

```math
B(s)\ge64r+63k\qquad(s\ge1).
```

The code computes `k` only when `64r` alone fails the exact comparison
`(64r)^6>=8^17`. Every qualifying system satisfies
`(64r+63k)^6>=8^17`, proving the stated bound `rho_bare>=17/24`.
An omitted order-one count means the zero-order certificate sufficed;
it does not denote zero pivots.

The smallest reported certificate floor is `448`, at coefficient code
`70`, first coordinate set `{0,1,3}`. In the recorded monomial order,
this is `Q=b_0+b_1+b_0b_1`. Its row zero-order rank is seven, so
`64*7=448` already suffices. Its later row profile and every physical
endpoint entry are independently recovered by the checker.

## 3. Reproduction

```sh
make verify-quadratic-register
python3 tools/check_quadratic_register_crossings.py \
  --output /tmp/quadratic-register-crossings.json
```

The standard-library checker independently evaluates every truth table
from its literal Boolean monomials, enumerates all 65,536 permutations,
tests which are full cycles, and checks every balanced partition.
It records the monomial order, a digest of every accepted function
code, the exact rank histogram, and the pressure comparisons. The
minimum case also has its full polynomial endpoint checked and its
complete row profile recovered by block-Toeplitz ranks.

Saved output: [`quadratic-register-crossings.json`](../evidence/quadratic-register-crossings.json).
The complete polynomial construction and transfer retain the preceding
work's attribution. This catalogue and its controls used substantial
OpenAI Codex assistance, without new Lean formalization.
