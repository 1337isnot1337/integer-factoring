# A degree-independent floor from the generic axis faces

## Question and result

Base: `5c2e66d`. The registered symmetric endpoint probe fails to supply
the desired prefix after setting all other variables to one. That
specialization can hide generic face rank. Instead combine actual face
rank with the exact complementary Smith pairing. This retains the coupled
axis variables and gives a general restriction without another kernel sweep.

Let `d>=2`, let m be positive and even, and let the d polynomial matrices
`B_a` over characteristic-two fields have order m and satisfy

```
B_a^2 = x_a^(2D) I,              D>0.
```

The other variables may belong to the coefficient field. Write `r_a`
for the rank of `B_a(0)` over that field, and put

```
c = (sum_a r_a)/(d*m),
rho = min(1/2, c+1/d).
```

For the common-scale Gaussian profile pressure in the
[shared-minor note](multiaxis-gaussian-profile-floor.md),

```
R_d(p,t) > rho       (0<p<1, 0<t<=1).                    (1)
```

Thus an average quarter-rank face gives `rho=min(1/2,1/4+1/d)`.
For **every d>=4** this is at least `1/sqrt(d)`. It excludes the strict
half-score profile criterion under these assumptions, independently of
degree, field size and coefficient choices. No complete factoring bound
changes. Non-Cayley kernels without this square identity, charged physical
types, other rank estimates and other algorithmic routes remain outside.

## Complementary orders and fractional retained channels

Let `s_(a,1)<=...<=s_(a,m)` be the nonnegative Smith orders. The inverse
identity `B_a^(-1)=x_a^(-2D) B_a` gives

```
s_(a,i)+s_(a,m+1-i)=2D.
```

There are exactly r_a zero orders, and the first m/2 orders are at most D.
In particular `0<=c<=1/2`. Pool these low halves across all d axes. Retain
every zero order, then enough of the remaining low-half orders to give
`N=d*m*rho` channels. If N is not integral, retain a fraction of one
channel. Every retention weight lies in `[0,1]`; this is an inequality
for the logarithmic partition sum and does not posit a fractional tensor.

The retained count is between the zero count and the available low-half
count. Its total order is at most

```
(N-d*m*c)D <= mD.                                      (2)
```

The function `f(s)=log(1+t^s)` is decreasing and convex for `0<t<1`.
Weighted Jensen and (2) therefore give

```
sum_(a,i) log(1+t^s_(a,i))
  > N log(1+t^(mD/N))
  = d*m*rho log(1+t^(D/(d*rho))).                       (3)
```

The strict inequality is valid because at least the entire high halves
are discarded and every finite order contributes a positive logarithm.
For `u=t^(D/(d*rho))`, subtracting the physical charge `mDp log t`
from (3) and using

```
log(1+u)-p log u >= H(p)
```

proves (1). At t=1 the pressure is `log(2)/H(p)>=1>rho` directly.
The floor is pointwise strict; its infimum need not exceed rho.
The proof needs neither equal face ranks nor equal profiles at the axes.

## A supported coefficient gives the required face minor

Now specialize to the subset-local Cayley model. Its crossblock is
multiplication on `2^(d-1)` exterior monomials by

```
F_a = sum_(A subset [d]\{a}) c_({a} union A) epsilon_A.
```

Here `c_T` is homogeneous of degree D, depends only on variables in T,
and has endpoint one. In particular `c_{a}=x_a^D`, so the displayed
square identity holds. Suppose a coefficient of size j in `F_a(0)`
is nonzero. Use columns indexed by all subsets I of the complement of
its support A, and rows `A union J` over the same complement. The entry
can be nonzero only when `I subset J`; every diagonal entry is that
coefficient. This is a nonzero triangular minor of size `m/2^j`.
This supported-monomial argument is reused from the
[bi-exterior degree budget](bi-exterior-degree-budget.md), with credit.

For a permutation-symmetric coefficient family, if some c_T with
`|T|=k` is not divisible by `product_(a in T) x_a`, permutation symmetry
supplies such a face at every axis with j=k-1. Consequently

```
R_d(p,t) > min(1/2, 2^(1-k)+1/d).                       (4)
```

For k=2 or 3 and d>=4, this is at least `1/sqrt(d)`. In the quarter-rank
case and d>4 the exact comparison is

```
d*(1/4+1/d)^2-1 = (d-4)^2/(16d) > 0.
```

At d=4 the floor is exactly 1/2 and (1) still excludes a strict crossing.
Thus any winning symmetric Cayley kernel in these dimensions must have
every pair and triple coefficient divisible by all its variables.
This is necessary, not sufficient. More generally k=4,5,6 give the same
exclusion for d>=47,223,959 respectively. The exact rational comparison
in (4), rather than a rounded square root, determines these thresholds.

## Why the all-one specialization missed it

The single registered endpoint fixture has d=9, D=8. For a nonempty
S of size k let `r=2^v2(k)` and set

```
c_S = elementary_r((x_i^(8/r))_(i in S)).
```

Lucas parity gives `binom(k,r)=1 mod 2`, so every coefficient has degree
8 and endpoint one. Set `c_empty=sum_(S nonempty)c_S` for the zero
diagonal Cayley kernel. At one axis, specializing the other eight
variables to one leaves the symmetric exterior coefficients

```
a_j(z)=z^(8/2^v2(j+1)),       j=0,...,8.
```

The credited [simultaneous exterior blocks](recursive-symmetric-exterior-blocks.md)
give the complete specialization profile

```
1^2,2^52,3^2,4^16,5^2,6^52,7^2,
9^2,10^52,11^2,12^16,13^2,14^52,15^2.
```

Its first 86 orders sum to 258, above `8*86/3`, so no prefix of more
than `256/3` channels meets that particular Jensen test. This disproves
the proposed specialization-to-prefix implication, not a universal
pressure floor.

Generically, however, `c_{aij}(x_a=0)=x_i^8+x_j^8` is nonzero.
The triangular minor gives rank at least 64 on every face, hence
**`R_9>13/36>1/3`** by (1). A literal F4 specialization with one other
variable equal to omega preserves that minor; the all-one specialization
annihilates every face coefficient. This is a concrete instance where
the generic obstruction is invisible at that symmetric specialization.
The fixture is closed by the general lemma; no neighboring search follows.

## Exact controls and limits

Run `make verify-multiaxis-face-rank`. The
[checker](../tools/check_multiaxis_face_rank_floor.py) verifies rational
fractional-selection budgets on unequal paired profiles, exact threshold
comparisons, the literal 64-by-64 F4 triangular minor and scalar-restricted
face rank. It computes all four small symmetric block profiles with an
untruncated fraction-free polynomial oracle and reconstructs all 256
specialized orders. Its [receipt](../evidence/multiaxis-face-rank-controls.json)
retains the profiles, selected minor, rational budgets and source hashes.
The ordinary argument, not a finite dimension screen, proves (1).

This adds a global face-rank restriction to the earlier endpoint and
quotient constraints. It supplies no winning crossing, stronger factoring
score, Lean formalization or independent mathematical review. Substantial
OpenAI assistant support; the Cayley, Gaussian transfer, supported-monomial
minor and simultaneous-block framework retain their preceding credits.
