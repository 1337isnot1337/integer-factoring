# An exact obstruction for the complete coupled two-place pencil

Every nonzero member of the eleven-parameter pencil from
[the zero-loop construction](joint-zero-rational-loops.md) has typed
Smith parameter at least `8/11`, at the inherited physical heights.
This is strictly above `1/sqrt(2)`. The result holds over every
characteristic-two field and permits every coefficient choice, rather
than just the finite `F4` catalogue.

The zero polynomial gives the original reference and is excluded from
the `8/11` floor. Other heights, actual coefficient decompositions and
other generator directions remain open. This is an obstruction to the
stated Smith-profile certificate, not a general tensor-rank lower bound
or an improved factoring bound. The score remains `1/2`.

## 1. A rank-two update controls the whole zero profile

Write the strict polynomial family as

```math
A_q=A+q(t)(K+t^4L),\qquad \deg q\le11,\quad q(1)=0,
\qquad U_q=I+\frac{q(t)}{t^8(t+1)^4}N.
```

Let `R=D_H A`, and normalize its zero matrix as `S=t^(-8)R`. The
checked commutator identity is

```math
[N,R]=t^{16}(t+1)^4J,
\qquad t^{-8}R_q=S+qJ,
```

with a constant binary rank-two matrix `J`. Take one minimum-valuation
minor at every prefix of the reference Smith profile. Its determinant
in `S+wJ` is a polynomial in `w` of degree at most two. All 25 selected
minors have no linear term. The only nonzero quadratic terms occur at
prefixes 9 through 17, with valuations

```text
prefix:     9  10  11  12  13  14  15  16  17
w^2 term:  26  35  44  54  64  75  86  98 116.
```

The constant terms attain the reference normalized prefixes

```text
42,51,60,70,80,91,102,114,126.
```

Put `a=ord_0(q)`. If `a>=8`, the gauge is a local invertible energy
unit at zero, so the whole zero profile equals the reference. For
`a<=7`, the quadratic term is uniquely smallest at prefixes 9 through
16 and reduces their ceilings by `16-2a`. At prefix 17 the only possible
tie is `a=5`. Convexity bounds that prefix by the average of the adjacent
ceilings: `(108+138)/2=123<126`. Every other prefix has a minor independent
of `q`. Consequently the reference weakly prefix-dominates the zero
profile for every polynomial `q`.

For a fixed nonzero degree `d<=7`, using `a<=d` gives the stronger safe
ceilings `reference_prefix-(16-2d)` at prefixes 9 through 16. All other
reference ceilings remain valid. These are the zero bounds used below.

## 2. Infinity minors and finite jet stability

Normalize the opposite-cut matrix at infinity by its original pole
eight. The pencil has the form

```math
C_g(z)=C(z)+g(z)\Delta(z),\qquad
g(z)=h(z^{-1}),\quad\operatorname{ord}_z g=12-d,
```

where `h=q/(t^8(t+1)^4)` and `d=deg(q)` lies in `1,...,11`.
The checker retains one minimum-valuation reference minor at every
prefix, expanded as a polynomial in a formal variable `g`. If its
coefficient valuations are `v_j`, a unique minimum of
`v_j+j(12-d)` gives an exact valuation for that minor after substitution,
over any field. Tied minima are omitted; no genericity assumption is
made about their cancellation.

There is an additional exact bound. Since `C_g=C mod z^(12-d)`, all
Smith orders strictly below `12-d` equal their reference values. This
follows by reducing the local cokernel modulo that power: each Smith
summand is truncated at the same order.

The lower convex hull of these prefix ceilings bounds the increasing
Smith-prefix polygon. The checker constructs such a hull separately
for each of the eleven degrees. If no high-prefix bound is available,
only the initial bounded part is used in the partition lower bound.

## 3. Integer pressure certificates

For increasing vectors whose partial sums are bounded above by a hull,
the exponential partition sum is bounded below by the sum from the
hull slopes. One direct proof uses

```math
\sum_i(c-a_i)_+=\max_k\left(kc-\sum_{i=1}^k a_i\right)
```

and integrates this hinge inequality against `tau^2 exp(-tau*c)`.
Rounding each hull slope upward further weakens the partition lower
bound and leaves only integer exponents.

For each degree `d`, let `P_d(x)` be the product of these two rounded
partition polynomials, with `x=exp(-tau/8)`. The two normalizing powers
from the pole and the common zero factor cancel. At `rho=8/11`, the
physical partition is

```math
F(x)=1+x^{11}+x^{33}+x^{44}+x^{55}.
```

The exact pressure comparison is `P_d(x)>F(x)^(32/11)`. Weighted
arithmetic-geometric mean gives the sufficient polynomial bound

```math
F^{32/11}=F^2 F^{10/11}\le\frac{F^2+10F^3}{11}.
```

Thus the checker proves positivity of the integer polynomial

```math
D_d(x)=11P_d(x)-F(x)^2-10F(x)^3.
```

Every `D_d` has a factor `x^8`. Its quotient has a strictly positive
Bernstein expansion on `[0,1]`. The degrees used are

```text
d:                 1   2   3   4   5   6   7   8   9  10  11
Bernstein degree: 328 240 360 376 392 440 424 440 456 472 488.
```

All scaled Bernstein coefficients are positive integers, with minimum
22 in every case. Degree six uses an elevation by 32; the other cases
need none. This proves the comparison for every finite nonnegative
scale and hence the uniform typed parameter floor `8/11`. Finally,
`2*8^2>11^2` proves that this exceeds `1/sqrt(2)`.

## 4. Exact replay and computational scope

```sh
make verify-joint-two-place-obstruction
```

The two-variable determinant computation substitutes the formal variable
by the 512th power of the local variable. Every determinant coefficient
has local degree less than 512, so its coefficients decode without
aliasing. The substitution is a ring homomorphism; fraction-free exact
division therefore computes the encoded genuine minors. The large
substitution order also keeps the selected minimum-valuation path equal
to the reference path.

The checker independently verifies all 50 decoded minors at sixteen
`F16` parameter pairs, using 800 direct Gaussian determinants. These
finite controls supplement the symbolic calculation. The pressure proof
uses only integer polynomials and binomial coefficients; numerical
optimization and field sampling are not assumptions of the theorem.

Sources:
[`check_joint_two_place_obstruction.py`](../tools/check_joint_two_place_obstruction.py)
and [`joint-two-place-obstruction.json`](../evidence/joint-two-place-obstruction.json).
The evidence retains both full minor expansions, all eleven integer
pressure polynomials, Bernstein coefficient counts/minima and their
fingerprints. This is ordinary proof with local author replay, without
new Lean formalization or independent mathematical acceptance.
Substantial OpenAI assistant support.
