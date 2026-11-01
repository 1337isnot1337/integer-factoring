# A Plucker change to the high-axis quadratic flag

## Research question and exact changed assumption

The [balanced monomial nine-axis kernel](high-axis-cayley-realization.md)
is closed by a rank-120 quadratic augmentation minor. Can a homogeneous
symmetric-polynomial coefficient make this quadratic flag decomposable,
while preserving the antipodal endpoint and almost balanced axis orders?
This changes the coefficient relations that made that minor nonzero.
It is not a refinement of the excluded monomial profile.

Fix nine axes and degree **64**. For each nonempty subset S of size k,
write `64=k*e+r`, `0<=r<k`, and set

```
c_S=(product_(a in S) x_a)^e * elementary_r((x_a)_(a in S)).
```

For k=1,...,9, the remainders are `0,0,1,0,4,4,1,0,1`.
Every binomial coefficient `binom(k,r)` is odd. All c_S are homogeneous
of degree 64 and evaluate to one at the all-one point. With
`c_empty=sum_(S nonempty)c_S`, the same Cayley/Pfaffian construction
therefore has zero diagonal and the exact antipodal endpoint.
Every axis crossblock still squares to `x_a^128 I`.

For a triple `{a,b,c}` the coefficient is
`(x_a*x_b*x_c)^21*(x_a+x_b+x_c)`. At the a-axis, its leading term is
`x_a^21*x_b^21*x_c^21*(x_b+x_c)`. Consequently the quadratic part of
that jet is the product of the two linear forms

```
(sum_(b!=a) x_b^22*epsilon_b) * (sum_(b!=a) x_b^21*epsilon_b).
```

These forms are independent over the other-variable rational function
field. Multiplication by their product has rank **64**, rather than 120.
Thus the previous nonzero quadratic-minor certificate does not apply.
This is only a changed flag; all higher channels and degree costs remain.

The registered first diagnostic computes complete specialized profiles
for exactly three nine-axis patterns of the other variables: all one;
one omega and seven ones; and the repeating sequence
`(1,omega,omega^2,1,omega,omega^2,1,omega)`. Keep every result, including
unresolved queries. Specialization gives optimistic generic orders,
never a generic positive certificate. If these profiles already fail,
seek an exact structural obstruction or pivot; do not extend the field
or sample catalogue without a reason. A favorable profile requires
symbolic determinantal divisibility/factorization and rational outward
pressure before any complete-bound claim. The target remains `rho<1/3`
with the full effective nine-axis transfer; score is still `1/2`.

This uses the credited homogeneous Cayley Gaussian source, multiaxis
transfer and Smith tools. The decomposable quadratic identity is ordinary
algebra. Directed by Cinna Davis / 1337isnot1337 with substantial OpenAI
assistant support; no new Lean or external novelty is claimed.


## Completed diagnostic and all-scale stopping certificate

All three [registered specializations](../evidence/high-axis-plucker-specializations.json)
resolved all 256 orders. Their numerical minima are `0.4217827364`,
`0.4250254489` and `0.4257603101`, respectively. These values are
optimistic specialized fits, with no generic positive conclusion.

The [independent checker](../tools/check_high_axis_plucker.py) confirms
that the quadratic multiplication rank is 64 at a nontrivial F4
specialization. However, the higher jets still exclude this kernel.
This follows without trusting the large native profiles.

Specialize the other eight variables to one. The matrix is then binary.
Let R(L) be the rank over F2 of its multiplication map modulo `x^L`,
a literal block Toeplitz matrix of order `256L`. If its complete Smith
orders are s_i, the exact identity is

```
R(L)=sum_i max(L-s_i,0).
```

The checker constructs the binary rows directly and eliminates them
using integer bit operations, independently of polynomial Smith code.
It finds **R(25)=742** and **R(26)=836**. Hence exactly **94** orders
are at most 25, and their sum is

```
26*(836-742)-836 = 1608.
```

The square identity `B^2=x^128 I` proves the generic rank needed here;
no claim that unresolved later pivots vanish is made. Specialization
cannot decrease the generic determinantal valuations, so every axis
of the original multivariate kernel has `Delta_94<=1608`.
The coefficients are binary, and the rank pivots persist over every
characteristic-two extension field.

Since `94>256/3` and `1608/94<64/3`, convexity and monotonicity give

```
log Z_a(t) >= 94*log(1+t^(1608/94))
           >= (256/3)*log(1+t^(64/3)),  0<t<1.
```

Summing the nine cuts and applying the binary entropy inequality with
`u=t^(64/3)` proves **`R_9(p,t)>=1/3` at every positive scale and
common occupied frequency**. This closes this kernel's degree-cap
pressure calculation. It does not exclude actual coefficient-rank
cancellation, other physical restrictions, or different high-axis kernels.

The [receipt](../evidence/high-axis-plucker-controls.json) preserves all
26 Toeplitz ranks, endpoint/degree checks, 6,561 independently expanded
nonzero-pattern entries and 56 independent untruncated small orders.
The changed quadratic flag was real; it was insufficient because the
rest of the low-order profile still supplies a stronger prefix bound.
The registered stopping condition is met. Do not continue this source by
varying degree or sampling more field points without a new mechanism.
The complete factoring score remains `1/2`.

```
make verify-high-axis-plucker
```
