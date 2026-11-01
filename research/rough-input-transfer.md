# Small-factor preprocessing improves the modular transfer

## Result

Trial division before a weighted divisor query gives a stronger support
bound than a worst-case estimate for the divisor function. This reduces
the nonzeroness and recovery overhead to `exp(O(sqrt(n)))`. With a
different choice of the same parameter, the coefficient rings have
polynomial dimension and the established half-rate evaluator can be
implemented in polynomial space.

The argument is elementary. It uses the weighted network and cyclic
nonzeroness test already proved in
[scale-one-modular-transfer.md](scale-one-modular-transfer.md) and
[sparse-divisor-nonzeroness.md](sparse-divisor-nonzeroness.md).
The arithmetic evaluator remains a separate hypothesis in the general
transfer theorem.

**Rough-input transfer theorem.** Fix a prime `p`. Suppose weighted
divisor-prefix polynomials for inputs of at most `n` bits have uniformly
generated arithmetic circuits over `F_p` of size and construction cost
at most `A(n)`, correct over every commutative `F_p`-algebra.
For every integer `1<=s<=n`, least-prime-factor recovery and complete
factorization have deterministic bit complexity

```math
\exp\!\left(O\!\left(
\log(A(n)+2)+s+\frac ns+\log(n+2)\right)\right).
\tag{1}
```

The constants do not depend on `s,n` or the input integer.
Taking `s=ceil(sqrt(n))` gives

```math
\exp\!\left(O\!\left(\log(A(n)+2)+\sqrt n\right)\right).
\tag{2}
```

Consequently, if `A(n)=exp(O(n/(log n)^alpha))` for any fixed
`alpha>0`, the factoring time has that **same exponent `alpha`**.
The previous additive `n/log n` term and its artificial cap at
`alpha=1` are not needed for this transfer.
This removes a transfer limitation; it does not itself construct an
evaluator with exponent above one half.

**Half-rate factoring in polynomial space.** There are fixed constants
`C,c>0` and one deterministic finite-alphabet machine which, for every
`n>=2` and every `2<=N<2^n`, computes its least prime factor and complete
prime factorization using time

```math
C\exp\!\left(\frac{Cn}{\sqrt{\log(n+2)}}\right)
```

and at most `C n^c` work-tape cells. This is an ordinary mathematical
consequence of the existing scale-one modular crossing, the streamed
local factors in [modular-factor-space.md](modular-factor-space.md),
and Sections 4–6 below. The half endpoint is proved explicitly in
Section 6. No new Lean formalization is used.

## 1. No small prime factor means few divisors

Let `2<=N<2^n`. Test the integers
`d=2,3,...,min(2^s,N)` in increasing order, returning the first
divisor found. The first divisor is prime and is the least prime factor.
This phase uses at most `2^s poly(n)` bit operations and polynomial
space: maintain the candidate and apply ordinary integer division.

Suppose this phase finds no divisor. Every prime factor of `N` is
strictly larger than `2^s`. Write

```math
N=\prod_p p^{a_p},\qquad
\Omega(N)=\sum_p a_p.
```

Then

```math
s\Omega(N)<\log_2N<n,\qquad
\tau(N)=\prod_p(a_p+1)\le2^{\Omega(N)}
\le 2^{\lceil n/s\rceil}.
\tag{3}
```

The first inequality is strict because each contributing prime exceeds
`2^s`. The divisor bound uses only `a+1<=2^a` for integers `a>=1`.
The algorithm does not know the prime factorization; it uses the
consequence (3) after certifying that there is no small divisor.

Put

```math
q=\left\lceil\frac ns\right\rceil,\qquad
K=2^q,\qquad M=n(K-1)+1,\qquad Q=4(M+1)^2.
\tag{4}
```

Every divisor-prefix polynomial, with the term for divisor one removed,
has at most `K` nonzero terms. Its exponents remain less than `2^n`.
The support bound is uniform over every prefix subsequently queried.

## 2. The same finite family of cyclic rings is sufficient

The sparse nonzeroness proof applies with this new `K`. For clarity,
fix one support exponent `d0` of a nonzero polynomial `G`. If
`G modulo (X^r-1)` vanishes for a prime `r`, some other support
exponent collides with `d0` modulo `r`. Thus `r` divides

```math
\prod_{d\ne d_0}|d-d_0|<2^{n(K-1)}.
```

At most `n(K-1)` distinct primes can divide this positive product.
Consequently any `M` distinct odd primes test nonzeroness exactly:

```math
G\ne0\quad\Longleftrightarrow\quad
\exists r\text{ in the list},\
G\bmod(X^r-1)\ne0.
\tag{5}
```

The elementary central-binomial argument in the companion note supplies
at least `M` odd primes at most `Q`. Enumerating candidates and testing
divisibility costs a fixed polynomial in `Q`. For `n>=2`,

```math
M+1=n(K-1)+2\le nK,\qquad
Q\le4n^2\,2^{2q}.
\tag{6}
```

Hence the number of queries, their maximum cyclic-ring dimension, and
the cost of generating their moduli are all
`exp(O(n/s+log(n+2)))`.

### A sharper elementary modulus bound

For the space calculation, the same prime-supply argument gives the
stronger bound

```math
Q_1=8M\left\lceil\log_2(M+2)\right\rceil.
\tag{6a}
```

There are at least `M` odd primes at most `Q1`. To verify this,
put `ell=ceil(log2(M+2))>=2`. For `M>=2`,

```math
\log_2Q_1\le2\ell+2,\qquad
Q_1-1-(M+2)(2\ell+2)
=(6M-4)\ell-2M-5
\ge10M-13>0.
```

The central-binomial estimate from the companion note is
`pi(v)>=(v-log2(v+1))/log2(v)` for even `v>=2`.
Using `log2(Q1+1)<=log2(Q1)+1` therefore gives
`pi(Q1)>M+1`. After discarding prime two there are enough odd
primes. For `M=1`, `Q1=16` and prime three suffices.

The algorithm may use `min(Q,Q1)` as its candidate limit. Both bounds
are proved independently, so taking the smaller one preserves prime
supply. In particular, the maximum needed ring dimension is

```math
O\!\left(n2^q(\log(n+2)+q)\right).
\tag{6b}
```

This does not require the quotient ring to be a field. Represent an
element of `F_p[X]/(X^r-1)` by its `r` coefficients. Addition,
multiplication by cyclic convolution, equality testing, and sequential
tape access have fixed polynomial cost in `r`. The prime `p` is
fixed, so coefficients have constant bit length. Schoolbook
multiplication is sufficient.

The weighted input for bit position `j` is the monomial
`X^(2^j mod r)`. Repeated doubling computes its exponent; no
polynomial of degree `2^n` is materialized.

## 3. Recover the least factor, including prime inputs

Use the original prefix convention: `k` low bits remain free and
the upper prefix is `x`. Define

```math
G_{N,k,x}(X)=
\sum_{\substack{d\mid N,\ d>1\\\lfloor d/2^k\rfloor=x}}X^d.
\tag{7}
```

The weighted network first computes the same sum including `d=1`.
Subtract `X` exactly when the prefix contains one.
The term `X^N` is retained, so the initial prefix `k=n,x=0`
always has a permitted divisor.

At each bit, test the zero child using (5). If it is nonempty, take
it; otherwise take the one child. The current prefix remains nonempty,
and choosing zero whenever possible gives its smallest permitted
integer at the leaf. That integer is the smallest divisor greater than
one, hence the least prime factor. If `N` is prime, the algorithm
returns `N`.

There are `n` prefix tests. Each test uses at most `M` cyclic
evaluations and one subtraction per modulus. Uniform circuit generation,
ring evaluation, prime generation, and ordinary indexing therefore cost
a fixed polynomial in `A(n),Q,n`. Combining this with the trial phase
proves (1).

For complete factorization, divide by the returned prime and repeat.
There are at most `n` calls. The same original threshold `2^s` may
be used in every call. A smaller remaining input either has a divisor
at most this threshold, or still satisfies the same `K` bound.
The definition of `A(n)` already covers all input lengths at most
`n`, so the same time bound applies.

For (2), put `s=ceil(sqrt(n))`. Then
`ceil(n/s)<=s`, so the two preprocessing terms are `O(sqrt(n))`.
For every fixed `alpha>0`,
`sqrt(n)=o(n/(log n)^alpha)`, which proves the claimed uncapped
logarithmic-power transfer.

## 4. An explicit time-space parameter

For space bounds, it is useful to choose an integer `1<=q<=n` first and
set

```math
s=\left\lceil\frac nq\right\rceil.
\tag{8}
```

Then `ceil(n/s)<=q`, so one may use `K=2^q` in (4).
Suppose, additionally, that the arithmetic evaluator can be executed
with polynomially many bits of control state and polynomially many
simultaneously live ring elements, while making at most `A(n)`
operations and queries. Suppose its scalar-factor and table-query
routines have the same polynomial-space property.

The same factoring algorithm then has

```math
\begin{aligned}
\text{time}
&\le
\exp\!\left(O\!\left(
\log(A(n)+2)+\frac nq+q+\log(n+2)\right)\right),\\
\text{space}
&\le n^{O(1)}\,2^{O(q)}.
\end{aligned}
\tag{9}
```

Prime candidates can be generated sequentially, and each cyclic
evaluation can be discarded before the next modulus is considered.
There is no need to store the entire prime list or all queried outputs.
Integer counters need only `O(n+q+log(n+2))` bits.
Equation (6) bounds the representation of one ring element by
`O(n²2^(2q))` field coefficients. This proves the space claim,
including a fixed polynomial number of live elements.

### Generic polynomial-space transfer

Choose `q=floor(log2(n))` for `n>=2`. Then `q=Theta(log n)`,
`K<=n`, and (6b) gives ring dimension `O(n² log(n+2))`.
Under the streamed evaluator hypothesis, least-factor recovery and
complete factorization therefore use polynomial space and time

```math
\exp\!\left(O\!\left(
\log(A(n)+2)+\frac n{\log(n+2)}\right)\right).
\tag{9a}
```

Thus every fixed logarithmic saving exponent `0<alpha<=1` is
preserved in polynomial space. The fastest abstract recovery choice
`s=ceil(sqrt(n))` in (2) is a different point of the tradeoff: it
permits larger coefficient rings in order to minimize the extra time.

### The half-rate choice

For `n>=2` choose, using exact integer arithmetic,

```math
q=\left\lfloor\sqrt{\lfloor\log_2n\rfloor}\right\rfloor,
\qquad s=\left\lceil\frac nq\right\rceil.
\tag{10}
```

Here `q>=1`, `s=O(n/sqrt(log n))`, and `q<=log2(n)`.
Consequently

```math
K=2^q\le n,\qquad M\le n^2,\qquad Q\le4n^4.
\tag{11}
```

The sharper bound (6b) improves the maximum ring dimension to
`n^(1+o(1))`: indeed `q=O(sqrt(log n))`,
`M<=n2^q=n^(1+o(1))`, and its additional logarithm is
`n^o(1)`. This choice keeps the trial phase within half-rate time
while using smaller rings than the generic polynomial-space choice.

If the evaluator has time `exp(O(n/sqrt(log n)))` and the
streamed-space property above, the result is therefore

```math
\text{time}=\exp\!\left(O\!\left(\frac n{\sqrt{\log n}}\right)\right),
\qquad
\text{space}=n^{O(1)}.
\tag{12}
```

Unlike a generic circuit-size assertion, the streamed-space hypothesis
must be proved for the evaluator. Writing the whole circuit or all of
its tables would not establish (12). The next section gives the needed
execution for the grid evaluator, subject to polynomial-space local
factor routines.

## 5. Execute the grid evaluator without storing its tables

Use the explicit schedule, stopping rule, and volume bounds from
[subexponential-piece-transfer.md](subexponential-piece-transfer.md).
At the top, `P=K_ell<=n^(1/40)` and the information volume is
`V_ell=P²lambda<=4717n`. Along every visited path,

```math
V_j\le D\rho^{\ell-j}V_\ell=O(n),\qquad
\ell=O(\log\log n).
\tag{13}
```

Store a table as a procedure with its integer descriptors. A query to a
child table enumerates its grouped local settings sequentially, forms
one product at a time by querying parent tables and scalar factors, and
adds the result to one ring accumulator. It does not construct a row of
the table. A base-case clique similarly enumerates one complete edge
labelling at a time, multiplies its node values, and adds to an
accumulator.

At one grid level the stored descriptors consist of:

- The current piece index at each crossing site, its row and column
  widths, and the chosen orientation.
- One grouping code, found by a sequential search through the source's
  finite coloring space.
- The current local-setting, edge-label, and product-loop counters.
- The table-query arguments and a constant number of ring accumulators
  for each active loop or subroutine.

All of these have polynomial length in `n`. In particular, there
are at most `K_j²/4` crossing sites, each piece index has
`O(lambda_j)` bits, and their total is `O(V_j)=O(n)`.
The integer widths, group maps, and legal local settings also have
polynomial descriptions. Grouping tests depend only on these integer
data, so they need no array of ring values.

The structural recursion has at most `ell` frames. A table query
calls successive parent-table procedures at most once per level before
reaching an original top-table query. Thus nested calls use only a
polynomial number of frames and stored ring elements. The local factor
routines must themselves use polynomial space; this is an additional
requirement beyond their exponential-time bound.

The top table is also evaluated by a sequential sum over original-gate
settings. The construction in
[half-top-constructor.md](half-top-constructor.md), Sections 3–5, gives
`O(n)` bits per setting and polynomial-time fixed-setting checks.
Its weighted version adds a cyclic monomial. A query therefore uses
polynomial space without writing the exponential table.

### Recalculation is already charged in the time bound

This execution does not rely on memoization. At a step with parent
volume `V_j`, a grouped-table query costs at most

```math
2^{C V_j/K_{j-1}}\,n^{C}
\tag{14}
```

parent queries and local operations, including the `2^O(lambda_j)`
factor routine. The volume analysis explicitly multiplies these
query costs along a path. It proves

```math
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
\le D V_\ell\rho^\ell
\sum_{j\ge1}\frac{\rho^{-j}}{K_{j-1}}
=O(n\rho^\ell).
\tag{15}
```

The convergent series in (15) is the reason repeated queries do not add
an extra exponential factor. Branching over piece assignments is also
already included by multiplying their counts along a worst path.
The original top-table query contributes `2^O(n/P)`, which is
smaller than the half-rate bound. A polynomial-time sequential-tape
implementation of the descriptors changes only the constant in that
exponent.

Thus the existing operation-count proof applies to the streamed
execution as well as to a stored circuit. Once the scalar factors have
`2^O(lambda)` time and `poly(lambda)` space algorithms, the exact
half-rate contraction has the space property required in (12).

## 6. Local factor interface for the unconditional half-rate consequence

### An exact half endpoint from the existing moment estimates

The modular construction already gives more quantitative control than
its fixed-epsilon corollary states. Write `a=A_2=pi²/24` and
`L=lambda ln(2)`. Equations (11)–(12) of
[modular-crossing.md](modular-crossing.md) give, for every sufficiently
small `t>0`, some integer weight `H` with

```math
\ln q_2(H)\ge\frac{2a}{t}
-C_0\bigl(t^{-3/4}+\ln(1/t)\bigr),
\qquad
H=\frac a{t^2}+O(t^{-7/4}),
\tag{16}
```

and the uniform cumulative bound

```math
\ln\sum_{j\le E}q_2(j)
\le2\sqrt{aE}+O(\ln(E+2)).
\tag{17}
```

Choose a sufficiently large fixed constant `C1` and set

```math
\frac{2a}{t}=L+C_1L^{3/4}.
```

For all sufficiently large `lambda` the positive error term on the
right dominates the loss in (16), so `q_2(H)>=2^lambda`. Moreover,

```math
H=\frac{L^2}{4a}+O(L^{7/4}).
```

The least permitted label energy `h` is at most this `H`. This
argument does not assume monotonicity of the individual counts.
Together with (17) and the explicit `h=O(lambda²)` bound, it gives

```math
L-O(\log(\lambda+2))
\le2\sqrt{ah}
\le L+O(\lambda^{3/4}).
\tag{18}
```

The two channel spaces have the source dimensions
`sum_(i=0)^(2h) q_2(s+i)` and
`sum_(i=0)^(2h) q_2(2h-s-i)`. By (17), their bit widths obey

```math
\begin{aligned}
r_s+c_s
&\le \frac{2\sqrt a}{\ln2}
\bigl(\sqrt{2h+s}+\sqrt{2h-s}\bigr)+O(\log(h+2))\\
&\le2\sqrt2\,\lambda+O(\lambda^{3/4}),
\qquad |s|\le2h.
\end{aligned}
\tag{19}
```

Cap each separate width at `2lambda` by the direct factorization if
needed. The additive error in (19) is
`O(lambda/loglog(lambda+4))`. The piece count is polynomial in
`lambda`, so its entropy satisfies the relaxed transfer with, for
example, `delta=1`. Thus that transfer applies with the **exact**
value `rho=1/sqrt(2)` and yields
`A(n)=exp(O(n/sqrt(log n)))`. No differentiated asymptotic for
`q_2(h)`, monotonicity of `q_2`, or limiting choice of `rho` is
needed.

### Polynomial-space factor values

The characteristic-two quotient crossing family in
[modular-crossing.md](modular-crossing.md) uses regular partitions:
their parts are odd and distinct. A regular partition of weight at most
`4h` has at most `2sqrt(h)` parts, because the first `m` odd
positive integers sum to `m²`. The label energy satisfies
`h=O(lambda²)`.

The memoized coefficient routine in that note is enough for time but
does not by itself give polynomial space. Its replacement is a streamed
normal-order coefficient formula: enumerate the common unaffected
subsets of input and output modes; for each subset, compute the
corresponding permanent of an `O(lambda)`-order matrix by a sequential
subset formula. In characteristic two the permanent also equals the
determinant. Angular coefficients use polynomial entries of degree
`O(h)`, followed by one coefficient extraction.

Channel indices and regular partitions can be enumerated sequentially
using partition-count dynamic programs through weight `4h`. Those
programs have polynomially many cells, each with `O(lambda)` count bits.
An exact piece value can be obtained by summing its original
factorization over the intermediate channel sequentially. Consequently
the direct `2lambda` width cap also preserves the time-space interface;
no large Gaussian-elimination table is required for that cap.

The complete local proof, including the formulas and exact checks, is
recorded in [modular-factor-space.md](modular-factor-space.md). With that
polynomial-space routine, Sections 4–5 yield the unconditional time-space
bound (12) from the existing exact half-rate modular family.

This conclusion improves the implementation resource bound. It makes
no claim of a logarithmic saving exponent above `1/2`.

## 7. Reproduction

The checker `tools/check_rough_input_transfer.py` tests the residual
support bound, the exact parameter inequalities, prefix recovery after
small-factor preprocessing, cancellation in individual cyclic rings,
and a control showing why the improved support bound cannot be used
before the small-factor test. These finite tests supplement the
uniform proof.
