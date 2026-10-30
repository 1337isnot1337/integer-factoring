# A score-preserving transfer for arbitrary vanishing width error

## Result and scope

The published score permits a logarithmic saving exponent with an `o(1)`
loss. That distinction weakens the quantitative crossing hypothesis.
For the existing grid recursion, an arbitrary `o(lambda)` error in the
total channel width preserves the score. A bound of
`O(lambda/loglog lambda)` is needed by the previous argument for its
stronger constant-factor runtime statement, but is unnecessary for the
score.

**Theorem.** Fix a prime `p`, a real number `1/2 <= rho < 1`, and
`delta > 0`. Suppose there is one uniform effective family of scale-one
crossing decompositions over `F_p`, with the precise interface below.
Suppose its widths and number of pieces satisfy

```math
\max_{i\in I_\lambda}(r_i+c_i)
 \le 4\rho\lambda+o(\lambda),
\qquad
\log_2\max(2,|I_\lambda|)
 \le \frac{\lambda}{(\log_2(2\lambda))^{\delta-o(1)}}.
\tag{1}
```

Then one deterministic finite-alphabet machine computes the least prime
factor of every binary integer `N >= 2`, and has score at least

```math
s=\min\{\log_2(1/\rho),\delta\}.
\tag{2}
```

Equivalently, for every fixed `epsilon > 0`, its worst-case time on
`n`-bit inputs is at most

```math
2^{\,n/(\log n)^{s-\epsilon}}
\tag{3}
```

for all sufficiently large `n`. The machine is the same for every
`epsilon`. Complete prime factorization has the same score.

This theorem does not supply a family with `rho < 1/sqrt(2)`, and does
not establish a factoring exponent above `1/2`. It removes a quantitative
error-rate obligation from future constructions. In particular,
`log |I_lambda| = O(lambda^beta)` for a fixed `beta < 1` places no
restriction on (2): every fixed `delta` is then available.

The score definition is pinned to the
[integer-factoring challenge data](https://github.com/bigospeedrun/site/blob/8ef5354ffd9cc448958c82cec54a1b5922399975/data/integer-factoring.json).
The underlying grid identities and grouping are Appendix A of
[Anthropic's evaluator at 16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The written factoring constructor is
[half-top-constructor.md](half-top-constructor.md), using Rudich and
Rousseau's source at
[c91cc3f](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
This note extends the ordinary
[volume-cutoff transfer at 1476525](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/subexponential-piece-transfer.md)
and uses its
[rough-input recovery theorem](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/rough-input-transfer.md).
Those pinned branch proofs, rather than a new formalization, are the
additional dependencies.

## 1. The full effective hypothesis

For every positive integer `lambda`, let

```math
E_\lambda(x,y,u,v)=[x=v][y=u],
\qquad x,y,u,v\in\{0,1\}^{\lambda}.
```

There is a finite index set `I_lambda` and tensors `A_i` such that
`sum_i A_i = E_lambda` over `F_p`. Each `A_i` has exact factorizations
through `2^r_i` and `2^c_i` coordinates for the row cut
`(x,y)|(u,v)` and column cut `(x,u)|(y,v)`, respectively. Assume:

1. `0 <= r_i,c_i <= 2lambda` for every index and every length.
2. A fixed program enumerates all indices, computes their integer widths,
   and evaluates any factor entry. Complete index enumeration and each
   individual query take at most `2^(C lambda)` bit operations for one
   fixed `C`. Index descriptions use at most `C lambda` bits.
3. The identities and algorithms are total and valid for every length.
   It is enough to handle finitely many exceptional lengths by any
   fixed effective exact family.
4. The two asymptotic conditions (1) hold uniformly over all pieces of
   each length. The entropy condition means that for every `eta > 0`,
   there is a constant `C_eta` and a threshold after which

   ```math
   f(\lambda):=\log_2\max(2,|I_\lambda|)
   \le C_\eta\frac{\lambda}
       {(\log_2(2\lambda))^{\delta-\eta}}.
   \tag{4}
   ```

The entry-query and enumeration costs are hypotheses, not consequences
of the ranks. They rule out advice containing uncomputed factorizations.
The factors need not have small workspace for this theorem.

The exact scalar identity extends to every commutative `F_p`-algebra.
This is required when evaluating the weighted divisor network in a
cyclic polynomial ring. No division by a possibly zero scalar occurs.

## 2. The algorithm does not use an asymptotic convergence rate

Fix a sufficiently large integer `x0 >= 18` once, depending only on
the family's fixed programs and constants. Use the integer schedule

```math
K_0=2^{x_0},\qquad
g_j=2^{3j+x_0}K_j,\qquad
K_{j+1}=g_j^2,\qquad
\zeta_j=2^{-j-5}.
\tag{5}
```

Choose the largest `ell` such that `P=K_ell` satisfies `P^40 <= n`.
Use trial division for the fixed finite initial range where no such
level exists. The schedule has the exact closed form

```math
\log_2 K_j=(3x_0+6)2^j-6j-2x_0-6.
\tag{6}
```

In particular,

```math
n^{1/400}<P\le n^{1/40},\qquad
\log n=\Theta(2^\ell),\qquad
\ell=O(\log\log n).
\tag{7}
```

The algorithm finds the level by integer multiplication and comparison.
The next schedule candidate has `O(log n)` bits because
`K_(j+1) <= K_j^10`. It never computes `rho`, `delta`, `epsilon`, or
a real-valued logarithm to choose a branch.

For each weighted divisor-prefix query, construct the same top clique as
in the cited transfer. Its information volume is
`V_ell=P^2 lambda_ell <= 4717n`; its complete generation and table-query
cost is `2^O(n/P)`. At a recursive call of level `j`, set
`V_j=K_j^2 lambda_j` and apply these exact rules:

- Stop and enumerate edge labels if `j=0` or `V_j <= ceil(sqrt(n))`.
- Otherwise enumerate crossing-piece assignments, choose the orientation
  with the smaller actual integer total, find a valid balanced grouping
  by the existing exhaustive search, and recurse at level `j-1`.

The square root is an integer square root. The cutoff uses the original
input bit length. Every call that performs a grid step consequently has

```math
\lambda_j=V_j/K_j^2>\sqrt n/P^2\ge n^{9/20}.
\tag{8}
```

All crossing identities, grouping procedures, and stopping evaluations
are exact at every input size. The asymptotic hypotheses affect their
time bound, not their correctness or choice of instructions.

## 3. A varying geometric envelope

Define the nonnegative relative width error

```math
a(\lambda)=\max\left(0,
 \frac{\max_i(r_i+c_i)}{4\lambda}-\rho\right).
\tag{9}
```

The separate caps imply `0 <= a(lambda) <= 1-rho`; (1) says
`a(lambda) -> 0`. For the analysis only, put

```math
u_n=\sup_{\lambda\ge n^{9/20}}a(\lambda),\qquad
\nu_n=u_n+n^{-19/40},\qquad
\beta_n=\rho+\nu_n.
\tag{10}
```

The supremum ranges over integer label lengths. Then `nu_n -> 0`.
For all sufficiently large `n`,

```math
\tfrac12\le\rho\le\beta_n
 \le\overline\rho:=(1+\rho)/2<1.
\tag{11}
```

At a grid step from level `j` to `j-1`, the source row/column identity
and balanced grouping give

```math
V_{j-1}\le
 (\rho+a(\lambda_j)+3\zeta_{j-1})V_j+K_{j-1}^2.
\tag{12}
```

This is the same grouping estimate as in the cited transfer, with the
actual total-width bound inserted. The caps `r_i,c_i <= 2lambda`
retain its degree and single-edge hypotheses. No new grouping theorem
is asserted.

Since `K_(j-1)^2 <= K_j <= P` and `V_j > sqrt(n)`,

```math
K_{j-1}^2/V_j\le n^{-19/40}.
```

Equations (8)--(12) therefore imply

```math
V_{j-1}\le(\beta_n+3\zeta_{j-1})V_j.
\tag{13}
```

For every visited path and every level reached on it,

```math
\begin{aligned}
V_j
&\le V_\ell\beta_n^{\ell-j}
 \prod_{t=j}^{\ell-1}\left(1+\frac{3\zeta_t}{\beta_n}\right)\\
&\le D\,V_\ell\beta_n^{\ell-j},
\qquad D=\exp\left(\frac3{16\rho}\right).
\end{aligned}
\tag{14}
```

Here `sum_(t>=0) zeta_t = 1/16`. The constant `D` is independent of
`n`, the branch, and the stopping level. From (11),

```math
\sum_{\text{visited }j}V_j
 \le \frac{D V_\ell}{1-\overline\rho}=O(n).
\tag{15}
```

Thus arbitrary vanishing width error does not introduce a depth factor
into the entropy charge. Its only relevant effect is to replace the
fixed contraction `rho` by `beta_n -> rho`.

## 4. Every source of time is charged

### Piece assignments

There are at most `K_j^2/4` crossing sites at a level-`j` call. The
logarithm of its number of assignments is at most

```math
\frac{K_j^2}{4}f(\lambda_j)
=\frac{V_j}{4}\frac{f(\lambda_j)}{\lambda_j}.
\tag{16}
```

Fix `0 < eta < delta`. Equations (4) and (8) bound (16) by
`O_eta(V_j/(log n)^(delta-eta))` at every expanded call, uniformly
over the path. Summing with (15) gives

```math
O_\eta\left(\frac n{(\log n)^{\delta-\eta}}\right).
\tag{17}
```

At a node, bound the sum of child costs by the number of assignments
times the maximum child cost. Taking logarithms and iterating bounds
the full recursion tree by charges along one worst path. This argument
does not identify the tree with a single branch.

### Local factors, index generation, and repeated table queries

The existing grouped-table routine takes at most

```math
2^{C_1 V_j/K_{j-1}} n^{C_2}
\tag{18}
```

parent queries and local operations. Its original local-state sum has
at most `2^(2V_j/K_(j-1))` terms. A factor call or complete index
enumeration costs at most `2^(C lambda_j)` by hypothesis; this fits
(18) because `lambda_j <= V_j/K_(j-1)`. The polynomial number of
sites, counters, and fixed-size coefficient manipulations also fit
(18). No table read or factor evaluation is treated as a free oracle.

Repeated parent queries multiply these overheads. By (14), their
total logarithmic contribution is at most

```math
\begin{aligned}
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
&\le D V_\ell\beta_n^\ell
 \sum_{j\ge1}\frac{\beta_n^{-j}}{K_{j-1}}\\
&\le \frac{D V_\ell\beta_n^\ell}{2^{17}-1}
=O(n\beta_n^\ell).
\end{aligned}
\tag{19}
```

Indeed `K_(j-1) >= 2^(18j)` and `beta_n >= 1/2`. The
`O(ell log n)` logarithmic indexing charge is lower order. The
original top-table query contributes `O(n/P)` to the exponent and
is included once in this multiplied query bound.

### Grouping and leaves

At level `j`, exhaustive grouping searches at most
`K_(j-1)^(g_(j-1))` colorings, where `g_(j-1)=sqrt(K_j)`.
Its logarithmic cost is at most `sqrt(P) log P`, apart from
polynomial indexing. Summing over the levels contributes at most
`O(n^(1/80)(log n)^2)`, which is `o(sqrt n)`.

A cutoff leaf costs `2^O(sqrt n)` arithmetic operations. A level-zero
leaf costs `2^O(V_0)=2^O(n beta_n^ell)` by (14). All inherited
table-query costs have already been included in (19).

Consequently the weighted arithmetic evaluator, including uniform
generation and structural work, has size and bit-operation cost over
the fixed base field at most

```math
A(n)\le
2^{\,O_\eta\left(
n\beta_n^\ell+
\frac n{(\log n)^{\delta-\eta}}+
\sqrt n+\frac nP
\right)}.
\tag{20}
```

This also counts an arithmetic execution over an arbitrary commutative
`F_p`-algebra before charging the representation and arithmetic of
one algebra element. A fixed finite-tape simulation changes the
constant in the exponent.

## 5. Bit complexity and factor recovery

For completeness, the rough-input reduction supplies the following
explicit final step. Trial-divide through `2^b`, where
`b=ceil(sqrt n)`, and return the first divisor if one exists. This
takes `2^O(sqrt n)` bit operations. Otherwise every prime factor is
larger than `2^b`, so

```math
\Omega(N)<n/b,\qquad \tau(N)\le2^{\Omega(N)}\le2^b.
\tag{21}
```

The weighted prefix polynomial, with divisor one removed, consequently
has at most `2^b` terms. Its exponents have at most `n` bits. The cited
sparse nonzeroness argument tests it in
`F_p[X]/(X^q-1)` for a deterministic family of primes `q`. Both the
number of rings and their largest dimension are `2^O(sqrt n)`.
Their generation, all coefficient arithmetic, and sequential tape
access have the same exponential bound up to fixed powers.

One nonzero term is isolated by at least one of the cyclic reductions:
for a chosen support exponent, every failing prime divides a difference
from another support exponent, and their product has fewer than
`n(2^b-1)` distinct prime divisors. This is an exact test and permits
ordinary lexicographic prefix recovery. The polynomial retains the term
for `N`, so prime inputs also terminate correctly. At most `n` prefix
tests return the smallest divisor greater than one, which is prime.

Evaluating (20) over these represented rings and performing all
recovery queries therefore takes

```math
T(n)\le 2^{O(\log_2(A(n)+2)+\sqrt n)}.
\tag{22}
```

This charges the coefficient-vector dimension outside the recursive
arithmetic-operation count. There is no unit-cost extension-field or
unbounded-integer assumption. Repeated division and least-factor calls
give complete factorization with a polynomial overhead.

## 6. Extracting the score

Put `alpha=log2(1/rho)` and
`alpha_n=log2(1/beta_n)`. Then `alpha_n -> alpha`. The fixed
schedule (7) gives, with a constant uniform in large `n`,

```math
\beta_n^\ell
 \le \frac{C}{(\log n)^{\alpha_n}}
 =\frac1{(\log n)^{\alpha-o(1)}}.
\tag{23}
```

Also `n/P <= n^(399/400)`, so `sqrt n+n/P` is smaller than
`n/(log n)^t` for every fixed `t`. Equations (20)--(23) show that,
for every sufficiently small fixed `eta > 0`,

```math
\log_2T(n)
 \le C_\eta\frac n{(\log n)^{s-\eta}},
\qquad s=\min(\alpha,\delta),
\tag{24}
```

eventually. Given the score tolerance `epsilon`, use
`eta < epsilon/2` in the analysis. The remaining factor
`(log n)^(epsilon-eta)` eventually exceeds `C_eta`; hence (24)
implies (3). Neither `eta` nor `epsilon` is an input or a parameter
of the machine.

## 7. Computability and the stronger constant-factor bound

No computable modulus of convergence for (1) is required. The machine
does not query the tail supremum in (10), use it to select a level, or
branch according to a threshold supplied by an existence proof. The
exact all-length programs, fixed schedule, actual-width orientation,
and volume cutoff already define one total machine. The tail supremum
is an analytical upper bound on that machine. Unknown finite
thresholds enter only the constants and eventual thresholds in (24).

If an explicit computable tail envelope `u(L)` is available, with
`a(lambda) <= u(L)` for every `lambda >= L`, one can replace `u_n`
by `u(ceil(n^(9/20)))`. This gives a quantitative version without
changing the algorithm. In particular,

```math
\beta_n^\ell
\le \rho^\ell\exp(\nu_n\ell/\rho),
\qquad
\alpha-\alpha_n
\le\frac{\nu_n}{\rho\ln2}.
\tag{25}
```

When `a(lambda)=O(1/loglog lambda)`, the exponent
`nu_n ell` is bounded. With the stronger entropy bound
`f(lambda)=O(lambda/(log lambda)^delta)`, this recovers
`T(n)=exp(O(n/(log n)^s))`.

For arbitrary `a(lambda)=o(1)`, that stronger conclusion does not
follow from this recurrence. For example, consider the scalar recurrence
of depth `ell=m^2` with contraction `rho+1/m`. The ratio to the
zero-error envelope is

```math
\left(1+\frac1{\rho m}\right)^{m^2}
 \ge 1+\frac m\rho\longrightarrow\infty,
```

while the loss in its logarithmic saving exponent is
`log2(1+1/(rho m)) -> 0`. This illustrates the exact distinction:
the score is preserved, but a bounded multiplicative error in the
runtime exponent need not be. It is a counterexample for an abstract
recurrence estimate, not a lower bound on a factoring algorithm.

## 8. Verification

Run the accompanying exact-arithmetic checker:

```sh
python3 tools/check_vanishing_slack_transfer.py
```

It checks the source schedule identity, integer cutoffs, rounded volume
recurrences with deliberately slow vanishing errors, a finite geometric
envelope that also bounds its infinite tail, the query-overhead series,
and the exact rational recurrence illustrating the score/runtime
distinction. It reports the scope and counts in JSON. These are finite
checks of the recurrence arguments, not new tensor constructions or an
execution of the asymptotic factoring machine. The all-input result is
the ordinary proof above. No Lean files are changed.
