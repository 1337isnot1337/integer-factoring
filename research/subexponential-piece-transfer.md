# A transfer theorem with more than polynomially many crossing pieces

## Result

The grid evaluator can use substantially more than polynomially many
crossing pieces. A direct recursion with an information-volume cutoff
controls their enumeration cost.

Fix a finite field `F_p` and constants `1/2<=rho<1`, `delta>0`, and
`C>=1`. Suppose that, for each positive integer `lambda`, a scale-one
crossing family has:

1. An exact crossing sum and both exact factorizations over `F_p`.
2. Separate channel widths at most `2lambda`.
3. Total channel width at most

   ```math
   r_i+c_i\le4\rho\lambda+
   \frac{C\lambda}{\log_2\log_2(\lambda+4)}.
   \tag{1}
   ```

4. Piece entropy bounded by

   ```math
   f(\lambda):=\log_2\max(2,|I_\lambda|)
   \le\frac{C\lambda}{(\log_2(2\lambda))^\delta}.
   \tag{2}
   ```

5. Uniform algorithms that enumerate the index set, compute its widths,
   and evaluate every factor in at most `2^(C lambda)` steps.
   Indices have at most `C lambda` bits. Complete enumeration, rather
   than an oracle for an unspecified index set, is included in this bound.

The bounds may instead hold for all sufficiently large `lambda`, with
uniform total algorithms at the finitely many remaining lengths.

Put `alpha=log2(1/rho)`. Then a deterministic finite-alphabet machine
computes the least prime factor, and the complete prime factorization,
of every `n`-bit integer in time

```math
\exp\!\left(
O\!\left(
\frac n{(\log n)^\alpha}
+\frac n{(\log n)^\delta}
+\frac n{\log n}
\right)\right).
\tag{3}
```

Thus the exponent is `min(alpha,delta,1)`.
There is no extra `log log n` factor when `delta=alpha`.

In particular, `log |I_lambda|=O(lambda^beta)` for any fixed
`beta<1` preserves the original `rho`-derived exponent.
The width slack in (1) can also be replaced by any smaller bound,
including `O(lambda^gamma)` with `gamma<1`.
Polynomially many pieces and logarithmic width slack are sufficient,
but are stronger hypotheses than this argument needs.

This is a conditional extension of the evaluator interface.
It does not construct a family with `rho<1/sqrt(2)`.
The exact grid identities, grouping, and scale-one base change are the
same as in [scale-one-modular-transfer.md](scale-one-modular-transfer.md).
The proof below changes the recursion analysis and adds an explicit cutoff.

## 1. Schedule and the exact cutoff

Use the source's schedule, with one fixed integer `x0>=18`:

```math
K_0=2^{x_0},\quad
\zeta_j=2^{-j-5},\quad
g_j=2^{3j+x_0}K_j,\quad
K_{j+1}=g_j^2.
\tag{4}
```

The source proves

```math
\log_2K_j=(3x_0+6)2^j-6j-2x_0-6,
\qquad K_{j+1}\le K_j^{10}.
\tag{5}
```

Choose the largest `ell` with `P:=K_ell` satisfying `P^40<=n`.
The finitely many input lengths without such a level can use trial division.
Then

```math
n^{1/400}<P\le n^{1/40},\qquad
\rho^\ell=\Theta((\log n)^{-\alpha}),\qquad
\ell=O(\log\log n).
\tag{6}
```

The [written top constructor](half-top-constructor.md) supplies a
`(P,lambda)` clique of information volume `P^2lambda<=4717n`.
Its construction cost is `2^O(n/P)`. The weighted version has the same
geometry, as proved in the scale-one transfer note.

At any recursive call with level `j`, label length `lambda_j`, and
information volume

```math
V_j=K_j^2\lambda_j,
```

apply this algorithmic rule:

- If `j=0` or `V_j<=ceil(sqrt(n))`, evaluate the clique by exhaustive
  enumeration of its edge labels.
- Otherwise perform one grid step, enumerate its piece assignments,
  choose the smaller actual row or column total, find the source's
  balanced grouping, and recurse at level `j-1`.

The cutoff uses the original input length `n`.
The comparison is exact: compute `ceil(sqrt(n))` with integer arithmetic.
The algorithm uses neither a numerical approximation to `rho` nor
an estimate of the runtime to choose its orientation.

Every call that performs a grid step has

```math
\lambda_j
=\frac{V_j}{K_j^2}
>\frac{\sqrt n}{P^2}
\ge n^{9/20}.
\tag{7}
```

This explicit lower bound is the reason the larger index sets are safe.
No asymptotic lower bound on the actual descendant volume is assumed.
A branch whose volume becomes small is handled by the cutoff.

## 2. One grid step with the new width slack

At level `j>=1`, put `K=K_(j-1)` and `g=g_(j-1)`, so
`K_j=g^2` and `V=g^4lambda`.
There are at most `g^4/4` crossing sites.

The source's row/column calculation gives, under (1),

```math
U\le
\left(\rho+\frac1g+
\frac{C}{4\log_2\log_2(\lambda+4)}\right)V.
```

The separate `2lambda` width caps retain the degree and individual-edge
bounds needed by its grouping lemma. Its case (ii), with
`zeta=zeta_(j-1)`, therefore gives a child width `L` satisfying

```math
K^2L\le
\left(\rho+3\zeta_{j-1}+
\frac{C}{4\log_2\log_2(\lambda+4)}\right)V+K^2.
\tag{8}
```

This is the source's display (42), with the new explicit slack substituted.
The rounding term is retained.
Since `K^2<=K_j<=P` and the parent has `V>sqrt(n)`,

```math
\frac{K^2}{V}\le n^{-19/40}.
\tag{9}
```

For sufficiently large `n`, (7) also gives

```math
\log_2\log_2(\lambda+4)
\ge\tfrac12\log_2\log_2 n.
```

Thus every visited grid step obeys

```math
V_{j-1}\le
(\rho+3\zeta_{j-1}+\varepsilon_n)V_j,
\qquad
\varepsilon_n=
\frac{C}{2\log_2\log_2 n}+n^{-19/40}.
\tag{10}
```

The estimates apply to every piece assignment.
They also apply to a child that immediately reaches the cutoff.

The primary source is
[Appendix A, Lemma 8 and displays (41)--(44)](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf)
at `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.
Only its exact grouping construction is used here.
The original polynomial-piece runtime theorem is not being applied
outside its hypotheses.

## 3. Volumes form a uniformly bounded geometric envelope

For any visited path from level `ell` to level `j`, (10) gives

```math
V_j
\le V_\ell\rho^{\ell-j}
\prod_{t=j}^{\ell-1}
\left(1+\frac{3\zeta_t+\varepsilon_n}{\rho}\right).
```

The product is bounded by

```math
\exp\!\left(
\frac3\rho\sum_{t\ge0}\zeta_t+
\frac{\ell\varepsilon_n}{\rho}\right).
```

The first sum is `1/16`.
Equations (5)--(6) give `ell<=log2log2 n` for large `n`,
so `ell epsilon_n=O(1)`.
There is consequently a fixed constant `D`, independent of the input,
the path, and the stopping level, such that

```math
V_j\le D\rho^{\ell-j}V_\ell,
\qquad
\sum_{\text{visited }j}V_j
\le\frac{D}{1-\rho}V_\ell=O(n).
\tag{11}
```

This argument permits a finite number of early-level factors in (10)
greater than one. It does not assume contraction at every single level.

The width slack in (1) is not charged additively as a runtime exponent.
Its effect is the bounded multiplicative product above.
This distinction permits slack as large as `lambda/loglog lambda`.

## 4. Enumerating pieces has no extra depth factor

At a visited node, the number of piece assignments is at most
`|I_lambda|^(g^4/4)`. Its base-two logarithm is at most

```math
\frac{g^4}{4}f(\lambda)
\le
\frac{C V}{4(\log_2(2\lambda))^\delta}
\le C'\frac{V}{(\log_2 n)^\delta},
\tag{12}
```

where (7) supplies a uniform constant `C'`.
Summing (12) along any path and using (11) gives

```math
O\!\left(\frac n{(\log n)^\delta}\right).
\tag{13}
```

There is no multiplication by `ell` in this estimate.
Such a multiplication would replace the geometric sum of the actual
volumes by an unnecessarily large bound.
This is why the endpoint `delta=alpha` loses no `log log n` factor.

For completeness, the recursive operation count can be bounded by a
maximum over paths. Each node's sum over piece assignments is at most
its number of assignments times the largest child cost.
Taking logarithms adds its local branching charge.
Repeatedly applying this inequality gives the sum of branching charges
along one worst path, together with its local query and termination costs.
It does not assume that the recursion has only one branch.

## 5. Generation, table queries, and stopping costs

The source grouping evaluates a child table entry by summing at most
`2^(2V/K)` terms. Each term makes polynomially many calls to parent
tables and factor values. The uniform generation and factor algorithms
add at most `2^(C lambda)` operations.
All integer index manipulation has size polynomial in `n`, because
the index bit lengths, volumes, and parameters along the path are `O(n)`.

Thus the multiplicative query overhead at a level has logarithm at most

```math
C_1\frac V{K_{j-1}}+C_2\log(n+2).
\tag{14}
```

Here `lambda=V/K_j^2<=V/K_(j-1)` absorbs factor generation.
The same estimate pays for constructing the complete piece list.

Using (11), its volume-dependent sum is bounded by

```math
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
\le D V_\ell\rho^\ell
\sum_{j\ge1}\frac{\rho^{-j}}{K_{j-1}}
=O(n\rho^\ell).
\tag{15}
```

The series converges uniformly for `1/2<=rho<1`.
For example, `K_(j-1)>=2^(18j)` from (4), while `rho^(-j)<=2^j`.
Its sum is at most `1/(2^17-1)`.

Finding the grouping enumerates `K_(j-1)^(g_(j-1))` colorings.
The logarithm of this work is at most `sqrt(P) log2 P` per level,
plus polynomial indexing work.
Summing it over the depth contributes
`O(n^(1/80)(log n)^2)`, which is smaller than `sqrt(n)`
for large `n`.

At level zero, brute force costs `2^O(V_0)`, with
`V_0=O(n rho^ell)`.
At an earlier cutoff, it costs `2^O(sqrt(n))`.
The accumulated parent-query overhead has already been charged in (14).
These statements remain valid over an arbitrary commutative
`F_p`-algebra: the count is of arithmetic circuit operations and
structural work, not unit-cost operations on an unbounded representation.

Combining (13)--(15), the stopping costs, and the top constructor gives
one uniformly generated weighted circuit of size

```math
A(n)=
\exp\!\left(
O\!\left(
n\rho^\ell+\frac n{(\log n)^\delta}
+\sqrt n+\frac nP
\right)\right).
\tag{16}
```

Equation (6) turns the first term into `O(n/(log n)^alpha)`.
The polynomial savings in the last two terms make them smaller than
that term. This proves the arithmetic-circuit part of (3).

Finally, the [sparse divisor transfer](sparse-divisor-nonzeroness.md)
charges every cyclic-ring representation and sequential tape operation,
adding `O(n/log n)` to the logarithmic runtime. That proves (3)
and the complete factoring claim.

## 6. What mere subexponential cardinality does not imply

The condition `log|I_lambda|=o(lambda)` alone does not make the
branching estimate (12) bounded by a fixed logarithmic saving.
For example, consider the entropy envelope

```math
f(\lambda)=
\left\lceil
\frac{\lambda}{\log_2\log_2(\lambda+4)}
\right\rceil.
\tag{17}
```

This is `o(lambda)`, but for every fixed `epsilon>0`,

```math
\frac{f(\lambda)}{\lambda/(\log\lambda)^\epsilon}
\longrightarrow\infty.
```

A literal grid expansion with that many pieces has a branching charge
of order `V/loglog lambda`, outside every fixed
`V/(log lambda)^epsilon` target.
This is a counterexample to a proposed inference from a cardinality
envelope. It is not a lower bound on every algorithm using such a family:
zero terms, shared structure, or a different summation method might avoid
literal enumeration.

For a precise integer subsequence, set `lambda_t=2^(2^t)`.
The entropy fraction is asymptotic to `1/t`, while its ratio to a target
with exponent `epsilon=1/k` grows as `2^(t/k)/t`.
Taking `t=km` gives the exact comparison `2^m/(km)`.

## 7. Reproduction and scope

Run:

```sh
python3 tools/check_subexponential_piece_transfer.py
```

The checker exercises the schedule, rounded volume recurrence, exact
label lower bound, early cutoffs, geometric branching sum, and local
table-query sum with rational arithmetic.
It also checks the explicit entropy counterexample.
It does not construct new crossing tensors or replace the proof above.

The grid and grouping algorithms are prior work from the pinned
Anthropic development. The top multiplication network is prior work
from the pinned factoring reference. This generalized interface,
cutoff analysis, and checker used substantial OpenAI assistant support.
No new Lean formalization was run.
