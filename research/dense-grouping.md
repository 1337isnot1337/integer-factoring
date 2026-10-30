# Grouping the dense width graphs produced by a grid step

## Result and status

The universal grouping theorem in `research/grouping-improvement.md` applies
to arbitrary weighted graphs and uses `g = O(D^2 K^2 log K)` vertices. The
row and column width graphs actually used by the recursion are dense. Every
off-diagonal weight lies between two thirds of a common upper bound and that
upper bound. For this restricted family, the following elementary theorem
uses only `g = O(D^2 K log(DK))` vertices.

This note independently derives a finite specialization of the bounded-edge
grouping route already present in Appendix A of `anthropics/formal-math`,
pinned at `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`. The existing formal
result is described below. This is not a claim that the grouping idea or
the resulting base-two schedule is new.

This note proves the finite grouping theorem, provides an exact deterministic
construction, and derives a parameter schedule. Its factoring rates remain
conditional on integrating that schedule and grouping construction into the
actual evaluator and proving its complete cost. No changed factoring machine
or new Lean theorem is claimed here.

The separate [certified k6 machine result](certified-k6.md) has exponent
`alpha = 1/30` and uses `G = 5K^6`, `K_next = 25K^12`. It does not use this
note's asymptotically degree-two schedule. The stronger rates derived in
Section 6 remain conditional; the [Appendix B route](appendix-b-factoring-route.md)
to `alpha = 1/2` also remains conditional on its stated formal bridge.

Source revision: `IsaacRudich/SubexponentialDeterministicFactoring` at
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. Repository base inspected:
`1337isnot1337/integer-factoring` at
`97521a3dfc89c9f36a7bd2a61a8f155eaa5314f0`.

## 1. Exact finite theorem

Let `K,D >= 2` be integers. For positive integers `x`, write

```math
 b(x)=\lfloor\log_2 x\rfloor+1.
```

Set

```math
 L=b(K)+b(D)+8,\qquad g=576D^2KL.
 \tag{1}
```

Let `B` be a positive integer and let `w_uv` be symmetric nonnegative integer
weights on a graph with `g` vertices, with zero diagonal. Suppose that

```math
 2B\le 3w_{uv}\le3B\qquad(u\ne v).
 \tag{2}
```

Define

```math
 d_u=\sum_v w_{uv},\qquad U=\sum_u d_u.
```

Let `P` be any integer with `d_u <= P` for every vertex. For a coloring
`chi:[g]->[K]`, let `s_a` be the size of class `a`. For distinct `a,b`, let
`Y_ab` be the weight of edges with one endpoint in each class, with each edge
counted once. Then there exists a coloring satisfying

```math
 Ks_a\le2g\quad\text{for every }a,
 \qquad
 DK^2Y_{ab}\le DU+gP\quad\text{for every }a\ne b.
 \tag{3}
```

All parameter computations and all tests in (2)--(3) use integer arithmetic.
No assumption about the values of the tensor-network entries is made. Only
the nonnegative integer channel widths are the weights in this theorem.

### A more general sufficient condition

The same proof works for any integer `g >= 2` satisfying

```math
 (g+2)K\exp\!\left(-\frac{g-1}{72D^2K}\right)<1.
 \tag{4}
```

Equation (1) is a deliberately conservative, explicit integer choice that
satisfies (4). There is no assertion that 576 is an optimal constant.

## 2. A weighted upper-tail inequality

Let `X_i` be independent Bernoulli random variables with common success
probability `p`. Let `0 <= a_i <= A`, where `A > 0`, and put

```math
 Z=\sum_i a_iX_i,\qquad \mu=p\sum_i a_i.
```

For `0 < eta <= 1/2` and `mu > 0`,

```math
 \Pr\{Z\ge(1+\eta)\mu\}
 \le\exp\!\left(-\frac{\eta^2\mu}{3A}\right).
 \tag{5}
```

Here is a proof with all constants. For `theta >= 0`, convexity of the
exponential on `[0,1]` gives

```math
 e^{\theta a_i/A}\le1+\frac{a_i}{A}(e^\theta-1).
```

Independence and `1+x <= exp(x)` therefore imply

```math
 \mathbb E e^{\theta Z/A}
 =\prod_i(1-p+pe^{\theta a_i/A})
 \le\exp\!\left(\frac\mu A(e^\theta-1)\right).
```

Apply Markov's inequality and take `theta=log(1+eta)`. The resulting exponent
is `-(mu/A)h(eta)`, where

```math
 h(x)=(1+x)\log(1+x)-x
 =\int_0^x\log(1+t)\,dt
 \ge\frac{x^2}{2(1+x)}.
```

For the last inequality, `log(1+t) >= t/(1+t) >= t/(1+x)` on `[0,x]`.
For `0 <= x <= 1/2`, `2(1+x) <= 3`. This proves (5).

## 3. Simultaneous vertex and degree balances

Color every vertex independently and uniformly from `[K]`. Set

```math
 \eta=\frac1{4D}\le\frac18.
```

For each vertex `u` and color `b`, define

```math
 X_{ub}=\sum_{v:\chi(v)=b}w_{uv}.
```

For fixed `u,b`, its summands are independent weighted Bernoulli variables.
The zero diagonal causes no problem. Its mean is `d_u/K`, each weight is at
most `B`, and (2) gives `d_u >= (2/3)(g-1)B`. Applying (5),

```math
 \Pr\left\{X_{ub}>(1+\eta)d_u/K\right\}
 \le\exp\!\left(-\frac{g-1}{72D^2K}\right).
 \tag{6}
```

No independence between different pairs `(u,b)` is used.

For each color `a`, define the degree sum

```math
 Z_a=\sum_{u:\chi(u)=a}d_u.
```

Its mean is `U/K`. Every weight is at most `(g-1)B`, while
`U >= (2/3)g(g-1)B`. Hence (5) also gives

```math
 \Pr\{Z_a>(1+\eta)U/K\}
 \le\exp\!\left(-\frac{g}{72D^2K}\right).
 \tag{7}
```

Finally, `s_a` is binomial with mean `g/K`. Markov's inequality applied to
`2^{s_a}` gives

```math
 \Pr\{s_a>2g/K\}
 \le(e/4)^{g/K}\le e^{-g/(3K)}.
 \tag{8}
```

The last inequality uses `2 log 2 - 1 > 1/3`. For example,
`log 2 = integral_1^2 dx/x > 2/3` follows by integrating the tangent line of
the convex function `1/x` at `x=3/2`; strictness holds away from that point.
The right side of (8) is no larger than the right side of (6).

There are `gK` events in (6), `K` events in (7), and `K` events in (8).
Their total failure probability is bounded by the left side of (4).

To check (4) for the integer choice (1), note that `L >= 12`, `g >= 2`, and

```math
 \frac{g-1}{72D^2K}\ge\frac{g}{144D^2K}=4L.
```

Also `g+2 <= 2g`, `1152 < 2^11`, `K < 2^{b(K)}`, `D < 2^{b(D)}`, and
`L <= 2^L`. Therefore

```math
 (g+2)K\le1152D^2K^2L
 <2^{11+2b(K)+2b(D)+L}=2^{3L-5}<e^{3L}.
```

The total failure probability is less than `e^{-L}<1`. Thus one coloring
satisfies all three sets of inequalities simultaneously.

For this coloring and distinct `a,b`, symmetry and the once-per-edge
convention give the exact identity

```math
 Y_{ab}=\sum_{u:\chi(u)=a}X_{ub}.
```

Consequently

```math
 Y_{ab}\le\frac{1+\eta}{K}Z_a
 \le\frac{(1+\eta)^2}{K^2}U.
```

Since `D >= 2`,

```math
 (1+1/(4D))^2
 =1+\frac1{2D}+\frac1{16D^2}
 \le1+\frac1D.
```

Finally `U <= gP`, so

```math
 DK^2Y_{ab}\le(D+1)U\le DU+gP.
```

Together with the size bounds, this proves (3).

## 4. Uniform deterministic construction and its cost

Enumerate all `K^g` colorings in lexicographic order, stopping at the first
coloring satisfying the integer tests (3). The existence proof guarantees
termination. The algorithm does not sample randomness, evaluate exponentials,
or receive an advice string.

For one coloring, class sizes and every `Y_ab` can be accumulated by a scan
of the `g(g-1)/2` off-diagonal entries. The work is
`poly(g,K,log(B+1),log(D+1),log(P+1))` bit operations. Thus a conservative
construction cost is

```math
 K^g\operatorname{poly}(g,K,\log(B+1),\log(D+1),\log(P+1)).
 \tag{9}
```

This is the same finite-search mechanism already allowed in the candidate
integration. The stronger intermediate inequalities (6)--(8) only prove
existence; the implementation can check (3) directly. A faster grouping
algorithm is not needed for the recurrence calculation.

## 5. Mapping to the actual grid graphs

The pinned `FactoringAux/Level.lean` defines

```math
 \texttt{rowNet.c}(a,a')
 =G\lambda+\sum_{q\in\binom{[G]}2}
       c_r(\lambda,i_{(a,a'),q}),
 \tag{10}
```

with the analogous formula for columns. Its `Net.cc` chooses the ordered
version of this width off the diagonal, makes it symmetric, and sets the
diagonal to zero. Thus `w=Net.cc` is the weight matrix in the theorem, and
`U=Net.total` has exactly the ordered-pair convention used above.

Both the original and compressed widths satisfy

```math
 \frac43\lambda\le c_r(\lambda,i),c_c(\lambda,i)\le2\lambda.
 \tag{11}
```

For the original widths this follows from their form
`min(2*lambda, 4*ceil(lambda/3)+3*i)` and its reversed-index version.
For the compressed widths, replace `3*i` by the nonnegative integer
`ceil(log_2(7^i))`, or by `ceil(73*i/26)`. Every argument of the minimum
is at least `4*lambda/3`.

Substituting (11) into (10) gives, for distinct vertices,

```math
 \frac{2G^2+G}{3}\lambda
 \le w_{uv}\le G^2\lambda.
 \tag{12}
```

For positive `lambda`, set `B=G^2*lambda`. Equation (12) implies (2).
The source already proves `Net.Fits(G^3*lambda)`, so use `P=G^3*lambda`.
The theorem produces exactly the group-size and pair-width conclusions of
`GoodAccuracy` in `formal/GroupingThreshold.lean`.

For `lambda=0`, all these widths are zero. The pair conditions hold for
every coloring. A balanced coloring, with class sizes at most `ceil(G/K)`,
satisfies the size bound because this schedule has `G>=K`. This case must
be handled explicitly rather than applying a tail bound with `B=0`.

The existing `Net.value_group` identity is independent of how the valid
coloring is found. A formal integration must still discharge its width and
size hypotheses for the actual new selector; the finite theorem does not
automatically alter that selector or the compiled evaluator.

### Why the path obstruction does not apply

The quadratic path obstruction in PR #9 concerns a guarantee for *every*
weighted graph. A nontrivial path has many zero off-diagonal weights and
cannot satisfy (2) for any `B>0`. No hypothesis of that obstruction is
silently weakened: the new theorem explicitly restricts the graph family.
It does not disprove the obstruction or imply a lower bound on factoring.

## 6. Parameter growth and conditional rate

For levels numbered from the bottom, set

```math
 \begin{aligned}
 K_0&=48,\\
 D_d&=(d+48)^2,\\
 L_d&=b(K_d)+b(D_d)+8,\\
 G_d&=576D_d^2K_dL_d,\\
 K_{d+1}&=G_d^2.
 \end{aligned}
 \tag{13}
```

The graph at level `d` has exactly `G_d` vertices, so the finite theorem
applies at every positive width. All of (13) is integer arithmetic.

These parameters obey

```math
 K_d\ge48\cdot2^d\ge d+48,\quad
 D_d\le K_d^2,\quad L_d\le K_d,\quad
 K_d\le G_d\le K_d^8.
 \tag{14}
```

For the bound on `L_d`, use `b(D_d)<=2b(K_d)`, and
`3b(K)+8<=K` for `K>=48`. The latter follows, for example, from
`4b(K)<=K` for `K>=32`: if `2^r<=K<2^{r+1}` and `r>=5`, then
`4(r+1)<=2^r`. Thus `3b(K)+8<=3K/4+8<=K`. Finally
`576D_d^2K_dL_d<=576K_d^6<=K_d^8`, since `K_d>=48`.

The coarse bounds `K_d^2<=K_{d+1}<=K_d^16` imply
`log K_d<=16^d log 48`. Writing `a_d=log K_d`, the exact recurrence is

```math
 a_{d+1}=2a_d+2\log576+8\log(d+48)+2\log L_d.
 \tag{15}
```

The coarse upper bound gives `b(K_d)=O(16^d)`, while
`b(D_d)=O(log(d+48))`. Hence `log L_d=O(d+1)`, and the nonnegative error
in (15) is `O(d+1)`. Divide (15) by `2^{d+1}` and sum. The series
`sum(d+1)/2^{d+1}` converges, proving

```math
 \log K_d=\Theta(2^d).
 \tag{16}
```

Choose the largest `ell` with `K_ell^40<=n`, and set `Ptop=K_ell`.
Then

```math
 \ell=\log_2\log n+O(1),\qquad
 n^{1/640}<P_{\rm top}\le n^{1/40}.
 \tag{17}
```

The lower inequality follows from `n<K_(ell+1)^40<=K_ell^640`.
Small input lengths can use the existing finite fallback. The original top
multiplication-network construction can retain its source grouping
`3*Ptop^8`; only recursive grid groupings change.

The intended recurrence is

```math
 N_d\le(\rho+1/G_d+1/D_d)N_{d+1}+cG_d^4,
 \tag{18}
```

where `(rho,c)=(11/12,2)` for the original crossing factors, or
`rho=(8+log_2 7)/12`, `c=3` for the compressed factors. The width must be
chosen identically for all piece branches, using the corresponding existing
common-width formula.

Both error sequences are summable. For instance,
`sum 1/D_d<=1/47` and `sum 1/G_d<=sum 1/K_d<=1/24`.
For either value of `rho>9/10`, the product of the factors in (18) over any
level interval is at most a fixed constant times the corresponding power
of `rho`. Also `G_d^4=K_(d+1)^2<=Ptop^2` on every recursive level.
The existing all-level volume argument therefore gives

```math
 N_d=O(\rho^{\ell-d}n+\ell P_{\rm top}^2),\qquad
 \max_{d<\ell}\frac{N_{d+1}+N_d}{K_d}
 =O(\rho^\ell n+\ell P_{\rm top}^2).
```

The second estimate uses `K_d>=48*2^d` and `2*rho>1`.

The search exponent is still at most
`G_d log_2 K_d <= sqrt(Ptop) log_2 Ptop`. The total logarithm of the
choice-path count is `O(ell*Ptop^2*log n)`. The table construction, internal
channel loops, positive multiplier, and small-modulus checks need exactly
the cost accounting described in `research/integration.md`; this theorem
does not eliminate any of those loops. The old top-level overhead remains
`O(n^(639/640)+n^(19/20))` under the same conservative cutoff.

**Conditional consequence.** If that evaluator construction and cost bridge
are established for (13), the main rate is

```math
 T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right),
 \qquad \alpha=\frac{\log(1/\rho)}{\log2}.
```

This is `alpha=log(12/11)/log 2` for the original factors, and
`alpha=log(12/(8+log_2 7))/log 2` for the compressed factors. Each is twice
the corresponding candidate rate with the old base-four schedule.

The finite grouping theorem above has an ordinary mathematical proof in this
note. Its Python checks cover finite certificates, not a Lean proof of that
theorem or of the asymptotic rate. The factoring statement in this section
is conditional on its named integration obligations. The separate certified
machine exponent remains `1/30`.

## Relationship to the existing formal development

At the pinned `anthropics/formal-math` revision,
`eth/proof/statements/AppRunningTimeStatement.lean` defines the second branch
of `AppRunningTime.Cases` by

```math
 K\mid g,\quad w_{uv}\le P/g\;(u\ne v),\quad
 K\log(2gK)\le2\zeta^2g.
```

`AppRunningTime.someGrouping`, assembled in
`eth/proof/rates/GeneralPieces/Main.lean`, proves the resulting grouping
existence statement. Its proof is in
`eth/proof/rates/GeneralPieces/Group/Lemma.lean` and
`eth/proof/rates/GeneralPieces/Group/CaseII.lean`.
The integer pairing width in that public statement is the ceiling of
`(U+2*zeta*g*P)/K^2`; the case-II proof proves a strict bound before the
ceiling. This rounding distinction must be respected when importing it into
the present repository's floor-based width schedule.

For the actual grid graphs, `P=G^3*lambda` and `B=G^2*lambda`, so the
bounded-edge hypothesis is exact. It does not require a positive lower
edge-weight bound. The accompanying `bounded-edge-grouping.md` gives an
elementary additive-concentration proof for this stronger grid-relevant
form, with another explicit integer schedule. The dense proof above allows
the smaller choice `P=max d_u`, which need not meet `B<=P/g`, but that extra
generality is not needed by the current factoring grid.

## Sources

- [Pinned grid construction and width sums](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Level.lean).
- [Pinned original piece widths](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Pieces.lean).
- [Existing grouping predicate and value interface](https://github.com/1337isnot1337/integer-factoring/blob/97521a3dfc89c9f36a7bd2a61a8f155eaa5314f0/formal/GroupingThreshold.lean).
- [Universal grouping obstruction and dense-graph suggestion](https://github.com/1337isnot1337/integer-factoring/pull/9).
- [Existing common-width and cost accounting](https://github.com/1337isnot1337/integer-factoring/blob/97521a3dfc89c9f36a7bd2a61a8f155eaa5314f0/research/integration.md).
- [Formal bounded-edge grouping statement](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).
- [Formal grouping implementation and value theorem](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/Lemma.lean).
