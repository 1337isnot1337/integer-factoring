# A sharper grouping lemma for the deterministic factoring engine

## Status and scope

This note gives an elementary proof of a stronger finite grouping lemma. The
lemma is independent of the factoring reduction. It improves the node-count
requirement in `FactoringAux/Grouping.lean` and `FactoringAux/Lemma1.lean` of
Isaac Rudich and Louis-Martin Rousseau's repository at commit
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.

The proof below is an ordinary mathematical proof. It has not been translated
into Lean or checked by the Lean kernel. The final section states exactly how
to use it in the counting engine and identifies the code changes still needed.

The principal improvement is to separate the allowed discrepancy from the
number of groups. With discrepancy parameter `D`, the number of vertices
needed is `O(D^2 K^2 log K)`, instead of the source's `3 K^8` at discrepancy
`1/K`.

## 1. Definitions and theorem

Let `g`, `K`, `D`, and `P` be integers with `K >= 2`, `D >= 2`, and `P >= 0`.
Let `w_uv = w_vu` be nonnegative integer weights on the unordered pairs of
distinct vertices of a graph on `g` vertices. Assume

\[
\sum_{v\ne u}w_{uv}\le P \qquad\text{for every vertex }u.
\]

Write the total over ordered pairs as

\[
U=\sum_u\sum_{v\ne u}w_{uv}=2\sum_{u<v}w_{uv}\le gP.
\]

For a coloring `chi : [g] -> [K]`, let `s_a` be the size of color class `a`,
and, for distinct colors `a,b`, let

\[
Y_{ab}=\sum_{u<v}w_{uv}
\mathbf 1_{\{\chi(u),\chi(v)\}=\{a,b\}}.
\]

Each edge between two different groups is counted once in `Y_ab`. Thus this
is the same quantity as the source's `(N.coarsen chi).c a b`; `U` retains
the source's convention of counting each original edge twice.

**Theorem.** If

\[
g\ge 12D^2K^2\log K,
\tag{1}
\]

where `log` is the natural logarithm, then there is a coloring with

\[
Ks_a\le 2g \quad\text{for all }a,
\qquad
DK^2Y_{ab}\le DU+gP \quad\text{for all }a\ne b.
\tag{2}
\]

In particular, every pair of groups shares at most

\[
\left\lfloor\frac{DU+gP}{DK^2}\right\rfloor
\tag{3}
\]

bits. All conclusions are deterministic existence statements. A finite
deterministic search, described in Section 4, finds such a coloring.

A convenient parameter that uses integer arithmetic only is

\[
g=12D^2K^2(\lfloor\log_2K\rfloor+1).
\tag{4}
\]

Since `floor(log_2 K)+1 > log_2 K >= log K`, (4) satisfies (1).

## 2. Weighted edges reduced to matchings

For `P > 0`, replace each edge `{u,v}` by `w_uv` parallel edges. The resulting
loopless multigraph has maximum degree at most `P`.

Color its edges greedily with `s=2P-1` colors. When an edge is processed,
at most `P-1` already-colored edges meet each endpoint, so at most `2P-2`
colors are forbidden. Therefore the greedy construction succeeds. Each
edge color class is a matching; empty classes are allowed.

Now give every vertex an independent uniform random color from `[K]`.
Fix distinct colors `a,b`, and let `Z_j` count the edges in matching `j`
whose endpoint colors are `a,b`. If that matching contains `m_j` edges,
then

\[
Z_j\sim\mathrm{Binomial}(m_j,p),\qquad p=2/K^2.
\]

Independence holds within each matching because its edges have disjoint
endpoints. No independence is asserted between different matchings.
Moreover,

\[
Y_{ab}=\sum_{j=1}^sZ_j,
\qquad
\mu:=\mathbb E Y_{ab}=p\sum_jm_j=U/K^2.
\]

For any `theta >= 0`, generalized Holder's inequality gives

\[
\begin{aligned}
\mathbb E e^{\theta Y_{ab}}
&\le\prod_{j=1}^s
   (\mathbb E e^{s\theta Z_j})^{1/s}\\
&=\prod_{j=1}^s
   (1-p+pe^{s\theta})^{m_j/s}\\
&\le\exp\left(\frac\mu s(e^{s\theta}-1)\right).
\end{aligned}
\tag{5}
\]

The last inequality is `1+x <= exp(x)` for `x >= 0`. Parallel copies of
one original edge may be perfectly dependent, but they belong to different
matchings, and (5) explicitly allows this dependence.

This argument uses the multigraph only to prove existence and the tail bound.
The deterministic search does not need to construct or store it.

## 3. Tail bounds and the simultaneous coloring

If `mu > 0`, Markov's inequality applied to (5), with

\[
\theta=s^{-1}\log(1+t/\mu),
\]

gives

\[
\Pr(Y_{ab}\ge\mu+t)
\le\exp\left[-\frac\mu s
  \big((1+t/\mu)\log(1+t/\mu)-t/\mu\big)\right].
\]

For `x >= 0`,

\[
(1+x)\log(1+x)-x
=\int_0^x\log(1+u)\,du
\ge\int_0^x\frac{u}{1+x}\,du
=\frac{x^2}{2(1+x)}.
\]

Consequently,

\[
\Pr(Y_{ab}\ge\mu+t)
\le\exp\left(-\frac{t^2}{2s(\mu+t)}\right).
\tag{6}
\]

When `mu=0`, every edge weight is zero and the desired pair condition holds
deterministically; no division by `mu` is needed.

Choose

\[
t=\frac{gP}{DK^2}.
\]

Using `s <= 2P`, `mu <= gP/K^2`, and `D >= 2`, we have

\[
2s(\mu+t)
\le 4P\frac{gP}{K^2}(1+1/D)
\le\frac{6gP^2}{K^2}.
\]

Thus, for `P > 0`,

\[
\Pr\left(Y_{ab}>\frac{U+gP/D}{K^2}\right)
\le\exp\left(-\frac{g}{6D^2K^2}\right)
\le K^{-2}.
\tag{7}
\]

The last inequality uses (1). The union bound over the `K(K-1)/2`
unordered color pairs is less than `1/2`.

For each color `a`, `s_a` is binomial with mean `g/K`. Its elementary
Chernoff bound at twice the mean is

\[
\Pr(s_a>2g/K)
\le(e/4)^{g/K}\le e^{-g/(3K)}.
\tag{8}
\]

For completeness, the first bound follows by applying Markov's inequality
to `2^{s_a}` and using
`E 2^{s_a}=(1+1/K)^g <= exp(g/K)`; the second uses
`2 log 2 - 1 > 1/3`.

By (1), the union bound over all group-size failures is at most

\[
K e^{-g/(3K)}\le K^{1-4D^2K}<\tfrac12.
\tag{9}
\]

Equations (7)--(9) show that the probability of any failure is less than
one. A coloring satisfying (2) therefore exists. For `P=0`, all pair
conditions hold automatically, and the same group-size argument applies.

This completes the proof of the theorem, including zero-weight and
zero-total cases.

## 4. Deterministic construction and bit complexity

Enumerate the `K^g` colorings in lexicographic order and stop at the first
one satisfying (2). Test the conditions with the exact integer comparisons

\[
Ks_a\le2g,
\qquad DK^2Y_{ab}\le DU+gP.
\]

No real arithmetic, random choices, or knowledge of the successful coloring
is used in this algorithm. The proof guarantees termination before all
`K^g` colorings have been rejected.

For each coloring, compute all class sizes and pair sums by scanning the
`g(g-1)/2` entries of the upper triangular weight matrix. This needs
`O(g^2+K^2)` integer additions and comparisons, plus polynomial work in
`g`, `K`, and the bit lengths of `P,D`. With schoolbook integer arithmetic,
a conservative bound is

\[
K^g\,\mathrm{poly}(g,K,\log(P+1),\log(D+1))
\tag{10}
\]

bit operations. Space is polynomial in the same quantities. A single
uniform program computes the parameters and performs the enumeration;
this construction introduces no advice string or family of unrelated
machines.

The proof does **not** claim a polynomial-time grouping algorithm. Replacing
exhaustive search by conditional expectation would require a separately
specified, efficiently evaluable pessimistic estimator. It is unnecessary
for the asymptotic route below.

## 5. Immediate corollary at the source's discrepancy

Taking `D=K` recovers the source's label-width target,

\[
\frac{U+gP/K}{K^2}=\frac{KU+gP}{K^3},
\]

but needs only `g=O(K^4 log K)` vertices. The original source requires
`g=3K^8`. This corollary alone lowers the effective recursion parameter
growth from `K_next=K^{16+o(1)}` to `K_next=K^{8+o(1)}`.

## 6. Using a free discrepancy parameter in the counting engine

The following algebra specifies a stronger modification than simply
taking `D=K`.

For one grid step, the source proves that one of the row and column
networks has `G` nodes, fits `P=G^3 lambda`, and has ordered total at most

\[
U_{\max}(G,\lambda)
=\left\lfloor
  \frac{12G^3\lambda+11G^4\lambda+24G^4}{12}
 \right\rfloor.
\tag{11}
\]

These are the parameter-generic statements `rowNet_fits`, `colNet_fits`,
and `min_total` in `FactoringAux/Level.lean`, assembled as
`total_le_Umax` in `FactoringAux/EngineGroup.lean`.

Apply the theorem of Section 1 with `g=G` and define the new width by

\[
L_D(K,G,\lambda)
=\left\lfloor
 \frac{D U_{\max}(G,\lambda)+G^4\lambda}{DK^2}
 \right\rfloor.
\tag{12}
\]

Writing `N=G^4 lambda` and `N'=K^2 L_D`, equations (11)--(12) give

\[
N'\le\left(\frac{11}{12}+\frac1G+\frac1D\right)N+2G^4.
\tag{13}
\]

There is no extra positive rounding error in (13), since both width
definitions round down. The valid grouping guarantees that these rounded
widths remain large enough to store every crossing edge label.

The number of nodes in each new group remains at most `2G/K`. Thus the
same direct table-construction procedure enumerates at most

\[
2^{(2G/K)P}=2^{2N/K}
\]

assignments for the incident bits of one group, before polynomial factors
and the other entry computations already present in the source are
included. Exhaustive grouping costs `K^G`; it is accounted for explicitly
by (10), and must not be omitted from a complete exponent calculation.

For levels counted from the bottom, a useful parameter schedule is

\[
\begin{aligned}
K_0&=48,\\
D_d&=(d+48)^2,\\
G_d&=12D_d^2K_d^2(\lfloor\log_2K_d\rfloor+1),\\
K_{d+1}&=G_d^2.
\end{aligned}
\tag{14}
\]

The relation `K_{d+1}=G_d^2` is exactly what the square-grid transformation
requires. It is not an optional simplification.

Here `log K_d = Theta(4^d)`, so the admissible number of levels under any
fixed-power cutoff `K_ell^A <= n` is

\[
\ell=\frac{\log\log n}{\log4}+O(1).
\tag{15}
\]

Also, `sum_d (1/G_d+1/D_d)` converges. The product of the multiplicative
coefficients in (13) is therefore at most a constant times `(11/12)^ell`.
Thus the main contraction term becomes

\[
n(11/12)^\ell
=O\left(\frac{n}{(\log n)^\alpha}\right),
\qquad \alpha=\frac{\log(12/11)}{\log4}\approx0.062765441.
\tag{16}
\]

The recurrence analysis, additive errors, top-level construction, all
choice-path enumeration, modular zero testing, and conversion to a
fixed finite-alphabet Turing machine must all be included before (16)
can be promoted to an unconditional end-to-end factoring theorem.
In particular, a change to (14) does not change the already compiled
upstream machine automatically.

The simplest top-level adaptation is to retain the source's original
top grouping with `g_top=3K_ell^8` and its proven initial table construction,
and use the new schedule only for the recursive contraction chain. The
old cutoff `K_ell^40 <= n` remains conservative. This avoids needing a
new multiplication-network construction or a new top-level grouping proof.

## 7. Source-level proof obligations

The finite grouping theorem has been proved above. A complete Lean
implementation of the improved algorithm still needs the following
specific changes:

1. Generalize `Net.lemma1` and `good_exists` to the assumptions of (1) and
   to the width in (3), replacing the exact hypothesis `g=3K^8`.
2. Add `D` to `Lnext`, `nextT`, and the good-code predicate. Replace
   `Lnext` by (12); keep the row/column value identities unchanged.
3. Replace the `Kb` and `Gs` recurrences by (14) in a new engine namespace.
   Do not silently modify an upstream theorem and continue citing its
   existing compiled proof.
4. Generalize the uniform expression evaluator to read the level index,
   compute `D_d`, and use the new width. All arithmetic remains integer.
5. Prove bounds for the actual new expression, including `K^G`, path
   counts, table writes, modulus enumeration, and bit lengths. Existing
   parameter-generic table and compiler lemmas can be reused only after
   their hypotheses have been checked against the modified definitions.
6. State a new final theorem and audit its axioms. A plain-text asymptotic
   derivation does not establish that the new Lean theorem compiles.

## Sources

- [Pinned source repository](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
- [Grouping lemma used by the source](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Lemma1.lean).
- [Grid transformation and its size bounds](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Level.lean).
- [Current width and recursion definitions](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean).
- [Parameter-generic exponent accounting](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostReq.lean).

The matching decomposition and probability estimates used here are proved
in full above; they are not unverified calls to an external concentration
theorem.
