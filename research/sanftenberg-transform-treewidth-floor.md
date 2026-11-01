# A treewidth floor for the public semiprime macro network

## Scope and source

This is an ordinary graph argument about the **unmodified** bit-variable
primal graph of the public `RiemannGaussian` macro network, at
[`dbsanfte/RiemannGaussian` commit `f5fd3ddb407972f9ca65298fbe9be3225fe06635`](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635).
The repository's current `main` still points to this source revision. The
claim concerns standard dense-table variable elimination (or a tree
decomposition of that primal graph). It is not a lower bound on factoring,
arithmetic circuits, structured tensor contraction, or representations that
compress the transform relations.

## A complete-bipartite line-graph minor

Write `ell` for the input bit length, `T = transformLength(ell)`, and let
`r` and `t` be the source's balanced transform axes, so `r*t=T` and
`r <= t <= 2r`. A word has `w=T+1` bit positions in this macro graph.
These are the literal source parameters in
[`SemiprimeMacroGraph.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroGraph.lean)
and
[`SemiprimeMacroParameters.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroParameters.lean).

Use one of the two forward transforms. In its bit-variable primal graph,
variables appearing together in one local table scope are adjacent. For each
cell `(j,k)` and bit position `s < w`, the first-transform output bit and
the corresponding scaled-word output bit are both in the scale factor's
scope, hence adjacent. Contract this edge. The contracted vertex is now in
both:

* the first-transform table at row `j`, whose scope contains all `t` input
  and all `t` output words; and
* the second-transform table at column `k`, whose scope contains all `r`
  input and all `r` output words.

The resulting selected vertices are indexed by `(j,(k,s))`. For fixed `j`
they form a clique, and for fixed `(k,s)` they form a clique. Therefore the
contracted primal graph contains the line graph `L(K_{r,t w})` as a subgraph:
the vertices of `K_{r,t w}` are row indices `j` and pairs `(k,s)`, and each
selected vertex represents their incident edge. Additional primal edges do
not affect this minor.

Harvey and Wood's Theorem 2 for the line graph of a complete multipartite
graph, specialized to `K_{a,b}`, gives

```math
\operatorname{tw}(L(K_{a,b})) \geq \frac{ab}{2}-a-b+\frac12.
```

With `a=r`, `b=t w`, we have `ab=rtw=T(T+1)` and `r+t w=O(T^{3/2})`.
Thus the source primal graph has treewidth at least

```math
\frac12T^2-O(T^{3/2})=\Omega(T^2)=\Omega(ell),
```

because the source proves `16*ell <= T^2 < 64*ell`. The external theorem is
in Harvey and Wood,
[*Treewidth of the Line Graph of Complete and Complete Multipartite Graphs*](https://arxiv.org/html/1210.8205v2),
Theorem 2.

## Consequence and next representation to test

Every elimination order therefore has width `Omega(ell)` and dense
intermediate tables of size `2^Omega(ell)`. Merely choosing a better
elimination order for the public network's ordinary dense local tables cannot
give a subexponential evaluation. This is a method-specific obstruction,
consistent with the separate acquired-boundary width estimate in
[`sanftenberg-all-input-prefix-adapter.md`](sanftenberg-all-input-prefix-adapter.md).
It does **not** rule out exploiting the Fourier/arithmetic structure inside
those tables, factoring each dense transform relation into a circuit, or
using a compressed tensor representation.

## Radix-2 butterfly tables reach the half-rate width floor

The pinned compiler's two axes `r,t` are powers of two, with `r*t=T`, and its
two-axis transform is a Cooley--Tukey factorization of the length-`T` DFT. This
is proved in
[`SemiprimeMacroFermatNetwork.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroFermatNetwork.lean#L109-L138)
and the axis powers are recorded in
[`SemiprimeMacroParameters.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroParameters.lean#L156-L171).
Consider the standard radix-2 expansion of those two axes, and encode each
butterfly as **one extensional Boolean table** on its two input and two output
words, each of width `w=T+1`. This section is specifically about that
four-word table model; it does not model a butterfly as a bit-level arithmetic
circuit.

There are

```math
d=\log_2 r+\log_2 t=\log_2 T
```

butterfly stages over `T=2^d` word positions. Index positions by
`x\in\{0,1\}^d`; stage `i` pairs positions `x` and `x\oplus e_i`. In the
primal graph of the four-word tables, take for each position `x` and bit
`s<w` the connected branch set formed by the wire variables carrying that bit
at that position through the stages; straight-through wires connect successive
stages. Every butterfly table contains all bits of both input and both output
words. Thus the branch sets form cliques of size `w` for each `x`, with
complete joins between the cliques at `x` and `x\oplus e_i`. The resulting
minor is

```math
Q_d\boxtimes K_w,
```

the strong product of the `d`-cube and a `w`-clique. Kaul's
[Theorem 1](https://arxiv.org/html/2607.16778) gives
`tw(G\boxtimes H) >= (tw(G)+1)(tw(H)+1)-1`; Chandran and Kavitha's
[hypercube result](https://doi.org/10.1016/j.disc.2005.12.011) gives
`tw(Q_d)=Theta(2^d/sqrt(d))`. Therefore

```math
\operatorname{tw}(Q_d\boxtimes K_w)
  \ge w(\operatorname{tw}(Q_d)+1)-1
  =\Omega\!\left(\frac{wT}{\sqrt{\log T}}\right)
  =\Omega\!\left(\frac{\ell}{\sqrt{\log \ell}}\right).
```

Here `w=T+1` and `16*ell <= T^2 < 64*ell`. Thus a dense-table variable
elimination on the radix-2-expanded network still has an intermediate table
of size `2^Omega(ell/sqrt(log ell))`. This reaches the same asymptotic width
scale as the current half-rate benchmark, so this table model cannot establish
a better logarithmic saving exponent merely through stagewise elimination.
This lower bound is not an upper bound for the butterfly evaluator, and it
does not rule out a compact bit-level circuit for each add, subtract, or
twiddle multiplication.

The root-four twiddles do simplify at bit level: they are signed cyclic
rotations with carry normalization, as derived in
[`sanftenberg-root-four-bit-circuits.md`](sanftenberg-root-four-bit-circuits.md).
The later generic multiplier theorem in
[`multiplier-treewidth-and-relation-cycles.md`](multiplier-treewidth-and-relation-cycles.md)
rules out sublinear treewidth for whole bit-variable multiplication networks
with explicit operand and product bits, even when the product is then pinned
to `N`. The source audit in
[`riemann-gaussian-semiprime-source-audit.md`](riemann-gaussian-semiprime-source-audit.md)
checks that the current D1 carry-target representation satisfies those
hypotheses. The remaining open evaluator question is whether structured
arithmetic contraction can exploit its local tables; this treewidth floor does
not bound such evaluators or networks whose operands exist only in transform
coordinates.

No upstream code was copied and no Lean rebuild was run. This ordinary
argument does not establish a stronger factoring exponent; the repository's
score remains `1/2`. Substantial OpenAI assistant support was used in this
source audit and graph analysis.
