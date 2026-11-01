# Fusion projectors supply middle channels and smaller interfaces

The height-three braid quotient has a useful nonmonomial specialization:
put a triplet projector on every central pair. The central rank is exactly
the Fibonacci label count, rather than the full bit-space dimension.
Its boundary-degree coefficients satisfy

```math
r(T_b)\le h_{n,b}M_n,\qquad
c(T_b)\le\binom nb v_b,                         \tag{1}
```

where `M_n=Fibonacci(n-1)` is the full n-strand vacuum dimension,
`v_b` is the vacuum dimension on 2b strands (`v_0=1`), and
`h_(n,b)=min(b+1,n-b+1)`.

The second bound uses a general vacuum-interface factorization. It
replaces the physical `4^b` interface by `v_b=O(phi^(2b))`.
For `b<=ceil(n/5)`, (1) gives the candidate rank profile
**`(L^(1+o(1)),L^(3/2+o(1)))`** on a common binary alphabet.
At eight sites the b=2 piece has exact ranks **(17,28)** on the full
13-label alphabet.

These channels provide a source of middle-rank mass to combine with
the [nilpotent seam's asymmetric pieces](nilpotent-seam-crossings.md).
A compatible mixed spanning identity is still open. No stronger
complete factoring score is established.

## 1. Vacuum interfaces instead of physical bit channels

Consider a TL word between two n-strand vacuum encoders and their
norm duals, in the same order-five height-three quotient as the
[path model](quantum-braid-path-quotient.md). Fix b occurrences of
`e_(n-1)` crossing the block boundary. Suppose every other operation
acts inside a block or is a scalar. The resulting four-tensor has

```math
c(T)\le v_b.                                    \tag{2}
```

To prove this, cut vertically between the two physical blocks. Each
selected e is a cup followed by a cap, so the cut leaves 2b interface
endpoints on either side. All other external strands on each side
are closed by the vacuum encoder and dual. Bending the endpoints
in contour order therefore gives two vacuum tangles on 2b strands.

The negligible-wall ideal is killed by every closed cup pairing,
including this gluing. Reduce both tangles to the height-three fusion
basis of the 2b-strand vacuum quotient. Its Gram norms are units.
Absorb the fixed cup/cap orientation and norm factors in one side.
The resulting factorization is

```math
T(P,Q,R,S)=\sum_{H\in\mathcal V_{2b}}
                 A_H(P,R) B_H(Q,S),
\qquad |\mathcal V_{2b}|=v_b.                    \tag{3}
```

This is the actual column cut. It proves (2) over the fixed field,
without lifting a small quotient rank to an ambient invariant subspace.
The argument also applies to linear combinations of block tangles.
For b=0 the empty vacuum has dimension one.

Selecting b terms from n marked boundary gates gives the bound
`binom(n,b)v_b`. Additional polynomial coefficient convolutions or
auxiliary channels multiply that count by their explicit channel
numbers. In particular the sixteen auxiliary channels of the
[all-degree envelope](nilpotent-seam-rank-envelope.md) can use v_b
in place of4^b; the interface proof is unchanged after fixing those
channels.

## 2. The central triplet specialization

Let `q=8 in F16`, `delta=q+q^-1=7`, and

```math
p_2=I+\delta^{-1}e.
```

Since `e^2=delta e` in characteristic two, p_2 is an idempotent
and kills e. It is the two-strand triplet projector. In a scalar
central gate `qI+x e`, the value `x=q/delta=5` gives `q p_2`.
Thus this is an effective linear combination of central-degree
coefficients, rather than a claim about one monomial coefficient.

The central seam contains the disjoint pairs with indices
`0,2,...,2n-2`. Let `P_n` be their product of p_2 projectors in the
2n-strand vacuum quotient. Its exact rank is

```math
\operatorname{rank}P_n=M_n=\operatorname{Fibonacci}(n-1).
                                                               \tag{4}
```

Here is the fusion count behind (4). The image assigns object2 to
each adjacent pair. In the negligible-wall quotient, its fusion rule
is `2 tensor2=0 direct-sum2`: the ordinary additional object4 is
the discarded f_4 channel. The local fusion resolutions have unit
Gram factors for the retained channels, so their paths are a basis,
not just a dimension upper bound. Starting at object0, fusion by2
therefore has the two-state transition matrix

```math
\begin{pmatrix}0&1\\1&1\end{pmatrix}.
```

The vacuum multiplicity after n objects is Fibonacci(n-1).
For comparison, this argument does not assign rank M_n to the
ordinary physical tensor product of triplet subspaces. It concerns
the vacuum quotient, and the earlier path comparison transfers its
matrix coefficients to the boundary encoders.

## 3. The polynomial candidate family

Use the anti-diagonal circuit order of the
[two-seam construction](nilpotent-seam-crossings.md). On the central
seam use p_2; on the n boundary-seam cells use `qI+y e_i`; elsewhere
use `qI+e_i`. Apply the scalar `q^(n-n^2/2)` and project to the
middle-zero boundaries. Write the resulting polynomial as

```math
T(y)=\sum_{b=0}^n T_b y^b.
```

The scalar convention agrees exactly with specializing central
`qI+x e_i` gates at x=5 and applying the usual block phase. It has
no effect on either rank. No strict endpoints are asserted for T(y).
Each constant T_b can be appended to an existing strict crossing by
an isolated term `(t^k+t^(k+1))T_b`, which vanishes at both endpoints.
Polynomially many such additions require only polynomial degree.

There are n/2 boundary gates before the central projector and n/2
after it. Taking degree b gives at most h_(n,b) pairs of outer
coefficients. Equation (4) bounds every resulting row factor by M_n,
proving the first half of (1). Applying (2) to each selected subset
of boundary gates proves the second half.

For binary width lambda choose the same least even n and common
exact binary alphabet as for the other braid families. Then
`M_n=Theta(L)`, n is `O(lambda)`, and restriction cannot increase
either rank. The path implementation and exact dense factors have
`2^O(lambda)` construction cost. All tensors are over the fixed F16.
The existing [fixed-field descent](fixed-field-descent.md), applied
after amplification, handles the extension without changing the
uniform exponent; no local projection of every physical gate is used.

## 4. A certified middle range and a concrete mixed target

For every b, `v_b<=(21/8)^b`. At the binomial tilt `t=2/21`,

```math
\binom nb v_b\le(5/4)^n(21/2)^b.
```

For `b<=ceil(n/5)` this is at most a fixed factor times
`((5/4)(21/2)^(1/5))^n`. The exact inequality

```math
21^2 5^{25}<2^{67}
```

proves that this exponential base is smaller than
`(8/5)^(3/2)<phi^(3/2)`. Hence all these coefficients have
`r<=L^(1+o(1))`, `c<=L^(3/2+o(1))` and rank cost
`sqrt(rc)<=L^(5/4+o(1))`. Cut reflection interchanging Q,R gives
the reversed profile and preserves the target SWAP.

An explicit positive target is to combine these middle pieces and
their reflections with a fixed-width near-leading band of the
nilpotent family, also in both orientations, and the two cheap
pairing tensors I and its cut reflection. The band has row bound
`L^(1/2+o(1))`, including when the top allowed coefficient vanishes;
the [all-degree counter](nilpotent-seam-rank-envelope.md) supplies
the corresponding defect-count proof. Its other cut is at most L^2.
For any fixed band width the whole proposed list is polynomial size
and has budget `L^(5/4+o(1))`, below `L^sqrt(2)`.

This is a **candidate list**, not a spanning theorem. Unlike the
edges alone, it includes pieces with possible middle-rank mass.
The scalar profiles `(2,1/2),(3/2,1),(1,3/2),(1/2,2)` are not
immediately excluded by the basic minimum/weighted budgets. These
are exponent upper profiles; actual ranks and compatible coefficients
must still be established. Higher necessary constraints also remain
in force.

If the crossing belongs to this explicit span at all sufficiently
large widths, exact linear solving finds its coefficients within
`2^O(lambda)`. Expose a solution for `SWAP+I` in separate positive
monomial degrees with constant coefficient I. This gives the strict
crossing needed by the existing post-half transfer. Proving that
mixed spanning statement, or changing the kernel while preserving
the rank bounds, is the open research task.

## 5. Finite replay

[The compiler](../tools/quantum_braid_fusion_seam.py) propagates F16
path coefficients. The [checker](../tools/check_quantum_braid_fusion_seam.py)
computes full4/6/8-site profiles, checks all selected boundary subsets
at those lengths against (2), and verifies the global triplet-projector
ranks on one through eight pairs. Independent literal scalar propagation
checks complete sampled polynomials at all F16 arguments. The
[source-pinned receipt](../evidence/fusion-projector-seams.json) records
exact counts and rank pairs. It also constructs the proposed mixed
list with central gaps0,1,2 and tests its complete span. The crossing
belongs to that span at four sites; it is outside at six and eight
sites (span dimensions32 and54). These finite failures do not prove
an all-length obstruction, and they provide no positive evidence for
the current list at larger widths. A modified kernel or a structural
large-width identity remains necessary.

```sh
make verify-braid-fusion-seam
```

NumPy is needed for the finite matrix replay. The complete factoring
score remains `1/2`; these are ordinary uniform rank factors and
finite controls, without a complete decomposition or new Lean proof.
