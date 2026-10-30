# A top-clique compiler for weighted deterministic gate networks

## Result and exact scope

The existing top constructor can accept more than one fixed-output
multiplication network. It can compile any deterministic-message gate
network satisfying the two weight bounds below. In particular, ports of
several multiplication networks can be joined before the top grouping.

Let `n>=1`. A gate `v` has a finite setting set `S_v`, with its cardinality
and an index-to-setting routine supplied, a weight in `{0,1}`, and one deterministic message for each
incident wire. An undirected wire `e` has a nonnegative integer bit width
`w_e`; its endpoints accept exactly when their messages agree. The messages
must be smaller than `2^w_e`. Set

```math
s_v=\lceil\log_2\max(1,|S_v|)\rceil,\qquad
a_v=s_v+\sum_{e\ni v}w_e,\qquad
W=\sum_v a_v,\qquad D=\max_v a_v.                 \tag{1}
```

An empty setting set makes the network value zero and is handled before
compilation. Suppose the network has polynomially many gates and wires,
its metadata, cardinalities, and local index/setting/message queries take polynomial bit time
in `n`, and, for a fixed constant `C>=1`,

```math
W\le Cn,\qquad D\le Cn^{3/4}.                    \tag{2}
```

For every integer `P>=2` with `P^40<=n`, a deterministic algorithm writes
a `P`-vertex clique of common positive label width `lambda`, preserving
the exact network value, with

```math
P^2\lambda\le(1+3C)n,\qquad
H\le 8Cn/P.                                     \tag{3}
```

Every table entry is a sum of at most `2^H` local setting contributions.
The complete construction takes `2^O(n/P)` bit operations. The constants
are independent of the particular network and of `P`.

The same construction over a fixed finite field sums the local
contributions in that field. Over a characteristic-two field it returns
the parity of accepted settings. No nonzeroness-to-existence assertion is
part of this compiler: a separate recovery argument must justify using
that parity.

This is an ordinary proof using the existing grouping theorem. It adds
no Lean theorem and does not establish a factoring exponent above `1/2`.

## 1. Deterministic clustering balances both relevant costs

Put `G=3P^8`. Start with `G` empty clusters. Visit the original gates in a
fixed order and place each gate in a cluster of smallest current total
`a`-weight, breaking ties by index. Immediately before placing a gate of
weight `a_v`, the smallest load is at most `(W-a_v)/G`. Therefore every
final cluster load is at most

```math
Q=\lfloor W/G\rfloor+D.                           \tag{4}
```

This also holds for empty clusters and zero-weight vertices. The
construction uses only the integer metadata, not the exponentially large
local tables.

The sum of setting bits in a cluster is at most its load. Its crossing
wire degree is at most its load as well: an internal wire is removed,
and every remaining incident width was included in (1). Thus the
cluster network satisfies `Fits Q`. Its ordered total wire width is at
most

```math
U=2\sum_e w_e\le W.                              \tag{5}
```

Here and below, ordered total width counts each undirected wire twice.
Confusing this convention with unordered width would invalidate the
direct application of the source grouping theorem.

## 2. Apply the existing grouping theorem once

The pinned declaration `Factoring.Aux.good_exists` applies to every
network on `G=3P^8` vertices with `Fits Q` and ordered total at most `W`.
It supplies a map from clusters to `P` final groups such that each group
contains at most `floor(2G/P)` clusters and each unordered group pair
has total crossing width at most

```math
L=\left\lfloor\frac{PW+GQ}{P^3}\right\rfloor.      \tag{6}
```

Enumerate the `P^G` maps and test these integer inequalities. The source
existence theorem makes the search total. This is a concrete search;
the proof does not use an uncomputed choice of grouping.

Choose `lambda=max(1,L)`. The total setting bits in a final group are
at most

```math
H=\lfloor2G/P\rfloor Q\le 2(W+GD)/P.             \tag{7}
```

Equations (4)--(6) also give

```math
P^2\lambda\le P^2+W+W/P+GD/P.                    \tag{8}
```

Since `P^40<=n`, we have `P^2<=n` and

```math
GD\le3CP^8n^{3/4}\le3Cn^{19/20}\le3Cn.          \tag{9}
```

Substituting (2) and (9) into (7)--(8), and using `P>=2`, proves (3).
The exponent `40` is a convenient sufficient choice inherited from the
existing composition; no optimality claim is made for it.

## 3. Exact tables and the cost of writing them

For each pair of final groups, list its crossing original wires in one
fixed order. Both endpoints use this same order. Concatenate their
deterministic messages and pad with zeros to `lambda` bits. Reversed
group order does not reverse the meaning of an original endpoint's
message.

For a table entry at a group, enumerate the original settings of gates
in that group. Keep a tuple exactly when all local gate weights are one,
all internal wire messages agree, every outgoing message matches its
prescribed label block, and all padding bits are zero. Sum the accepted
tuples over the chosen coefficient ring.

There is no further independent sum over internal wire labels. An
accepted tuple determines every such label uniquely. Expanding the
product of the group tables and then summing the shared external labels
therefore counts every accepted original setting exactly once. This
proves value preservation over the integers and every commutative ring.

One entry visits at most `2^H` tuples, each with polynomial work. There
are `P*2^((P-1)lambda)` entries. By (3),

```math
H+(P-1)\lambda=O(n/P).                           \tag{10}
```

For integer tables every entry uses at most `H+1` bits. A fixed finite
field uses a constant number of bits per entry. Both fit the claimed
bound. The grouping search also fits, because

```math
\log_2(P^G)=3P^8\log_2P\le3P^9\le3n/P.          \tag{11}
```

All wire ordering, integer capacities, table output, and grouping tests
are included in `2^O(n/P)`. A large implicit table is not treated as a
unit-cost input to the evaluator.

## 4. Exposing multiplication outputs without new witnesses

The pinned source uses `L=2^(2s+1)`, `A=2^(s+1)`, and
`M=2^(4L)+1`. Its `faP` and `faQ` gates store the operand digits.
The transform layers compute the ordinary integer convolution, since
every convolution coefficient is strictly below `M`.

A carry gate stores a convolution residue `c_j` and incoming carry
`k_j`. Its deterministic output digit and next carry are

```math
z_j=(c_j+k_j)\bmod2^L,\qquad
k_{j+1}=\lfloor(c_j+k_j)/2^L\rfloor.              \tag{12}
```

The original local weight compares `z_j` to a fixed digit of `N`.
Remove just that comparison to expose `z_j` as an output port. Retain
the first incoming carry zero and final outgoing carry zero. For fixed
operands the transform values, convolution residues, all carries, and
all output digits are still unique. Conversely, the source transform
and carry correctness identify those output digits with the product.
No new free digit variable is introduced.

An output port can be copied to an input port of another multiplication
gadget. If their digit widths differ, split the equality into the
overlapping global bit intervals. A copied bit belongs to exactly one
source output digit and one target input digit. Concatenate multiple
intervals joining the same pair of gates into one wire. The added
unordered width equals the number of copied bits, including each
requested fanout copy. Pinned zero bits and fixed final output bits are
local weight conditions and need no wires.

This explains the interface used by an arithmetic-chain application.
It does not assert that every arbitrary collection of multiplication
gadgets satisfies (2). Each application must prove its total bit budget,
maximum gate load, and uniqueness of its arithmetic witnesses.

## Sources and reproduction

The grouping and deterministic-message representation are due to Isaac
Rudich and Louis-Martin Rousseau, at commit
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1):

- [`EngineGroup.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean), `good_exists`, gives exactly (6), with ordered total width.
- [`Gates.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Gates.lean), `Fac.eval_numTest` and `GNet.value_toNet`, identify accepted message settings.
- [`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean), `mS`, `mV`, and `mW`, give the local data and the output comparison changed in Section 4.
- [`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean), `ib_mid_fa`, `cv_lt`, and `val_cv`, prove the integer convolution computed through the auxiliary residue ring.
- [`MulCount.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean) supplies the source correctness argument for carry settings.

The direct table method is the one in the repository's
[half-rate top constructor](half-top-constructor.md). The change here
is a joint weighted clustering argument for an arbitrary admissible
network and the explicit variable-output port interface. The source's
Apache-2.0 attribution is retained. The new proof and checker used
substantial OpenAI assistant support.

```sh
python3 tools/check_weighted_network_ports.py \
  --output evidence/weighted-network-ports.json
```

The checker replays the actual transform formulas and carry ports,
including a composite auxiliary modulus, compares them to independent
integer multiplication, and checks joined gadgets and rejection controls.
These finite checks support the port implementation; the all-size
compiler follows from Sections 1--4 and the cited source theorems.
