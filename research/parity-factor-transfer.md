# Scalar parity suffices for the factoring transfer

## Result and scope

The [fixed-prime successor](fixed-prime-factor-transfer.md) extends this scalar reduction to every fixed prime characteristic. The characteristic-two construction below has smaller constants and avoids comparison gadgets by using lower-half divisor intervals.

Integer factorization reduces deterministically to a polynomial number
of **scalar characteristic-two counts** of explicitly specified
multiplication predicates. The predicates include repeated squaring,
but their multiplication networks retain total bit width `O(n)`.
They can therefore be converted to written cliques at the same
asymptotic scale as the existing divisor-prefix construction.

This removes the cyclic-ring family and the requirement that an
evaluator be correct over every algebra over its coefficient field.
The evaluator may be an arbitrary algorithm that returns one bit:
the contraction of a written clique over `F_2`.
The result does not supply a faster evaluator or raise the currently
achieved factoring exponent above one half.

**Scalar transfer theorem.** There is an absolute constant `C` with
the following property. Let `E(P,lambda)` bound the running time of
one uniform deterministic finite-alphabet algorithm that evaluates
every explicitly written `P`-vertex clique over `F_2`, with `lambda`
bits on each edge and arbitrary `0/1` vertex tables. Its running time
includes reading the complete input.

For an input `2<=N<2^n`, choose an integer `P>=2` with `P^40<=n`.
There is a deterministic algorithm for complete prime factorization,
and hence for the smallest prime factor, with time at most

```math
\operatorname{poly}(n)\left(
  2^{Cn/P}+E\!\left(P,\left\lfloor\frac{32000n}{P^2}\right\rfloor\right)
\right).
\tag{1}
```

Every clique passed to the evaluator is written out in full. Its
tables contain single bits. The reduction uses no arithmetic over
polynomial rings, no reconstruction of integer counts, and no random
isolation step. A choice of `P` computable in polynomial time gives
one uniform algorithm;
any fixed finite initial range can be handled by trial division.

More generally, the arithmetic argument in Sections 1--3 gives a
polynomial-time oracle reduction from factoring to the parity queries
defined there. Sections 4--7 prove the quantitative network closure
needed for (1), rather than treating those queries as uncharged oracles.

The prior multiplication network and the generic grouping theorem are
from Isaac Rudich and Louis-Martin Rousseau's reference at revision
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The parity witness family, free-output composition, and quantitative
mass argument below are ordinary arguments built on those sources.
No new Lean formalization is claimed.

## 1. Every composite has an odd witness family

For a positive integer `N`, let

```math
N=\prod_q q^{e_q},\qquad
J=\lfloor\log_2(n-1)\rfloor.
```

For `1<=j<=J`, define the nontrivial power-root divisors

```math
R_j(N)=\{a\ge2:a^{2^j}\mid N\}.
\tag{2}
```

Also define the lower-half divisors

```math
L(N)=\{d:2\le d<\sqrt N,\ d\mid N\}.
\tag{3}
```

Every element of either set is a proper divisor of `N`.
For (2), `a^(2^j)<=N` and `j>=1` imply `a<=sqrt(N)<N`.

**Odd-witness lemma.** If `N` is composite, at least one of the
sets `R_1(N),...,R_J(N),L(N)` has odd cardinality. If `N` is prime,
all of these sets are empty.

First suppose some prime exponent is at least two. Let

```math
e_{\max}=\max_q e_q,\qquad j=\lfloor\log_2e_{\max}\rfloor.
```

Then `j>=1` and `2^j<=e_max<2^(j+1)`. An integer `a` satisfies
`a^(2^j)|N` exactly when each exponent in `a` is between zero and
`floor(e_q/2^j)`. These upper bounds are all zero or one, with at
least one equal to one. Consequently

```math
|R_j(N)|
=\prod_q\left(1+\left\lfloor\frac{e_q}{2^j}\right\rfloor\right)-1
=2^t-1
\tag{4}
```

for some positive integer `t`. This count is odd. The level is in
the advertised range: `2^(e_max)<=N<2^n` gives `e_max<=n-1`,
so `j<=J`.

If every positive exponent is one, `N` is squarefree. Write `t` for
its number of distinct prime factors. For composite `N`, `t>=2`.
The involution `d -> N/d` has no fixed point and pairs the `2^t`
divisors. Exactly half lie below `sqrt(N)`. Removing one gives

```math
|L(N)|=2^{t-1}-1,
\tag{5}
```

which is odd. This proves completeness in both cases. For prime `N`,
neither a proper lower-half divisor nor a nontrivial power-root
divisor exists.

The algorithm does not know the factorization, `e_max`, or which case
applies. It tests the stated list of `O(log n)` parity queries.
The exponent formula is used only to prove that the list is complete.

## 2. Odd parity gives a witness, not necessarily the least one

For any of the sets `S` in Section 1, define the prefix query

```math
\pi_S(k,x)
=\#\{a\in S:\lfloor a/2^k\rfloor=x\}\pmod2,
\qquad 0\le k\le n.
\tag{6}
```

After finding a set with `pi_S(n,0)=1`, start with `(k,x)=(n,0)`.
At each step, evaluate the lower child's parity `pi_S(k-1,2x)`.
If it is one, take that child. Otherwise take the upper child
`(k-1,2x+1)`. Since the parent parity is one, the selected child
also has parity one. After `n` steps, the selected single integer
belongs to `S`, and hence is a proper divisor of `N`.

An even child may contain witnesses. Therefore this procedure does
not generally return the least witness. For example,

```math
R_1(36)=\{2,3,6\},
```

and this parity descent returns six: the branch containing two and
three has even parity. Six is a valid proper divisor. It must not be
reported directly as the smallest prime factor.

To factor `N` completely, apply the procedure recursively to both
the returned divisor and its quotient. If every root parity query
is zero, the odd-witness lemma certifies that the current integer
is prime. All descendants can be padded to the original `n` bits;
the same range of levels remains valid.

The leaves are prime factors with multiplicity. Their number is
`Omega(N)<=n-1`, because `2^Omega(N)<=N<2^n`. A full binary split
tree with that many leaves has at most `2n-3` nodes. Thus there are
only `O(n)` calls, each with `O(log n)` root queries and at most
`n` parity-prefix steps. Sort the resulting primes, or retain their
minimum if only the smallest factor is required.

## 3. Lower-half queries use the existing prefix predicate

Compute

```math
B=\lfloor\sqrt{N-1}\rfloor
```

by integer arithmetic. The strict condition `d<sqrt(N)` is exactly
`d<=B`, including when `N` is a square. Intersect `[2,B]` with the
interval specified by the current prefix. Any such integer interval
can be decomposed into at most `2n` disjoint binary-prefix intervals.
For example, repeatedly take the largest aligned power-of-two block
starting at the left endpoint that still fits in the interval.

Sum the ordinary divisor-prefix counts for those blocks in `F_2`.
Disjointness makes the result exactly (6). Empty intersections have
parity zero. This adds only a polynomial number of queries and needs
no new comparison network.

For a power-root query, the multiplication construction below counts
all `a>=1` satisfying `a^(2^j)|N` and the prefix. Subtract the one
known solution `a=1` when its prefix matches. This is one final bit
operation, not an extra network computation. Its cofactor is uniquely
`N`, so exactly one setting is removed.

## 4. Expose the multiplication network's input and output bits

Here are the precise source features used in the composition.
For an intended input width `w>=1`, choose the least power of two
`S>=1` such that `w<=4S^4`, and set

```math
A=2S,\qquad \ell=2S^2,\qquad f=8S^2+1,\qquad C=\ell^2=4S^4.
\tag{7}
```

Thus `C` is the source's input bit capacity, and

```math
w\le C\le16w,\qquad S^4\le4w.
\tag{8}
```

The source has four families of `A` gates and one carry family of
`A^2=2ell` gates. Its two inputs have `ell` digits, each of `ell`
bits. An input gate stores a subset of those digits. The other
gates store received transform residues or an incoming carry.
Messages are deterministic functions of the local setting.

At carry position `t`, the outgoing product digit is

```math
\bigl(\operatorname{carry}_t+c_t\bigr)\bmod2^\ell,
\tag{9}
```

and the next carry is its quotient on division by `2^ell`.
The original fixed-output predicate compares (9) to the corresponding
digit of `N`. To expose a free output, remove that comparison and
retain the checks that the first and final carries are zero. Attach
one-bit ports to all `2C` output bits in (9). Attach one-bit ports
to the `C` bits of each input as well. Pin input bits above the
intended width `w` to zero.

For each pair of input values, all transform messages and carries
are uniquely determined. The outputs are their full product, with
`2C` output bits. The final carry is zero because the product of
two `C`-bit inputs is less than `2^(2C)`. Conversely, any accepted
setting has those forced messages, carries, and output bits.
There is exactly one extension of each valid input pair.

This follows from the actual source definitions and proof:

- [`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
  defines `mS`, `mV`, `mW`, and the digit comparison at each carry gate.
- [`MulCount.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean)
  proves `good_pinned`, identifies the received coefficients through
  `val_fC`, and proves the carry-chain multiplication identity.

The transform's internal arithmetic modulo `2^(8S^2)+1` specifies
finite gate functions. The outer network counts valid settings over
`F_2`. These are different roles; reducing the outer count to parity
does not replace the internal arithmetic used to define a gate.

Connect two bit ports by an equality edge to identify their values.
For squaring, feed the same bit value to both input copies. An output
bit used by the next square has two equality edges, one to each input
copy. Each exposed bit vertex therefore has at most three incident
one-bit edges: one to its arithmetic gate and at most two copy edges.
All bit settings are uniquely forced once the arithmetic inputs are
chosen. Padding bits are pinned to zero, so they create no multiplicity.

## 5. Repeated squaring has total structural size `O(n)`

Fix `r=2^j<=n-1` and put `m=ceil(n/r)`. A root divisor `a` can be
stored in `m` bits: `a^r<=N<2^n` implies `a<2^(n/r)<=2^m`.
Build the sequence

```math
a_0=a,\quad a_{i+1}=a_i^2\quad(0\le i<j),\qquad a_jb=N.
\tag{10}
```

The square at stage `i` has intended input width `2^i m`.
Its full output has width `2^(i+1)m`. The final multiplication
uses input width `max(n,rm)`, padding `b` and `a_j` as needed.
Pin all output bits of that final multiplication to `N`, including
the leading zero bits. Do not truncate an intermediate product.

The sum of the intended multiplication input widths is

```math
\sum_{i=0}^{j-1}2^im+\max(n,rm)
=rm-m+\max(n,rm)<4n,
\tag{11}
```

since `rm<n+r<=2n`. Every individual intended width is less than
`2n`. Equation (11), not merely the number `O(log n)` of
multiplications, is the size estimate used below.

The prefix constraints are checks on the original bits of `a`.
If the requested padded prefix has a nonzero bit above its `m`
available bits, its count is zero immediately. Otherwise pin the
corresponding input bit ports. These checks add no edges or settings.

For every `a>=1` with `a^r|N`, the powers in (10), the cofactor
`b=N/a^r`, and every internal setting are uniquely determined.
Conversely, an accepted setting satisfies the exact equalities (10).
Thus the network counts the desired root divisors once each.

Keeping high bits is essential. For example, `3^4=81` reduces to
17 modulo `2^5`. A five-bit truncation would invent a witness for
the prime 17. The construction retains all eight allocated power
bits for that example and rejects it.

## 6. A quantitative mass bound survives composition

For a vertex `v`, let `d_v` be its incident bit width and let `b_v`
be an explicit integer with at most `2^(b_v)` local settings.
Define its mass as `d_v+b_v`. This pays both for interface width
and for enumerating settings when tables are written.

For the source gadget in (7), the ordered total internal bit width
is at most `10A^2 f`. Use `b_v=2Af` for each non-carry gate and
`b_v=2f` for each carry gate; their sum is at most `10A^2 f`.
These bounds follow from the displayed setting types and the per-gate
width bounds in
[`Cluster.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Cluster.lean).

There are `4C` exposed bit vertices. Each has one setting bit and
at most three incident bits, hence mass at most four. Attaching
their arithmetic edges also adds `4C` to the masses of the original
gates. Charging all possible copy edges through these port-degree
bounds, one composed gadget has total mass at most

```math
20A^2f+20C\le200C\le3200w.
\tag{12}
```

The inequality uses `A^2f<=9C`. It remains valid when some ports
are pinned or unused.

For individual vertices, a non-carry gate has mass at most `6Af`;
a carry gate has mass at most `6f`; and a bit port has mass at
most four. Since `Af<=18S^3` and `S^4<=4w`, every vertex has
mass at most `512w^(3/4)`.

Summing (12) over (11) gives, for the complete power predicate,

```math
\sum_v(d_v+b_v)\le12800n,
\qquad
\max_v(d_v+b_v)\le1024n^{3/4}.
\tag{13}
```

The ordinary multiplication predicate for a lower-half prefix has
the same bounds. The network has polynomially many vertices and
edges, and every local setting and message can be generated and
checked in `poly(n)` bit time. The source's arithmetic is explicit,
including its modular powers; no gate table is assumed to be
available without construction.

## 7. Convert those networks to fully written parity cliques

Set `G=3P^8` and define the following computable integer bounds:

```math
U_0=12800n,\quad
h=\lceil(n^3)^{1/4}\rceil,\quad
Q=\left\lceil\frac{U_0}{G}\right\rceil+1024h.
\tag{14}
```

First assign vertices to `G` clusters greedily, always placing the
next vertex in a cluster of least current total mass. The load of
the selected cluster before placement is at most the average.
Consequently each cluster's final mass is at most total mass
divided by `G`, plus the largest vertex mass. Equation (13) makes
`Q` a valid integer upper bound.

After coarsening, each cluster has incident width at most `Q`, the
ordered total width is at most `U0`, and its setting-count exponent
is at most `Q`. Internal equality edges only remove possibilities;
they do not require independent channel assignments.

Apply the reference's generic
[`EngineGroup.good_exists`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean)
to this `G`-cluster network. It supplies a grouping into `P` groups
such that each group has at most `floor(2G/P)` clusters and each
pair shares at most

```math
\left\lfloor\frac{P U_0+GQ}{P^3}\right\rfloor
\tag{15}
```

bits. Find a grouping by enumerating all `P^G` colorings and
testing these integer width and cardinality conditions. The
theorem is generic for the coefficient ring and network; its use
does not depend on the original multiplication graph's shape.

Let `lambda` be the larger of one and (15). Let
`H=floor(2G/P)Q`, so each final group has at most `2^H` local
setting tuples to enumerate. The cutoff `P^40<=n` gives

```math
G\le3n^{1/5},\qquad h\le2n^{3/4},\qquad Gh\le6n.
```

Using (14),

```math
GQ\le12800n+G+1024Gh\le18947n.
```

Since `P>=2` and `P^2<=n`, these imply the explicit bounds

```math
P^2\lambda\le32000n,\qquad PH\le38000n.
\tag{16}
```

Use the common width `lambda_star=floor(32000n/P^2)` for every
query. It is at least `lambda`. Append zero padding to each
pair's actual message blocks, with the same block order at both
ends. Padding bits are required to be zero in both vertex tables.

For one vertex-table entry, enumerate the original local settings
in its group. Check its local `0/1` weights, every internal message
equality, and every prescribed outgoing bit. Toggle a parity bit
for each accepted tuple. The result is the exact `F_2` table entry.
The product and sum of these tables count each global good setting
once, since every edge's shared label is determined by its equal
messages. This is the same direct-coarsening codec justified in
[the existing top constructor](half-top-constructor.md).

There are `P*2^((P-1)lambda_star)` entries, each constructed by at
most `2^H` polynomial-time checks. By (16), writing all tables costs
`2^O(n/P)` bit operations. The coloring search fits as well:

```math
\log_2(P^G)=3P^8\log_2P\le3P^9\le3n/P.
```

The initial network, all mass bounds, greedy clustering, and each
coloring test use polynomial time and polynomial bit lengths.
There is no real-number operation in the construction: `h` is
computed by integer fourth-root search and every other displayed
algorithmic parameter uses integer arithmetic.

## 8. Complete scalar-evaluator accounting

Each power-prefix query uses one clique evaluation and one known-bit
correction for the root divisor one. Each lower-half prefix query
uses at most `2n` ordinary prefix evaluations. Section 2 uses only
polynomially many such queries, including recursive factorization.
All instances are constructed at the original padded width `n`,
with the same `P` and `lambda_star`.

Their construction cost is `2^O(n/P)` and their scalar evaluation
cost is at most `E(P,lambda_star)`. This proves (1). The usual
sequence of clique sizes and shrinking evaluation coefficients can
therefore be used exactly as in the earlier transfer arguments,
without a separate `exp(O(n/log n))` or `exp(O(sqrt(n)))` term
for nonzeroness testing.

For example, if a proposed scalar evaluator has the required
written-input cost leading to `exp(O(n/(log n)^alpha))` on these
instances, for some fixed `alpha>0`, this reduction retains that
same exponent. The quantitative condition is the actual bound in
(1), including table construction; a field-operation count with
unwritten input tables would not establish it.

The currently available half-rate evaluator can still be used by
lifting the written `0/1` tables to integers and reducing its exact
answer modulo two. This reproduces the current exponent one half.
The new interface also permits a future evaluator designed directly
for scalar parity, even if it has no valid interpretation over
larger coefficient algebras. It does not reduce the internally
added precision in the present Cauchy-deformation evaluator.

## 9. Reproduction and exact controls

Run

```sh
python3 tools/check_parity_factor_transfer.py \
  --output evidence/parity-factor-transfer.json
```

The checker uses explicit integer enumeration as the parity reference
oracle. It does not implement the hypothetical fast clique evaluator.
The recorded run includes:

- all 65,534 integers from 2 through 65,535 in an independent
  prime-exponent check of the odd-witness lemma;
- complete factorization of every integer from 2 through 16,384,
  using padded parity descent and recursive splitting;
- 48,879 direct power-root cardinality checks;
- 483,304 parity-prefix queries;
- 113,618 repeated-squaring width and mass controls;
- 96 deterministic greedy-clustering controls and 189 exact
  large-cutoff constructor inequalities;
- 1,384 independent free-output carry-chain checks;
- explicit negative controls for mistaking parity descent for
  least-witness recovery and for truncating intermediate powers.

The finite evidence verifies the arithmetic and implementation
details of the reduction. The uniform network bound and its
all-input correctness are supplied by the argument above and the
identified source identities. The evidence file records the exact
checker hash and source revision.
