# A written top-clique constructor for the half-rate factoring route

## Statement and scope

Let N >= 2 have binary length n. Let P >= 2 satisfy P^40 <= n. There is one
deterministic algorithm which, given N, P and a divisor-prefix query
0 <= k <= n, 0 <= x < 2^(n-k), writes an Appendix B instance I with P rows.
It also computes positive integers lambda and B such that

~~~math
\begin{aligned}
\operatorname{WFe}(P,\lambda,B,I),\qquad
\operatorname{zval}(P,I)
  &=\#\{d:d\mid N,\ \lfloor d/2^k\rfloor=x\},\\
P^2\lambda&\le4717n,\\
B&\le6769n/P.
\end{aligned}
~~~

The complete construction, including the search for a grouping and writing
every table entry in binary, takes 2^O(n/P) bit operations. All constants
are absolute. They do not depend on P, the recursion level, or the prefix.

This is an ordinary mathematical proof of the top constructor and data
adapter. It uses the pinned factoring reference's multiplication-network
correctness and grouping theorem. It does not call that reference's recursive
counting engine. The separate Appendix B evaluator is needed to turn the
instance into a faster factoring algorithm. No new Lean proof is asserted.

The reference is Isaac Rudich and Louis-Martin Rousseau's
[SubexponentialDeterministicFactoring at c91cc3f](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The target data interface is Anthropic's
[formal-math/eth at 16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Both projects use Apache-2.0. The new argument below concerns their
composition as ordinary algorithms.

## 1. Retain the source's top grouping

Choose the least s >= 0 with n <= 2^(4s+2), and put

~~~math
S=2^s,\quad A=2S,\quad f=8S^2+1,\quad G=3P^8.
~~~

The least-choice rule gives

~~~math
n\le4S^4,\qquad S^4\le4n.
\tag{1}
~~~

For s=0 the second inequality follows from n>=1. For s>0, minimality gives
2^(4(s-1)+2)<n, hence S^4<4n. This slightly sharpens the source's convenient
all-input inequality S^4<=4n+4.

Use the source's multiplication network and its clustering map into G
clusters. Its original gate count is 4A+A^2. The clusters split each of the
four A-element gate families and the A^2-element carry family by the rule
floor(rG/m) for the corresponding family size m.

Write

~~~math
Q=7(\lfloor A/G\rfloor+2)Af
  +3(\lfloor A^2/G\rfloor+2)f,\qquad U_0=10A^2f.
~~~

The source proves that the cluster network has weighted degree at most Q
and ordered total width at most U_0. The theorem good_exists, applied with
G=3P^8, supplies a grouping chi of these clusters into P groups satisfying

~~~math
|\chi^{-1}(a)|\le\lfloor2G/P\rfloor,\qquad
Y_{ab}\le L:=\left\lfloor\frac{P U_0+GQ}{P^3}\right\rfloor
\quad(a<b).
\tag{2}
~~~

Find such a grouping by enumerating all P^G colorings in radix-P order and
testing (2) with integer arithmetic. The test uses only the integer width
matrix and group sizes. It evaluates no table entries. Stop at the first
valid coloring. Existence makes this a total deterministic procedure.

Set lambda=max(1,L). The increase preserves every capacity bound.

Sources: [Cluster.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Cluster.lean)
(clMG, clP, clusterNet_fits, clusterNet_total_le),
[EngineGroup.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean)
(good_exists, firstGood_spec), and
[EngineCostRate.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostRate.lean)
(Ntop_S).

## 2. Build the tables from original gate settings

Every original multiplication gate has a finite set of settings, a 0/1
weight, and one deterministically computed message for each neighbor. An
edge accepts precisely when its two messages agree. This follows directly
from Fac.numTest and eval_numTest, not from a numerical bound on an
unevaluated channel sum.

Let eta send each original gate to its cluster and then to its chi-group.
Use the single coarsening by eta. For each unordered pair of final groups,
list all crossing original edges in increasing order of their original
endpoint pair. Concatenate their messages in this same order at both ends.
An original edge whose orientation is reversed by eta still uses its
original endpoint's message; group order does not select a different
message function.

The width of this concatenation is Y_ab. It equals the width obtained by
first coarsening to clusters and then coarsening by chi: both sums include
each crossing original edge once. The lists may use different block orders,
but the direct coarsening fixes a common order without requiring identical
intermediate table bytes. Append lambda-Y_ab zero bits.

For a final group a and a tuple theta of prescribed labels to other groups,
define its table entry as the number of local original-gate settings for
which:

1. Every local gate has weight 1.
2. Every edge internal to the final group has equal messages at its ends.
3. Every outgoing edge message equals its prescribed block of theta.
4. Every appended padding bit is zero.

This is exactly the source's Net.toClique formula for the directly
coarsened network. To see why, expand its node weight and half-edge factors.
The node weight is the product of the original 0/1 weights and internal
edge equality tests. Each half-edge factor is the product of its original
message tests. A full factor on an internal edge is already an equality
test, because its unique possible channel is the message itself.

There is consequently no additional enumeration over internal channels.
Enumerating each local original-gate setting once and doing these checks
computes the entry exactly.

Expanding the product of all group tables and summing over one shared
label per final group pair now counts each globally valid original setting
once. Its external labels are uniquely determined by its messages and
zero padding. Conversely, every accepted local-setting and shared-label
tuple gives one globally valid original setting. The source's
mulNet_value theorem identifies these global settings with the divisors
having the specified prefix.

Sources: [Gates.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Gates.lean)
(Fac.numTest, eval_numTest), [MulNet.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
(mS, mW, mV), [MulCount.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean)
(mulNet_value), [Net.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Net.lean)
(coarsen, orient, value_coarsen), and
[Clique.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Clique.lean)
(toClique, value_group).

## 3. Local setting counts and entry sizes

The source bounds each non-carry gate's setting count by 2^(2Af), and each
carry gate's by 2^(2f). Its elementary cluster-size argument therefore gives
at most 2^W settings in one cluster, where

~~~math
W=8(\lfloor A/G\rfloor+2)Af
  +2(\lfloor A^2/G\rfloor+2)f.
~~~

A final group contains at most floor(2G/P) clusters. Put

~~~math
H=W\lfloor2G/P\rfloor,\qquad B=H+1.
\tag{3}
~~~

Each entry is a sum of at most 2^H terms in {0,1}. Thus it is a nonnegative
integer at most 2^H, strictly below 2^B. Its binary representation uses at
most B bits. This is a local entry bound; a global O(n)-bit modulus is not
used in the construction.

The exact source inequalities give

~~~math
\begin{aligned}
P^2L&\le720S^4+306GS^3,\\
PH&\le720S^4+648GS^3.
\end{aligned}
\tag{4}
~~~

For the second inequality, use P floor(2G/P)<=2G and

~~~math
GW\le10A^2f+G(16Af+4f).
~~~

Then A^2f<=36S^4, Af<=18S^3, and f<=9S^3 prove (4). The first line is the
source theorem Ntop_S.

Since P^40<=n<=4S^4 and P^4>=4,

~~~math
P^{36}\le S^4,\qquad P^9\le S,\qquad GS^3=3P^8S^3\le3S^4/P.
\tag{5}
~~~

Combining (1), (4), (5), and P>=2 gives

~~~math
\begin{aligned}
P^2\lambda
 &\le P^2+720S^4+306GS^3\\
 &\le n+(720+918/P)4n\le4717n,\\
PH&\le(720+1944/P)4n\le6768n.
\end{aligned}
~~~

Finally, n/P>=1 absorbs the extra 1 in B=H+1. This proves the bounds in
the statement.

Source: [EngineCostTop.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostTop.lean)
(gsz_le_w, Wc, cluster_s_le), with Wc_mul_le in EngineCostRate.lean.

## 4. The exact off-diagonal codec

At a source clique vertex a, a query has P blocks of lambda bits. Its
self block is zero. Appendix B instead passes P blocks whose self block
is the empty word, with all other blocks of length lambda.

For each vertex, order the other vertices increasingly. If their numeric
labels are theta_b, use the row index

~~~math
j_a=\sum_{b<a}\theta_b\,2^{\lambda b}
    +\sum_{b>a}\theta_b\,2^{\lambda(b-1)}.
\tag{6}
~~~

The individual labels and their concatenation use lowest-bit-first order.
For every 0<=j<2^((P-1)lambda), the inverse map takes these successive
lambda-bit blocks and inserts an empty block at a. To obtain the source
query, replace that empty block by lambda zero bits.

Both maps are bijections on legal off-diagonal queries. In a whole
contraction, use the same lambda-bit label on the unordered edge {a,b}
at both endpoints. The source's symmetric zero-diagonal arrays and
Appendix B's edge labellings are then in bijection. Every corresponding
product of table entries is equal. Thus the target zval equals the
source clique value, including over the integers.

The padding test in Section 2 is required at both endpoints. Padding
adds one permitted extension, the all-zero extension, rather than
2^(lambda-Y_ab) choices. Empty final groups have one empty local setting;
their row is 1 when every incident label is zero and 0 otherwise.

The instance is the pair (lambda in unary, the P written rows). Its rows
have exactly 2^((P-1)lambda) entries, and (3) bounds their absolute values.
The target predicate WFe follows directly.

Sources: [PaddedDescent/F/Data.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/PaddedDescent/F/Data.lean)
(PInst, look, unflat, zval, WF),
[ClauseBound/Alg/Interfaces.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Alg/Interfaces.lean)
(WFe), and
[StepMachines/Prog/Chunks.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/StepMachines/Prog/Chunks.lean)
(toWord, ofWord, block).

## 5. Complete bit-operation cost

The following explicit algorithm avoids reliance on an old full-recursion
cost theorem:

1. Compute the parameters, original gate list, cluster map, and width matrix.
2. Enumerate and test the P^G grouping codes until a valid code is found.
3. For each final group and each row index in (6), enumerate all its local
   original-gate settings, apply the checks in Section 2, and increment a
   binary counter for each accepted setting.
4. Append the binary counter to its row, and write the resulting lists in
   the target encoding.

All quantities in step 1 have polynomial bit length in n. The number of
original gates and all index ranges in one gate's message or weight
evaluation are polynomial in n. The messages are computed by finite sums
and modular products in the ring with modulus 2^(8S^2)+1. Repeated squaring
computes its powers in polynomial bit time. In particular the source's
tInv is an explicit power, not a run-time nonconstructive choice. The
message width is O(sqrt(n)). A local setting and its counters use O(n)
bits by Section 3.

Consequently each grouping test and each fixed-setting table check costs
n^c for a fixed c. Using sequential tape arithmetic and appending output,
the complete cost is at most

~~~math
n^c\left(P^G+P\,2^{(P-1)\lambda+H}+1\right).
\tag{7}
~~~

This deliberately pays for every local setting separately for every
written table entry. It requires neither free table access nor unit-cost
arithmetic. No exponent in the polynomial depends on the recursion depth.

Since log_2 P<=P and P^10<=n,

~~~math
\log_2(P^G)=3P^8\log_2P\le3P^9\le3n/P.
~~~

Also, (P-1)lambda+H <= (4717+6768)n/P. Thus (7) is 2^O(n/P).
Its polynomial factors are absorbed because P<=n^(1/40).

The target's tree encoding has four symbols. A list adds two bracket
symbols per cell; a binary integer has its sign and absolute-value bits.
For E=2^((P-1)lambda) and entries of absolute value at most 2^B, its exact
definitions give the coarse bound

~~~math
\operatorname{len}(I)\le
6+4\lambda+P\bigl(4+E(3B+10)\bigr)=2^{O(n/P)}.
~~~

The binary output is therefore charged in (7). Conversion between this
four-symbol encoding and an ordinary finite-tape binary implementation
has a fixed polynomial cost. This finishes the claimed constructor lemma.

Sources: [Transform.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean)
(faOut, midOut, ibOut, tInv),
[ToTuring/Enc.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/Enc.lean), and
[ToTuring/TM/Tree.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/TM/Tree.lean).

## 6. Finite implementation checks

[half_clique_codec.py](../tools/half_clique_codec.py) implements a small
reference for the adapter. Its grouped-table constructor can scatter each
valid local setting into its unique boundary row entry. This is the same
finite sum as step 3, with a different loop order.

Run:

~~~sh
python3 -m unittest tests.test_half_clique_codec -v
~~~

All seven tests passed. They check the off-diagonal index bijection,
arbitrary signed source tables, internal message constraints, reversed
group orientation, empty groups, zero padding, and all 256 binary message
maps of a two-edge path in both grouping orientations. The path contractions
are also checked modulo 4, 6, 8, and 9. A one-endpoint channel-order
corruption is a negative control and changes the expected value.

These checks exercise finite value and codec errors. They do not run the
large factoring algorithm or replace the mathematical cost proof.
