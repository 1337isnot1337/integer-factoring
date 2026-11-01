# A complete crossing resolution by central fusion occupancy

The earlier [graded prefix projectors](invertible-fusion-graded-seams.md)
provide tunable ranks but miss the crossing in their complete small spans.
This note resolves every central pair into its vacuum and triplet channels,
then groups all placements by occupancy. The resulting polynomially many
pieces have an exact crossing identity at every length. A favorable complete
rank sum, needed for a stronger factoring bound, remains unproved.
The positive open question is cancellation between these complete sectors,
using the [compatible spectral moves](fused-spectral-transport.md).

## 1. Orthogonal occupancy sectors

Work in the same height-three vacuum quotient over F16, with q=8,
delta=q+q^-1. Put

```math
p_0=\delta^{-1}e,\qquad p_2=I+p_0.
```

The TL relation gives p0^2=p0, p2^2=p2, p0 p2=0 and p0+p2=I.
The central seam consists of the n disjoint adjacent pairs of 2n
strands. For a subset S of these pairs define P_S to use p2 at S
and p0 elsewhere. All P_S are mutually orthogonal idempotents and
sum to I. Define

```math
Q_k=\sum_{|S|=k}P_S,
\qquad \prod_{j=0}^{n-1}(p_{0,j}+u p_{2,j})
       =\sum_{k=0}^n u^k Q_k.                           \tag{1}
```

Thus the Q_k are orthogonal idempotents giving a complete resolution,
using n+1 sectors instead of listing 2^n placements. This grouping is
valid in characteristic two: different projector images form a direct
sum. An even number of equal-dimensional images does not cancel.

Let d_0=1 and d_k=Fibonacci(k-1) for k>=1; in particular d_1=0.
Then

```math
\boxed{\operatorname{rank}Q_k={n\choose k}d_k.}         \tag{2}
```

Indeed each p0 pair carries the tensor unit, and each p2 pair carries
object2. The [retained local fusion resolutions](fusion-projector-seams.md)
give 2 tensor2=0 direct-sum2, with unit Gram factors. The rank of P_S
is therefore the vacuum multiplicity d_k, independently of placement.
Orthogonality adds those ranks as integers. This is a vacuum-quotient
rank assertion, not the rank of an ordinary ambient tensor product.

Completeness gives the exact integer identity

```math
\sum_{k=0}^n {n\choose k}d_k=Fibonacci(2n-1).            \tag{3}
```

The right side counts all 2n-strand return paths. Equivalently, the
binomial transform of the two-state fusion matrix has eigenvalues
phi^2 and phi^-2. No characteristic-zero category classification is
used in the rank argument.

## 2. Complete two-parameter pieces and their endpoint

Use the anti-diagonal block-braid order of the
[nilpotent seam compiler](nilpotent-seam-crossings.md), at even n.
Replace its whole central layer by (1), its n marked boundary gates
by qI+y e, and keep qI+e at the other cells. Include the scalar
q^(-n-n^2/2), and take vacuum boundary matrix coefficients. Write

```math
T(u,y)=\sum_{k,b=0}^n T_{k,b}u^k y^b.
```

The two seams are disjoint at even n. Since

```math
qI+e=q^{-1}p_0+q p_2
     =q^{-1}(p_0+q^2p_2),
```

the exact block-braid endpoint gives

```math
\boxed{\sum_{k,b=0}^n q^{2k}T_{k,b}=E,\qquad
E(P,Q,R,S)=[P=S][Q=R].}                               \tag{4}
```

The factor q^-n precisely restores the n central crossings; q^(-n^2/2)
removes the vacuum block phase. All coefficients and weights are explicit,
and the right side has scale one. This supplies a uniform complete
crossing resolution, rather than an incomplete corner candidate list.

A strict univariate family can also be obtained without losing these
pieces. Use u=q^2 t, y=t^(n+1), and then
G(t)=(1+t)I+t T(q^2 t,t^(n+1)). Its endpoints are I and E.
The b-windows for the T coefficients are disjoint, and the extra I
is explicit; the degree is O(n^2). The direct resolution (4) is simpler
when only the pieces are needed.

The highest occupancy k=n agrees, up to the scalar q^(-2n), with the
[all-triplet seam](fusion-projector-seams.md). At k=0 the central
rank is one; k=1 vanishes identically. All intermediate occupancies
are retained in (4).

## 3. Both cut bounds, with every convolution channel

Let h_(n,b)=min(b+1,n-b+1) and
s_(n,k)=min(n/2,k)-max(0,k-n/2)+1. The same-piece bounds are

```math
r(T_{k,b})\le h_{n,b}{n\choose k}d_k,
\qquad
c(T_{k,b})\le s_{n,k}{n\choose b}d_{2b}.               \tag{5}
```

Both bounds may be capped by L^2, where L=Fibonacci(n-1) is the
boundary label count. To prove the row bound, convolve the two outer
boundary-degree halves around Q_k; there are at most h_(n,b) terms,
each factoring through its rank (2).

For the column bound, split Q_k across the two physical blocks:
Q_k=sum_i Q_i^(left) tensor Q_(k-i)^(right). Exactly s_(n,k) occupancy
splits are feasible. For each, fix one of the binom(n,b) marked
boundary subsets. The [vacuum-interface lemma](fusion-projector-seams.md)
factors that tangle through its d_(2b) vacuum interface. Summing all
those channels proves (5). A single triplet in a block cannot simply
be dropped: the open interface can carry its nonzero charge.

The row-count entropy at k/n=p is H(p)+p log(phi), maximizing at
p=1/phi. The peak rank has exponent two in L, up to polynomial factors.
Thus orthogonal occupancy completeness does not itself make the middle
channels cheap. The desired square-root rank sum remains open; (4)
should not be presented as an improvement beyond the complete score1/2.

For exactly 2^lambda labels choose the same least even n and common
boundary subset as in the [path compiler](quantum-braid-path-quotient.md).
Then n=O(lambda), L=Theta(2^lambda), restriction preserves (4) and
cannot increase either rank. The explicit path computation and dense
rank factors cost 2^O(lambda), with only O(n^2) pieces. The fixed-field
[whole-tensor descent](fixed-field-descent.md) remains available. A
favorable complete rank-sum estimate is still required before the
factoring transfer yields any stronger score.

## 4. Evidence and next step

The [implementation](../tools/quantum_braid_fusion_occupancy.py) constructs
all sectors and all complete boundary coefficients. Its
[checker](../tools/check_quantum_braid_fusion_occupancy.py) compares
literal placement projectors, exact ranks, full crossing endpoints,
complete two-cut profiles, and independent scalar path evaluations.
The [receipt](../evidence/fusion-occupancy-resolution.json) pins its sources.

```sh
make verify-braid-fusion-occupancy
```

Next: use spectral transport to couple occupancy sectors to boundary
channels, or construct a genuinely different complete family. Filtering
out the costly middle pieces from (4) is not a proof of completeness.
No new Lean formalization or stronger full factoring bound is claimed.
