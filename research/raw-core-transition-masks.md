# Constant transition masks on an extreme Gaussian core

## Construction

Let \(G_M(x,y,u,v;w)\) be the raw binary Gaussian on a
power-of-two interval of length \(M\). Retain its complete physical
alphabet and let \(c(x)\in\{0,1,2,3\}\) record the occupations
of modes \(0\) and \(M-1\). Choose a constant mask

\[
m:\{0,1,2,3\}^4\longrightarrow\mathbb F_4,
\qquad m(a,b,b,a)=1\quad\text{for all }a,b.
\tag{1}
\]

Define the polynomial tensor

\[
T_m(x,y,u,v;w)
=m(c(x),c(y),c(u),c(v))G_M(x,y,u,v;w).
\tag{2}
\]

Its maximum entry degree is no larger than the Gaussian's. Binary
masks give binary tensors. General masks use only the fixed extension
\(\mathbb F_4\). The endpoint, degree accounting, and finite
alphabet are unchanged by replacing a charge indicator with (1).
This is a search family; no improved factoring exponent is asserted.

## Exact typed endpoint

The raw Gaussian has the typed crossing endpoint. In a tensor power,
each coefficient of (2) is the corresponding Gaussian coefficient
multiplied by the product of its local masks. Therefore an off-diagonal
typed endpoint coefficient, which was already zero, remains zero. On
a diagonal typed coefficient, each local tuple is \((x,y,y,x)\),
so every factor of the mask equals one by (1).

Both local flattenings at the endpoint are triangular after fixed
coordinate permutations. Their diagonal entries are the retained
tuples in (1). The masks preserve their triangular support and their
unit diagonal, hence both polynomial flattenings have full generic
rank. This rank statement does not rely on truncated series detecting
the absence of later pivots.

## Bounded exact profiles and arbitrary core gauges

Use the regular creation and annihilation support reductions from
the [raw extreme charge construction](raw-extreme-charge-projections.md).
For \(M\ge8\), both cuts have ten active Gaussian modes and
explicit Gaussian spectators. The full physical core algebra acts
within that active space, independently of the interval length.

The regular physical core coordinate bases intertwine every core
occupation projector. Consequently an arbitrary mask in (1) acts
on the cached Gaussian coordinate matrix by an exact entry mask:
the entry's two row labels and two column labels specify the four
arguments of \(m\). No new phase decomposition or matrix-series
construction is necessary for each mask.

These coordinate bases also permit arbitrary diagonal gauges on
each of the four physical core alphabets, at both places. All
noncore Gaussian modes and all physical alphabet factors remain in
the pressure calculation. The active matrix has order at most
\(2^{10}=1024\); it does not replace the complete physical alphabet.

Every matrix entry still conserves the integer Gaussian charge.
Thus the matrix splits into exact square charge blocks before
Smith elimination. For ten active modes the largest such block
has order \(\binom{10}{5}=252\). Every alleged off-block entry
is checked to vanish exactly. The generic full-rank statement then
gives full rank for each square block, and the complete Smith
profile is the sorted union of their profiles. Arbitrary diagonal
physical gauges preserve this decomposition.

## Two finite mask subfamilies

There are ten unordered pairs \(\{a,b\}\) from the four core
states, allowing repetitions. Any partition \(\Pi\) of these
ten pairs gives the admissible binary mask

\[
m_\Pi(x,y,u,v)
=\mathbf1[\{x,u\}\text{ and }\{y,v\}
\text{ lie in the same part of }\Pi].
\tag{3}
\]

Charge projections are a special case: partition pairs by
\(f(a)+f(b)\) in a finite abelian group. Finite enumeration over
the prime moduli \(3,5,7,11,13\) produces 104 distinct pair
partitions. This number is a finite enumeration result, not a
completeness claim for all groups. Binary partitions of the ten
pairs supply another small family that does not require additive
charges on individual states.

The fully separated partition keeps precisely the transitions with
\(\{x,u\}=\{y,v\}\). The single-part partition recovers the
Gaussian. These provide useful controls at both ends of the family.

## Implementation and replay

The constructor is
[`raw_core_mask_profiles.py`](../tools/raw_core_mask_profiles.py).
The separate exact block wrapper is
[`blocked_core_charge_pressure.py`](../tools/blocked_core_charge_pressure.py).
It leaves the earlier constructor and saved charge evidence unchanged.

The checker
[`check_raw_core_mask_profiles.py`](../tools/check_raw_core_mask_profiles.py)
compares both original polynomial flattenings with the bounded
active construction under quarter-integer physical gauges. It
includes binary pair masks and a general \(\mathbb F_4\) mask,
plus literal two-copy typed endpoint coefficients. Saved results
are in
[`raw-core-mask-profiles.json`](../evidence/raw-core-mask-profiles.json).

```sh
python tools/check_raw_core_mask_profiles.py --output evidence/raw-core-mask-profiles.json
```
