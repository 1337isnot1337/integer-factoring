# Implicit residue ranges sharpen the inverse precision budget

## Constructive result

A union of binary-encoded odd progressions admits an exact residue-count
range oracle. At each dyadic level it finds the minimum and maximum
occupancy of all odd classes by sweeping only the progression endpoints.
It does not enumerate the residue classes or the matrix nodes.

These ranges give an explicit inverse precision certificate for the
pole Cauchy matrix. They extend the
[joint analytic contraction algorithm](pole-cauchy-analytic-kernels.md)
from single progressions to arbitrary explicit progression unions,
including the balanced padding of every factoring rectangle. For a
balanced square domain of size `m`, the inverse loss improves from
`k-1` to

\[
 \boxed{L=\max(0,k-1-v_2(m)).}                    \tag{1}
\]

The two sides need not be the same set. Balance through the displacement
level is sufficient. For `c_A` and `c_B` progression pieces on the
two sides, even the unbalanced precision loss is at most
`(k-1)(2c_A+2c_B-3)`. Thus the entire regular contraction has polynomial
bit cost in the explicit input representation, independently of the
number of nodes.
The singular factor-count contraction is still unresolved, and the
factoring score remains one half.

## 1. Count all residue classes with a short endpoint sweep

An odd progression is given by `(a_0,m_0)` and contains
`a_0,a_0+2,...,a_0+2(m_0-1)`. At dyadic level `ell>=1`, replace odd
integers by their indices `(a-1)/2` modulo `Q=2^(ell-1)`. Write

\[
 m_0=qQ+r,\qquad 0\le r<Q.                       \tag{2}
\]

This progression contributes `q` to every residue count and one extra
to a cyclic interval of length `r` beginning at `(a_0-1)/2 mod Q`.
Split a wrapping interval at zero. For a union of `c` disjoint
progressions, sum the constant contributions and sweep the at most
`4c+2` interval endpoint positions. The count is constant between
successive integer endpoints, so its exact minimum and maximum are
attained on the swept segments. This includes unoccupied residues.

The work is `O(c log(c+1))` integer operations per level, on integers
whose lengths are polynomial in the input lengths and `ell`. The
number `Q` is not an enumeration cost. The implementation first checks
that the positive odd progressions are ordered and disjoint.

## 2. A computable cofactor bound

Use `M_(a,b)=2^k w/(N+2^k w-ab)` with the usual positive odd inputs,
square domains, and displacement bound. At each level, let

\[
 a_\ell^-\le r_{\ell,c}\le a_\ell^+,
 \qquad b_\ell^-\le s_{\ell,c}\le b_\ell^+         \tag{3}
\]

be the count ranges for the row domain and the reciprocal column
domain. The map `b -> N/b` permutes odd residue classes, so the latter
range equals the range of the original column domain. Both ranges
therefore come from the short sweeps in section 1. Put

\[
 \epsilon_\ell=\max(a_\ell^+-b_\ell^-,
                          b_\ell^+-a_\ell^-),
 \qquad
 L=\sum_{\ell=1}^k\max(0,2\epsilon_\ell-1).
 \tag{4}
\]

Every residue imbalance `Delta_(ell,c)=r_(ell,c)-s_(ell,c)` has
absolute value at most `epsilon_ell`. The
[inverse-entry cofactor identity](pole-cauchy-conditioning.md) gives

\[
 v_2((M^{-1})_{b,a})
 =\sum_{\ell=d+1}^k
   (1-\Delta_{\ell,a}+\Delta_{\ell,N/b}),
 \qquad d=\min(k,v_2(N-ab)).                      \tag{5}
\]

Each summand is at least `1-2 epsilon_ell`. Omitting positive lower
bounds and extending the negative part to all levels proves

\[
 \boxed{2^L M^{-1}\text{ is integral}.}           \tag{6}
\]

Equation (4) is a sufficient certificate. It can overestimate the
actual loss, because it discards which residues realize the extremes
and which ancestor levels are shared by a given entry. Its validity
does not assume that the two residue-count patterns are aligned.

### Every explicit progression union has a polynomial loss bound

For a domain of `c>=1` pieces and total size `m`, (2) shows that
every count lies strictly between `m/Q-c` and `m/Q+c`. Indeed its
constant part is a sum of floors, differing from `m/Q` by less
than `c`, and its extra part lies between zero and `c`.
Consequently the difference between any row count and any column
count has absolute value strictly below `c_A+c_B`. Because these
counts are integers,

\[
 \epsilon_\ell\le c_A+c_B-1.
\]

At level one both counts equal `m`, so the imbalance is zero. Hence

\[
 \boxed{L\le(k-1)(2c_A+2c_B-3).}                 \tag{6a}
\]

The numerical loss is therefore polynomial in the explicit piece
counts and `k` for every accepted domain, without a balance assumption.
The sharper swept value (4) is still useful for reducing the actual
compiler degree. This representation bound was identified during an
independent audit of the residue-range argument.

## 3. Balanced domains save v2(m) further bits

Suppose both domains of common size `m` are balanced at every level
through `k`: each odd class contains either `floor(m/Q)` or
`ceil(m/Q)` points. If `Q` divides `m`, both distributions are
constant, so `epsilon_ell=0`. Otherwise both ranges contain those
two consecutive values, so `epsilon_ell=1`. Therefore the nonzero
terms in (4) are exactly the levels satisfying

\[
 \ell-1>v_2(m).
\]

This proves (1). The sweep certifies balance by checking that each
range has width at most one. Its average is `m/Q`, so this width test
is equivalent to the required floor/ceiling condition.

In particular, if `m=2^s`, the loss is at most `max(0,k-s-1)`.
The existing [balanced padding construction](pole-cauchy-conditioning.md)
produces at most two odd progressions on either side without adding
true factor pairs. These domains now pass directly through the same
effective certificate interface.

## 4. Separate compilation precision from moment precision

Let a regular joint circuit be multiplied by an external scalar
`2^-r`, and request `p>=1` absolute output bits. Once (6) is certified,
compile the circuit pointwise at precision

\[
 Q_{\rm compile}=p+r+L.                           \tag{7}
\]

The compiled integer polynomial differs from the circuit by a
multiple of `2^(p+r+L)` at every odd pair. Contracting the difference
against `M^-1` and multiplying by `2^-r` leaves error divisible by
`2^p`, by (6).

The polynomial moment engine returns the normalization `2^k M^-1`.
Consequently compute its two jets and scalar residue at precision

\[
 Q_{\rm moment}=p+r+k,                            \tag{8}
\]

and return scale `2^(-k-r)`. Equations (7) and (8) have different
purposes and need not be equal. Exact chosen polynomial lifts are
valid even when the compiler has fewer coefficient bits than the
moment engine; the preceding pointwise argument controls their error.

By (6a), all required costs are polynomial in the explicit encoded
domains and the requested precision, including for unbalanced unions.
The algorithm charges the actual value of `L` through (7). An optional
caller-supplied precision cap can reject a computation before allocating
the compiler workspace; it is a resource limit, not an extra theorem
hypothesis. The number of progression pieces is charged as ordinary
explicit input, rather than supplied by a compressed description that
implicitly lists exponentially many pieces.

## 5. Exact controls

The [implementation](../tools/pole_cauchy_domain_precision.py) is
checked by [check_pole_cauchy_domain_precision.py](../tools/check_pole_cauchy_domain_precision.py).
The checker compares 1,080 implicit residue ranges with complete
finite residue enumeration. It checks forty exact dense matrix
inverses, including both balanced padding and unbalanced unions,
and compares forty complete joint contractions against those inverses.
It also checks invalid domains and an exceeded precision cap.

The large example has `N=2^401+45`. It pads an original row progression
of size `2^180+7` against a column progression of size `2^200`, giving
a square domain of size `m=2^200` without creating factor pairs.
The resulting displacement is `k=615`. The certified loss is `414`,
so a twelve-bit result uses 426 compilation bits and 627 moment bits.
The joint circuit from the preceding note compiles to total degree
twelve. An independent sixteen-bit run agrees at the requested
precision. Neither run enumerates a matrix node.

The [receipt](../evidence/pole-cauchy-domain-precision.json) retains
every level's count ranges, the precision budgets, source hashes,
and both large residues.

```sh
python3 tools/check_pole_cauchy_domain_precision.py \
  --output evidence/pole-cauchy-domain-precision.json
```
