# Factorial normalization leaves an inverse-window count

## The concrete one-variable to two-variable bridge

The [complete one-variable pole trace](circuit-box-pole-trace.md) removes
singular normalization for an ordinary or globally regular circuit in
one variable. Applying it to a product over the other factor interval
suggests removing a large known factorial content first.

That baseline removal is valid. The remaining valuation is only
`O(n)` per row on an `n`-bit factor rectangle without a hit. Its exact
formula, however, is a sum of modular-inverse window indicators. An
increment of the capped valuation recovers the exact factor count.
Thus a fast **valuation-only** aggregation would suffice, but a small
bound on the residual valuation does not itself provide that aggregation.

This note also gives a worst-case width certificate for the direct
dyadic continuation. Even empty rectangles for semiprimes with no small
factor can force `Omega(2^(n/2))` retained states under the established
empty-box, corner and singleton pruning rules. The obstruction concerns
that strategy, not factoring algorithms in general.

## 1. The exact factorial baseline

Let `N,a,b_0` be positive and odd, and let the column progression be

\[
 b_j=b_0+2j,\qquad 0\le j<M.
\]

In the dyadic integers put

\[
 z(a)=\frac{N/a-b_0}{2}.
\]

Then

\[
 \prod_{j=0}^{M-1}(N-ab_j)
 =2^M a^M M!\binom{z(a)}M.
\tag{1}
\]

The known dyadic content is

\[
 B_M=M+v_2(M!).
\tag{2}
\]

It can be exponentially large as a numerical exponent, but its binary
description is short. The remaining binomial is dyadically integral.
The factors `a^M` and the odd part of `M!` are units and do not affect
the valuation question.

For `r>=1`, let `rho_r(a)` be the least nonnegative residue of `z(a)`
modulo `2^r`. The number of indices in the progression with
`z(a)=j mod 2^r` is

\[
 n_r(a)
 =\left\lfloor\frac M{2^r}\right\rfloor
  +\mathbf 1\{\rho_r(a)<M\bmod2^r\}.
\tag{3}
\]

Summing valuations and using Legendre's factorial formula gives, when
the row product is nonzero,

\[
 v_2\binom{z(a)}M
 =\sum_{r\ge1}\mathbf 1\{\rho_r(a)<M\bmod2^r\}.
\tag{4}
\]

This is an exact carry formula, including arbitrary binary lengths `M`.
It uses neither a factorial expansion nor an enormous product residue.

## 2. The nonzero residual valuation has a short height bound

Assume every nonzero row difference satisfies

\[
 |N-ab_j|<2^h.
\]

If the row has no zero, its maximum normalized distance valuation is

\[
 s=\max_j v_2(z(a)-j)<h-1.
\]

For `r>s`, no index satisfies the congruence in (3). Every summand in
(4) is zero or one. Hence

\[
 0\le v_2\binom{z(a)}M\le s<h-1.
\tag{5}
\]

The enormous baseline has therefore been removed exactly, leaving only
a polynomial-size valuation threshold. For `M=2^d` there is a sharper
description. The first `d` tail lengths in (4) vanish. All higher-level
hits belong to the same nearest index modulo `2^d`, and

\[
 v_2\binom{z(a)}{2^d}
 =\max\{0,\max_j v_2(z(a)-j)-d\}.
\]

This does not make the nearest index an ordinary regular function of
the row parameter: it is the principal residue of an inverse.

## 3. Capped valuations recover exactly the factor-count slope

Define

\[
 V_T(a)=\sum_{j=0}^{M-1}\min\{v_2(N-ab_j),T\},
\]

with `min(infinity,T)=T`. It is also the exact content of the displaced
row product at `N+2^T w`, where `w^2+w+1=0`. If

\[
 T\ge1+\lfloor\log_2M\rfloor,
\]

the factorial baseline has entirely appeared, and

\[
 R_T(a):=V_T(a)-B_M
 =\sum_{r=1}^{T-1}\mathbf 1\{\rho_r(a)<M\bmod2^r\}.
\tag{6}
\]

If also `T>=h-1`, then

\[
 R_{T+1}(a)-R_T(a)
 =\#\{j:N=ab_j\}.
\tag{7}
\]

Indeed, the increment counts differences divisible by `2^(T+1)`.
The height bound excludes every nonzero difference. A row contains at
most one exact factor pair.

For any row progression `A`, summing (7) yields

\[
 \sum_{a\in A}R_{T+1}(a)-\sum_{a\in A}R_T(a)
 =\#\{(a,b)\in A\times B:ab=N\}.
\tag{8}
\]

This is the valuation-only version of the global pole-count target.
Computing the two sums modulo a power of two exceeding
`min(|A|,|B|)` suffices to recover the integer count. Exact scalar
valuations have only polynomial binary output length even when their
numerical values are enormous.

For a single supplied row, all residues in (6) are computed by unit
inversion and binary arithmetic in polynomial time. The unresolved
step is summing them over a long row interval.

To see its precise form, write `a=a_0+2i`. Then

\[
 z(a_0+2i)
 =\frac{(N-b_0a_0)/2-b_0i}{a_0+2i}.
\tag{9}
\]

This is a globally regular rational circuit. The operation missing from
the normalized-average theorem is the comparison of its **principal
residue** with `M mod 2^r`. Each such comparison is a modular-inverse
window. It is not an affine floor sum. An ordinary polynomial or
globally odd-unit circuit for the residue selector would encounter the
existing [sharp selector precision bound](regular-analytic-selector-precision.md).

Equations (1)--(9) leave a narrower positive target than evaluating the
entire normalized binomial: aggregate its capped valuation or its final
increment. They do not supply a fast algorithm for that target.

## 4. Exponential width persists on empty semiprime rectangles

The [bilinear pole-tree reference](interval-pole-trace.md) retains the
two parity children that can satisfy the product congruence. It prunes
empty intersections, rectangles whose product corners have the same
strict sign relative to `N`, exact corner hits, and singleton axes.

Let `H=2^s`, with `s>=4`, and take the odd integers in

\[
 I=[H,3H),\qquad J=[3H,9H).
\]

Choose primes

\[
 H/2<p<H,\qquad 8H<q<16H,
 \qquad N=pq.
\tag{10}
\]

Such primes exist for every `H` by Bertrand's postulate; a primary proof
is [Erdős, 1932](https://www.renyi.hu/~p_erdos/1932-01.pdf). In particular,

\[
 4H^2<N<16H^2.
\]

Neither prime factor belongs to `I`, and neither do 1 or `N`. Therefore
the rectangle contains no factor pair. Both prime factors have size
comparable to `H`, so this is not a small-prime or perfect-power example.

At depth `r`, the retained residue pairs satisfy

\[
 ab\equiv N\pmod{2^r},\qquad a,b\text{ odd}.
\]

There are exactly `2^(r-1)` pairs: each odd residue of `a` determines one
odd residue of `b`. For every `1<=r<=s-3`, both interval intersections
are nonempty and contain more than one point.

Put `d=2^r<=H/8`. In every such state, its smallest represented product
is less than

\[
 (H+d)(3H+d)
 \le\frac{225}{64}H^2<4H^2<N.
\]

Its largest represented product is at least

\[
 (3H-d)(9H-d)
 \ge\frac{1633}{64}H^2>16H^2>N.
\]

No corner-sign test can prune any of these states. There is no corner
hit because the rectangle has no factor pair. Thus depth `s-3` alone
contains

\[
 2^{s-4}=H/16
\]

live states, and the complete retained tree through that depth has
`H/8-1` states. Since `N` has at most `2s+4` bits, this is
`Omega(2^(n/2))` for the stated strategy.

The same argument stops earlier if a leaf routine enumerates a fixed
or polynomial number of rows before calling the one-variable oracle.
It still leaves exponentially many retained states, up to that
polynomial factor. Small final output size does not bound the work of
this search: every rectangle in (10) has output zero.

This proof does not rule out a shared state contraction, a different
representation, arithmetic pruning beyond the stated rules, or a
faster factoring algorithm. It rules out claiming that the new
one-variable normalization alone controls the width of the existing
two-variable recursion.

## 5. Exact checks

The [checker](../tools/check_circuit_box_inverse_windows.py) and
[receipt](../evidence/circuit-box-inverse-windows.json) record:

- 400 exact rational identities for the binomial baseline;
- 3,171 capped-content comparisons against literal row products;
- 389 nonzero-row residual height bounds;
- 400 single-row and 400 full-rectangle factor-count slope identities;
- two nonenumerative rows of length `2^120+17`, one containing a factor
  and one containing none;
- 8,178 literal retained states across twelve semiprime rectangles.

The two long rows have the same known baseline, about `2^121`. At caps
122 and 123, the factor row's residual changes from 116 to 117. The
nonfactor row's residual stays at six. These are exact binary-integer
calculations, not evaluations at the enormous baseline precision.

The largest finite width control uses

\[
 H=32768,\qquad N=4302094417=16411\cdot262147.
\]

The checker proves both displayed factors prime by exact trial division.
The specified rectangle has no factor pair, but all 2,048 states at
depth twelve survive the stated pruning tests. It checks all 4,095
retained states through that depth.

```sh
python3 tools/check_circuit_box_inverse_windows.py \
  --output evidence/circuit-box-inverse-windows.json
```

The exact baseline and slope formulas identify the remaining global
operation. The width family closes one direct continuation strategy.
Neither result establishes a stronger complete factoring exponent.
