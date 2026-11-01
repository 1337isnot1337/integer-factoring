# Positive ratio lookahead gives a complete adaptive search

## Constructive improvement inside the inverse-window representation

The [weighted Kloosterman blocks](circuit-box-kloosterman-blocks.md)
evaluate a prescribed fine ratio class in polynomial bit time. Searching
all such classes costs `2^(n/2)*poly(n)` in the balanced factoring
regime. A positive lookahead count permits a smaller adaptive search.

For an odd `N`, an odd ratio `c mod 2^h`, and two ordinary intervals,
define

\[
 T_h(c)=\#\{(a,b)\in I\times J:
       ab\equiv N\pmod{2^{2h}},\quad
       a\equiv cb\pmod{2^h}\}.
 \tag{1}
\]

The count in (1) is computed exactly by at most four square roots and
four affine interval-congruence counts. Each step of the search adds
one ratio bit while checking two product bits. A class with zero count
has no surviving descendant and is discarded.

Suppose `I,J` contain positive integers, their endpoints have `O(n)`
bits, `n=bits(N)`, and

\[
 \max_{a\in I,b\in J}ab\le C N
 \tag{2}
\]

for a fixed constant `C`. The complete rectangle factor count is then
computable in

\[
 \boxed{2^{\,n/3+O(n/\log(n+2))}}
 \tag{3}
\]

finite-alphabet bit time and polynomial bit space. A full smallest-prime-
factor algorithm follows with the same bound by charged rectangle and
prefix queries. All counts, roots, endpoints, and intermediate integers
have polynomial bit length.

This is a faster complete algorithm for this inverse-window search
representation. It is weaker than the repository's established
`exp(O(n/sqrt(log n)))` factoring theorem, so its score remains `1/2`.
It is also weaker than older general factoring algorithms. In particular,
a cube-root elementary-operation bound already appears in Lehman's
[*Factoring Large Integers* (1974)](https://doi.org/10.1090/S0025-5718-1974-0340163-2),
and Harvey–Hittmeir's
[*A log-log speedup for exponent one-fifth deterministic integer
factorisation*](https://arxiv.org/abs/2105.11105) has a stronger
general bit bound. No historical improvement to general factoring is
claimed here. The result supplies a positive, fully evaluated and
analyzed search within the present arithmetic programme.

## 1. Two product bits per ratio bit

Put `R=2^h`. Any pair counted by (1) is odd. Its low row residue
`r=a mod R` satisfies

\[
 r^2\equiv cN\pmod R.
 \tag{4}
\]

There are at most four such odd residues, and successive binary lifting
constructs their complete list in polynomial bit time. For each root,
put

\[
 s=Nr^{-1}\bmod R,
 \quad
 U_r=\{x:r+Rx\in I\},
 \quad
 V_s=\{y:s+Ry\in J\}.
 \tag{5}
\]

These are ordinary intervals obtained by ceiling division. The original
pair has the unique form `a=r+Rx`, `b=s+Ry`. Expanding the product gives

\[
 ab=rs+R(sx+ry)+R^2xy.
\]

Hence `ab=N mod R^2` is equivalent to

\[
 sx+ry\equiv\frac{N-rs}{R}\pmod R.
\]

Because `s` is odd, divide by it. Equation (4) gives `r/s=c mod R`, so
with

\[
 t_r=\frac{N-rs}{R}s^{-1}\pmod R,
\]

the exact count is

\[
 \boxed{T_h(c)=\sum_{r^2\equiv cN\ (R)\atop0\le r<R}
       L_R(c,t_r;U_r,V_s),}
 \tag{6}
\]

where `L_R` counts `x+cy=t_r mod R` between the two displayed intervals.
The roots are odd automatically. The interval and Euclidean floor-sum
formula for `L_R` is proved in Section 4 of the preceding note. It
handles signed interval endpoints and complete periods without
enumerating their elements. The implementation of (6) uses that exact
routine.

This counts a positive predicate. The signed frequency quantity
`2 C_m-C_(m-1)` in the preceding note cannot generally be used for
early pruning: its vanishing can result from cancellation. Equation
(6) removes that issue by computing the stronger product congruence
itself.

## 2. A nested positive search

At level one, the single odd ratio class is `c=1 mod 2`. A live class
`c mod 2^h` has two possible children:

\[
 c\pmod{2^{h+1}},\qquad c+2^h\pmod{2^{h+1}}.
\]

Their counts satisfy

\[
 T_{h+1}(c)+T_{h+1}(c+2^h)
 =\#\{(a,b)\in I\times J:
     ab\equiv N\pmod{2^{2h+2}},\ a\equiv cb\pmod{2^h}\}
 \le T_h(c).
 \tag{7}
\]

In particular, a zero parent has two zero children. Test a node using
(6); discard it if its count is zero, otherwise visit both children.

Let

\[
 B=\max_{a\in I,b\in J}|ab-N|,
 \qquad H=\max\{1,\lceil\mathrm{bits}(B)/2\rceil\}.
 \tag{8}
\]

The four rectangle corners determine `B`; for positive intervals it
suffices to use the minimum and maximum products. The convention is
`bits(0)=0`. At depth `H`, `2^(2H)>B`. Therefore the product congruence
in (1) is exact equality, and summing the terminal counts gives

\[
 \sum_{c\bmod2^H\atop c\text{ odd}}T_H(c)
       =\#\{(a,b)\in I\times J:ab=N\}.
 \tag{9}
\]

All true factor pairs survive every test. Every terminal pair is a
true factor pair. Distinct ratio classes are disjoint. These facts
prove correctness without a heuristic distribution assumption.

An empty interval gives zero directly. If `N` is outside the minimum
and maximum positive rectangle products, the count is also zero
directly. The complexity analysis below concerns the remaining case.

## 3. A uniform live-class bound

Let `P` be the largest positive product in the rectangle, and put

\[
 \Delta(P)=\max_{1\le v\le P}\tau(v).
 \tag{10}
\]

At level `h`, all the counted products belong to the progression

\[
 v=N+z4^h,\qquad 1\le v\le P.
 \tag{11}
\]

There are at most `P/4^h+1` such **positive** integers. No `tau(0)`
occurs. The offset `z=0`, when in range, means `v=N` and is included
with the ordinary finite bound `tau(N)`. Positivity of both rectangle
coordinates is essential to this step. Arbitrary signed boxes from the
local identity (6) do not satisfy this particular global proof.

Each positive `v` has exactly `tau(v)` ordered positive factor pairs
before the rectangle restrictions. Thus the number `A_h` of live ratio
classes satisfies

\[
 \boxed{A_h\le
   \min\left\{2^{h-1},
        \left(\frac{P}{4^h}+1\right)\Delta(P)\right\}.}
 \tag{12}
\]

The second bound is on actual integer products, not an estimate for
equidistribution modulo a power of two. Every live class contains a
counted pair, and different classes contain disjoint pairs.

### The divisor bound used here

The standard elementary uniform bound is

\[
 \log_2\Delta(P)=O\left(\frac{\log_2 P}{\log\log(P+2)}\right).
 \tag{13}
\]

For completeness, let `L=bits(P)`, and split the prime factors of
`v<=P` at `Y=2^u`, where `u` is of order `(log_2 L)/2`. There are at
most `Y` primes below that threshold, and each exponent is at most
`L`. Their divisor contribution is at most `(L+1)^Y`. The total
number of prime factors above `Y`, counted with multiplicity, is at
most `L/u`; using `e+1<=2^e` bounds their divisor contribution by
`2^(L/u)`. Therefore

\[
 \log_2\tau(v)\le Y\log_2(L+1)+L/u
       =O(L/\log(L+2)).
\]

The first term is `O(sqrt(L) log(L+2))`, which is smaller than the
second asymptotically; bounded small `L` can be absorbed into the
constant. This proves (13) without a prime-number estimate or a
factorization supplied to the algorithm. The repository's
[optimized divisor bound](optimized-divisor-count-bound.md) supplies
sharper finite bounds when those constants are useful.

### Summing the visited nodes

Every tested nonroot node has a live parent, so the total number of
tested nodes is at most

\[
 1+2\sum_{h=1}^{H-1}A_h.
 \tag{14}
\]

Set `d=log_2 Delta(P)` and split the sum near
`h0=ceil((log_2 P+d)/3)`. Below this level the first part of (12)
forms a geometric sum. Above it, use the second part. This yields

\[
 \sum_{h=1}^{H}A_h
 =O\left((P\Delta(P))^{1/3}+H\Delta(P)\right).
 \tag{15}
\]

The same bound holds if the splitting point lies beyond the terminal
depth, by using only the first geometric sum. Under (2),
`log_2 P=n+O(1)` and `H=O(n)`. Combining (13)–(15) gives

\[
 \#\text{tested nodes}=2^{n/3+O(n/\log(n+2))}.
 \tag{16}
\]

Every tested node uses `O(n)` binary square-root lifting steps, a
constant number of inversions on `O(n)`-bit residues, and at most four
floor counts. Each floor count performs `O(n)` Euclidean divisions
on integers of `O(n)` bits. The interval counts and every intermediate
coefficient also have `O(n)` bits. School arithmetic therefore gives
polynomial bit cost per node; multiplying (16) by that polynomial
preserves (3). In particular, the floor-count arithmetic is not free.

A depth-first traversal stores only one pending sibling per level.
There are `O(n)` levels, and a node label has `O(n)` bits. The
reference implementation retains `O(n)` level totals as well. These
structures and the arithmetic scratch space require `O(n^2)` bits,
up to a conventional polynomial overhead for the chosen elementary
integer routines. No list of all live nodes is needed.

## 4. An arbitrary factor and the smallest prime factor

A nonzero terminal class permits recovery of some factor pair by
bisecting the row interval and reusing (6) at the fixed terminal
precision. The partner is obtained by exact division `b=N/a`. Such
a divisor need not be the smallest prime factor: an arbitrary
successful rectangle or ratio class can contain a larger or composite
divisor.

The complete reference algorithm obtains the smallest prime factor
directly. Return two for an even input. For odd `N`, visit powers of
two `A=2,4,8,...` in increasing order, with

\[
 I_A=[A,\min(2A,\lfloor\sqrt N\rfloor+1)),\qquad
 J_A=[\lfloor N/(2A)\rfloor+1,\lfloor N/A\rfloor+1).
 \tag{17}
\]

Every divisor from two through `sqrt(N)` lies in exactly one first
interval and has its partner in the associated second interval.
Every positive product in these rectangles is strictly below `2N`,
so (2) holds with the fixed constant `C=2` for every query. All
endpoints and their square-root construction have polynomial bit cost.

Apply the adaptive exact count to the rectangles in this order. If
all counts vanish, `N` is prime and is its own smallest prime factor.
Otherwise, the first nonzero rectangle contains the smallest divisor.
Binary-search its row interval using the complete adaptive count on
each prefix. These restricted rectangles still satisfy `ab<2N`.
The resulting least divisor is prime: a composite least divisor would
have a smaller prime divisor of `N`.

There are `O(n)` rectangle queries before the first hit and `O(n)`
prefix queries afterward. Multiplying (3) by these polynomial factors
does not change its exponent. All queries are charged in the checker.
This route avoids an implicit assumption that an arbitrary recovered
divisor is already prime; no uncharged factoring recursion is used.

For a complete prime factorization one may repeatedly apply the
smallest-prime-factor routine to the remaining cofactor. There are at
most `n` nontrivial output steps, each on at most `n` bits, so the same
upper bound holds after another polynomial factor. That optional
extension is separate from the implemented smallest-factor endpoint.

## 5. Why the earlier empty-rectangle example is no longer fatal

The [inverse-window audit](circuit-box-inverse-windows.md) proved a
large width for a tree that retained the current product congruence and
used empty intervals, corner inequalities, and singleton tests to
prune. Its explicit rough-semiprime family was not a lower bound on
other pruning rules.

The present algorithm checks the extra `h` product bits in (1) by a
whole affine count at each depth `h`. No individual high-bit
completion is enumerated. This stronger exact test removes classes
that have many geometric candidates but no pair satisfying the
stronger congruence. Bound (12) charges the survivors by nearby
positive integer products and their divisor counts.

The checker uses the same family: `S=2^s`, primes
`S/2<p<S` and `8S<q<16S`, `N=pq`,

\[
 I=[S,3S),\qquad J=[3S,9S).
\]

Its only positive divisor rows are `1,p,q,N`, all outside `I`, so the
rectangle count is independently zero. Its product corners straddle
`N`, and the maximum product is below `7N`, satisfying (2) with a
fixed constant. The primes in the controls are certified by complete
trial division through their square roots.

| Input bits | Adaptive tested nodes | Fixed terminal ratio classes | Earlier tree states certified through depth `s-3` |
| --- | ---: | ---: | ---: |
| 23 | 169 | 4,096 | 127 |
| 27 | 513 | 16,384 | 511 |
| 33 | 1,917 | 131,072 | 4,095 |
| 39 | 7,161 | 1,048,576 | 32,767 |
| 45 | 29,049 | 8,388,608 | 262,143 |
| 51 | 115,929 | 67,108,864 | 2,097,151 |

The last column is a proved partial-tree lower count for the earlier
rules, not its complete measured running time. The adaptive column
counts every visited ratio node through termination. These finite
comparisons support the representation change; the uniform theorem
is (12)–(16).

## 6. Exact validation and next target

The [implementation](../tools/circuit_box_ratio_lookahead.py),
[checker](../tools/check_circuit_box_ratio_lookahead.py), and
[receipt](../evidence/circuit-box-ratio-lookahead.json) record:

- 800 literal ratio counts and 800 positive child-nesting identities;
- 400 literal adaptive rectangle counts;
- 859 complete level profiles, including exact live-class counts,
  exact congruence-pair totals, and the parent-to-child query count;
- 859 positive-product divisor-sum upper bounds, with zero explicitly
  excluded from the integer product list;
- smallest-prime-factor agreement with independent trial division for
  every input from 2 through 1,199;
- four huge interval controls through ratio precision 512, using
  complete residue periods and four affine counts;
- six independently empty semiprime rectangles through 51 bits;
- complete smallest-prime-factor runs at 21, 27, and 33 bits, with all
  rectangle queries, prefix queries, ratio nodes, and floor-count calls
  retained.

The 33-bit factoring control is
`4,297,850,923 = 65,537 * 65,579`. It returns `65,537` after 20 complete
rectangle/prefix queries, 10,017 visited ratio nodes, and 23,545 affine
count calls. This is a finite correctness control, not a claim of
practical competitiveness with general factoring software.

```sh
python3 tools/check_circuit_box_ratio_lookahead.py \
  --output evidence/circuit-box-ratio-lookahead.json
```

Use `--small-only` to omit the 39-, 45-, and 51-bit empty-rectangle
searches. The full receipt includes them. No Lean work is involved.

The current lookahead checks product precision exactly twice the ratio
precision because the remaining bilinear term then vanishes. A useful
further improvement must either check substantially more product bits
per class without paying for their high-bit completions, or combine
many of the remaining positive classes with a charged collective
operation. Replacing the extra bilinear congruence by a supplied oracle
would not establish either improvement.
