# Classical small roots give a quarter-bit terminal

## Result and attribution

The [lattice-line terminal](circuit-box-ratio-lattice.md) recovers one
ratio class when its modulus has about one third of the input bits.
Classical small-root recovery reduces this to one quarter. This note
gives explicit integer parameters, handles every factor-size range,
and supplies an exact reference implementation with no external lattice
or polynomial package.

For an odd `n`-bit integer `N`, set

\[
 h=\left\lceil n/4\right\rceil+2,\qquad R=2^h.
 \tag{1}
\]

Given an odd residue `r mod R` and a dyadic band `[2^k,2^(k+1))`,
all divisors of `N` in that band and residue class can be found in
polynomial bit time and space. A ratio class `a=cb mod R`, `ab=N`,
requires at most four such residue queries per band. Streaming all
residues and all relevant bands gives a complete smallest-prime-factor
algorithm in

\[
 \boxed{2^{n/4}\mathrm{poly}(n)\text{ bit time and polynomial bit space}.}
 \tag{2}
\]

This is a classical bound, not a new general factoring result. It is
weaker than the repository's established `1/2` theorem and does not
increase the score. The contribution to this lane is a fully evaluated
terminal with explicit parameter, integer-root, and outer-enumeration
costs.

The underlying theorem is due to Coppersmith, Howgrave-Graham, and May.
[Cohn–Heninger, Theorem 1.1](https://arxiv.org/pdf/1008.1284) states its
monic-polynomial form: small integers `x` with
`gcd(f(x),N)>=N^beta` can be recovered in polynomial time when
`|x|<=N^(beta^2/deg(f))`. Their Theorem 2.1 states the LLL norm bound
used below. The residue-class application is also treated directly by
[Coppersmith–Howgrave-Graham–Nagaraj,
*Divisors in Residue Classes, Constructively*](https://doi.org/10.1090/S0025-5718-07-02007-8).
The present parameter calculation specializes these established results;
it does not assume a new multivariate small-root theorem.

## 1. Every factor size has a valid linear equation

Fix `0<=k<n`, put `A=2^k`, and retain only positive candidates at most
`N`. A divisor in the prescribed residue class has the unique form

\[
 a=r+Rx,\qquad A\le a<2A,
 \qquad 0\le x< X:=2^{k+1-h}.
 \tag{3}
\]

The actual interval for `x` is obtained by ceiling division and can be
shorter than the displayed bound. Empty intervals are discarded. An
interval with one integer is handled by direct division. In
particular, no lattice operation is needed when `k<h`.

Since `N` is odd, `R` is invertible modulo `N`. Let

\[
 f(x)=x+c,\qquad c\equiv rR^{-1}\pmod N.
 \tag{4}
\]

For every genuine divisor `a` in (3), `a` divides `f(x)`: multiply
(4) by `R` and use `Rx+r=a`. Conversely, a recovered root is retained
only after checking that `r+Rx` belongs to the requested band and
divides `N`. Thus the method neither assumes that the divisor is prime
nor mistakes a large common divisor of `f(x),N` for the desired row.

For the standard small-root statement, take `beta=k/n`. Then
`N^beta<=2^k<=a`. Moreover,

\[
 \begin{aligned}
 \log_2 N^{\beta^2}-\log_2 X
 &\ge \frac{(n-1)k^2}{n^2}-(k+1-h)\\
 &\ge \frac{(k-n/2)^2}{n}+1-\frac{k^2}{n^2}>0.
 \end{aligned}
 \tag{5}
\]

The final inequality holds for `k<n`. This checks the small-root
condition across the whole range, including highly unbalanced and
composite divisors. For clarity, the next section gives a concrete
integer lattice and a sufficient dimension bound, so the implementation
does not evaluate irrational powers or use a hidden epsilon choice.

## 2. Integer parameters and a short-polynomial certificate

For integers `1<=m<d`, use the degree-ordered polynomials

\[
 \begin{array}{ll}
 N^{m-i}f(x)^i,&0\le i<m,\\
 x^j f(x)^m,&0\le j<d-m.
 \end{array}
 \tag{6}
\]

Every value at a genuine root is divisible by `a^m`. The coefficient
vectors after substituting `x` by `Xx` form a triangular integer basis
with determinant

\[
 D=N^{m(m+1)/2}X^{d(d-1)/2}.
 \tag{7}
\]

LLL finds a nonzero lattice vector `v` with
`||v||_2<=2^((d-1)/4) D^(1/d)`. The corresponding integer polynomial
`g` obeys

\[
 |g(x)|\le\|v\|_1\le\sqrt d\,2^{(d-1)/4}D^{1/d}
 \qquad(0\le x<X).
 \tag{8}
\]

If the final bound is less than `A^m`, divisibility by `a^m>=A^m`
forces `g(x)=0`. This is the usual short-polynomial step, specialized
here to a linear polynomial and an explicit dyadic band.

Write `L=k+1-h`. Using `N<2^n` and `X=2^L`, a sufficient condition
for (8) is the following integer comparison:

\[
 d^{2d}<2^E,
 \quad
 E=4kmd-d(d-1)-2nm(m+1)-2Ld(d-1).
 \tag{9}
\]

No real logarithm is needed to check (9). The code compares `E` with
the bit length of `d^(2d)`.

There is always a valid choice

\[
 d=2n,\qquad m=2k
 \tag{10}
\]

when a nontrivial interval needs processing. To prove this, subtract
the logarithm of the upper bound in (8) from `km`. With (10), that
margin is at least

\[
 (k-n/2)^2+n(h-n/4-3/2)-h/2+3/4
       -\tfrac12\log_2(2n).
 \tag{11}
\]

Since `n/4+2<=h<=n/4+11/4`, the quantity in (11) is at least

\[
 (k-n/2)^2+3n/8-5/8-\tfrac12\log_2(2n).
 \tag{12}
\]

For `n>=8`, use `log_2(2n)<=n/2`; the margin is at least
`(k-n/2)^2+(n-5)/8>0`. For `n=2,3` there is no nontrivial band.
For `4<=n<=7`, direct substitution checks the seven possible pairs
with `h<=k<n`. Their integer budgets in (9) are
`72,290,276,420,238,350,574`, respectively, and all satisfy (9).

The reference implementation tries dimensions in increasing order and
uses the first pair `(d,m)` satisfying (9). This finite parameter
search performs polynomially many operations and stops by `d=2n`.
It often produces a much smaller lattice. Larger constant guard values
in (1) are allowed and can further reduce the local lattice dimension.
The complete search uses the two guard bits displayed in (1).

The basis has `O(n)` rows and `O(n^2)` coefficient bits. Exact rational
LLL therefore takes polynomial bit time. The implementation recomputes
Gram–Schmidt data after elementary row changes and records the integer
combination producing its first vector. It may stop early when

\[
 \sum_j |g_j|X^j<A^m.
 \tag{13}
\]

Equation (13), nonzeroness, and the recorded lattice combination are
checked exactly. These three checks certify the root implication even
when the full reduced basis has not yet been completed. If early
stopping does not occur, the proved LLL bound guarantees (13) at the
end. There is no heuristic success condition or fixed iteration cap.

## 3. Integer roots are isolated in polynomial time

The auxiliary polynomial has degree below `d<=2n`. Its coefficients
have polynomial bit length. Finding its integer roots by testing every
integer below `X` would lose the bound; the implementation uses a Sturm
sequence and interval bisection instead.

Let `g_d` be the leading coefficient and set `Q=|g_d|+1`. For an
integer endpoint `t`, evaluate the Sturm variation at

\[
 t-1/Q.
 \tag{14}
\]

The fraction in (14) has reduced denominator `Q`, which is larger
than `|g_d|`. By the rational-root theorem, it cannot be a root of
`g`. Thus all interval boundaries are unambiguous, even when `g` has
repeated roots or rational noninteger roots. The number of distinct
real roots between two such endpoints is the difference of their
Sturm variations.

Bisect an integer interval only when that difference is positive. At
a one-integer leaf, evaluate `g` at that integer and retain it exactly
when the value is zero. A real root can pass through at most `O(log X)`
levels, and there are at most `d-1` real roots. Hence only
`O(d log X)` nonempty intervals and their immediate empty siblings are
visited. Euclidean polynomial remainders, positive-content
normalization, and rational evaluations have polynomial bit cost.
The output contains at most `d-1` integer candidates.

All divisor candidates are then checked in the original equation.
The search uses polynomial bit space, including its lattice
transformation certificate and its depth-first root-isolation stack.

## 4. Ratio queries and complete factoring

For a prescribed odd ratio `c mod R`, every factor pair satisfies

\[
 a^2\equiv cN\pmod R.
 \tag{15}
\]

At most four odd row residues solve (15). The complete square-root
lifting routine in the preceding notes computes them in polynomial
time. Apply the present local endpoint to each row residue, then test
the complementary column and the original interval conditions. An
arbitrary positive row interval can be intersected with `[1,N+1)` and
split among at most `n` dyadic bands. Therefore a quarter-bit ratio
query gives an exact positive factor count and a complete local factor
list in polynomial time. No assumption about equal-sized factors is
used.

For the complete smallest-prime-factor search, the implementation
streams odd row residues directly. This has the same asymptotic cost
as grouping them by the ratios in (15). For each increasing band
`[2^k,2^(k+1))` that meets `[2,isqrt(N)+1)`, query the `R/2` odd
residues and retain the least recovered row at most `isqrt(N)`. At the
first band with a candidate, return that least row. All smaller
possible divisors were covered in earlier bands, so it is the smallest
nontrivial divisor and therefore prime. If no band contains a divisor,
return `N`; even inputs return two immediately.

There are `O(n)` bands and `O(2^(n/4))` residue classes. Every local
query, including its output list, has polynomial cost. This proves
(2) for composite, prime, square, and unbalanced inputs. The minimum
and operation counters are streamed, so the complete factor list is
never stored. Repeatedly removing the smallest factor gives a complete
factorization with at most `n` such calls, preserving (2).

The same local implementation also gives an exact rectangle-count
algorithm: partition the row interval into bands, stream residues, and
test whether each complementary divisor belongs to the column
interval. Unlike the cube-root lattice-line proof, this corollary does
not require every product in the original positive rectangle to be
`O(N)`; all genuine rows are first restricted to `[1,N]`.

## 5. Checks and exact certificates

Run:

```sh
python3 tools/check_circuit_box_ratio_small_roots.py \
  --output evidence/circuit-box-ratio-small-roots.json
```

The [reference code](../tools/circuit_box_ratio_small_roots.py),
[checker](../tools/check_circuit_box_ratio_small_roots.py), and
[receipt](../evidence/circuit-box-ratio-small-roots.json) separate:

- exact universal and selected parameter inequalities;
- polynomial integer-root lists checked by literal enumeration,
  including repeated roots and rational noninteger roots;
- a degree-eight root-isolation problem on `[-2^300,2^300)`, containing
  the two integer roots `1` and `2^256`;
- local divisor lists checked by direct enumeration;
- all smallest-prime-factor answers from `2` through `1999`;
- independently reconstructed short-polynomial lattice certificates
  for four known-factorization inputs through 255 bits;
- complete reference factoring examples that actually invoke LLL.

The large local controls use integers `3^u 5^v`. Their supplied prime
factorization independently generates every divisor, so the checker
compares the complete recovered list, including a neighboring residue
class. The largest interval has exponentially many integers, but its
root isolation visits only a polynomial number of intervals. Larger
constant guard values used in some local controls are recorded in the
receipt. Those examples are local endpoint tests, not complete
255-bit factorization benchmarks.

## 6. Remaining global issue

The quarter-bit endpoint solves each prescribed class completely.
The exponential cost in (2) is the explicitly enumerated residue grid.
Ordinary improvements to the LLL approximation constant change the
local overhead; they do not by themselves provide a shared query over
many residue classes. A stronger result for this programme must reduce
that collective search cost or connect it to another proved source of
arithmetic structure. The present theorem does not assert a new
small-root exponent, a new general factoring record, or a better
repository score.
