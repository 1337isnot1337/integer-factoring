# Fixed-field count values and a boundary parity provider

## Two results and their scope

Evaluating the globally cleared inverse at **zero** gives a common
odd unit times the number of proper factor rows. That unit is the
same for every threshold. Differences of two values therefore give
an interval count up to a common unit. Using the prime-only windows
from the [single-characteristic transfer](one-prime-bezout-recovery.md),
only `p=ceil(log_2(n+1))` bits are needed over the **fixed** quadratic
2-adic field, with `O(n)` scalar queries and polynomial overhead.

A uniform provider with cost `poly(n) 2^O(p)` would consequently
give polynomial-time deterministic factoring. Even the weaker cost
`poly(n) 2^O(2^p)` would give `exp(O(sqrt(n)))` factoring after
roughness preprocessing. These are explicit sufficient provider
targets, not new complete bounds. The general scalar evaluator
remains open. A simpler direct suffix test below also records the
initial `p=ceil(n/s)+1` precision tradeoff.

There is also an **implemented restricted provider** at precision one.
For two ordinary odd intervals whose symmetric difference has `k`
nodes, it computes the true normalized scalar, and the complete
cleared pair modulo two, in `O(k poly(n))` bit operations. It visits
only the boundary nodes. Identical intervals require only a square
test. This supplies genuine singular examples on implicit large
domains, while leaving the general suffix-count problem open.

Files: [implementation](../tools/bezout_count_values.py),
[checker](../tools/check_bezout_count_values.py), and
[evidence](../evidence/bezout-count-values.json).

## 1. The fixed-field value contract

Let `N>=3` be odd with bit length `n`. Fix an odd threshold
`3<=a_0<=N` and an integer precision `p>=1`. Put

\[
 u=n+p,\quad m=2^u,\quad
 A=a_0+2[0,m),\quad B=3+2[0,m),\quad E=2m-2-u.
 \tag{1}
\]

As in the [fixed-scale theorem](regularized-bezout-factor-moments.md),
choose `H` from the product corners, put `K=H+p`, and work in
`Z_2[omega]`, where `omega^2+omega+1=0`. Define

\[
 z=N+2^K\omega,\qquad
 P(X)=\prod_{a\in A}(X-a),\qquad
 Q_z(X)=\prod_{b\in B}(z-bX),
 \tag{2}
\]

and let `V_z Q_z=1 mod P`, `deg V_z<m`. The requested scalar is

\[
 v(N,a_0;p)=\widehat V(0)\pmod {2^p},\qquad
 \widehat V=\omega2^{K+2E}V_z.
 \tag{3}
\]

This uses the same global normalization as the previous coefficient
interface. The large exponent `K+2E` is encoded in binary. A provider
must evaluate the normalized quantity before reduction.

Let `F_(a_0)` be the proper divisors `a>=a_0` with `N/a>=3`,
and let `D_(a_0)` be their number. The residue theorem gives

\[
 \widehat V(X)\equiv
     P(X)\sum_{a\in F_{a_0}}\frac{c_a}{X-a}\pmod {2^p}.
 \tag{4}
\]

The padding in (1) makes the two derivative units equal to the same
odd unit `alpha_0` on every actual factor pair. Repeated squaring
also gives `a^m=1 mod2^p`. The proof in the leading-coefficient
note therefore applies at this arbitrary positive precision:

\[
 c_a\equiv\alpha_0^{-2}a\pmod {2^p}.
 \tag{5}
\]

Now zero is outside both odd root sets, and `P(0)` is odd. Evaluating
(4) at zero proves the count identity

\[
 \boxed{v(N,a_0;p)
   \equiv-P(0)\alpha_0^{-2}D_{a_0}\pmod {2^p}.}
 \tag{6}
\]

The multiplier is a unit. Thus whenever `D_(a_0)<2^p`, the
normalized value is zero exactly when the selected support is empty.
The decision procedure need not compute `P(0)` or `alpha_0`.
Their residues are nevertheless available from the existing unit
interval-product and normalized factorial algorithms.

There is a stronger uniformity fact. Modulo `2^p`, the progression
in (1) visits every odd residue `2^(u-p+1)` times. This repetition
count is even. The product of all elements in a finite abelian
group has square one, by pairing elements with their inverses.
Also `m` is even. Therefore

\[
 P(0)\equiv1\pmod {2^p},\qquad
 \boxed{v(N,a_0;p)\equiv-\alpha_0^{-2}D_{a_0}\pmod {2^p}.}
 \tag{6a}
\]

The derivative unit `alpha_0` depends only on the common interval
length, not on its first node. Thus (6a) has exactly the same odd
multiplier at every threshold for fixed `N,p`. This includes
precisions one and two. The height bound and displacement may vary
with the threshold; the proved normalization removes that variation.

### A short supplied jet also gives this value

Write `Y=X-1`, and let `v_d` be the coefficient `d` positions below
the top degree `m-1` in `widehat V`. The fixed-scale theorem gives
`v_2(v_d)>=d`, so only `min(m,p)` entries matter. Consequently

\[
 \widehat V(0)
 \equiv\sum_{d=0}^{\min(m,p)-1}(-1)^{m-1-d}v_d\pmod {2^p}.
 \tag{7}
\]

The implementation provides this linear-time jet postprocessor. It
requires the actual globally correct coefficients. Its divisibility
checks are necessary consistency conditions, not a certificate of
the supplied values' global correctness.

## 2. Roughness pays for the needed precision

Choose an integer `s>=1` and trial-divide `N` by integers below
`2^s`. An encountered divisor yields the smallest prime factor by
increasing-order search. If the search reaches the square root,
it also certifies primality when no divisor occurs.

Otherwise every prime factor of the remaining input is at least
`2^s`. Write `Omega(N)` for the number of prime factors with
multiplicity. Since `N<2^n`,

\[
 \Omega(N)<n/s,\qquad
 \tau(N)=\prod_{q^e\parallel N}(e+1)
       \le2^{\sum e}=2^{\Omega(N)}<2^{n/s}.
 \tag{8}
\]

Every suffix count is smaller than `tau(N)`. Therefore

\[
 p=\lceil n/s\rceil+1
 \tag{9}
\]

suffices in (6) for every threshold query. The nonempty suffix
predicate is monotone. Its largest true odd threshold is the largest
proper divisor, and division by that divisor gives the smallest
prime factor. At most `O(n)` scalar queries are needed.

If `C(n,p)` bounds the provider's total cost uniformly for the
contracts (1)--(3), the resulting algorithm has bound

\[
 T(n)\le2^s\operatorname{poly}(n)
       +O(n)C(n,\lceil n/s\rceil+1)
       +\operatorname{poly}(n).
 \tag{10}
\]

In particular,

\[
 C(n,p)\le\operatorname{poly}(n)2^{O(p)}
 \quad\Longrightarrow\quad
 T(n)\le\operatorname{poly}(n)2^{O(s+n/s)}.
 \tag{11}
\]

Taking `s` near the integer square root of `n` yields the conditional
`exp(O(sqrt(n)))` bound. An exponential dependence on precision is
therefore compatible with a useful provider here. A bound polynomial
in the unexpanded interval length would not meet (11).

### The primitive-scale alternative

The exact residue polynomial `W` in the
[displaced extraction theorem](displaced-bezout-extraction.md) has
coefficient content `mu`. If the support is nonempty, (6) and its
strict count bound show

\[
 \mu\le v_2(W(0))=v_2(D_{a_0})<p.
 \tag{12}
\]

The unit factor in (6) cannot alter this valuation. Since `K=H+p`,
the exact stabilized primitive exponent is `r_K=K+2E-mu`.
Hence it exceeds `2E+H-1` exactly when the support is nonempty;
the empty-support upper bound is the one already proved in the
displaced extraction note. A uniform evaluator for this scale test
with the corresponding cost in (11) would give the same transfer.
The transfer does not implement that evaluator.

## 3. Prime-only windows reduce precision to `O(log n)`

The [single-characteristic transfer](one-prime-bezout-recovery.md)
uses windows that contain only prime divisors. The same argument
applies to differences in (6a), and needs no change of coefficient
field.

After optional trial division below `2^s`, maintain a certified lower
bound `L>=2^s` for every prime factor. Every divisor in

\[
 [L,J),\qquad J=\min\{L^2,\lfloor\sqrt N\rfloor+1\},
 \tag{12a}
\]

is prime: a composite divisor has at least two prime factors, each
at least `L`, so it is at least `L^2`. The number of prime divisors
in any subinterval of this window is at most `Omega(N)<n/s`.
Take

\[
 R=\lceil n/s\rceil,\qquad
 p=\lceil\log_2(R+1)\rceil.
 \tag{12b}
\]

For odd thresholds `a<b`, subtracting (6a) gives the same unit times
the divisor count in `[a,b)`. Round each endpoint of (12a) upward
to the next odd integer; this preserves the odd divisors in the
window. The interval count is less than `2^p`, so its two scalar
residues are equal exactly when the interval is empty.

An empty window permits `L<-J`. If `L` passes the square root,
the input is prime. A nonempty window can be split in two, always
retaining its first nonempty half. Every tested subinterval still
contains only prime divisors. Binary search returns the smallest
one, and the certified lower bound excludes every earlier prime.

There are `O(log n)` successive squaring windows and `O(n)` total
scalar queries, including the final binary search. In particular,
with `s=1` no odd trial division is required and

\[
 \boxed{p=\lceil\log_2(n+1)\rceil}
 \tag{12c}
\]

suffices for deterministic factoring by supplied normalized values
over the fixed field. No factorial unit is needed by the decoder:
it compares the two residues directly.

The more precise cost transfer is

\[
 T(n)\le2^s\operatorname{poly}(n)
       +O(n)C\!\left(n,\left\lceil
          \log_2(\lceil n/s\rceil+1)\right\rceil\right)
       +\operatorname{poly}(n).
 \tag{12d}
\]

A bound `C(n,p)<=poly(n)2^O(p)` makes (12d) polynomial at `s=1`.
Even `C(n,p)<=poly(n)2^O(2^p)` gives
`T(n)<=poly(n)2^O(s+n/s)`, since `2^p<=2R`. The latter implies
the conditional square-root exponential bound by taking `s` near
`sqrt(n)`. No algorithm meeting either provider bound is claimed.

The implementation includes both the direct suffix wrapper in
section 2 and this stronger window wrapper. A complete suffix
count may be a nonzero multiple of `2^p`; the window method never
uses an isolated zero suffix value as an empty-support certificate.
The decoder requires the actual residues for subtraction. Individual
zero tests or primitive-scale indicators do not supply those differences.

## 4. An actual provider from a small symmetric difference

Here the domains may have independent power-of-two lengths:

\[
 A=a_0+2[0,r),\qquad B=b_0+2[0,s),\qquad
 r=2^u,\quad s=2^v.
 \tag{13}
\]

Use the fixed regularization at precision one, with clearing exponent
`K+E_A+E_B`. Every surviving factor weight is one modulo two.
Writing `D(A,B)` for the ordered factor-pair count, we have

\[
 [X^{r-1}]\widehat V=\widehat V(0)
       =D(A,B)\pmod2.
 \tag{14}
\]

Let `C=A intersect B`. The pair set has the disjoint partition

\[
 F(A,B)=F(C,C)\sqcup F(A\setminus B,B)
                       \sqcup F(C,B\setminus A).
 \tag{15}
\]

On `F(C,C)`, the map `(a,b) -> (b,a)` pairs every off-diagonal
solution. The only possible fixed point is `(sqrt(N),sqrt(N))`.
Thus

\[
 \begin{aligned}
 D(A,B)\equiv{}&[N\text{ is a square and }\sqrt N\in C]\\
 &+\sum_{a\in A\setminus B}[a\mid N\text{ and }N/a\in B]\\
 &+\sum_{b\in B\setminus A}[b\mid N\text{ and }N/b\in C]
 \pmod2.
 \end{aligned}
 \tag{16}
\]

Each displayed sum needs only one divisibility and interval-membership
test per boundary node. The interval intersection and its at most
four boundary segments are found from binary endpoints. If
`k=|A symmetric_difference B|`, the complete cost is
`O((k+1) poly(input length))`. No interior node is visited.

The implementation also discards boundary nodes excluded by positive
product bounds before it applies the work limit. For the first sum,
only `ceil(N/max B)<=a<=floor(N/min B)` can contribute. For the
second sum, use the corresponding bounds from `C`, and skip it if
`C` is empty. Thus the actual visit count can be smaller than `k`.

All roots reduce to one modulo two. The complete normalized pair is
therefore also determined by this parity:

\[
 \boxed{\widehat V\equiv D(A,B)(X-1)^{r-1},\qquad
        \widehat U\equiv D(A,B)(X-1)^{s-1}\pmod2.}
 \tag{17}
\]

The implementation returns these single-power representations with
binary degrees. It checks the boundary size against a supplied work
limit before it starts enumeration.

When `A=B`, formula (16) is exactly a square test. For equal-length
overlapping intervals shifted by `d` node positions, at most `2|d|`
boundary nodes are examined. These are genuine global coefficient
computations for the stated families, including singular factor
incidence. They do not compute the general suffix oracle in (3):
its two interval starts can have an exponentially large separation,
and the resulting boundary work is charged in full.

## 5. Evidence and remaining obligation

The checker compares the value postprocessor and the boundary provider
against 48 complete small rational inverse computations. It also
exercises the roughness wrapper on 393 input-and-trial-bound pairs,
using 795 explicitly known-factor value queries and 683 exact
derivative-weight controls. Those known-factor queries test the
reduction, not a general value provider.
The prime-window wrapper independently handles the same 393 controls
using 675 known-factor value queries. A further case has 32 proper
divisors and a zero initial suffix value modulo 32; the difference
of its first two prime-window values still recovers the prime three.

Three further controls use the implemented boundary provider on
implicit large domains. Two 201-bit domain counts require zero
boundary visits and distinguish a square from a nonsquare. A third
uses a 303-bit input and two 301-bit domain counts; three feasible
boundary visits, from six raw boundary nodes, return the exact
global normalized parity. These three
scalars are computed, rather than merely specified by contracts.
An excessive-boundary case is rejected before enumeration.

A separate 402-bit general contract has precision 22 and size
exponent 424. That contract example does not evaluate the general
scalar. The distinction is retained in the evidence receipt.

```sh
python3 tools/check_bezout_count_values.py --output /tmp/bezout-count-values.json
cmp /tmp/bezout-count-values.json evidence/bezout-count-values.json
```

The remaining target is a uniform evaluator for the scalar in (3),
or its short jet in (7), with a bound strong enough for (12d).
The boundary provider supplies one limited singular family. Its
boundary parameter is not bounded by `poly(n)2^O(p)` on general
factor-search queries, so the complete proved score remains `1/2`.
