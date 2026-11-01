# A quadratic trace and a two-depth primitive Bezout testbed

## Result and scope

The [deep matched construction](matched-progression-bezout.md) has a
useful nonuniform-depth testbed. The matching denominators now have
two different valuations, and their cancellation is a weighted
quadratic reciprocal sum. That sum still has a polynomial binary-size
unit evaluator.

Let `m=2^t`, `t>=1`, `D=2^h`, and let `a>D(m-1)` be odd. Put

\[
 A_i=a+Di,\quad B_i=a-Di\quad(0\le i<m),\qquad
 N=a^2+D^2w>0,\qquad w\equiv1\pmod4.             \tag{1}
\]

For `h>=2t-1`, the exact least primitive Bezout exponent is

\[
 R=2hm+2m-3t.                                   \tag{2}
\]

Define `n=m-1`, `F=v2(n!)=m-1-t`, `n!_odd=2^(-F)n!`, and

\[
 T_t(w)=2^{2-t}\sum_{j=0}^{n}
                \frac{\binom nj^2}{j^2+w}.       \tag{3}
\]

The number `T_t(w)` is an odd 2-adic integer. Its residue modulo
`2^p` is computable in polynomial bit time using at most `4p+8`
initial Mahler coefficients, independent of `m`.

Let `P=prod_i(X-A_i)`, `Q=prod_j(N-B_jX)`, and let `U,V` be
their unique degree-bounded primitive pair, with
`UP+VQ=2^R`. In the explicit precision range

\[
 1\le p\le h-2t+2,                              \tag{4}
\]

their **complete** residues are

\[
 \begin{aligned}
 U(X)&\equiv a\,T_t(w)(n!_{\rm odd})^{-2}(X-a)^n,\\
 V(X)&\equiv-a^{1-m}T_t(w)(n!_{\rm odd})^{-2}(X-a)^n
                    \pmod{2^p}.                 \tag{5}
 \end{aligned}
\]

This is a constructive testbed, not an unknown-matching factoring
algorithm. The restriction `w=1 mod4` already excludes a common root.
For unrestricted `w`, incidence in these reflected domains is exactly
the elementary condition `w=-i^2` with `0<=i<m`. The checker includes
both incident and nonincident examples in that wider family. The
complete factoring score remains `1/2`.

Files: [constructor](../tools/quadratic_coalescence_bezout.py),
[checker](../tools/check_quadratic_coalescence_bezout.py), and
[evidence](../evidence/quadratic-coalescence-bezout.json).

## 1. A normalized weighted-moment identity

For `0<=k<=n`, the Vandermonde identity gives

\[
 M_k:=\sum_{j=0}^{n}\binom nj^2\binom jk
       =\binom nk\binom{2n-k}{n-k}.              \tag{6}
\]

For example, extract the coefficient of degree `n-k` from
`(1+z)^(n-k)(1+z)^n` after taking out `binom(n,k)`.
Since `n=2^t-1`, binary digit sums give the exact order

\[
 v_2(M_k)=t-v_2(k+1).                            \tag{7}
\]

Indeed `binom(n,k)` is odd and
`pop(n-k)=t-pop(k)`,
`pop(2n-k)=t+1-pop(k+1)`. Substitute these in the
factorial valuation formula. The consecutive odd units update through

\[
 \frac{M_{k+1}}{M_k}
       =\frac{(n-k)^2}{(k+1)(2n-k)}.             \tag{8}
\]

Thus all needed units follow from the initial central-binomial unit
and short modular updates. The initial unit uses three normalized
factorial queries, not the exponentially long factorial integer.

## 2. The quadratic Mahler coefficients have a linear tail bound

Set `f(z)=2/(z^2+w)` and `b_k=Delta^k f(0)`. Every `f(j)` is
2-adically integral: its denominator has valuation zero for even `j`
and one for odd `j`. For all `k`,

\[
 v_2(b_k)\ge v_2(k!)-\left\lfloor\frac{k+1}{4}\right\rfloor.
                                                               \tag{9}
\]

To prove (9), take a square root `alpha` of `-w` in a quadratic
extension. It is a unit. For even integers `j`, `j+alpha` and
`j-alpha` are units; for odd `j`, both have valuation `1/2`.
The latter follows from `v2(j^2+w)=1` and the difference `2alpha`.
Partial fractions and the exact forward-difference formula give

\[
 b_k=\frac{(-1)^k k!}{\alpha}
 \left(\frac1{(-\alpha)_{k+1}}-
       \frac1{(\alpha)_{k+1}}\right).            \tag{10}
\]

Each product has exactly `floor((k+1)/2)` odd-index factors.
Taking the valuation lower bound in (10), and rounding up because
`b_k` is rational, proves (9).

The first coefficients have two useful extra cancellations:

\[
 v_2(b_0)=1,\quad v_2(b_1)=0,\quad
 v_2(b_3)\ge3,\quad v_2(b_7)\ge5.               \tag{11}
\]

For `b_3`, the exact numerator over
`prod_(j=0)^3(w+j^2)` is `72(w-1)`.
For `w=1+4z`, the numerator of `b_7` over
`prod_(j=0)^7(w+j^2)` is

\[
 512(121275+732060z-590940z^2+35280z^3).
\]

These denominators have orders two and four, respectively.
Together, (9)--(11) imply

\[
 v_2(b_k)\ge v_2(k+1)\qquad(k\ne1).            \tag{12}
\]

For `k=0,2,4,5,6` this is immediate; indices three and seven use
(11). For `k>=8`, use
`v2(k!)-floor((k+1)/4)>=floor(k/2)-1` and the elementary
comparison with `v2(k+1)`. Checking eight through eleven and then
using `(k-3)/2>=log2(k+1)` for `k>=12` suffices.

Newton interpolation on `0,...,n` is a finite exact identity, so

\[
 T_t(w)=2^{1-t}\sum_{k=0}^{n}b_k M_k.           \tag{13}
\]

The `k=1` term is odd. Every other term is even by (7) and (12).
This proves both the integrality and oddness asserted in (3).

Moreover the order of the `k`th normalized term is at least

\[
 \left\lfloor\frac{k}{2}\right\rfloor-v_2(k+1).
                                                               \tag{14}
\]

For all `k>=4p+8`, this is at least `p`. One may use
`v2(k+1)<=log2(k+1)` and `2^(p+3)>=4p+9` to verify the endpoint;
the continuous lower bound `(k-1)/2-log2(k+1)` then increases.
Consequently only the first `4p+8` coefficients in (13) are needed.
They come from the same number of literal values of `f`, followed
by an ordinary difference table. This lists `O(p)` interpolation
samples, not the `m` nodes of the original domains.

The normalized moment shift is `1-v2(k+1)`. Guard precision
`max(0,max_k(v2(k+1)-1))=O(log p)` makes every division by its
power-of-two denominator exact at the requested output precision.
The code includes this guard before forming the difference table.

## 3. Nonuniform paired depths and the exact primitive exponent

In the coordinate `Z=(X-a)/D`, the reciprocal roots are

\[
 g(j)=\frac{j+\lambda w}{1-\lambda j},\qquad
 \lambda=D/a,\qquad
 g(j)-j=\frac{\lambda(j^2+w)}{1-\lambda j}.       \tag{15}
\]

Their paired depths are `h` for even `j` and `h+1` for odd `j`.
For distinct indices, all differences have the valuation of the
index difference, which is at most `t-1`. Indeed

\[
 \begin{aligned}
 g(j)-i&=\frac{j-i+\lambda(w+ij)}{1-\lambda j},\\
 g(j)-g(l)&=\frac{(j-l)(1+\lambda^2w)}
                   {(1-\lambda j)(1-\lambda l)}.
 \end{aligned}                                                   \tag{16}
\]

Let `d_r(lambda)` be the divided difference of the indicator that
is zero on `0,...,m-1` and one on `g(0),...,g(r-1)`. In its
Lagrange sum, the `j`th term is

\[
 \frac1{\prod_{i<m}(g(j)-i)
             \prod_{l<r,l\ne j}(g(j)-g(l))}.     \tag{17}
\]

For `r=m`, replacing each nonmatching factor by its constant term
and retaining the matching first-order factor gives

\[
 d_m(\lambda)=\frac1{\lambda(n!)^2}
       \sum_{j=0}^{n}\frac{\binom nj^2}{j^2+w}
       +\mathcal E.                             \tag{18}
\]

This is a convergent 2-adic expansion at the actual `lambda`.
Every nonmatching relative correction in (16) has order at least
`h-t+1`; all other unit factors have at least that order as well.
Because `binom(n,j)` is odd, each leading summand in (18) has
order at least `-h-2F-1`. Therefore

\[
 v_2(\mathcal E)\ge-2F-t.                       \tag{19}
\]

By Section 2, the displayed sum in (18) has order `t-2`.
The relative error is consequently at least `h-2t+2>=1`.
This proves the exact leading divided-difference order

\[
 v_2(d_m)=-h-2F+t-2.                            \tag{20}
\]

For the earlier coefficients, (17) directly gives

\[
 v_2(d_r)\ge-h-1-F-v_2((r-1)!).                 \tag{21}
\]

There is no cancellation assumption in this bound. Its two internal
factorial denominators have orders `F` and at most `v2((r-1)!)`;
the matching denominator contributes at most `h+1`.

Write the interpolating indicator as `P U_0`, where `U_0` is the
rational Bezout polynomial before normalization. Its Newton form is

\[
 U_0(X)=\sum_{r=1}^{m}D^{-(m+r-1)}d_r(\lambda)
                    \prod_{j<r-1}(X-N/B_j).      \tag{22}
\]

Equations (20) and (22) put the leading coefficient at valuation
`-R`, with
`R=2hm+2F-t+2=2hm+2m-3t`.
For `s=m-r>=1`, the normalized earlier coefficient has order at least

\[
 hs+F-v_2((r-1)!)-t+1\ge h-t+1.                \tag{23}
\]

All roots in (22) are integral, so these inequalities prove that
`2^R U_0` is integral and has an odd leading coefficient. Since
`Q` has an odd leading coefficient, exact polynomial division of
`2^R-P(2^R U_0)` by `Q` gives an integral `V`. Its leading
coefficient is also odd. This proves optimality of (2).

## 4. Complete coefficients in the precision window

Under (4), the relative error in (18) vanishes at output precision.
The normalized leading coefficient from (18) is exactly

\[
 a T_t(w)(n!_{\rm odd})^{-2}\pmod{2^p}.          \tag{24}
\]

Every earlier Newton term vanishes by (23). Every reciprocal root
shift `N/B_j-a=Dg(j)` is divisible by `2^h`; expanding the remaining
top Newton product in `X-a` therefore leaves only its highest power
modulo `2^p`. This proves the complete first formula in (5).

The reversed input polynomials also have only constant coefficients
at that precision. Their constants are one and
`prod_j B_j=a^m mod2^p`, respectively; `m` is even. The reversed
Bezout identity gives the second formula in (5). The reversed
scalar correction lies beyond degree `m-1`. This establishes the
entire output, including every omitted coefficient.

The parameters and `R` remain binary integers. Formula (13), the
short odd-unit updates, and normalized factorial queries have bit
cost polynomial in `t,p,log(1+|w|)`. Formula (5) adds binary modular
powering in `a,m`. No literal integer `2^R` is constructed.

## 5. Exact controls and the unknown-matching boundary

The checker compares (3) with literal rational sums, checks (6)
and (7), and compares every coefficient in (5) with complete
rational Bezout inversion. It includes the endpoint `h=2t-1` and
the largest allowed precision. Huge inputs compare independent
precision runs and retain zero enumerated original domain nodes.

For arbitrary `w` and `h>=t`, a common root in the reflected domains
would require

\[
 a(i-j)=D(w+ij).
\]

Since `a` is odd and `|i-j|<m<=D`, reduction modulo `D` forces
`i=j`, after which `w=-i^2`. The broader incidence controls use
this exact equivalence and compare both signs of the outcome with
literal products. They do not infer incidence from a modular
Bezout identity whose scalar right-hand side has vanished modulo
the working precision.

The new contribution is the quadratic weighted trace and its use
in an exact primitive normalization with nonuniform paired depths.
An arbitrary ordinary interval has an unknown reciprocal matching
permutation and does not reduce to the same square test. Assembling
its normalized coefficients still requires a separate interface.

Run:

```sh
python3 tools/check_quadratic_coalescence_bezout.py \
  --output evidence/quadratic-coalescence-bezout.json
```

The proof, implementation, and exact controls used substantial
OpenAI assistant support. No new Lean formalization is claimed.
