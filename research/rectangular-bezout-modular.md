# Complete primitive Bezout polynomials on unequal arithmetic domains

## Constructive result

The modular Bezout constructor now allows arbitrary independent domain
sizes. It computes the complete optimal primitive pair modulo `2^p`
in polynomial bit time, without padding either domain or expanding
its polynomial degree.

Let `m,s>=1`, let

\[
 A_i=a+D_Ai\ (0\le i<m),\qquad
 B_j=b+D_Bj\ (0\le j<s),
\]

where `a,b>0` are one modulo four and the positive steps have
`v2(D_A)=v2(D_B)=2`. Let `N>0` be three modulo four. Set

\[
 P(X)=\prod_i(X-A_i),\qquad Q(X)=\prod_j(N-B_jX).
                                                               \tag{1}
\]

The [rectangular denominator theorem](rectangular-bezout-denominator.md)
gives the exact least exponent

\[
 R=\max_{0\le k<m}\left(s+k-v_2\binom{s+k-1}{k}\right),
 \qquad \Delta=m+s-1-R\le\mathrm{bits}(m+s-1).   \tag{2}
\]

Let `U_0,V_0` be the unique rational Bezout pair of degrees below
`s,m`, respectively. Given numerical precision `p>=1`, the algorithm
returns every coefficient of `U=2^R U_0`, `V=2^R V_0` modulo `2^p`
as at most `p+Delta` possible highest powers of `X-1` in each
polynomial. All omitted coefficients are zero modulo `2^p`.
Each polynomial has an odd coefficient, although its leading
coefficient need not be odd.

The cost is polynomial in the complete binary input length and `p`.
The domain sizes and normalization exponent remain binary integers.
The construction extends the [equal-size compiler](balanced-bezout-modular.md)
using two independent low-surplus binomial expansions. Its unit-product
routine also permits odd higher jet coefficients, a case not covered
by the older coefficientwise-even step rule.

The separation modulo four already proves coprimality. This computes
global normalized coefficient data for that family, not a general
common-root discriminator. The complete factoring score remains `1/2`.

Files: [constructor](../tools/rectangular_bezout_modular.py),
[checker](../tools/check_rectangular_bezout_modular.py),
[binary scale helper](../tools/rectangular_bezout_scale.py), and
[evidence](../evidence/rectangular-bezout-modular.json).

## 1. Normalized differences from low symmetric degrees

The two child-coordinate prefixes are

\[
 x_i=(A_i-1)/4,\qquad y_j=(N/B_j-3)/4.           \tag{3}
\]

Both are 2-adic isometries of their consecutive index sets. For the
first, the step `D_A/4` is odd. For the second, inversion and
multiplication by odd units preserve each difference valuation.

For any such prefix `z_0,...,z_a`, define

\[
 C_{a,a+i}=a![z_0,\ldots,z_a]{Z\choose a+i}.
                                                               \tag{4}
\]

The companion theorem proves `C_(a,a+i)` is integral and `C_(a,a)=1`.
The coefficient identity from the earlier balanced proof gives

\[
 C_{a,a+i}=
 \frac{[T^i]\left(\prod_{j=0}^{a+i-1}(1-jT)\right)
                    \left(\prod_{j=0}^{a}(1-z_jT)\right)^{-1}}
      {(a+1)(a+2)\cdots(a+i)}.                  \tag{5}
\]

To compute all `0<=i<p`, retain additional precision

\[
 g=v_2((a+1)\cdots(a+p-1))
       =O(p\log(a+p)).                          \tag{6}
\]

The numerator in (5) is computed modulo `2^(p+g)`. Divide by the
known power-of-two part of the short denominator, then invert its
odd part modulo `2^p`. Integrality from (4) justifies that division;
the guard (6) preserves every requested bit.

The first product in (5) uses small-index elementary functions of
`0,...,a+i-1`, computed by exact power sums and Newton identities.
The second product uses the unit progression construction below.
Only degrees below `p` and polynomial guard precision are retained,
even when `a` is exponential.

## 2. Unit progressions with arbitrary higher jet coefficients

Work at precision `P` and truncated degree `d`. Let `u` be a unit
jet, let `W=v/u`, and suppose only that the constant coefficient
of `W` is even. Its other coefficients may be odd. The finite identity

\[
 \prod_{i=0}^{L-1}(u+iv)
 =u^L\sum_{j=0}^{L}e_j(0,\ldots,L-1)W^j        \tag{7}
\]

still has a short exact truncation. A retained coefficient of degree
at most `d` in `W^j` uses at most `d` nonconstant factors. If the
constant coefficient has valuation `c>=1`, that coefficient has
valuation at least `c(j-d)`. Therefore only

\[
 j\le d+\lfloor(P-1)/c\rfloor                  \tag{8}
\]

can survive. If the constant coefficient is zero, `W^(d+1)=0`
already in the truncated polynomial ring.

The number of elementary symmetric functions in (7) is polynomial
in `P,d`. Their integer bit lengths are polynomial in those values
and `log L`. Modular binary powering computes `u^L`. Thus this
extension has a fully charged polynomial bit bound.

For the `x` prefix in (3), a factor of its reverse node product is

\[
 1-\left((a-1)/4+(D_A/4)i\right)T.
\]

Its progression step has zero constant coefficient, so the exact
degree cutoff applies. For the reciprocal prefix,

\[
 1-y_jT=\frac{B_j+(3B_j-N)T/4}{B_j}.            \tag{9}
\]

The numerator progression step is `D_B+(3D_B/4)T`: its constant
is even and its linear coefficient is odd. Formula (8) is the
needed extension. The odd denominator progression is evaluated at
degree zero. These formulas compute every node product in (5)
without listing its nodes.

## 3. Two surplus indices give each exact Newton unit

Fix a Newton index `k`. Put `alpha=k`, `beta=s-1`, and let `D_(alpha,beta)`
be the parity divided difference on the normalized row and reciprocal
column prefixes from the companion theorem. Define

\[
 \Omega_{\alpha,\beta}=\prod_{v=-\alpha}^{\beta}(1+2v),\qquad
 B_{\alpha,\beta}=
 \frac{(-1)^\beta\binom{\alpha+\beta}{\alpha}}
      {\Omega_{\alpha,\beta}}.
\]

Its [double-Mahler formula](rectangular-bezout-denominator.md) gives
the odd normalized unit

\[
 J_{\alpha,\beta}:=D_{\alpha,\beta}/B_{\alpha,\beta}
 \equiv\sum_{i+j<p}
 \frac{(-1)^j2^{i+j}(\alpha+\beta+1)_{i+j}
       C_{\alpha,\alpha+i}(x)C_{\beta,\beta+j}(y)}
 {\prod_{u=1}^i(1-2\alpha-2u)
  \prod_{u=1}^j(1+2\beta+2u)}\pmod{2^p}.        \tag{10}
\]

There are `p(p+1)/2` terms. Every denominator displayed in (10)
is odd, and each rising product has fewer than `p` factors.
All omitted terms have valuation at least `p`; that conclusion
comes from the double-Mahler theorem, not an empirical truncation.

The single odd product `Omega` is a long arithmetic progression
product and is computed by (7) at degree zero. The odd part of the
binomial is computed from three normalized factorial units, while
its valuation comes directly from binary digit sums. No factorial
or binomial integer of exponential bit length is formed.

Let `gamma=(-1)^s prod_j B_j`, the odd leading coefficient of `Q`,
and write `binom(alpha+beta,alpha)=2^kappa*b_odd`. The scaled Newton
coefficient of `V` is

\[
 c_k=-(-1)^\beta\,
 2^{R-s-k+\kappa}\,
 b_{\rm odd}\Omega_{\alpha,\beta}^{-1}
 J_{\alpha,\beta}\gamma^{-1}.                  \tag{11}
\]

This records both the affine divided-difference scale and the
leading coefficient of the reciprocal root polynomial. The signs
are checked against complete rational inversion, not only against
valuation formulas.

The column surplus list in (10) is common to every row prefix and
is computed once. For descending `k`, the odd binomial unit updates
by multiplication by `odd(k)/odd(k+s-1)`. Also
`Omega_(k-1,beta)^(-1)=Omega_(k,beta)^(-1)(1-2k)`.
These updates require only binary integers and modular unit operations.

## 4. Optimal scaling and nonmonotone coefficient windows

The [binary scale helper](../tools/rectangular_bezout_scale.py)
computes (2) by inspecting a final binary-length window. Its
binomial valuation routine is the exact identity

\[
 \kappa_k=\operatorname{pop}(k)+\operatorname{pop}(s-1)
              -\operatorname{pop}(s+k-1).        \tag{12}
\]

Binary carries bound `kappa_k` by `bits(m+s-1)`. The last index
therefore already attains a value within that many units of
`m+s-1`, and earlier indices outside the final window cannot improve
the maximum. This charges the entire scale computation.

The exact normalized gap at index `k` is

\[
 \delta_k=R-s-k+\kappa_k.                       \tag{13}
\]

Every index below `max(0,R-s-p+1)` has gap at least `p` and can be
omitted. The complete remaining interval contains at most `p+Delta`
indices. Within this interval the gaps need not be monotone, so the
implementation inspects every index before omitting a coefficient.
For example, with `m=5,s=2,p=3`, indices four and two survive while
the intervening index three vanishes. This case is retained in the
checker.

For a surviving coefficient, (10) needs only precision `p-delta_k`.
The separate short-factorial guard in (6) is still included. The
large binary exponent `R` is used to compute the small gap; only
gaps below `p` become literal powers of two.

## 5. The output contains all coefficients

Write `Y=X-1`. The Newton form of the complete second polynomial is

\[
 V(X)=\sum_{k=0}^{m-1}c_k\prod_{i<k}(X-A_i).    \tag{14}
\]

A term `ell=m-1-k` places below the largest Newton degree has
coefficient valuation at least `ell-Delta`. Every root shift
`A_i-1` has valuation at least two. Its contribution to the
monomial coefficient `d` places below the top degree thus has
valuation at least

\[
 \ell-\Delta+2(d-\ell)\ge d-\Delta.
\]

The symmetric rectangular theorem gives the same coefficient bound
for `U`: the reciprocal root shifts `N/B_j-1` have valuation one,
which still gives `ell-Delta+(d-ell)=d-Delta`. Therefore every
coefficient beyond reverse degree `p+Delta-1` is zero modulo `2^p`
in both complete polynomials.

The implementation accumulates (14) through that depth using
short reverse products of its affine prefixes. Let stars denote
reversal at the degrees `m,s,s-1,m-1` of `P,Q,U,V`, respectively.
Then

\[
 U^*=-V^*Q^*/P^*\pmod{T^s}.                    \tag{15}
\]

The reversed scalar correction `2^R T^(m+s-1)` lies outside the
degree range in (15), including when one domain has size one.
The divisor `P^*` has constant coefficient one, so this is a unit
series division with no precision loss.

Every degree-`j` coefficient of `Q^*/P^*` has valuation at least
`j`. Thus this ratio itself needs only degrees below `p`, even
when the output defect is much larger than `p`. Combining that
short ratio with `V^*` gives (15) through the complete output
window. The independent tail theorem above, rather than a formal
division heuristic, certifies all omitted coefficients.

There need not be a unit in the highest degree. For `m=4,s=2`,
the scale maximizer precedes the final index; modulo two the
highest `V` coefficient vanishes and the next one is odd. The
constructor preserves this case and verifies that each complete
output still has at least one odd coefficient.

## 6. Bit cost and exact controls

Let `L` be the binary length of the complete input, including both
domain sizes. The scale search and coefficient window have
`O(L+p)` indices. Each retained coefficient uses `O(p^2)`
low-surplus terms and unit progression products of polynomial
degree and precision. The largest short-factorial guard is
`O(p(L+log p))`. Exact small-index power sums, Newton identities,
binary modular powering, and modular inversion all have polynomial
bit cost in these quantities. Thus the complete algorithm has a
polynomial bit bound in `L+p`.

No loop ranges over either domain. No intermediate precision is
`R+p`, and no literal integer `2^R` is formed. The output exponent
and all polynomial exponents are binary encoded. The construction
charges the potentially large factorial-ratio guard separately
from the output precision; it does not assume these are equal.

The exact checker verifies:

- 64 complete rational rectangles and all 612 coefficients of
  their normalized pairs;
- 194 rational low-surplus divided differences;
- 96 literal progression products with arbitrary higher jet
  coefficients;
- nonmonotone Newton gaps and primitive pairs whose leading
  coefficient is even;
- five invalid input contracts.

Three additional cases have independent domain-size bit lengths
`(201,174)`, `(131,201)`, and `(1,257)`. Each is evaluated at two
precisions. The checker compares reductions, verifies the complete
short-support polynomial identity, and checks zero enumerated nodes.
The middle case has defect 200 and needs 205 working bits for four
output bits. Its complete outputs have only three nonzero terms
each; the large guard and binary support offset remain explicit.

Run the constructor controls with:

```sh
python3 tools/check_rectangular_bezout_modular.py \
  --output evidence/rectangular-bezout-modular.json
```

The evidence contains the exact parameters, coefficient outputs,
control counts, and hashes of the implementation and its dependencies.
This is a constructive extension of the companion rectangular
denominator theorem. It supplies global primitive coefficient
data for a residue-separated family; it does not change the
unconditional factoring exponent or score.
