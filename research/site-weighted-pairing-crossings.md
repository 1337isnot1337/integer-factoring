# Bivariate site-weighted crossing identities

## Results and scope

This note studies a constant-auxiliary-dimension construction that changes
the local polynomial at every position. It gives an exact crossing identity
with polynomially many coefficients, an exact rank computation, a general
rank bound for coefficient refinement, and a complete asymptotic exclusion
for centrally symmetric weight multisets. Asymmetric site weights remain
a different search family. Exact asymmetric profiles through sixteen
positions are reported below; none certifies an improved factoring score.

Work over `F_2`, or any extension field. Let `S=E+C`, so `S+C=E`.
Choose nonnegative integer weights `a_1,...,a_t`, and put `W=sum_i a_i`.
Define

\[
 F_t(z,w)=\bigotimes_{i=1}^t(S_i+wz^{a_i}C_i)
        =\sum_{k=0}^t\sum_{u=0}^W w^kz^u A_{k,u}.        \tag{1}
\]

Then

\[
 F_t(1,1)=E_t.                                          \tag{2}
\]

There are at most `(t+1)(W+1)` coefficient pieces. In particular, if
`W=t^{O(1)}` and the weights are uniformly computable in polynomial time,
this is a uniform polynomial-piece family. Each tensor entry is computed
by a bivariate polynomial recurrence of polynomial size. Full coefficient
construction and ordinary exact rank factorization cost `2^{O(t)}` bit
operations. The construction therefore fits the uniform transfer interface
in [vanishing-slack-score-transfer.md](vanishing-slack-score-transfer.md).
It still needs a suitable simultaneous rank bound to improve the score.

The auxiliary dimension in (1) is one. Order dependence comes from the
different site weights. Changing the weights alone is not the claimed
escape from other fixed local profiles: the local construction itself is
the two-pair tensor `S+wz^{a_i}C`.

## 1. Exact cut ranks from homogeneous square-zero polynomials

Use the rank-one square-zero local map `N` and the representation

\[
 V=J_2\oplus2J_1
\]

from [two-state-run-rank-compression.md](two-state-run-rank-compression.md).
Define

\[
 f_{k,u}(\epsilon)=
 \sum_{\substack{A\subseteq[t]\\|A|=k,\ \sum_{i\in A}a_i=u}}
       \prod_{i\in A}\epsilon_i,
 \qquad
 R_{k,u}=\operatorname{rank}_{V^{\otimes t}} f_{k,u}(N).
                                                               \tag{3}
\]

The column cut is normalized by an invertible local matrix so that `S`
acts as the identity and `C` as `N`. In the row cut, `C` is already the
identity and `S` is `N`. Complementing the selected sites therefore gives
the exact identities

\[
 \boxed{r(A_{k,u})=R_{t-k,W-u},\qquad c(A_{k,u})=R_{k,u}.} \tag{4}
\]

Every polynomial in (3) is homogeneous of degree `k`. Its rank is computed
using

\[
 R_{k,u}=\sum_{J\subseteq[t]}2^{t-|J|}
 \operatorname{rank}_{K[\epsilon_J]/(\epsilon_i^2)}
          m_{f_{k,u}(\epsilon_J,0)}.                    \tag{5}
\]

The native engine enumerates the active subsets `J`, groups their monomials
by `(k,u)`, and uses exact homogeneous multiplication ranks. It applies the
unit, linear, common-factor, quadratic-pivot, and graded-block reductions
proved in the preceding note. No approximate matrix rank is used.

The different tensors `A_{k,u}` with a nonempty subset class are linearly
independent: `S,C` are locally linearly independent, so their tensor words
are linearly independent, and the different subset classes are disjoint.
In particular, projection to one Hamming count is an intrinsic linear map
on the coefficient span, including arbitrary coefficient recombinations.

## 2. A general refinement inequality

Let

\[
 Q_{t,k}=\sum_{|A|=k}\prod_{i\in A}N_i,
 \qquad U_{t,k}=\operatorname{rank}Q_{t,k}.               \tag{6}
\]

Choose a local Jordan basis in which `N` sends one basis vector of degree
zero to a vector of degree one and annihilates the other two trivial
vectors. Give the degree-one vector at position `i` weight `a_i`.
The total vector space is graded by weights `0,...,W`.

If `D(z)` scales a basis vector of total weight `v` by `z^v`, then over
`K(z)`

\[
 Q_{t,k}(z)=D(z)Q_{t,k}D(z)^{-1}
          =\sum_u z^u f_{k,u}(N).                       \tag{7}
\]

The coefficient at `u` is the direct sum of the blocks of `Q_{t,k}` from
weight `v` to weight `v+u`. Each block is a submatrix of `Q_{t,k}`, and
there are at most `W+1` blocks. Hence

\[
 R_{k,u}\le(W+1)U_{t,k}.                                \tag{8}
\]

Also `Q_{t,k}=sum_u f_{k,u}(N)`, so rank subadditivity gives

\[
 U_{t,k}\le\sum_uR_{k,u}.                               \tag{9}
\]

These statements hold for arbitrary nonnegative weights. They imply

\[
 \max_{k,u}r(A_{k,u})c(A_{k,u})
 \le(W+1)^2\max_kU_{t,k}U_{t,t-k}.                      \tag{10}
\]

The right side is an upper bound, not a matching lower bound for
asymmetric weights. A better simultaneous contraction would require an
improved tradeoff between the complementary levels in (4). Knowing that
some weighted level has large rank in each cut separately does not prove
that the same piece has both large ranks.

## 3. Centrally symmetric weights: a full-span exclusion

Assume the weight multiset is centrally symmetric: there is a permutation
`pi` and an integer `c` such that

\[
 a_i+a_{\pi(i)}=c\quad\text{for every }i.                \tag{11}
\]

Arithmetic progressions are one example. Repeated weights are allowed.

### 3.1 The exact finite rank-sum inequality

First suppose `t` is even and put `k=t/2`. Then `W=kc`.
For any linear combination of central coefficients, the polynomial in
the column cut has the form

\[
 f=\sum_u\alpha_u f_{k,u}.
\]

Complementation replaces the coefficient of a subset `A` by the
coefficient of its complement. But (11) gives
`weight(pi(A))=kc-weight(A)=W-weight(A)` for every `k`-subset.
Thus complementation on this central span is exactly the permutation
`pi` of algebra variables. By (4), the row and column ranks of every
central-span tensor are equal, including arbitrary recombinations.

Let `Pi` project a general tensor in the full coefficient span to Hamming
count `t/2`. In each cut, this is homogeneous-degree extraction for the
local nilpotents. The same argument as (7)--(8), now with every site
weight equal to one, shows

\[
 r(\Pi T)\le(t+1)r(T),\qquad
 c(\Pi T)\le(t+1)c(T).                                 \tag{12}
\]

The two left-hand ranks agree. For any decomposition `E_t=sum_j T_j`
with every `T_j` in the full bivariate coefficient span, the row matrix of
`Pi(E_t)` is `Q_{t,t/2}`. Therefore

\[
 \begin{aligned}
 U_{t,t/2}
 &\le\sum_j r(\Pi T_j)\\
 &\le(t+1)\sum_j\sqrt{r(T_j)c(T_j)}.
 \end{aligned}
\]

We have proved the exact finite inequality

\[
 \boxed{\sum_j\sqrt{r(T_j)c(T_j)}
       \ge\frac{U_{t,t/2}}{t+1}.}                       \tag{13}
\]

There is no restriction on the number of summands, and coefficients may
lie in any characteristic-two extension field.

If `t` is odd, central symmetry supplies a site of weight `c/2`. Fix its
four local labels to `(0,1,1,0)`, where `S=1,C=0,E=1`. This restriction
does not increase either cut rank. It deletes that site and leaves a
centrally symmetric even-length family. Consequently (13) holds with
`t'=2 floor(t/2)` in place of `t` on its right side. This also prevents
odd lengths from avoiding the leading obstruction.

### 3.2 The correct ordinary two-pair reference exponent

For clarity, the ordinary two-pair coefficient family has leading
contraction approximately **0.81127812**, rather than `3/4`.
The `3/4` theorem for a larger position-invariant tensor space is a
different lower bound; it is not a matching upper bound for this pencil.

Here is the standard representation-theory calculation needed to evaluate
(13). It uses the genuine characteristic-two tilting direct sums summarized
and sourced in the
[pinned SL2 coefficient-rank note](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/sl2-coefficient-ranks.md).
The relevant primary inputs are the tilting property of tensor powers in
[Doty, Lemma 13](https://arxiv.org/pdf/math/0610591v3), and the binary tensor
description of `T(n)` in
[Sutton--Tubbenhauer--Wedrich--Zhu, Proposition 3.4(d) and Corollary 3.6](https://arxiv.org/pdf/2105.07724v2).
No additive-rank inference from a Weyl filtration is used.

The module `2·1 plus V` has character `2+x+x^{-1}`. Write its `t`th tensor
power as a direct sum of `mu_n` copies of `T(n)`, where `0<=n<=t`.
The binary tensor description gives `dim T(n)<=(t+1)^2`, and its unipotent
coefficient of degree `k` is nonzero exactly when `k<=n`: the twist weights
have subset sums covering every integer from zero to `n`.

The integer weight multiplicities are

\[
 c_k=\operatorname{coeff}_{x^k}(2+x+x^{-1})^t
     =\operatorname{coeff}_{x^{t+k}}(1+x)^{2t}
     =\binom{2t}{t+k}.                                  \tag{14}
\]

The character and block dimensions give

\[
 \frac{c_k}{(t+1)^2}\le U_{t,k}
 \le(t+1)^2\sum_{n\ge k}c_n
 \le(t+1)^3c_k.                                        \tag{15}
\]

For the lower inequality, the weight-`k` multiplicity is at most
`(t+1)^2 sum_{n>=k}mu_n`, and each such block contributes at least one to
the coefficient rank. For the upper inequality, use `mu_n<=c_n` and the
decrease of `c_n` for `n>=0`.

Binomial log-concavity places the largest leading product
`c_k c_{t-k}` at the middle. The factorial estimate for the middle value
then gives

\[
 \log_2 U_{t,t/2}
 =t\log_2\frac{16}{3\sqrt3}+O(\log(t+2))
 =2tH_2(3/4)+O(\log(t+2)).                              \tag{16}
\]

Thus, for polynomially bounded centrally symmetric weights, (10),
(13), and (15) show that the optimal rank-sum exponent in the entire
coefficient span is

\[
 \frac{\log_2 S_{\rm opt}}{2t}
 =H_2(3/4)+O\!\left(\frac{\log(t+2)}t\right)
 =0.811278124459\ldots+O\!\left(\frac{\log(t+2)}t\right). \tag{17}
\]

The same leading exponent holds for the optimized maximum rank product,
even when arbitrarily many summands are permitted. To see the lower bound,
select a basis from any proposed summand list and re-express `E_t` using
that basis. At most `(t+1)(W+1)` original summands are needed; their ranks
are unchanged by nonzero scalar multiplication. Combine this polynomial
bound on their number with (13). The original coefficient pieces give the
matching upper bound by (10) and (15).

This theorem concerns centrally symmetric weight multisets. It does not
establish (13) or (17) for asymmetric weights, and its asymptotic statement
is not a claim that every finite width is separately excluded by the
numerical lower bound (13).

Related work in [PR 27](https://github.com/1337isnot1337/integer-factoring/pull/27),
at [commit 4bb2b3b](https://github.com/1337isnot1337/integer-factoring/commit/4bb2b3b17548d2eb26b1ad8c81b1a30d842fd83f),
gives an affine-selector frame obstruction with the same entropy constant.
The argument here uses central symmetry of site weights and homogeneous
degree projection. No inclusion between the two model classes is asserted.

## 4. Site-weight finite profiles

The squared weights `a_i=i^2` have `W=O(t^3)` and at most `O(t^4)`
coefficient pieces. The integer-computable Beatty weights
`a_i=i+floor(sqrt(2i^2))` have `W=O(t^2)` and at most `O(t^3)` pieces.
These are explicit uniform families. The engine also tests a finite
irregular weight list without assigning it an asymptotic continuation.

As before, define `P=max r c`, `S=sum sqrt(r c)`,
`rho_max=log2(P)/(4t)`, and `rho_S=log2(S)/(2t)`.

| Weights | `t` | Nonzero pieces | `rho_max` | `rho_S` |
|---|---:|---:|---:|---:|
| Squares | 10 | 720 | 0.60961757 | 0.99837673 |
| Squares | 12 | 1,859 | 0.62422080 | 0.99670258 |
| Squares | 14 | 4,032 | 0.64853139 | 0.99470476 |
| Squares | 16 | 7,667 | 0.66669946 | 0.99281627 |
| Beatty | 12 | 553 | 0.68876820 | 0.97995943 |
| Beatty | 16 | 1,435 | 0.71636756 | 0.96617760 |
| Irregular list | 12 | 558 | 0.67293828 | 0.98447194 |

The irregular list is `0,1,2,4,5,7,8,12,14,18,21,23`.
Some individual Beatty prefixes are centrally symmetric: at `t=16`,
every pair of weights in opposite order sums to 40. The theorem in
Section 3 applies to that profile. The displayed square-weight profiles
and the irregular list are not centrally symmetric.
For square weights, the normalized maximum-product growth increments
from lengths 10 to 12, 12 to 14, and 14 to 16 are respectively about
`0.697237`, `0.794395`, and `0.793876`. Their difference from the
cumulative values shows why the early low values do not certify an
asymptotic contraction. At length sixteen, the rank-sum exponent is still
about `0.992816`, well above the strict finite-seed target
`1/sqrt(2)`.

## 5. Verification and reproduction

[`site_weight_crossings.cpp`](../tools/site_weight_crossings.cpp) implements
(3)--(5) with exact binary ranks and 128-bit rank products.
[`check_site_weight_crossings.py`](../tools/check_site_weight_crossings.py)
checks 24 full-tensor examples independently, matches the eight-position
uniform profile with the earlier kernel engine, and checks arbitrary
central projections of 24 symmetric-weight coefficient recombinations.
It includes an asymmetric negative control whose two central cut ranks
are unequal, so the symmetry hypothesis is exercised explicitly.

From the repository root:

```sh
g++ -O3 -std=c++17 tools/site_weight_crossings.cpp -o /tmp/site-weight-crossings
/tmp/site-weight-crossings profile 16 /tmp/site-square-t16.json 1 4 9 16 25 36 49 64 81 100 121 144 169 196 225 256
/tmp/site-weight-crossings profile 16 /tmp/site-beatty-t16.json 2 4 7 9 12 14 16 19 21 24 26 28 31 33 36 38
python3 tools/check_site_weight_crossings.py --output /tmp/site-weight-verification.json
```

A profile command takes the length, output path, and one nonnegative
integer weight per site. Every nonzero coefficient is recorded as
`[Hamming_count,weighted_sum,row_rank,column_rank]`. Exact integer square
root bounds can certify finite rank-sum comparisons; the displayed decimal
exponents describe the computed profiles and are not asymptotic evidence
on their own.

The complete seven-profile dataset, including all 16,824 nonzero coefficient
rank pairs, is retained in `evidence/site-weight-profiles.json`. Its auditor
checks every reachable bidegree, the complement rank identities, exact
maximum products, and integer lower and upper bounds for each rank sum.
It also flags the centrally symmetric Beatty prefix. To reproduce every
profile and the audit outside the repository, use:

```sh
python3 tools/reproduce_site_weight_profiles.py --directory /tmp/site-weight-profiles --output /tmp/site-weight-profiles/site-weight-profiles.json
```
