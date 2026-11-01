# Exact parity interpolation and balanced Bezout denominators

## Result

For the arithmetic domains

\[
 A=B=\{1,5,9,\ldots,1+4(m-1)\},\qquad N=3\pmod4,
\]

let `P(X)=prod_a(X-a)` and `Q(X)=prod_b(N-bX)`. The least exponent `r`
for a polynomial Bezout identity over the rationals with odd denominators is

\[
 UP+VQ=2^r,\qquad
 \boxed{r=m+v_2((m-1)!)=1+v_2((2m-1)!).}          \tag{1}
\]

This holds for every `m>=1`. More generally it holds for independently
deformed node sets whose enumerations satisfy

\[
 a_i,b_i=1\pmod4,\qquad
 v_2(a_i-a_j)=v_2(b_i-b_j)=2+v_2(i-j).            \tag{2}
\]

The [earlier residue-ball theorem](interval-bezout-precision.md) gives the
sharp all-node bound `r=2m-2` when `m=2^s+1`. Equation (1) determines the
finer denominator for every size under (2). Irregular node sets that
violate (2) need not obey (1).

The proof is an exact interpolation statement: parity divided differences
have the same valuation on every ordered dyadic isometry. It does not
supply a faster general rational gcd or factor-count algorithm. The complete
factoring exponent remains `1/2`.

## 1. Parity has a uniformly convergent binomial expansion

Let `epsilon(z)` be zero for even `z in Z_2` and one for odd `z`. Then

\[
 \epsilon(z)=\sum_{n=1}^{\infty}
        (-1)^{n+1}2^{n-1}{z\choose n}.           \tag{3}
\]

Here is a direct proof including the integrality needed below. For fixed
`n`, the rational polynomial `binom(z,n)` is continuous on `Z_2` and is an
integer at every nonnegative integer. These integers are dense in `Z_2`,
which is closed in `Q_2`; hence `binom(z,n)` belongs to `Z_2` for every
`z in Z_2`. The `n`th summand in (3) has valuation at least `n-1`
uniformly in `z`, so the series converges uniformly.

At a nonnegative integer `z`, the series terminates and the ordinary
binomial identity gives

\[
 \frac{1-(-1)^z}{2}
 =-\frac12\sum_{n=1}^{z}(-2)^n{z\choose n}.
\]

Thus (3) holds on a dense set. Both sides are continuous, proving it
throughout `Z_2`. No integrality assertion about binomial polynomials on
larger unramified coefficient rings is used.

## 2. Every divided difference has an exact valuation

Let `z_0,...,z_k in Z_2` satisfy

\[
 v_2(z_i-z_j)=v_2(i-j)\quad(i\ne j).             \tag{4}
\]

For `k>=1`, write `d_k=epsilon[z_0,...,z_k]` for the divided difference.
Then

\[
 \boxed{v_2(d_k)=k-1-v_2(k!).}                   \tag{5}
\]

Indeed, its Lagrange expression is a fixed finite linear combination of
the function values, with denominator valuations

\[
 v_2\left(\prod_{j\ne i}(z_i-z_j)\right)
  =v_2(i!)+v_2((k-i)!)\le v_2(k!).               \tag{6}
\]

For every `n`, binomial integrality therefore gives

\[
 v_2\left({Z\choose n}[z_0,...,z_k]\right)
       \ge-v_2(k!).                             \tag{7}
\]

Apply this fixed linear functional to the uniformly convergent series
(3). Terms with `n<k` vanish. At `n=k`, the divided difference is the
leading coefficient `1/k!`, so that summand has valuation
`k-1-v2(k!)`. Every remaining summand, with `n>=k+1`, has strictly larger
valuation by (7). The tail cannot cancel the first term. This proves (5).

The equality includes every prefix of a longer isometric node sequence.
It uses the full ordering condition (4), not merely pairwise distinctness
or the presence of the two parities.

## 3. The affine interpolation denominator

Let `F` be the polynomial of degree below `L`, for `L>=2`, interpolating
parity on `z_0,...,z_(L-1)`. Its Newton expansion is

\[
 F(Z)=d_0+\sum_{k=1}^{L-1}d_k
                \prod_{j=0}^{k-1}(Z-z_j),        \tag{8}
\]

where `d_0=epsilon(z_0)`. Fix an integral center `c in Z_2` and an integer
`h>=1`, and substitute `Z=(X-c)/2^h`. The resulting Newton basis is

\[
 B_k(X)=\prod_{j=0}^{k-1}(X-(c+2^hz_j)).
\]

Each `B_k` is monic with integral coefficients. The change of basis between
these polynomials and `1,X,...,X^(L-1)` is triangular with diagonal one,
and both it and its inverse have integral entries. Therefore the minimum
coefficient valuation is the same in either basis. The coefficient of
`B_k` for `k>=1` has valuation

\[
 v_2(2^{-hk}d_k)=(1-h)k-1-v_2(k!).               \tag{9}
\]

This is nonincreasing with `k`. The constant coefficient in the Newton
basis is integral, so the exact dyadic denominator exponent is

\[
 \boxed{\rho_L(h)=(h-1)(L-1)+1+v_2((L-1)!).}    \tag{10}
\]

For rational nodes and center, (10) is also the least power of two that
makes every coefficient rational with odd denominator. The integral
basis argument avoids assuming that the leading coefficient alone
controls an arbitrary polynomial's content.

## 4. Interleave the integer and reciprocal roots

Assume (2), index both node sets by `0<=i<m`, and put

\[
 z_{2i}=\frac{a_i-1}{2},\qquad
 z_{2i+1}=\frac{N/b_i-1}{2}.                     \tag{11}
\]

The first nodes are even and the second are odd. Differences within the
first set have valuation `1+v2(i-j)`. For the second set, `N,b_i,b_j` are
odd, so

\[
 v_2(z_{2i+1}-z_{2j+1})
   =v_2\left(\frac{N(b_j-b_i)}{2b_ib_j}\right)
   =1+v_2(i-j).
\]

Every cross-parity difference is a unit. Thus the interleaved sequence
satisfies (4) with `L=2m`.

Let `F` be its odd-parity interpolant. Then `F((X-1)/2)` is zero at
the roots of `P` and one at the roots of `Q`. If `V0` is the unique
degree-below-`m` inverse of `Q` modulo `P`, uniqueness of interpolation gives

\[
 1-F((X-1)/2)=V_0(X)Q(X).                        \tag{12}
\]

Both sides have degree at most `2m-1` and have the same values at all
`2m` distinct roots. The orientation in (12) matters: `V0Q` is the even
indicator, not the odd indicator.

By (10), the minimum coefficient valuation of the odd interpolant is
`-1-v2((2m-1)!)<0`. Subtracting it from one preserves that minimum.
The polynomial `Q` is integral and has odd leading coefficient. Gauss's
content identity therefore says that multiplication by `Q` preserves
the minimum coefficient valuation of `V0`. The minimum-denominator
argument in the preceding note proves that this is exactly the exponent
required by a Bezout identity of any degree. Hence

\[
 r=1+v_2((2m-1)!)=m+v_2((m-1)!),
\]

where the last equality follows by separating the even factors in
`(2m-1)!`. This proves (1).

## 5. Exact consequences and limits

For consecutive step-four domains, the exponent can be evaluated from
the binary encoding of `m`:

\[
 r=2m-1-s_2(m-1),                               \tag{13}
\]

where `s_2` counts one bits. This computes the exact required scalar
normalization without expanding a polynomial. It does not construct the
large-degree Bezout coefficient arrays by itself.

The theorem extends beyond arithmetic progressions. Independent
integer isometries of the two index sets preserve (2); examples include
`i -> i+2i^2`, XOR with a fixed integer, and sufficiently
divisible perturbations of individual nodes. Each application must
check the distance condition. The explicit irregular set
`B={1,5,9,17}` with `A={1,5,9,13}`, `N=3`, violates it and has `r=6`
instead of the inapplicable value five.

There is a useful constructive corollary for the scaled interpolant.
For `h=1`, scale (8) after substitution by `2^r`, with
`r=1+v2((L-1)!)`. Its `k`th Newton coefficient, for `k>=1`, has valuation

\[
 v_2((L-1)!)-v_2(k!).                            \tag{14}
\]

Modulo `2^p`, every term with `L-1-k>=2p` vanishes: the product of those
consecutive integers contains at least `p` even factors. At most `2p`
terminal Newton terms remain. This is a representation statement; it
does not by itself compute their huge-index coefficient units. The
low-surplus moment formula below makes that remaining task explicit.

## 6. Huge-index coefficient units reduce to low moments

Let `e_j` and `h_j` denote elementary and complete homogeneous symmetric
functions, respectively. For `k>=1,j>=0`, put

\[
 C_{k,j}=\sum_{t=0}^j(-1)^{j-t}
 e_{j-t}(0,1,\ldots,k+j-1)\,h_t(z_0,\ldots,z_k).
\]

Expanding the falling factorial and using
`Z^(k+t)[z_0,...,z_k]=h_t(z_0,...,z_k)` gives

\[
 {Z\choose k+j}[z_0,\ldots,z_k]
       =\frac{C_{k,j}}{(k+j)!}.                 \tag{15}
\]

The monomial identity follows, for example, by taking coefficients in
the divided-difference identity for `1/(1-ZT)`, whose value is
`T^k/prod_i(1-z_i T)`.

Define the normalized unit

\[
 D_k=(-1)^{k+1}k!\,2^{1-k}d_k.
\]

The proof of (5) shows that `D_k=1 mod2`. More precisely, for every
`p>=1`, the uniformly convergent expansion (3) and (15) give

\[
 \boxed{D_k\equiv\sum_{j=0}^{p-1}(-1)^j2^j
    \frac{C_{k,j}}{(k+1)(k+2)\cdots(k+j)}\pmod{2^p}.} \tag{16}
\]

The empty denominator is one. By (7), each displayed ratio belongs to
`Z_2`; it is not a free division by an even integer. To evaluate it,
retain an additional

\[
 v_2((k+1)\cdots(k+j))=O(j\log(k+j))
\]

bits before its dyadic division, and invert only the remaining odd part.
Terms with `j>=p` have valuation at least `p` and are absent from (16).

Thus the units of enormous-index divided differences require only
symmetric moments through degree `p-1`, and a polynomial guard precision
in `p` and the binary length of `k`. The interval elementary functions in
(16) have the existing exact power-sum construction. For the interleaved
integer/reciprocal nodes, the other moments are short jets of explicit
arithmetic products. This identifies a constructive continuation without
assuming that the node list or a factorial is expanded.

## 7. Independent exact checks

The [checker](../tools/check_balanced_bezout_denominator.py) uses exact
rational divided differences on ordinary integer nodes and nonlinear
rational isometries. It constructs the whole interpolation polynomial
on these finite controls, checks its values, and compares its coefficient
denominators after three different affine scalings with (10).

For the factoring specialization, it constructs reciprocal interleavings
from two independent index isometries. Scaled Lagrange inversion of `Q`
modulo `P` and a separate Hensel inverse in the normalized quotient agree
with (1). The irregular-node counterexample is retained to check the
scope of the isometry hypothesis.

The [evidence](../evidence/balanced-bezout-denominator.json) records 1,932
divided-difference valuations, 598 interpolation values, 138 affine
coefficient norms, and 52 independent pairs of Bezout inversions. It also
checks 132 low-surplus symmetric identities and 132 normalized-unit
truncations against direct exact rational divided differences.
Representative node sets and source hashes are retained. A binary-described
progression with `m=2^200+13` also receives its exact exponent from (13);
that record applies the theorem and does not expand the polynomial.

```sh
python3 tools/check_balanced_bezout_denominator.py \
  --output evidence/balanced-bezout-denominator.json
```

The general factor-count contraction and any stronger complete factoring
exponent remain unresolved.
