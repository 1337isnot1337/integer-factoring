# Every coordinate code of two reflection sites has pressure at least 5/7

## Result and scope

Take two copies of the four-letter reflection crossing. Retain any common
coordinate subset of their sixteen physical words on all four legs, with
at least two words retained. Keep the inherited natural endpoint gauges,
and allow every common-type distribution on the retained alphabet. Then
the local Smith-pressure criterion satisfies

```math
\boxed{\rho\ge\frac57>\frac1{\sqrt2}.}             \tag{1}
```

This is a complete finite result over all `65,519` nontrivial coordinate
codes. Exchanging the two sites reduces the exact computation to `33,269`
representatives. The proof uses integer coefficient inequalities valid at
every scale; a sampled numerical minimum is not used.

The computation is over `F_2`. It also applies after extending the
coefficient field to any field of characteristic two, since the ranks
and Smith orders of these binary polynomial matrices do not change under
field extension. The conclusion concerns these two sites, coordinate
subsets, and the stated gauges. It does not exclude general linear
retractions, different gauges, longer blocks, or arbitrary coefficient
recombinations. It does not improve the factoring exponent.

## 1. The two cuts share one normalized profile

Use the degree-twelve reflection seed `T(z;x,y,u,v)` from
[the auxiliary row-fiber calculation](triangular-auxiliary-row-fiber.md),
with physical height vector

```math
h=(0,3,1,4).
```

The two regular polynomial endpoint tensors can be written as

```math
\begin{aligned}
Z(z;x,y,u,v)&=z^{h(x)+h(y)-4}T(z;x,y,u,v),\\
W(w;x,y,u,v)&=w^{h(u)+h(v)+4}T(w^{-1};x,y,u,v).
\end{aligned}                                                \tag{2}
```

Every entry of `Z` is polynomial despite the displayed shift. The
explicit seed has the coefficient identity

```math
W(w;x,y,u,v)=Z(w;x,u,y,v).                         \tag{3}
```

The checker verifies (2)--(3) for all 256 local entries. Therefore the
row cut of `Z` and the column cut of `W` are the same polynomial matrix.
The identity persists under tensor powers and restriction to the same
coordinate code on all four legs.

Index a two-site word by `a+4b`, and set `H(a+4b)=h(a)+h(b)`.
For a code `C` of size `q`, let `sigma_C` be the complete increasing
Smith profile of the corresponding row-cut restriction of `Z tensor Z`.
It has `q^2` finite orders: evaluation at one is a permutation matrix.
Define

```math
g_C(r)=\sum_{\sigma\in\sigma_C}r^\sigma,
\qquad A_C(r)=\sum_{c\in C}r^{H(c)}.              \tag{4}
```

The zero profile before removal of the common scalar is `8+sigma_C`.
The regular infinity profile is `sigma_C`. Removing the extra infinity
regularization multiplies their pressure product by `s^8`, cancelling
the zero scalar. Thus the restored product and physical energies are

```math
B_C(s)=g_C(s^{-1})^2,
\qquad G(c)=4H(c).                               \tag{5}
```

This restoration is necessary: keeping the common scalar in one place
while omitting its charge at the other would change the pressure test.

## 2. The target becomes an integer polynomial inequality

The [common-type criterion](type-optimized-crossing-pressure.md) at a
trial exponent `rho` compares (5) with

```math
\left(\sum_{c\in C}s^{-G(c)/(4\rho)}\right)^{4\rho}.
```

This partition function already optimizes the common-type distribution.
At `rho=5/7`, put `x=s^{-1/5}`. For every finite `s>1`, the required
lower comparison is equivalent to

```math
P_C(x):=g_C(x^5)^7-A_C(x^7)^{10}\ge0,
\qquad 0<x<1.                                   \tag{6}
```

At `x=1`, its value is `q^14-q^10>0`. The exact replay certifies (6)
for every code by the following elementary coefficient test.

**Two cumulative sums suffice.** Remove the common monomial factor of
the two positive summands in (6), and write the remaining polynomial as
`P(x)=sum_(j=0)^D a_j x^j`. Put

```math
S_j=\sum_{i=0}^j a_i,
\qquad U_j=\sum_{i=0}^j S_i.                    \tag{7}
```

For every enumerated code, either all the `S_j` are nonnegative, or
all the `U_j` are nonnegative. In both cases
`S_D=P(1)=q^14-q^10>0`.

Indeed, the coefficients of `P(x)/(1-x)` through degree `D` are
the `S_j`; its remaining coefficients equal `S_D`. Similarly, the
coefficients of `P(x)/(1-x)^2` through degree `D` are the `U_j`,
and thereafter they form the sequence `U_D+r S_D`, for `r>=1`.
The relevant power series consequently has nonnegative coefficients
and is not zero. It converges and is strictly positive for `0<x<1`.
Multiplication by `(1-x)` or `(1-x)^2` proves (6).

Hence every finite-scale type root exceeds `5/7`, and its infimum over
all scales is at least `5/7`. Finally `2*5^2>7^2` proves the strict
comparison with the desired threshold in (1).

## 3. Exact enumeration and independent controls

There are `sum_(q=2)^16 binom(16,q)=65,519` codes. The involution
`a+4b -> b+4a` preserves both profiles and heights. The checker takes
one representative of each orbit and restores its multiplicity, then
verifies the binomial count separately at every alphabet size.

Each physical site preserves total occupation parity. The two-site
row cut therefore splits into four exact blocks, indexed by the two
parity differences. The checker verifies the ambient support identity
before using this decomposition. Every Smith order is resolved by
minimum-valuation elimination with an enlarged series window when
needed. No unresolved pivot is interpreted as a zero channel.

Polynomial powers in (6) are computed by exact Kronecker substitution.
All coefficients of either positive summand are nonnegative and have
total at most `q^14`. A binary slot wider than this bound prevents
every carry between coefficient positions. Decoding therefore gives
the exact integer coefficients, without modular reconstruction or
floating arithmetic.

The independent controls include all 65,536 two-site endpoint entries
and twenty binary block-Toeplitz ranks for five retained codes. The full
sixteen-letter profile is also checked against the subset-sum polynomial

```math
\bigl((1+r)(1+r^3)(1+r^5)(1+r^7)\bigr)^2.        \tag{8}
```

The evidence records the number of certificates of each cumulative
depth, the counts at every alphabet size, and a digest of every ordered
profile and certificate. Repeating the full checker reconstructs all
profiles from the seed; it does not trust a table of saved orders.

The completed enumeration has 27,770 certificates that pass after one
cumulative sum and 5,499 that require two. The largest polynomial degree
after removal of the common monomial is 2,240. The full replay took
approximately 566 seconds in the recorded environment.

## 4. Replay

```sh
python3 tools/check_reflection_coordinate_codes.py --full --progress \
  --output evidence/reflection-coordinate-codes.json
```

Without `--full`, the command checks the five explicit controls and
labels its output as a partial replay. Only the completed full output
supports the exhaustive claim.

See the [checker](../tools/check_reflection_coordinate_codes.py) and
[evidence](../evidence/reflection-coordinate-codes.json).
