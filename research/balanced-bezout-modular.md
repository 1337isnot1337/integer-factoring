# Complete modular Bezout polynomials with a binary optimal scale

## Constructive result

The [balanced-node denominator theorem](balanced-bezout-denominator.md)
identifies families whose optimal global Bezout denominator has an
exponential numerical exponent. This note computes their **complete
normalized Bezout polynomials modulo `2^p`** in polynomial bit time.
The output has at most `4p` possible powers of `X-1` in each polynomial,
with the degrees and optimal scale written in binary.

Fix positive integers `m,h,N,a,b,D_A,D_B`, with `h>=1`, satisfying

\[
 a\equiv b\equiv1\pmod{2^{h+1}},\qquad
 v_2(D_A)=v_2(D_B)=h+1,\qquad
 N\equiv1+2^h\pmod{2^{h+1}}.                       \tag{1}
\]

Put `A_i=a+D_A i`, `B_i=b+D_B i`, for `0<=i<m`, and define

\[
 P(X)=\prod_{i<m}(X-A_i),\qquad
 Q(X)=\prod_{i<m}(N-B_iX).                         \tag{2}
\]

There is a unique pair of rational polynomials `U_0,V_0`, each of
degree below `m`, for which `U_0P+V_0Q=1`. The exact denominator theorem
gives the smallest common dyadic clearing exponent

\[
 r=h(2m-1)-s_2(m-1),                              \tag{3}
\]

where `s_2` is binary digit sum. Set `U=2^r U_0`, `V=2^r V_0`.
Both belong to `Z_(2)[X]`, both have degree `m-1`, and both leading
coefficients are odd. In particular they are primitive at two.

Given the data in binary and numerical precision `p>=1`, the algorithm
returns every coefficient of `U,V modulo 2^p`, as short lists

\[
 U(X)\equiv\sum_{j=0}^{J_U}u_j(X-1)^{m-1-j},\qquad
 V(X)\equiv\sum_{j=0}^{J_V}v_j(X-1)^{m-1-j},\quad
 J_U,J_V<4p.                                    \tag{4}
\]

All omitted coefficients are zero modulo `2^p`. Its bit cost is
polynomial in the encoded input length and `p`. It enumerates neither
progression, expands no factorial, and never uses the numerical value
of `r` as a precision or allocation size.

This is stronger than a local Taylor evaluation: (4) represents the
entire polynomials over `Z/2^p`. The domain separation in (1) already
excludes a factor pair. The construction does not find an unknown
common root, and the complete factoring score remains `1/2`.

Implementation: [balanced_bezout_modular.py](../tools/balanced_bezout_modular.py).
Checker: [check_balanced_bezout_modular.py](../tools/check_balanced_bezout_modular.py).
Evidence: [balanced-bezout-modular.json](../evidence/balanced-bezout-modular.json).

## 1. Exact Newton coefficients and a terminal window

Interleave the roots by

\[
 x_{2i}=A_i,\quad x_{2i+1}=N/B_i,\qquad
 z_i=(x_i-1)/2^h.                               \tag{5}
\]

Every `z_i` lies in `Z_(2)`. The even-indexed nodes are even and the
odd-indexed nodes are odd. Equation (1) gives

\[
 v_2(z_i-z_j)=v_2(i-j).                           \tag{6}
\]

Let `f` be parity on the 2-adic integers and put `K=2m-1`.
The interpolating polynomial `F` of degree at most `K` has values
zero at every `A_i` and one at every `N/B_i`. Its Newton expansion is

\[
 F(X)=\sum_{k=1}^{K}c_k\prod_{i<k}(X-x_i),\qquad
 c_k=2^{-hk}f[z_0,\ldots,z_k].                   \tag{7}
\]

The constant Newton term is zero because `z_0` is even. The exact
divided-difference theorem in the companion note states, for `k>=1`,

\[
 v_2(f[z_0,\ldots,z_k])=k-1-v_2(k!).             \tag{8}
\]

It follows that the normalized coefficient `2^r c_k` has exact valuation

\[
 \delta_k=(h-1)(K-k)+v_2(K!)-v_2(k!).             \tag{9}
\]

These valuations are nondecreasing when `k` decreases. Among any
`K-k` consecutive positive integers, at least `floor((K-k)/2)` are
even. Therefore

\[
 K-k\ge2p\quad\Longrightarrow\quad\delta_k\ge p. \tag{10}
\]

Only the last `2p` possible Newton coefficients can survive modulo
`2^p`. This is independent of the binary digit sum of `m-1`.
It uses the exact coefficient valuations, not only the weaker fact
that a locally normalized inverse has integral coefficients.

## 2. Huge-order divided differences use only low symmetric degrees

Define the odd unit

\[
 D_k=(-1)^{k+1}\frac{k!}{2^{k-1}}f[z_0,\ldots,z_k]. \tag{11}
\]

The parity Mahler series is

\[
 f(z)=\sum_{n\ge1}(-1)^{n+1}2^{n-1}{z\choose n}. \tag{12}
\]

For fixed `k`, the isometry denominators in the Lagrange formula have
valuation at most `v2(k!)`. Since every binomial value on `Z_2` is
integral, `k!` times each divided difference of a binomial polynomial
is integral. Terms with `n<k` vanish. Thus, after the normalization
in (11), a term with `n=k+j` has valuation at least `j`.
Only `0<=j<p` is needed modulo `2^p`.

Write

\[
 E_n(T)=\prod_{0\le s<n}(1-sT),\qquad
 G_k(T)=\prod_{i=0}^{k}(1-z_iT),
\]

\[
 H_{k,j}=[T^j]\frac{E_{k+j}(T)}{G_k(T)},\qquad
 B_{k,j}=\prod_{s=1}^{j}(k+s),\quad B_{k,0}=1.   \tag{13}
\]

The divided difference of a monomial of degree `k+t` is the complete
homogeneous function `h_t(z_0,...,z_k)`. Expanding the falling factorial
and using this identity proves

\[
 {z\choose k+j}[z_0,\ldots,z_k]
    =\frac{H_{k,j}}{(k+j)!}.
\]

Consequently the exact finite formula for the required unit is

\[
 \boxed{\quad
 D_k\equiv\sum_{j=0}^{p-1}(-1)^j2^j
            \frac{H_{k,j}}{B_{k,j}}\pmod{2^p}.\quad} \tag{14}
\]

The Lagrange argument above also proves that every quotient
`H_(k,j)/B_(k,j)` is 2-adically integral. This justifies division by
the power-of-two part; it is not a formal division by a nonunit in
the residue ring.

### Every precision loss has polynomial size

The exact integer `B_(k,j)` is a product of at most `p-1` integers
with `O(log(k+p))` bits. Let

\[
 g=\max_{j<p}v_2(B_{k,j})
   =v_2(B_{k,p-1})\le(p-1)\log_2(k+p).           \tag{15}
\]

Compute all coefficients in (13) at precision `p+g`. Divide a
coefficient by the known power `2^v2(B_(k,j))`, then invert the odd
part of `B_(k,j)` modulo `2^p`. Equation (15) guarantees that all
`p` requested bits remain. The literal shifts and allocated precision
are polynomial in `p` and the binary input length. A factorial of
size `k!` is never formed.

The first `p` coefficients of `E_n` come from exact small-index
elementary symmetric functions of `0,...,n-1`. The existing
[power-sum and Newton recurrences](interval-unit-jets.md) compute them
in polynomial bit time. To compute `G_k`, split its interleaved prefix
into its even and odd nodes. An even-node factor is

\[
 1-\left(\frac{a-1}{2^h}+\frac{D_A}{2^h}i\right)T.
                                                               \tag{16}
\]

An odd-node factor is

\[
 1-\frac{N/B_i-1}{2^h}T
 =\frac{B_i+(B_i-N)T/2^h}{B_i}.                 \tag{17}
\]

Each numerator product in (16)--(17) is an affine progression of
jets with unit constant term and an even step coefficientwise.
The denominator in (17) is an odd integer progression product.
The same unit-jet engine therefore computes `G_k` and its truncated
inverse at precision `p+g`, without enumerating the nodes. The number
of retained coefficients is only `p`, even when `k` is exponential.

Finally, write `k!=2^v2(k!)*u_k`. Its odd part `u_k modulo 2^p` is
computed by the existing normalized affine-product algorithm. Equivalently,
split `1,...,k` by exact valuation; each of `O(log k)` buckets is an
odd arithmetic progression. All bucket endpoints are binary integers.
Equation (11) now gives the actual normalized Newton coefficient:

\[
 2^r c_k=(-1)^{k+1}2^{\delta_k}D_k u_k^{-1}.
                                                               \tag{18}
\]

For successive retained indices, the odd factorial part updates by
division by the odd part of `k`. Only one full factorial-unit query
is needed. The implementation also uses the smaller local precision
`p-delta_k` in (14).

## 3. Entire polynomial support is short in the shifted basis

Put `Y=X-1` and let `H=2^r F`. Every `x_i-1` is divisible by `2^h`.
Thus the coefficient at distance `j` from the leading term of

\[
 N_k(Y)=\prod_{i<k}(Y-(x_i-1))
\]

has valuation at least `hj`. Put `q=floor((p-1)/h)`. In the reverse
polynomial `T^K H(1/T+1)`, a retained Newton term contributes only
at degrees

\[
 (K-k)+j\le(2p-1)+q.                            \tag{19}
\]

Equation (10) discards every earlier Newton term. Equation (19)
discards every further coefficient of every retained term. Therefore
all coefficients of this reverse polynomial beyond

\[
 J_H=\min(K,2p-1+q)                             \tag{20}
\]

are zero modulo `2^p`. The implementation uses the actual smallest
retained `k`, which can reduce this bound.

Let `P^*(T)=T^m P(1+1/T)` and `Q^*(T)=T^m Q(1+1/T)`.
Their constant coefficients are units. For `j>=0`, their degree-`j`
coefficients, and those of their power-series inverses, satisfy

```math
 v_2([T^j](P^*)^{\pm1})\ge(h+1)j,\qquad
 v_2([T^j](Q^*)^{\pm1})\ge hj.                  \tag{21}
```

For the products this follows from the valuations of their roots
after translation: `A_i-1` is divisible by `2^(h+1)`, while
`N/B_i-1` has valuation exactly `h`. Multiplication by the odd leading
coefficient of `Q` changes no bound. The usual inverse-series recurrence
proves the same coefficient bounds for the inverses.

The interpolation conditions give the exact polynomial identities

\[
 H=UP,\qquad 2^r-H=VQ.                          \tag{22}
\]

Their degrees imply `deg U=deg V=m-1`. Reverse them, and retain
degrees at most `m-1`. The term `2^r T^K` in the second reversed
numerator lies beyond this range. Thus

\[
 U^* = H^*/P^*\pmod{T^m},\qquad
 V^* =-H^*/Q^*\pmod{T^m}.                      \tag{23}
\]

Equations (20)--(21) show that all omitted coefficients vanish modulo
`2^p` once the reverse index exceeds

\[
 J_U=\min(m-1,J_H+\lfloor(p-1)/(h+1)\rfloor),
 \qquad
 J_V=\min(m-1,J_H+q).                           \tag{24}
\]

Both lists have fewer than `4p` possible indices. Hence (23), computed
with truncated series of length `O(p)`, gives **every** coefficient
of the modular Bezout polynomials, not only selected high coefficients.
All exponent arithmetic and all degree comparisons use binary integers.

The two leading coefficients are odd because `delta_K=0`, `D_K` is
odd, and the leading coefficients of `P,Q` are units. The output
therefore retains nonzero normalized polynomial data even when the
scalar `2^r` is zero at the requested modulus.

## 4. Complete algorithm and bit costs

The implementation performs these finite operations.

1. Compute `K,r` and the nondecreasing gaps (9) by binary integer
   arithmetic. Stop at the first gap at least `p`; (10) bounds the
   number of retained indices by `2p`.
2. Compute the odd part of `K!` modulo `2^p` using valuation buckets.
   For each retained `k`, evaluate (14) at local precision
   `p-delta_k` and form (18).
3. Compute a reversed terminal prefix product to degree
   `floor((p-1)/h)` by (16)--(17), with the root shifts multiplied
   by `2^h`. Remove one final linear factor as `k` decreases. Its
   reversed factor has constant one, so the series division is by
   a unit. Accumulate the short reverse polynomial of `H`.
4. Compute the short reversed input products and apply (23).
   Equations (20)--(24) certify that the discarded coefficient
   tails are zero modulo `2^p`.

Every ring operation uses a polynomial number of coefficients. The
largest allocated precision is bounded by `O(p(log m+log p))`, in
addition to the encoded input coefficient lengths. Small-index
integer symmetric functions have polynomial bit length by the
power-sum bounds in the cited unit-jet note. Modular powering of
an interval unit uses `O(log m)` multiplication steps.

The implementation reduces the elementary-sum work in the guarded
jets by also using their nilpotent degree. Suppose an even increment
jet `W` has degree cap `d`, every coefficient has valuation at least
`v>=1`, and its constant has valuation `c>=v`. A retained coefficient
of degree at most `d` in `W^j` uses at most `d` nonconstant factors.
Its valuation is at least

\[
 cj-(c-v)d.                                      \tag{25}
\]

Thus terms beyond `floor((P-1+(c-v)d)/c)` vanish at precision `P`.
If the constant is zero, the power is zero already for `j>d`.
This is a valid extra cutoff in the same finite elementary-product
identity, and avoids computing coefficients that cannot survive.
The construction does not rely on a numerical timing observation
for its polynomial bound.

Because (4) is a complete polynomial representation, it also gives
evaluation and Taylor jets at any specified modular argument by
binary powering the displayed powers of `X-1`. There is no need
to recover the unnormalized rational polynomials or to allocate
absolute precision `r+p`.

## 5. Meaning for the extraction problem

The denominator family has a short **normalized global** representation,
including the optimal scale and all its modular coefficients. An
exponential value of `r` therefore does not by itself prevent this
kind of compressed global algebra. The raw scalar in a Bezout
identity and its normalized polynomial coefficients have different
precision requirements.

The arithmetic separation in (1) is essential to the proof. It
gives the isometry ordering, the exact parity interpolant, and the
known coprimality. Applying the construction to general overlapping
factor domains would require new information about their common
roots. No such extraction algorithm is inferred from this special
family.

The output is modular at requested precision. It does not encode
all exact odd rational denominators of `U_0,V_0` in one finite list.
It does give their complete primitive reductions at every requested
precision, consistently with the unique exact rational Bezout pair.

## 6. Independent controls

Run:

```sh
python3 tools/check_balanced_bezout_modular.py \
  --output evidence/balanced-bezout-modular.json
```

The small checker constructs the rational inverse of `Q modulo P`
by ordinary Lagrange interpolation at the explicit `A` nodes. It
then performs exact polynomial division to get the other Bezout
polynomial. This calculation uses neither parity interpolation nor
the compressed Newton formula.

Thirty-six varied progression cases compare all 308 normalized
coefficients and verify every complete polynomial identity. Another
171 controls form the parity divided differences directly by the
rational difference table. Five invalid domain or precision contracts
are rejected. Heights one through three, unequal odd step multiples,
and the case `m=1` are included.

The retained implicit controls use domain sizes `2^200+13`, `2^200`,
and `2^128+9`, with 402-bit `N`. They compare two requested precisions
and check the complete modular polynomial identity using the proved
short reverse representations. The exact power-of-two count tests
a large factorial-ratio guard, so the evidence does not hide that
precision cost. No input nodes are listed in these cases.

The new coefficient construction combines the companion parity
divided-difference theorem with the existing unit interval products.
It was developed with substantial OpenAI assistant support. It has
ordinary proofs and exact computational controls, without new Lean
formalization or a stronger factoring score.
