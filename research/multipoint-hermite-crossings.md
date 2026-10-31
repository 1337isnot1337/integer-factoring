# Multipoint Hermite interpolation for effective crossing pieces

## Result and scope

A polynomial crossing can use jets at several distinct points, including
infinity. The jet multiplicities share one degree budget. Each point uses
both of its local cut profiles, and the two cut estimates may use different
positive scales. This gives an explicit sufficient criterion for an
effective crossing family with polynomially many pieces.

The construction differs from the usual two-place coefficient estimate.
That estimate bounds the two ranks of one monomial coefficient using
opposite places. Here, a piece is one local jet; both rank bounds come from
that jet's place. A Hermite reconstruction identity then combines the
pieces from all selected places. Neither scheme is asserted to dominate
the other.

An explicit three-letter polynomial over `F_3` satisfies the multipoint
criterion with `rho < 19/20`. This is a complete ordinary construction,
including exact endpoint and rational pressure checks, but it is much
weaker than the existing `rho = 1/sqrt(2)` family. The achieved factoring
exponent remains one half. Section 6 proves a general limitation: every
certificate from this unweighted, same-place criterion has `rho > 3/4`.
The transfer theorem remains useful as a precise interface for extensions
that change its rank estimates or introduce physical type weights.

## 1. A finite polynomial seed and its local profiles

Let `F` be a fixed finite field, let `q >= 2`, and let

```math
K(z;x,y,u,v)\in F[z],\qquad 0\le x,y,u,v<q,
```

have degree at most the fixed integer `D`. Assume the exact endpoint

```math
K(1;x,y,u,v)=[x=v][y=u].                                  \tag{1}
```

Choose finitely many distinct points `a_j` of `P^1(F)`, none equal to one.
At a finite point use `K_j(t)=K(a_j+t)`. At infinity use

```math
K_\infty(t)=t^D K(t^{-1}).                               \tag{2}
```

Both are polynomial tensors. Flatten in the two cuts

```math
(x,y)\mid(u,v),\qquad (x,u)\mid(y,v).
```

Their generic ranks are `q^2`, because (1) is invertible in both cuts.
Write their nonnegative Smith orders over `F[[t]]` as
`nu_(j,1),...,nu_(j,q^2)` and `mu_(j,1),...,mu_(j,q^2)`.
For real `s >= 1`, define

```math
U_j(s)=\sum_i s^{-\nu_{j,i}},\qquad
V_j(s)=\sum_i s^{-\mu_{j,i}}.                            \tag{3}
```

The normalization in (2) is essential: infinity orders are those of the
degree-`D` homogeneous section. They are not the possibly negative orders
of `K(t^-1)` alone.

## 2. The multipoint criterion

Choose nonnegative rational numbers `d_j` with

```math
\sum_j d_j=D.                                           \tag{4}
```

Choose scales `s_j,r_j >= 1`. Put

```math
F_j=\log U_j(s_j)+\log V_j(r_j)
       +d_j\log(s_jr_j).                               \tag{5}
```

**Theorem.** If a rational `rho` satisfies `1/2 <= rho < 1` and

```math
\max_j F_j<4\rho\log q,                                \tag{6}
```

then the seed supplies effective scale-one crossing pieces over `F` with
parameter `rho`. There are `O(lambda)` pieces at binary label width
`lambda`. Their complete tables and both rank factorizations can be
constructed in `2^O(lambda)` bit operations, with fixed constants depending
only on the displayed seed and certificate.

The existing [scale-one finite-field evaluator](scale-one-modular-transfer.md)
and [fixed-prime recovery interface](fixed-prime-transfer-source-audit.md)
then apply. If `F` is a fixed extension of `F_p`, embed the original
`F_p` tables in `F` and perform the evaluator's operations there. The
result remains the embedded `F_p` contraction, so scalar recovery still
uses the same prime. A certificate with `rho < 1/sqrt(2)` would therefore
improve the factoring exponent. Section 6 shows that the unweighted
criterion (6) cannot supply such a certificate.

### A rational verification form

When all scales are rational, choose an integer `L > 0` for which
`Ld_j` and `4Lrho` are integers. Condition (6) is exactly the finite list
of rational comparisons

```math
U_j(s_j)^L V_j(r_j)^L(s_jr_j)^{Ld_j}<q^{4L\rho}.          \tag{7}
```

Thus floating logarithms are unnecessary for a final certificate.

## 3. A complete reconstruction identity

Let `m` be a positive multiple of the denominators in (4), and set
`N=Dm`. Pick one selected point `j_0` and use multiplicities

```math
h_j=md_j+[j=j_0],\qquad \sum_j h_j=N+1.                 \tag{8}
```

For a polynomial `P` of degree at most `N`, retain coefficients
`[t^k]P(a_j+t)`, `0 <= k < h_j`, at finite points. At infinity retain
`[t^k]t^N P(t^-1)` in the same range. These are Hasse jets: ordinary
derivatives divided by factorials are not used in positive characteristic.

The resulting linear map from polynomials of degree at most `N` to
`N+1` jet coordinates is invertible. If all coordinates vanish, the
finite factors `(z-a_j)^h_j` divide `P`. The infinity conditions also
give `deg P <= N-h_infinity`. The required finite vanishing degree is
`N+1-h_infinity`, which is larger. Therefore `P=0`. Equal finite
dimensions now give invertibility.

It follows that explicitly computable field elements `c_(j,k)` satisfy

```math
P(1)=\sum_j\sum_{k=0}^{h_j-1}c_{j,k}J_{j,k}(P)          \tag{9}
```

for every such `P`. Gaussian elimination on the `N+1` dimensional Hermite
matrix computes these coefficients in `poly(m)` field operations.

Apply (9) entrywise to `K(z)^(tensor m)`. Equation (1) gives the crossing
on `q^m` letters. Its pieces are exactly

```math
A_{j,k}=c_{j,k}[t^k]K_j(t)^{\otimes m}.                 \tag{10}
```

At infinity, tensoring (2) gives the required factor `t^(Dm)`. Thus there
is no lost degree shift in (10). Zero interpolation coefficients simply
give zero pieces.

## 4. Rank and coefficient costs

For one cut, a local Smith factorization has the form

```math
M(t)=U(t)\operatorname{diag}(t^{\nu_1},\ldots,t^{\nu_{q^2}})V(t),
```

where `U,V` are invertible over `F[[t]]`. In a tensor power, a channel
word `i_1,...,i_m` starts at order `nu_(i_1)+...+nu_(i_m)`. If this sum
is at most `k`, its contribution to coefficient `k` has rank at most
`k+1`: expand the corresponding outer product of two power-series
vectors and sum over their two coefficient orders. Consequently

```math
\begin{aligned}
\operatorname{rank}[t^k]M(t)^{\otimes m}
&\le (k+1)\#\{(i_1,\ldots,i_m):\textstyle\sum_a\nu_{i_a}\le k\}\\
&\le (k+1)s^k\left(\sum_i s^{-\nu_i}\right)^m,
\qquad s\ge1.                                         \tag{11}
\end{aligned}
```

The second inequality also holds at `s=1`. Apply (11) to both cuts in
(10). Since every retained `k` is at most `md_j`,

```math
r(A_{j,k})c(A_{j,k})
\le (Dm+1)^2\exp(mF_j).                               \tag{12}
```

There are exactly `Dm+1` listed pieces. Thus (6) gives the required
per-piece bound on the sum of the two binary factor widths, with a
logarithmic overhead. It also bounds the sum of square-root rank costs
by `poly(m) q^(2rho m)`.

For an arbitrary binary width `lambda`, choose the least admissible `m`
with `q^m >= 2^lambda`, then restrict every leg to the first `2^lambda`
words in a fixed order. This changes `m log_2 q` by only a fixed constant,
preserves the crossing endpoint, and cannot increase either rank.
Ordinary rank factorization after restriction caps each factor width at
`2lambda`. The strict margin in (6) absorbs rounding and the polynomial
factors.

For the finitely many widths before the strict margin absorbs the fixed
overheads, use the coordinate rank-one decomposition of the crossing.
Enlarging the fixed constants covers its finitely many pieces, factor
widths, and construction costs. Thus the effective family is defined at
every positive label width, not only eventually.

Every entry of `K(z)^(tensor m)` is a product of `m` fixed polynomials.
Writing all `q^(4m)` entries and their `O(m)` coefficients costs
`exp(O(m))` field operations. Taylor shifts and infinity reversal have
polynomial overhead per entry. Exact elimination on each flattened piece
also costs `exp(O(m))`. The field is fixed, so each field operation has
constant bit cost. These procedures construct complete tables and both
factorizations; the theorem does not assume a free tensor-entry oracle.

## 5. A concrete ternary certificate

Index the physical letters by `F_3`, and define

```math
K(z;x,y,u,v)
=[x+y=u+v]\bigl(1-(x-u-z(v-u))^2\bigr).                \tag{13}
```

All arithmetic in (13) is in `F_3`. It has degree at most two. At `z=1`,
the second factor is `[x=v]`; the conserved sum then forces `y=u`.
Hence (1) holds without an auxiliary boundary or an endpoint repair.

The complete local profiles are:

| Point | First cut orders | Second cut orders |
|---|---|---|
| `0` | `0` nine times | `0,1,1,1,1,2,4,4,4` |
| `2` | `0,0,0,1,1,1,2,2,2` | `0` nine times |
| infinity | `0,0,0,1,1,1,2,2,2` | `0` nine times |

These small finite calculations can be certified by block-Toeplitz ranks.
For the nontrivial second-cut profile at zero, the ranks at jet lengths
`1,...,5` are `1,6,12,18,27`. Their consecutive increments are
`1,5,6,6,9`, recovering all nine orders. For either other nontrivial
profile, the ranks at lengths `1,2,3` are `3,9,18`, recovering three
orders each at zero, one, and two. The endpoint proves generic rank nine,
so a later unseen pivot is impossible once all nine orders are recovered.
The checker constructs these Toeplitz matrices directly from (13).

Use

```math
d_0=15/14,\qquad d_2=d_\infty=13/28.
```

At zero choose scales `(s_0,r_0)=(1,9/5)`. At the other two places use
`(s_j,r_j)=(5/2,1)`. Define the rational numbers

```math
\begin{aligned}
B_0&=9\left(1+4/(9/5)+(9/5)^{-2}+3(9/5)^{-4}\right),\\
B_1&=27\left(1+(5/2)^{-1}+(5/2)^{-2}\right).
\end{aligned}
```

Exact integer comparisons after clearing denominators give

```math
B_0^{140}(9/5)^{150}<3^{532},\qquad
B_1^{140}(5/2)^{65}<3^{532}.                            \tag{14}
```

Here `532=4*140*(19/20)`. Equations (4) and (7) therefore hold with
`rho=19/20`. For `m=28r`, one explicit reconstruction uses multiplicities
`30r+1,13r,13r`, totaling `2m+1`.

## 6. A universal floor for the unweighted criterion

The flexibility to move interpolation points does not remove the global
degree budget. The following bound includes scalar factors that vanish
at selected points, including content at infinity from an oversized
degree bound.

**Theorem.** Put `n=q^2`, and define

```math
f_n=\log(2n)-\frac{n+1}{2n}\log(n+1).
```

For every scalar polynomial crossing in Section 1, every choice of
distinct rational places over the working field, every `d_j >= 0` with
sum `D`, and every pair of scales `s_j,r_j >= 1`,

```math
\max_j F_j\ \ge\ \log n+f_n.                           \tag{15}
```

Consequently every successful strict certificate (6) satisfies

```math
\rho>\frac12+\frac{f_n}{2\log n}>\frac34.               \tag{16}
```

The first lower bound in (16) is approximately `0.8871987352` at `q=2`,
`0.8666347614` at `q=3`, `0.8535669207` at `q=4`, and `0.8445011436`
at `q=5`. It tends to `3/4` as the alphabet grows. Thus the
unweighted same-place Hermite criterion cannot reach the known
`1/sqrt(2)` crossing parameter.

The proof also applies to a tensor with fixed auxiliary indices if both
of its local flattenings have **exact generic rank `n`**, their entries
are rearrangements of the same tensor entries, and the criterion uses
their complete lists of `n` Smith channels. For example, putting both
auxiliary indices on the row side of either cut leaves exactly `n`
columns. A scalar endpoint contraction equal to the crossing then proves
rank `n`. An endpoint contraction alone only proves rank at least `n`;
it does not justify this extension when an auxiliary index on each side
allows larger generic rank. No auxiliary transfer theorem is assumed
here: the assertion is a floor for a criterion with the stated channel
and common-content hypotheses.

### Remove homogeneous content before charging local orders

Homogenize every tensor entry to degree `D`. Let `G` be their common
homogeneous gcd, let `g=deg G`, and divide it out. The remaining tensor
has homogeneous degree `D'=D-g`. Its entries have no common zero on
the projective line. In particular, at every selected place one entry
is a unit, so each local cut has a zero Smith order.

Write `c_j=ord_(a_j) G`. At point `j`, both original profiles equal
their primitive profiles plus `c_j`. Write `u_j` and `v_j` for the two
primitive mean orders. For each cut, fix one nonzero `n` by `n` minor.
It is a homogeneous polynomial of degree `nD'`. At each selected place
its vanishing order is at least the sum of the local Smith orders.
The selected places are distinct and rational, so

```math
\sum_j u_j\le D',\qquad
\sum_j v_j\le D',\qquad
\sum_j c_j\le g.                                      \tag{17}
```

The same minor works at all places; a separate minor is not chosen to
spend the degree budget repeatedly. For non-rational closed points,
residue degrees would have to be included throughout. Those points are
outside the present statement.

Let `J_+={j:d_j>=c_j}` and, for these indices, put
`delta_j=d_j-c_j`. This set is nonempty, and

```math
\sum_{j\in J_+}\delta_j
\ge\sum_j(d_j-c_j)\ge D-g=D'.                         \tag{18}
```

Equations (17)--(18) imply that at least one `j` in `J_+` has
`u_j+v_j <= 2delta_j`. Otherwise the strict opposite inequalities
would contradict the summed bounds. This argument also covers
`D'=0` and `delta_j=0`.

At this chosen place, put `A=log s_j` and `B=log r_j`. Subtracting
the scalar content from the two profiles changes (5) exactly into

```math
F_j=\log\sum_i e^{-A\nu'_i}
    +\log\sum_i e^{-B\mu'_i}
    +\delta_j(A+B).                                  \tag{19}
```

Thus scalar content moves between the orders and the interpolation
allowance without being charged twice.

### One low mean and one zero channel force the floor

If `delta_j=0`, both primitive means vanish. All primitive orders are
nonnegative, so (19) is `2log n`, which implies (15).

Otherwise one of the two means is at most `delta_j`, and the other
is at most `2delta_j`. For the smaller mean, Jensen's inequality gives

```math
\log\sum_i e^{-A\nu'_i}+\delta_j A\ge\log n.           \tag{20}
```

For the other profile, retain one zero order and apply Jensen to the
remaining `n-1` orders. Their sum is at most `2n delta_j`. With
`x=delta_j B >= 0`, this gives

```math
\log\sum_i e^{-B\mu'_i}+\delta_j B
\ge \log\left(1+(n-1)e^{-2nx/(n-1)}\right)+x
\ge f_n.                                              \tag{21}
```

The last minimum occurs at
`x=(n-1)log(n+1)/(2n)`. Swapping the two cuts if necessary makes
(20)--(21) apply in either order. Their sum proves (15).

For completeness, `f_n > (log n)/2` for every `n>1`. Multiplying
the difference by `2n` gives

```math
h(n)=2n\log2+n\log n-(n+1)\log(n+1).
```

Here `h(1)=0` and `h'(n)=log(4n/(n+1))>0` for `n>1`, proving
the strict second inequality in (16).

The one-profile Jensen expression in (21) is the same generic budget
function used in earlier restricted-family floors. The new accounting
here is the shared degree budget over all Hermite places, including
homogeneous scalar content. This is a limitation of the displayed
sufficient criterion, not a lower bound for arbitrary crossing
decompositions, actual jet ranks, factoring algorithms, weighted physical
types, or the opposite-place coefficient method.

## 7. Reproduction and search scope

Run:

```sh
python3 tools/check_multipoint_hermite_crossings.py \
  --output evidence/multipoint-hermite-crossings.json
```

The standard-library checker verifies complete reconstructed crossings
over characteristics two, three, and five, including infinity jets and
nontrivial Hermite coefficients. It checks `84,819` complete tensor
entries, `72` exact jet-rank inequalities, the ternary seed's full local
profiles by independent Toeplitz ranks, and both rational comparisons
in (14). It also checks the global degree ledger and the cancellation of
nontrivial scalar content at both finite places and infinity. These are
finite controls for the separate all-length proofs.

Exploratory tests also applied the criterion to rational fused `sl_2`
polynomial matrices and to higher-alphabet versions of (13). They did
not approach the post-half threshold. Those numerical searches are not
used in the theorem or in the exact certificate above.

Files: [checker](../tools/check_multipoint_hermite_crossings.py) and
[exact record](../evidence/multipoint-hermite-crossings.json).
