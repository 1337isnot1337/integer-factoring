# Walsh twirling and principal Cauchy permanents

## Status

This note gives an exact diagonalization of the Boolean translation
average, all-size valuation formulas for principal Cauchy compounds,
and finite scalar decoders on a restricted subspace.

The permanent and determinant have different valuation behavior.
For every contiguous subset, their valuations agree and have a simple
dyadic formula. A circulant comparison proves the permanent valuation
for every evenly spaced subset of power-of-two size and even spacing.
The [arithmetic Cauchy theorem](arithmetic-cauchy-factorial-law.md) now
extends that conclusion to every subset size and proves the individual
squared-Cauchy Smith orders previously observed here.
On the even positions it differs from the determinant valuation by
a quadratic amount. This gives a concrete family of scalar modes
outside the exterior-algebra calculation.

There are also two limits. The translation average has an exact
`2^L` denominator at every nonzero local center, and the source's
grouping cut splits every row/column port pair. The whole principal
permanent span admits a high-height eight-wire collision hiding
112 bits. The restricted interval span has much better finite
conditioning: explicit dual certificates recover parity at precision
`a+3` for `L=2^a<=64`. The same six measurements do not suffice at
`L=128`. No improved general contraction or factoring exponent is
claimed.

## 1. The translation average is diagonal in the Walsh basis

For an arbitrary physical matrix `W`, let

\[
 K(W)_{Y,X}=\prod_{x=0}^{L-1}l_{Y_x}((WX)_x),
 \qquad l_0(v)=1-v,\quad l_1(v)=v.
\]

Let `P_g` translate Boolean states by bitwise addition of `g`, and put

\[
 \mathcal T(W)=2^{-L}\sum_{g\in\{0,1\}^L}P_gK(W)P_g.
 \tag{1}
\]

For a subset `S`, write `chi_S(X)=(-1)^(|S intersect X|)`. The Walsh
eigenvalue of (1) on `chi_S` is

\[
 \boxed{\quad \lambda_S(W)=\mathrm{per}(W_{S,S}).\quad}
 \tag{2}
\]

The empty permanent is one. This identity does not require a Cauchy
matrix.

**Proof.** Translation averaging removes the off-diagonal Walsh
coefficients. The diagonal coefficient before normalization is

\[
 \sum_{X,Y}\chi_S(X)\chi_S(Y)K(W)_{Y,X}
   =\sum_X\chi_S(X)\prod_{x\in S}(1-2(WX)_x).
\]

The product has Boolean degree at most `|S|`. Summing against
`chi_S` annihilates all its monomials except `X_S`. That coefficient
is `(-2)^(|S|) per(W_(S,S))`. The remaining character sum is
`(-1)^(|S|)2^(L-|S|)`, giving `2^L lambda_S`. Division by `2^L`
proves (2). Equivalently,

\[
 2^L\mathcal T(W)_{Y,X}
    =\sum_S\lambda_S(W)\chi_S(X\mathbin{\oplus}Y).
 \tag{3}
\]

For the Cauchy packet, `W(0)=I`, so `mathcal T(W(0))=I`.
At a nonzero grid point `t`, `W(t)` is a strictly triangular partial
shift. Every nonempty principal permanent is zero. Hence

\[
 \mathcal T(W(t))=2^{-L}J,
 \tag{4}
\]

where `J` is the all-ones matrix.

On the whole local disk `T=t+2Lz`, reduction modulo two is the same
partial shift, as a polynomial identity. Thus every nonempty
`lambda_S(t+2Lz)` is an even polynomial. The right side of (3) is
odd in every entry at every integral local evaluation point. The
denominator `2^L` in (1) is therefore exact; it cannot be discarded
as a formal averaging convention.

### 1.1 The source grouping cut splits the port pairs

For coordinate `i=v_2(t)` with unequal packet symbols, the source
assigns the row and column port of the same position to opposite
symbol groups. Position bit `i` specifies which orientation is used.
This is the packet convention in the source's
[balanced packet construction][source].

Equation (3) depends only on `X xor Y`. Swapping the row and column
bit at any chosen position therefore leaves this tensor unchanged.
For each position, name the port in the first group `U_x` and the
port in the second group `V_x`. Its grouped matrix is exactly

\[
 Q_{U,V}=\sum_S\lambda_S(W)\chi_S(U\mathbin{\oplus}V),
 \qquad Q=2^L\mathcal T(W),
\]

up to index permutations. Thus this cut has the same Smith orders
as the full row/column matrix `Q`. A cut that keeps both ports of
each position together has a different geometry and cannot be used
for this source step.

Over the fraction field, the grouped rank is the number of nonzero
principal permanents. In particular the polynomial matrix has full
rank `2^L`, since all its eigenvalues are one at `T=0`. Modulo two,
however, `Q` has rank one at every local center. Its higher Smith
orders, and the required division by `2^L` per full packet, both
matter for a finite-precision contraction argument. Full generic
rank alone does not settle that question.

### 1.2 Exponentially many modes survive the existing precision schedule

Choose one element from each antipodal pair
`{r,r+L/2}`, `0<=r<L/2`, and let `S` be the resulting subset.
There are `2^(L/2)` such subsets. For every modulus `2^j<=L/2`,
their residue counts are constant. At the remaining levels
`j=a,a+1` the counts are zero or one. Consequently all their
discrepancies in (5) are balanced, at every retained center.
Equations (6)--(8) give

\[
 v_2(\lambda_S(t+2L\zeta))=A_t(S)\le L
 \quad\text{whenever }\zeta\text{ is a unit}.
\]

Let `H` be the unnormalized Walsh matrix, so
`H^2=2^L I`. From (3),

\[
 H Q H=2^{2L}\mathrm{diag}_S(\lambda_S).
\]

Thus at least `2^(L/2)` Smith orders of `Q` are at most `3L`.
One can see this without diagonalizing `Q`: multiplication by
integral matrices cannot lower a Smith order, while the displayed
diagonal matrix has that many entries of valuation at most `3L`.

For `P` full packets on a fixed grouping cut, the same argument
uses the tensor Walsh matrix. It proves that their unnormalized
tensor product has at least

\[
 2^{LP/2}
\]

surviving Smith directions modulo `2^(3LP+1)`. In particular any
integral rank-one grouped decomposition at that precision requires
at least that many summands. The normalization of the corresponding
twirled product costs an additional `LP` bits.

This already lies below the source's existing Hermite precision
schedule for full packets: if their `LP` positions contribute that
many real edges, its `8M(a+1)` term alone exceeds `3LP`.
Therefore direct replacement by twirled forward packets, together
with general grouped rank decompositions at the existing precision,
does not yield a subexponential improvement. This is a whole-operator
statement for full packets. It leaves open lower-precision scalar
recovery, restrictions imposed by surrounding tensors, and reductions
that treat padded unused ports separately.

The [group-profile replay](../evidence/walsh-twirl-group-profile.json)
checks the literal source-oriented reshuffling and records finite
Smith profiles independently of this lower bound.

## 2. Two dyadic discrepancy quantities

Take `L=2^a`, a nonzero retained center `|t|<L`, and a unit local
point `zeta`. Use the Cauchy matrix from
[the Boolean Smith note](boolean-cauchy-smith-profile.md):

\[
 W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},\qquad T=t+2L\zeta.
\]

For a subset `S`, let

\[
 n_j(r)=|\{x\in S:x\equiv r\pmod{2^j}\}|,
 \qquad
 \delta_j(r)=n_j(r)-n_j(r-t).
\]

Define its linear and quadratic discrepancies by

\[
 A_t(S)=\frac12\sum_{j=1}^{a+1}\sum_{r\bmod2^j}|\delta_j(r)|,
 \qquad
 D_t(S)=\frac12\sum_{j=1}^{a+1}\sum_{r\bmod2^j}\delta_j(r)^2.
 \tag{5}
\]

Then

\[
 \boxed{\quad
  v_2(\mathrm{per}\,W_{S,S})\ge A_t(S),\qquad
  v_2(\det W_{S,S})=D_t(S).
 \quad}
 \tag{6}
\]

The empty subset has both quantities zero.

### 2.1 The assignment bound is exact before summing permutations

The common Cauchy numerator has valuation `a+1`. Its denominator
has valuation `min(v_2(t+x-y),a+1)`, with the minimum interpreted as
`a+1` when `t+x-y=0`. Thus an entry has valuation

\[
 a+1-\min(v_2(t+x-y),a+1)
   =\sum_{j=1}^{a+1}\mathbf1[x+t\not\equiv y\pmod{2^j}].
 \tag{7}
\]

At residue level `j`, any matching between `S+t` and `S` leaves
at least
`|S|-sum_r min(n_j(r-t),n_j(r))` pairs outside common residue
classes. This number is `(1/2)sum_r |delta_j(r)|`.

All these lower bounds can be attained simultaneously. Match within
the finest residue classes first, then match their unmatched elements
inside each parent class, and continue up the binary residue tree.
This produces the minimum-weight assignment and proves that its
weight is exactly `A_t(S)`. Each term of the permanent has at least
that valuation. Cancellation can raise the permanent's valuation but
cannot lower it, proving the first statement in (6).

### 2.2 The determinant has quadratic discrepancy

Cauchy's determinant formula gives

\[
 v_2(\det W_{S,S})=
 |S|(a+1)+2\sum_{\substack{x<y\\x,y\in S}}v_2(y-x)
   -\sum_{x,y\in S}\min(v_2(t+x-y),a+1).
\]

At each residue level, the first two terms count
`sum_r n_j(r)^2` and the last term counts
`sum_r n_j(r)n_j(r-t)`. Translation preserves the sum of squares,
so their difference is `(1/2)sum_r delta_j(r)^2`. Summing the
levels proves the second statement in (6).

This proof works over the unramified coefficient rings used by the
source, since the local point is a unit and the integer differences
have the displayed valuations.

## 3. Exact balanced cases, including every interval

Suppose all the discrepancies `delta_j(r)` belong to `{-1,0,1}`.
Then `A_t(S)=D_t(S)`. Every term of the determinant and permanent
is divisible by `2^(A_t(S))`. After that division their reductions
modulo two agree, because the signs disappear. The determinant is
a unit after division, by (6). Therefore

\[
 \boxed{\quad
 v_2(\mathrm{per}\,W_{S,S})=A_t(S)=D_t(S)
 \quad}
 \tag{8}
\]

in every such balanced case.

If `S` is a contiguous interval, its counts in each residue system
are a common integer plus zero or one. Translating those counts
leaves their differences in `{-1,0,1}`. Thus (8) holds for every
contiguous interval and every nonzero retained center.

In particular, at `t=1`, an interval of length `k>0` has

\[
 \boxed{\quad
 v_2(\lambda_{[c,c+k)}(1+2L\zeta))
       =a+1-v_2(k).
 \quad}
 \tag{9}
\]

At a residue level the two interval endpoints contribute zero when
`2^j` divides `k` and contribute one otherwise. This proves (9)
directly from (5). The location `c` does not matter: translation of
both principal index sets leaves the Cauchy submatrix unchanged.

The valuation statement concerns unit local points. It is also the
content of the full local polynomial in this case. Each physical
entry polynomial has content equal to its entry weight in (7).
Every permanent term is therefore coefficientwise divisible by
`2^(A_t(S))`. Equality at a unit point proves that its content
cannot be larger.

### A caution about runs

At `t=1`, the number of consecutive runs in `S` is a lower bound
on the ordinary order of vanishing in `z`. It need not be the exact
order. Two runs of lengths two and three separated by one missing
position have leading reduced permanent

\[
 \frac1{2\cdot3}-\frac1{1\cdot6}=0.
\]

Consequently the leading two-break coefficient cancels. The
assignment bound (6) remains valid and does not assume that a
minimum-weight permutation survives summation.

## 4. Evenly spaced sets distinguish permanents from determinants

Let `S={0,2,...,L-2}` and `n=L/2`. At any odd center the
discrepancy quantities are

\[
 A_t(S)=(a+1)n,\qquad D_t(S)=2n^2.
 \tag{10}
\]

For the quadratic quantity, the contributions at levels `1,...,a`
are `L^2/2^(j+1)` and the final level contributes `n`. The linear
contribution is `n` at every level.

The permanent reduces to the Cauchy block

\[
 Q_{ij}(T)=\frac1{T+2i-2j},\qquad 0\le i,j<n.
\]

### 4.1 A circulant permanent congruence

Let `n=2^b`. If an `n` by `n` circulant matrix `A` is congruent
to `J` modulo two over `mathcal O`, then

\[
 \boxed{\quad \mathrm{per}(A)\equiv n!\pmod{2^n}.\quad}
 \tag{11}
\]

Write `A=J+2B`. Expanding the permanent gives

\[
 \mathrm{per}(J+2B)
   =\sum_{k=0}^n2^k(n-k)!\,R_k(B),
\]

where `R_k(B)` is the weighted sum of all `k`-edge partial
matchings. Simultaneously translating the row and column indices
acts on those matchings and preserves every weight. For `0<k<n`,
the stabilizer size divides `k`: a selected row set is a union of
orbits of the stabilizer. Hence each matching orbit has size
divisible by `2^(b-v_2(k))`. It follows that
`v_2(R_k(B))>=b-v_2(k)`.

The corresponding summand has valuation at least

\[
 k+v_2((n-k)!)+b-v_2(k)
   =n-s_2(n-k)+b-v_2(k)\ge n.
\]

The last inequality uses
`s_2(n-k)=b-s_2(k-1)<=b-v_2(k)`. The `k=n` term is also
divisible by `2^n`, while the `k=0` term is `n!`. This proves
(11), including `n=1`.

### 4.2 Even-spacing Cauchy blocks of power-of-two size

Let `d` be even and `T` be a unit. For `0<=i,j<n=2^b`, put

\[
 C_{ij}=\frac1{T+d(i-j)}.
\]

Then

\[
 \boxed{\quad
  \mathrm{per}(C)\equiv T^{-n}n!\pmod{2^n},
  \qquad v_2(\mathrm{per}(C))=n-1.
 \quad}
 \tag{12}
\]

To prove this, set `A=TC` and compare it with the circulant matrix
`A^0_ij=T/(T+d((i-j) mod n))`. Their entrywise difference is
divisible by `2^(b+1)`, because the two denominators differ by a
multiple of `dn`. All entries of both matrices are odd.

The cofactor bound needed for row replacement is specifically for
size `m=n-1=2^b-1`. Every odd `m` by `m` matrix `J+2B` has
permanent divisible by `2^(v_2(m!))`: its expansion terms have
valuation at least

\[
 k+v_2((m-k)!)=m-s_2(m-k)\ge m-b=v_2(m!).
\]

This uses the bound `s_2(m-k)<=b`. It is not an assertion that an
arbitrary odd matrix of every size has factorial divisibility.

Replace the rows of `A^0` by those of `A` one at a time.
Each permanent difference is divisible by

\[
 2^{b+1+v_2((n-1)!)}=2^n.
\]

Equation (11) now gives `per(A)=n! mod 2^n`. Scaling back by
`T^n` proves (12).

### 4.3 Exact local content of evenly spaced packet modes

The same proof applies over `mathcal O[[z]]`. Let
`S={c,c+d,...,c+(n-1)d}` be contained in the packet, with `d`
positive and even and `n` a power of two. At every odd retained
center `t`,

\[
 \boxed{\quad
 \lambda_S(t+2Lz)
  \equiv 2^{\,n(a+2)-1}z^n
       \pmod{2^{\,n(a+2)}}.
 \quad}
 \tag{13}
\]

Indeed the common Cauchy numerator divided by `C'(0)` is
`2^(a+1)z u(z)` with `u(z)=1 mod 2`. Apply (12) to the
remaining block, whose `T=t+2Lz` is a formal power-series unit.
After dividing by `2^(n-1)`, its permanent is one modulo two.
This proves the coefficientwise congruence (13), and therefore
exact content `n(a+2)-1`. At every unit local point its valuation
has that same value.

For the even-position subset of (10), this gives the all-size
comparison

\[
 v_2(\lambda_S)=n(a+2)-1,\qquad
 v_2(\det W_{S,S})=2n^2,\qquad n=L/2.
\]

The resulting exact comparisons include:

| Packet length `L` | Subset size `n` | Permanent valuation | Determinant valuation |
|---:|---:|---:|---:|
| 4 | 2 | 7 | 8 |
| 8 | 4 | 19 | 32 |
| 16 | 8 | 47 | 128 |
| 32 | 16 | 111 | 512 |
| 64 | 32 | 255 | 2,048 |
| 128 | 64 | 575 | 8,192 |
| 256 | 128 | 1,279 | 32,768 |

There is a low-creation-precision scalar consequence for the twirled
operator. Take `A=delta_0` and the two binary columns
`B_+=(1+chi_S)/2` and `B_-=(1-chi_S)/2`. Their original
contractions are one and zero. Their twirled values differ by
`lambda_S`, so (13) gives an explicit binary collision at content
`n(a+2)-1` on every odd disk. The individual twirled values are
`(1+lambda_S)/2` and `(1-lambda_S)/2`; their half-integral
normalization must remain explicit. This does not assert a
collision for the untwirled source kernel.

### 4.4 The individual Smith profile is now proved at every size

The checker also evaluates the blocks by Borchardt's identity

\[
 \det C\,\mathrm{per}\,C=\det(C\circ C).
\]

Here the circle denotes entrywise multiplication. A primary proof
is [Singer (2004), Theorem 1.1][borchardt]. Small cases also use
a separate direct permanent computation. For consecutive indices,
the squared-Cauchy Smith orders are

\[
 2r+2v_2(r!)+v_2(r+1),\qquad 0\le r<n.
\]

The original tests through `n=128` gave
`v_2(per C)=v_2(n!)=n-s_2(n)`, including non-power-of-two sizes
through sixteen. Both statements now have ordinary all-size proofs in
[the arithmetic Cauchy theorem](arithmetic-cauchy-factorial-law.md).
Its sampling-moment lemma removes the size restriction from the permanent
formula; its integral orthogonal-polynomial basis proves the individual
Smith orders. The proof above remains a separate power-of-two argument.
Nonconsecutive index sets can have larger valuations.

Reversing columns makes `C circle C` a Hankel matrix
`1/(A+2i+2j)^2`. Two forward-difference transforms give entries

\[
 \frac{(-2)^{r+s}(r+s)!}
      {\prod_{h=0}^{r+s}(A+2h)}
       \sum_{h=0}^{r+s}\frac1{A+2h}.
\]

Every denominator is odd. This exposes factorial weights and
odd-progression harmonic sums. The new theorem supplies the additional
integral basis needed to prove the individual Smith formula.

## 5. Exact restricted scalar decoders

Put `f_k(T)=lambda_[0,k)(T)`, including `f_0=1`, and consider

\[
 F(T)=\sum_{k=0}^L c_k f_k(T),\qquad c_k\in\mathcal O.
 \tag{14}
\]

At zero, `F(0)=sum_k c_k`. The functions `f_k` are independent
over the rational coefficient field: their degrees
`k(2L-1)` are distinct and their leading coefficients are
\(k!/C'(0)^k\ne0\).

Let `F_(t,j)` denote the coefficient of `z^j` in `F(t+2Lz)`.
The [recorded exact duals](../evidence/interval-permanent-decoders.json)
recover `F(0) mod 2` from at most these six measurements:

\[
 F_{-1,0},\quad F_{-1,1},\quad F_{-1,2},\quad
 F_{-3,3},\quad F_{-1,3},\quad F_{1,1}.
 \tag{15}
\]

The two-wire case uses only centers `-1,+1`. Each positive record
contains integers `d_(t,j)` such that, for every `k=0,...,L`,

\[
 \sum_{t,j}d_{t,j}[z^j]f_k(t+2Lz)
       \equiv 2^{a+2}\pmod{2^{a+3}}.
 \tag{16}
\]

Linearity makes (16) an explicit decoder for every coefficient vector
in (14), with no bound on its height. Apply the same combination to
the measurements, divide by `2^(a+2)`, and reduce modulo two.

| `L` | Required measurement precision | Certified for the whole interval span |
|---:|---:|:---:|
| 2 | 4 bits | Yes |
| 4 | 5 bits | Yes |
| 8 | 6 bits | Yes |
| 16 | 7 bits | Yes |
| 32 | 8 bits | Yes |
| 64 | 9 bits | Yes |

At `L=128` the same six measurements have an explicit collision
modulo `2^10` with `sum_k c_k` odd. That negative record rules out
continuing these six measurements unchanged at the proposed
precision. It does not rule out a larger collection of jets or a
different restricted decoder.

### 5.1 Computing the required jets without a full Boolean matrix

The first three jets in (15) can be computed in polynomial time in
`L` and the precision. The following explicit formula also supplies
an independent replay of the large finite certificates.

At `t=1` put `u(z)=W_1(1+2Lz)` and `r_d=2L/d`. For a principal
interval of length `k`, every permutation has the last shift edge
broken. A term with `h` broken shift edges contains `z^h`.
Through degree three there are only one-, two-, and three-break
terms. Write

\[
 f_k(1+2Lz)=u(z)^k
       \bigl(p_{k,1}z+p_{k,2}z^2+p_{k,3}z^3\bigr)
         \pmod{z^4}.
\]

Their coefficients are

\[
 p_{k,1}=r_k,\qquad
 p_{k,2}=-r_k^2+\sum_{\ell=1}^{k-1}r_\ell r_{k-\ell},
 \tag{17}
\]

\[
\begin{aligned}
 p_{k,3}={}&r_k^3
 -\sum_{\ell=1}^{k-1}
      (r_\ell^2r_{k-\ell}+r_\ell r_{k-\ell}^2)\\
 &+\sum_{\substack{\ell,m,n>0\\\ell+m+n=k}}
 \left(r_\ell r_m r_n
       -r_m r_{\ell+m}r_{m+n}-r_m^2r_k\right).
\end{aligned}
 \tag{18}
\]

For three broken edges the residual three by three matrix has three
allowed permutations. They give the three terms in parentheses.
All other permutations have at least four broken edges and do not
affect these jets.

At `t=3`, for `k>=3`, the only contribution to its third jet
uses the three mandatory boundary rows. With
`w_3=C'(3)/C'(0)` this gives

\[
 [z^3]f_k(3+2Lz)
     =w_3^k\,\mathrm{per}
          \left(\frac{2L}{k+i-j}\right)_{0\le i,j<3}.
 \tag{19}
\]

The cases `k=1,2` are computed directly. Negative centers follow
from the exact reflection identity

\[
 f_k(-T)=\left(\frac{L-T}{L+T}\right)^k f_k(T).
 \tag{20}
\]

The implementation checks (17)--(20) against complete local
polynomials at the small lengths. It then verifies (16) directly
for every `k` in each larger certificate. No floating-point rank
or numerical tolerance enters those checks.

### 5.2 What a general packet representation would still need

Contiguous subsets of the same size give the same polynomial, so
this family has only `L+1` scalar directions. General principal
permanents have more directions. A low-content mode by itself does
not make an arbitrary linear combination easy to reconstruct.

To use (14) in the factoring evaluator, one would need a packet
representation whose surrounding tensors force the principal
Fourier coefficients into this interval span, or an efficient
decomposition into such spans with controlled precision and group
interaction. No such representation is supplied here. Enumerating
all subsets and choosing a separate permutation that makes each one
contiguous would restore an exponential number of choices per packet.

One concrete sufficient condition is a row tensor with an integral
Walsh expansion `A(X)=sum_I a_I chi_I(X)` supported on intervals and
the empty set. For an arbitrary integral column tensor `B`, the
twirled scalar then has (14) with
`c_k=sum_(|I|=k) a_I B_hat(I)`. There is no residual `2^L` scalar
denominator, and `F(0)=sum_X A(X)B(X)`. A large enough nonnegative
constant term makes bounded integer choices of `a_I` into a
canonical nonnegative row tensor without changing this support.

This sufficient tensor class already has a short sum of product
representations, one per interval. It therefore does not by itself
create a faster general evaluator. A useful new representation would
need the scalar-span condition to follow from more general surrounding
tensors or from a class closed under the recursive construction.

## 6. An explicit high-height collision in the whole principal span

The [eight-wire witness](../evidence/walsh-permanent-witness.json)
contains 57 nonzero integral coefficients `c_S`, each in
`[0,2^112)`. It satisfies

\[
 \sum_S c_S\equiv1\pmod2,\qquad
 \sum_S c_S\lambda_S(t+16z)\equiv0\pmod{2^{112}}
 \tag{21}
\]

as complete polynomials at all eight odd centers. At least one
nonconstant coefficient has valuation exactly 112. The stored
vector is replayed directly; no claim that a finite Smith search
proved an optimal threshold is needed.

There is also an explicit realization by canonical packet tensors
for the twirled operator. Let `H` be the unnormalized Walsh matrix,
take `A=delta_0`, and set

\[
 B(X)=(Hc)_X\bmod2^{120},
 \qquad 0\le B(X)<2^{120}.
 \tag{22}
\]

Because `H^2=256I`, the coefficients of this scalar twirled packet
are `c'_S=(HB)_S/256`. They are integers and satisfy
`c'_S=c_S mod2^112`. Thus (21) remains true. The original
contraction `B(0)` is odd. The saved tensor has creation precision
120 bits, and those representatives are fixed.

The canonical lift cannot cancel every valuation-112 nonconstant
coefficient: its coefficient changes are multiples of `2^112`,
and every nonempty `lambda_S` is an even local polynomial.
The empty mode affects only the constant coefficient.

This is a high-height collision. It is not evidence against an
`O(L log L)` loss for tensors created at binary precision, and it
does not transfer directly to the untwirled scalar packet. It does
show that diagonalizing the operator alone does not remove the
general scalar reconstruction problem.

## 7. Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/check_walsh_cauchy_permanents.py \
  --output /tmp/walsh-principal-permanents.json
```

The checker uses exact integer residues, rational constants, direct
permanents, minimum-weight assignments, and finite Smith elimination.
It independently compares the local-root and global-polynomial
kernel constructions. It replays every positive dual and the
limited-query collision, then checks all complete polynomials of
the 112-bit witness.

The finite squared-Cauchy profile observations are recorded
separately from the proved discrepancy identities and the explicit
decoder certificates in
[walsh-principal-permanents.json](../evidence/walsh-principal-permanents.json).

[borchardt]: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v11i1r48/pdf
[source]: https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex
