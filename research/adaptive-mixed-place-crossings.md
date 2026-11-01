# Adaptive opposite flags for crossing reconstruction

## Result and scope

The opposite-flag construction does not require a fixed proportion of
jets from each place. It can choose the cheapest next jet at every step,
using exact local channel counts. This produces one uniform crossing
family, and a finite list of rational interval inequalities can certify
its asymptotic rank bound.

The result strengthens the reconstruction interface in
[mixed-place-hermite-crossings.md](mixed-place-hermite-crossings.md).
It does not by itself supply a crossing parameter below `1/sqrt(2)`,
or improve the established factoring exponent `1/2`. Finite-power
rank costs and floating rate curves are exploratory until the uniform
criterion below has an exact certificate.

All arguments are ordinary mathematics. The exact implementation uses
integer arithmetic and finite-field linear algebra, without new Lean
formalization.

## 1. Seed, places, and exact jet costs

Let `F` be a fixed finite field. Let a fixed polynomial tensor
`K(z;x,y,u,v)` have physical alphabet size `q>=2`, degree at most `D`,
and the strict endpoint

```math
K(1;x,y,u,v)=[x=v][y=u].
\tag{1}
```

Choose nonempty disjoint rational place sets `A,B` on the projective
line over `F`, excluding the endpoint one. At finite `a` use
`K_a(t)=K(a+t)`, and at infinity use
`K_infinity(t)=t^D K(t^{-1})`. The complete homogeneous degree is
charged, including scalar content and unused degree at infinity.

At each place in `A`, retain the complete Smith orders of the first
cut `(x,y)|(u,v)`. At each place in `B`, retain those of the second
cut `(x,u)|(y,v)`. Every list has `n=q^2` nonnegative integers, since
both generic ranks are `q^2`. Write a local list as
`nu=(nu_1,...,nu_n)`.

For tensor length `m`, define

```math
N_\nu(m,j)=
\#\{(i_1,\ldots,i_m):\nu_{i_1}+\cdots+\nu_{i_m}\le j\},
\qquad
c_\nu(m,j)=(j+1)N_\nu(m,j).
\tag{2}
```

The coefficient-rank argument in the predecessor gives

```math
\mathrm{rank}[t^j]M(t)^{\otimes m}\le c_\nu(m,j).
\tag{3}
```

Indeed, a channel word of total order at most `j` contributes at most
`j+1` rank-one outer products to coefficient `j`. Words of larger
total order contribute zero. The cost in (2) is nondecreasing in `j`,
including its initial zero segment when the smallest local order is
positive.

These counts are computable in time polynomial in `m`. The coefficients
of

```math
\left(\sum_{i=1}^n X^{\nu_i}\right)^m
\tag{4}
```

give the exact order counts; cumulative sums give (2). The seed and
orders are fixed, the degree needed is `Dm=O(m)`, and every coefficient
is at most `n^m`, with `O(m)` bits. Ordinary dynamic programming suffices.
The implementation also uses exact carry-free integer packing of (4),
checked against independent word enumeration and convolution.

## 2. Cheapest-jet flags and their exact prefix bounds

Set `N=Dm`. For the row places, begin with the order-zero jet at every
place. Repeatedly select a currently available jet of smallest cost (2),
break ties by a fixed place order, and make the next jet at that place
available. Stop after `N+1` selections. Do the same independently at
the column places.

Every prefix contains an initial segment of jets at each place. Since
the individual sequences are nondecreasing, the merge selects the
globally cheapest jets. In particular, among every selection of `r`
jets that consists of initial segments, its prefix has the least sum
of the displayed costs. Ties do not affect that minimum.

Let `S_R(r)` be the sum of the first `r` row costs, capped by `n^m`,
and define `S_C(c)` similarly. Set `S_R(0)=S_C(0)=0`. The cap is the
full matrix dimension. Rank subadditivity gives the exact statements

```math
\begin{aligned}
\mathrm{rank}(L(K^{\otimes m}))&\le S_R(r)
 &&\text{if }L\text{ belongs to the first }r\text{ row jets' span},\\
\mathrm{rank}(L(K^{\otimes m}))&\le S_C(c)
 &&\text{in the second cut for the corresponding column span}.
\end{aligned}
\tag{5}
```

The two ranks in (5) use different flattenings. A zero prefix cost
means every tensor entry is killed by that entire functional span.

The flags remain opposite even though their order depends on `m` and
the exact costs. On the dual of polynomials of degree at most `N`,
write their prefixes as `R_r,C_c`. Disjoint Hasse-jet conditions give

```math
\dim(R_r\cap C_c)=\max(0,r+c-N-1).
\tag{6}
```

Vanishing at finite places imposes the stated divisibility conditions;
jets at infinity impose the complementary degree cap. The two sets
are disjoint, so their numbers of conditions add. This proves (6) in
every characteristic without ordinary derivative or factorial division.

Choose a nonzero functional

```math
L_k\in R_{k+1}\cap C_{N-k+1},\qquad 0\le k\le N.
\tag{7}
```

The intersection is one-dimensional, and `L_k` is not in `R_k`.
Thus the `L_k` form a dual basis. Solving one linear system expresses
evaluation at one as `sum_k a_k L_k`. The resulting pieces

```math
T_{m,k}=a_kL_k(K^{\otimes m})
\tag{8}
```

sum to the crossing and satisfy

```math
r(T_{m,k})c(T_{m,k})
\le S_R(k+1)S_C(Dm-k+1).
\tag{9}
```

This is an exact finite-power upper bound. It is not an assertion that
the displayed bounds equal the actual ranks.

## 3. A finite rational certificate for every tensor length

The following test proves a uniform asymptotic bound without inferring
one from finitely many values of (9).

Fix a rational grid

```math
0=x_0<x_1<\cdots<x_J=D.
\tag{10}
```

For each interval `j`, choose nonnegative rational row allocations
`a_(j,v)` for `v in A` and column allocations `b_(j,w)` for `w in B`
such that

```math
\sum_{v\in A}a_{j,v}\ge x_{j+1},\qquad
\sum_{w\in B}b_{j,w}\ge D-x_j.
\tag{11}
```

For each positive allocation choose a rational scale at least one,
denoted `s_(j,v)` or `t_(j,w)`. Put

```math
U_v(s)=\sum_i s^{-\nu_{v,i}},\qquad
V_w(t)=\sum_i t^{-\mu_{w,i}}.
```

**Adaptive interval criterion.** Suppose a fixed rational
`1/2<=rho<1` satisfies, for every interval and every pair of positive
row and column allocations,

```math
\boxed{
 U_v(s_{j,v})V_w(t_{j,w})
 s_{j,v}^{a_{j,v}}t_{j,w}^{b_{j,w}}<q^{4\rho}.
}
\tag{12}
```

Then the exact greedy flags supply effective scale-one crossing pieces
with parameter `rho` at every binary label width. There are `O(lambda)`
pieces, and their full tables and factorizations have construction bit
cost `2^O(lambda)`.

The allocations at different intervals **need not be nested**. They
are comparison sets used to bound the same exact greedy algorithm;
they are not prescribed prefixes of a second flag.

### Proof of the criterion

Consider a piece with `x_j<=k/m<=x_(j+1)`. At row place `v` with
positive allocation, the first

```math
\lceil m a_{j,v}\rceil+1
```

jets form a comparison set. Their union contains at least `k+1` jets
by (11). A zero allocation contributes no jets. The total allocation
is positive because `x_(j+1)>0`, so at least one place is present.
Every selected order is at most `m a_(j,v)+1`.

For every `s>=1`, the elementary generating-function bound is

```math
N_\nu(m,l)\le s^l U(s)^m.
\tag{13}
```

Therefore the sum of all costs in this comparison set is at most a
fixed polynomial in `m` times

```math
\left[\max_{v:a_{j,v}>0}
 U_v(s_{j,v})s_{j,v}^{a_{j,v}}\right]^m.
\tag{14}
```

The constants absorb the one extra jet, fixed rational scales, and
number of places. If a comparison list at one place has more than
`Dm+1` jets, truncate it there. Its union still contains the required
`k+1<=Dm+1` jets, and its cost can only decrease. All retained jets
are therefore defined on the correct polynomial space. The exact
greedy prefix is no more expensive than this comparison set.
The identical argument for column allocations gives a bound for
the first `Dm-k+1` column jets, using `D-x_j` in (11).

Multiplying the two bounds proves (9) is at most `poly(m) times C^m`,
where `C` is the largest left side of (12) over the finite certificate.
Its strict margin gives `C<q^(4rho)`. Taking binary logarithms supplies
the required per-piece sum of factor widths, with `O(log m)` overhead.

For any `lambda`, choose the least `m` with `q^m>=2^lambda`, then
restrict all legs to the first `2^lambda` words in a fixed order.
This preserves the crossing and cannot increase either rank. Exact
rank factorization caps either width at `2lambda`. The strict margin
absorbs rounding and polynomial factors. A fixed coordinate
decomposition covers the finitely many smaller widths if needed.

The exact merge, flag linear systems, and basis reconstruction have
polynomial bit cost in `m`. Writing every entry of the tensor power
and applying its `O(m)` functionals costs `exp(O(m))`; exact finite-field
elimination to factor each piece also costs `exp(O(m))`. The field and
seed are fixed, and `m=O(lambda)`, proving the claimed effectiveness.

If `F` is a degree-`d` extension of its prime field, choose an explicit
prime-field basis `e_1=1,e_2,...,e_d` and the coefficient functional
`ell` with `ell(1)=1`. Apply `ell` entrywise to every reconstructed
piece. Their sum is still the zero-or-one crossing. A rank-one product
`a b` satisfies

```math
\ell(ab)=\sum_{i=1}^d a_i\ell(e_i b),
\qquad a=\sum_{i=1}^d a_i e_i.
```

Thus either cut rank increases by at most the fixed factor `d`.
This adds only constant factor widths and construction overhead.
Using a coordinate functional avoids any assumption that the field
trace of one is nonzero. The existing
[scale-one transfer](scale-one-modular-transfer.md) then applies over
the prime field with all of its stated hypotheses.

To verify (12) exactly, choose one positive integer `L` clearing the
denominators of `4rho` and all allocations. Raise (12) to `L`. Every
quantity is then rational, so integer arithmetic after clearing
denominators decides the strict comparison. No numerical logarithm
is part of the certificate.

## 4. Scalar content gives free jets and cancels its degree

Suppose `Q(z)` is a nonzero scalar polynomial with `Q(1)=1`, considered
as a homogeneous section of degree `e>=deg Q`. Assume all of its
homogeneous zeros lie at the selected places `A union B`. Let `g_v`
be its order at a selected place, including `e-deg Q` at infinity.
Then `sum_v g_v=e`.

Multiplication of the seed by `Q` adds `g_v` to every local Smith
order. At tensor length `m`, it consequently adds `m g_v` initial
zero-cost jets at that place. The remaining exact costs satisfy

```math
\begin{aligned}
c'_v(m,m g_v+j)
 &=(m g_v+j+1)N_v(m,j),\\
c_v(m,j)&\le c'_v(m,m g_v+j)
       \le(m g_v+1)c_v(m,j).
\end{aligned}
\tag{15}
```

Let `g_A=sum_(v in A)g_v` and `g_B=sum_(v in B)g_v`. After the
`m g_A` additional free row jets and `m g_B` additional free column
jets are deleted, cheapest-prefix costs differ from the original
ones by at most fixed polynomial factors in `m`. This follows from
the optimal-prefix comparison in Section 2 applied in both directions
to the pointwise inequalities (15). Capping by the unchanged full
matrix dimension preserves the same inequalities.

The new degree is `D+e`. With `k'=k+m g_A`, the nontrivial
complementary prefix lengths become

```math
\begin{aligned}
k'+1-m g_A&=k+1,\\
(D+e)m-k'+1-m g_B&=Dm-k+1.
\end{aligned}
\tag{16}
```

Indices outside `0<=k<=Dm` have an entirely free prefix in one cut,
so their pieces vanish. Thus scalar content preserves the exponential
rate of this exact adaptive certificate. Its added degree has been
charged and canceled by actual zero jets; it has not been dropped
from the reconstruction.

This statement includes surplus degree at infinity. It assumes all
scalar zeros are available in the selected rational place sets; no
unaccounted zero over an extension field is permitted.

## 5. Rate curves and physical types

For exploration, define a local lower-tail entropy bound

```math
h_v(x)=\inf_{s\ge1}\{\log U_v(s)+x\log s\}.
\tag{17}
```

The adaptive row-prefix envelope can be estimated by distributing a
prefix budget among places and taking the maximum active local rate.
The column envelope is treated independently, then evaluated at the
complementary budget `D-x`. This is the motivation for the numerical
capacity routines. They are neither outward-rounded certificates nor
proofs of a new exponent. Initial zero jets and the first nonzero jet
need separate care at capacity boundaries.

Independent paired physical types and local leg gauges can be included
as in [the weighted mixed-place criterion](mixed-place-hermite-crossings.md#5-independent-paired-types-and-local-physical-gauges).
For fixed rational types, add each place's fixed mean-height shift to
its allocation exponent in (12). The right side becomes `exp(2rho H)`,
where `H` is the sum of the two type entropies. One common pair of
types must serve every place and every interval. Admissible tensor
lengths clear the fixed type denominators; rational heights can first
be cleared by fixed local ramification. Pair-swapping and retaining
the two coefficient indices separately gives equal external alphabets.

The unweighted exact implementation and checks below do not assert a
weighted improving witness. Those require their own complete rational
certificate, including every physical shift.

## 6. A complete certificate for a fused four-letter seed

The [complete fused-place catalogue](fused-all-place-profiles.md)
constructs the seed with `m=3`, encoded `theta=4`, physical alphabet
`q=4`, and degree `D=3`. Its singular places split in a fixed extension
of degree four over `F16`. Use the six recorded first-cut places as
row places and the nine disjoint second-cut places, including infinity,
as column places. The nonlinear closed points supply their distinct
conjugate root labels in this one fixed field.

The certificate
[fused-q4-adaptive-certificate.json](../evidence/fused-q4-adaptive-certificate.json)
gives a grid of 256 intervals, complete rational allocations and scales,
and

```math
\rho=\frac{31}{32}.
\tag{18}
```

Every interval satisfies (11)--(12) strictly by exact integer-power
comparisons. There are 2,832 local powered bounds. Comparing the row
and column maxima certifies every active pair in (12). The independent
checker also matches every selected profile and the seed parameters
to the complete catalogue, and verifies its recorded source hash.

Consequently this is an effective crossing family for **every** label
width with parameter `31/32`. The construction first uses `F_(2^16)`;
the descent above yields a family over `F2`, with at most a factor of
16 in either rank before the full-dimension cap. The uniform statement
follows from the rational interval proof, not the number of tested
tensor lengths.

The exploratory adaptive rate is approximately `0.966272`. That number
is not a certified optimum. The exact value in (18) remains above
`1/sqrt(2)`, so this seed does not improve the established factoring
exponent. No separation from an optimal fixed-allocation method is
asserted.

## 7. Exact checks and limits

Run

```sh
python3 tools/check_adaptive_mixed_place.py \
  --certificate evidence/fused-q4-adaptive-certificate.json \
  --output evidence/adaptive-mixed-place.json
```

The checker compares exact packed counts with independent word
enumeration and convolution; exhaustively compares greedy prefixes
with all small initial-segment allocations; constructs opposite
functional bases; and reconstructs complete original-tensor crossings.
It checks actual cut ranks against (9), scalar-content zero jets,
surplus degree at infinity, and a rational interval certificate whose
comparison allocations are explicitly nonnested.

The record is [adaptive-mixed-place.json](../evidence/adaptive-mixed-place.json).
These finite controls verify the algorithmic transcription. The uniform
crossing theorem is supplied by (10)--(14), not an extrapolation of the
finite rank results.

The seed, local Smith estimate, and original opposite-flag transfer are
prior work recorded in the linked notes. The adaptive merge, uniform
interval criterion, and scalar-content argument were developed with
OpenAI assistant support.
