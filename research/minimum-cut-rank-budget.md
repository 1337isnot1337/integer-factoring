# A minimum-cut rank budget and a twenty-two-label cutoff

## Results and scope

For every field and every exact decomposition of the `d`-label crossing,

```math
E_d(x,y,u,v)=[x=v][y=u]=\sum_i A_i(x,y,u,v),
```

write `r_i,c_i` for the ranks in the cuts `(x,y)|(u,v)` and
`(x,u)|(y,v)`. In addition to the usual two rank budgets, there is
a third universal constraint:

```math
\boxed{\qquad \sum_i\min(r_i,c_i)\ge d.\qquad}          \tag{1}
```

Combining (1) with the
[all-field swap/product bound](swap-product-rank-sum.md) gives an
exact finite certificate that

```math
\boxed{\qquad
\sum_i\sqrt{r_ic_i}>d^{\sqrt2}
\quad\text{for every }2\le d\le22.
\qquad}                                                \tag{2}
```

Thus a finite seed passing the strict sufficient criterion in
[rank-sum amplification](finite-seed-amplification.md) needs at
least **twenty-three labels**. This strengthens the separately verified
twelve-label cutoff. It does not assert that a winning twenty-three-label
seed exists, and it does not exclude other amplification methods or
better asymptotic crossing families.

## 1. A quotient bound for a sum of product operators

Let `V` have dimension `d`, let `S` be the swap on `V tensor V`,
and let

```math
B=\sum_{j=1}^k P_j\otimes Q_j.
```

For any `s`-dimensional subspace `U` of `V`, put
`L=sum_j P_j U`. Then `dim L<=ks`. Restrict the input to
`U tensor V` and quotient the first output factor by `L`.
Every product term vanishes under this quotient. Swap maps
`U tensor V` onto `V tensor U`, whose quotient has dimension
`s(d-dim L)`. Therefore

```math
\operatorname{rank}(S-B)
\ge s(d-ks),\qquad 0\le s\le d.                        \tag{3}
```

The right side may be negative; taking its maximum with the `s=0`
value gives the useful bound. When a four-tensor `B` has column-cut
rank `k`, its row matrix has a representation by `k` product
operators. Hence

```math
r(E_d-B)\ge\max_{0\le s\le d}s(d-k s),
\qquad k=c(B).                                         \tag{4}
```

The statement with the cuts exchanged also holds. This is an
all-field argument and assumes neither positivity nor invertibility.
For `k=1`, the sharper Fitting bound in the linked note is used
below. Equation (4) also handles higher cut ranks; for example,
at `d=13,k=2` it gives complement rank at least 21.

## 2. Proof of the minimum-cut budget

Discard zero pieces. Let `A` be the sum of pieces with `r_i<=c_i`,
and let `B` be the sum of all remaining pieces. Since `E_d=A+B`,
equation (4) with `s=1` gives

```math
r(A)\ge d-c(B).
```

Rank subadditivity now yields

```math
d\le r(A)+c(B)
\le\sum_{r_i\le c_i}r_i+\sum_{r_i>c_i}c_i
=\sum_i\min(r_i,c_i),
```

proving (1). Empty sums cause no problem: a zero matrix has rank
zero, and (3) remains valid when `k=0`.

## 3. The rank-one complement budgets

Put `N=d^2`. Let `a` be the largest row rank among pieces of
column rank one, or zero if there is no such piece. Define `b`
with the cuts exchanged. As proved in
[swap-product-rank-sum.md](swap-product-rank-sum.md),

```math
a,b\in\mathcal P_d:=\{0\}\cup\{uv:1\le u,v\le d\}.
```

For positive `a`, define

```math
h_d(a)=\max_{\substack{uv=a\\1\le u,v\le d}}\min(u,v),
\qquad
D_d(a)=N-\frac{h_d(a)(h_d(a)+1)}2.
```

The actual maximum piece `(a,1)` forces the row total to be at
least `a+D_d(a)`. Thus the valid budgets are

```math
R=\begin{cases}N,&a=0,\\\max(N,a+D_d(a)),&a>0,\end{cases}
\qquad
C=\begin{cases}N,&b=0,\\\max(N,b+D_d(b)),&b>0.\end{cases}        \tag{5}
```

They are lower bounds on the actual totals, not claims that every
profile with these totals is realizable. The profile `(a,b)=(0,1)`
is impossible because a piece `(1,1)` would make `a>=1`.

## 4. A three-variable dual certificate

Choose `u,v,w>=0` satisfying

```math
\begin{aligned}
Nu+2v+2w&\le\sqrt{2N},\\
2u+Nv+2w&\le\sqrt{2N},\\
u+v+w&\le1,\\
au+v+w&\le\sqrt a &&(a>0),\\
u+bv+w&\le\sqrt b &&(b>0).
\end{aligned}                                                \tag{6}
```

Then every permitted rank pair obeys

```math
ur+vc+w\min(r,c)\le\sqrt{rc}.                         \tag{7}
```

To check all ranks at least two, consider the triangle
`2<=r<=c<=N`. There the left side minus the right side is
`(u+w)r+vc-sqrt(rc)`. For fixed `c` it is convex in `r`, so
its maximum is at `r=2` or `r=c`. Along `r=2` it is convex
in `c`; along the diagonal it is linear. The three vertices
`(2,2),(2,N),(N,N)` therefore suffice. The opposite triangle
is identical. The first two constraints in (6) handle the mixed
vertices, and `u+v+w<=1` handles both diagonal vertices.

On a possible rank-one interval `(r,1)`, `1<=r<=a`, the
minimum rank is constant. Convexity reduces the check to `(1,1)`
and `(a,1)`. The other interval is symmetric. This proves (7).

Summing (7) and using (1) and (5) gives

```math
\sum_i\sqrt{r_ic_i}\ge Ru+Cv+dw.                       \tag{8}
```

### Retain the cost of forced maximum pieces

If `a>0`, remove one actual piece of ranks `(a,1)` before
applying (7) to the rest. Its minimum rank is one, so the remaining
minimum-rank budget is `d-1`. This proves the additional bound

```math
\sqrt a+(R-a)u+(C-1)v+(d-1)w.                         \tag{9}
```

There is a symmetric expression for `b>0`. If `a,b>0` and
`(a,b)!=(1,1)`, the forced pieces `(a,1)` and `(1,b)` are
distinct. Removing both gives

```math
\sqrt a+\sqrt b+(R-a-1)u+(C-b-1)v+(d-2)w.             \tag{10}
```

No number-of-pieces hypothesis is used. In particular, a shared
`(1,1)` piece is never removed twice. The constraints in (6)
remain valid for all remaining pieces, even if one of their
actual maxima becomes smaller after removal.

### A nonlinear surplus bound

Suppose every rank-one-cut piece, if any, has ranks `(1,1)`.
This includes the profiles `(a,b)=(0,0)` and `(1,1)`.
Put `m_i=min(r_i,c_i)`, `t_i=|r_i-c_i|`, and `M=sum_i m_i`.
Every non-diagonal pair has `m_i>=2`. Its cost is

```math
m_i+g_{m_i}(t_i),\qquad
g_m(t)=\sqrt{m(m+t)}-m
      =\frac{t}{\sqrt{1+t/m}+1}.
```

For fixed `t>=0`, the function `g_m(t)` increases with `m`.
Also `g_2` is increasing and concave, with `g_2(0)=0`, so
`sum_j g_2(t_j)>=g_2(sum_j t_j)`. The diagonal `(1,1)`
pieces have zero surplus and cause no exception to this argument.

Separate the row-heavy and column-heavy pieces. If `M<N`, the
total surplus in each orientation is at least `N-M`, by the two
total-rank budgets. Consequently

```math
\sum_i\sqrt{r_ic_i}
\ge M+2g_2(N-M)
=M-4+2\sqrt{2(N-M+2)}.
```

The last expression is nondecreasing for `M<=N`: its derivative
is `1-sqrt(2)/sqrt(N-M+2)>=0`. Use `M>=d` from (1).
If `M>=N`, the simpler bound `sum sqrt(r_i c_i)>=M>=N`
is sufficient; the expression below is at most `N` because
`x+2-2sqrt(2x)=(sqrt(x)-sqrt(2))^2>=0` with `x=N-d+2`.
Thus in either case

```math
\boxed{\quad
\sum_i\sqrt{r_ic_i}\ge d-4+2\sqrt{2(N-d+2)}.
\quad}                                                   \tag{11}
```

This is a general nonlinear consequence of the rank budgets.
It is not a dimension-specific fit. At `d=21` it gives
`17+2sqrt(844)>75.103`, closing both profiles that the linear
dual alone leaves below the target.

## 5. Exact certification through twenty-one labels

For each `2<=d<=21`, enumerate every unordered pair `a,b` in
`P_d`, omitting `(0,1)`. Replace each square root in (6), (9),
and (10) by the rational lower bound

```math
\ell(z)=\frac{\operatorname{isqrt}(10^{14}z)}{10^7}.
```

Intersect triples of boundary planes from (6) and the three
coordinate planes. Keep only nonnegative points satisfying every
rounded constraint. These are feasible for the original constraints
as well. At each point evaluate the applicable bounds (8)--(10),
and retain the largest one. Optimality is unnecessary: each retained
value is already a proved lower bound for its profile. In the two
applicable profiles also use the downward-rounded bound (11).

The checker uses integer determinants with a common denominator,
then exact rational comparisons. For every profile it verifies

```math
L^{408}>d^{577},\qquad 577^2-2\cdot408^2=1.             \tag{12}
```

Consequently `L>d^(577/408)>d^sqrt(2)`, proving (2).
The final twenty-one-label lower bound is

```math
L_{21}=\frac{65335133917}{878000000}
=74.413592160\ldots .                                  \tag{13}
```

All square-root rounding is downward. Floating-point values appear
only as readable summaries of exact rational certificates.

## 6. The twenty-two-label cutoff

### The first linear profile failure is not a tensor certificate

At `d=22`, the no-rank-one **rank-budget relaxation** receives the rounded bound
`98657731/1250000` from (11), below the target. The nonlinear
rank-budget relaxation is attained by the numerical rank pairs

```math
(18,18),\qquad (2,464),\qquad (464,2).
```

Their row and column totals are both `484`, their minimum-rank
total is `22`, and their cost is `18+2sqrt(928)<22^sqrt(2)`.
But they **cannot be the ranks of a crossing decomposition**.
For the piece with row rank two, (4) with the cuts exchanged
requires its column complement to have rank at least

```math
\max_s s(22-2s)=60.
```

The other two pieces have total column rank only `18+2=20`,
contradicting rank subadditivity. Thus this profile identifies a
missing higher-cut constraint in the present finite certificate;
it is not a passing seed or evidence of realizability.  The quotient
constraint closes the diagonal and rank-one-free cases, and Section 6.3
closes the remaining rank-one cases.

### The `d=22` profile with only diagonal rank-one pieces

The quotient constraint above does not by itself give a `d=22` cutoff, but
it does eliminate every decomposition whose rank-one pieces are all `(1,1)`.
Put
`M=sum_i min(r_i,c_i)`.  The minimum-cut budget gives `M>=d=22`.

For `M>=23`, the nonlinear surplus bound in the form

```math
\sum_i\sqrt{r_ic_i}\ge M-4+2\sqrt{2(N-M+2)}          \tag{14}
```

is increasing in `M` for `M<=N`, and at `M=23` it gives

```math
19+2\sqrt{926}>22^{\sqrt2}.                           \tag{15}
```

For `M>=N` the weaker bound `sum_i sqrt(r_ic_i)>=M>=N` is already
above the target.  Thus only `M=22` can survive.

Assume now that every piece either has both cut ranks at least two or is
the diagonal rank-one piece `(1,1)`.  Put
`A=sum_{r_i>c_i}(r_i-c_i)`, `B=sum_{c_i>r_i}(c_i-r_i)`.  The two rank
budgets give `A,B>=N-M=462`; the diagonal rank-one pieces have zero
surplus and are already counted in `M`.  If one orientation has all unbalanced
pieces with minimum rank at least three, then concavity of
`g_m(t)=sqrt(m(m+t))-m` and `g_m(t)>=g_3(t)` give the lower bound

```math
22+g_3(462)+g_2(462)
=17+\sqrt{1395}+\sqrt{928}>22^{\sqrt2}.               \tag{16}
```

If one orientation has at least two unbalanced pieces, then the same
concavity argument, applied with one surplus unit separated, gives

```math
22+g_2(1)+g_2(461)+g_2(462)
=16+\sqrt{928}+\sqrt6+\sqrt{926}>22^{\sqrt2}.        \tag{17}
```

The exact inequalities (15)--(17) are checked by raising both sides to
the power `408` and using `577^2>2*408^2`; the square roots are rounded
downward.

The only remaining possibility is one unbalanced piece in each
orientation, both with minimum rank two.  The column-heavy piece then
has row rank `2` and column rank `2+t` for some `t>=462`.  The other
pieces are the row-heavy piece, whose column rank is `2`, and diagonal
pieces (possibly including `(1,1)`) whose minimum ranks sum to
`22-2-2=18`.  Hence the column
complement of the column-heavy piece has rank `20`.  The quotient bound
(4) for row rank two requires that complement to have rank at least

```math
\max_s s(22-2s)=60,                                    \tag{18}
```

a contradiction.  The surviving diagonal pieces may include `(1,1)`
pieces; they contribute their minimum rank one to the complement, so the
complement is still `20`.  This excludes every `d=22` decomposition whose
rank-one pieces are all diagonal.

### Non-diagonal rank-one pieces

The remaining case has a rank-one piece `(1,b)` or `(a,1)` with
`a,b>=2`.  At `d=22` the checker enumerates all unordered pairs in
`P_22` and finds exactly fifteen pairs that are not already above the
target by the linear bounds (8)--(10):

```math
(0,0),(0,2),(0,3),(0,4),(0,5),(0,6),\quad
(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),\quad
(2,2),(2,3),(2,4).
```

Every one has both maxima at most six.  Thus in a possibly surviving
profile every rank-one piece has other cut rank at most six.

Assume that at least one rank-one piece is non-diagonal.  Remove all
rank-one pieces, and let `r>=1` be their number and `T` the sum of
their other cut ranks.  Put `s=T-r`.  At least one removed piece is
non-diagonal, so `s>=1`; since every other rank is at most six,
`T<=6r`, hence `s<=5r`.  A rank-one piece with other rank `z` costs
`sqrt z` and contributes at most `z` to either cut total.  Therefore
the total reduction in either cut is at most `T`, and concavity of the
square root gives

```math
\sum_i\sqrt{z_i}
\ge\sqrt{T-r+1}+(r-1)=\sqrt{s+1}+(r-1).          \tag{19}
```

Let `M'` be the minimum rank of the remaining rank-one-free pieces.
If a remainder exists, `M'>=max(2,22-r)` and both remaining cut totals
are at least `N-T`.  The generalized form of (14) gives

```math
\text{remainder cost}
\ge M'-4+2\sqrt{2(N-T-M'+2)}.                     \tag{20}
```

For fixed `r` the resulting lower bound

```math
L(r,s)=\sqrt{s+1}+(r-1)+M'-4
      +2\sqrt{2(N-(r+s)-M'+2)}                     \tag{21}
```

is concave in `s`, so its minimum on
`1<=s<=min(5r,N-M'-r)` is attained at an endpoint.  The checker
evaluates every endpoint for `1<=r<=N`; the minimum is at
`(r,s)=(1,1)`, `M'=21`, and equals

```math
\boxed{\quad 17+\sqrt2+2\sqrt{926}
=79.2747097\ldots>22^{\sqrt2}.\quad}             \tag{22}
```

If no remainder exists, every piece is rank-one.  In a failing profile
each piece has `row rank + column rank <=7`, while the two cut budgets
force the sum to be at least `2N=968`; hence `r>=139`.  The cost is at
least `r>22^{\sqrt2}`.  If `r>N`, the same cost bound is immediate.
Thus no non-diagonal rank-one piece can occur in a `d=22`
decomposition below the target.

Combining this with the diagonal case excludes every `d=22`
decomposition.  Therefore (2) holds for every `2<=d<=22`, and a finite
seed passing the strict rank-sum criterion needs at least
**twenty-three labels**.  The finite profile and endpoint comparisons
are verified by the checker; the field-independent inputs are (1), the
quotient bound (4), and the nonlinear surplus lemma (14).

## Verification

The checker independently verifies the unrounded endpoint inequalities
for every selected dual, the higher-rank quotient bound on small
finite-field examples, and the minimum-cut budget on exact crossing
decompositions. It also checks the `d=22` profile enumeration, the
rank-one endpoint reduction, and records `d=23` as the first dimension
not covered, without treating any relaxation as a positive seed.

```sh
python tools/check_minimum_cut_rank_budget.py \
  --output evidence/minimum-cut-rank-budget.json
```

The field-independent proof is Sections 1--4. The finite arithmetic
certificate in Section 5 tests every permitted pair of scalar maxima;
Section 6 adds the finite `d=22` profile and rank-one endpoint checks.
Neither section infers a universal theorem from random tensor samples.
