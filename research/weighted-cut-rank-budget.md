# Weighted cut budgets and a forty-six-label cutoff

## Results and scope

Let `F` be any field and let

```math
E_d(x,y,u,v)=[x=v][y=u]=\sum_i A_i(x,y,u,v),\qquad d\ge2.
```

Write `r_i,c_i` for the ranks of a nonzero piece in the cuts
`(x,y)|(u,v)` and `(x,u)|(y,v)`, respectively. Put
`S=\sum_i\sqrt{r_ic_i}`. The following inequalities hold for every
integer `1<=s<=d`:

```math
\sum_i\min(r_i,s^2c_i)\ge sd,\qquad
\sum_i\min(s^2r_i,c_i)\ge sd.                         \tag{1}
```

They strengthen the unweighted minimum-cut budget. In particular,

```math
\boxed{\quad
S\ge d\bigl(2H_{2d}-H_d-1\bigr),\qquad
H_m=\sum_{j=1}^{m}\frac1j.
\quad}                                                \tag{2}
```

The right side is `d(log d+2log 2+gamma-1+o(1))`. The exact
rational bound (2) alone excludes the strict finite rank-sum target
`S<d^sqrt(2)` through forty-one labels. Combining the same weighted
budgets with an actual piece of largest rank imbalance excludes the
target through **forty-six labels**, over every field. Therefore a
finite seed passing that particular sufficient criterion needs at
least forty-seven labels.

These are restrictions on the finite rank-sum criterion from
[finite-seed amplification](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/finite-seed-amplification.md). They
do not prove a lower bound for factoring, exclude different ways of
amplifying a seed, or establish a new factoring exponent. The previous
[minimum-cut certificate](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/minimum-cut-rank-budget.md) excluded seeds
through twenty-one labels. The new argument uses all weighted cut
budgets and a rational certificate; it does not enumerate tensors.

## 1. The quotient inequality and all subset choices

Let `V=F^d`, and let `S_V` swap the two factors of `V tensor V`.
If a four-tensor `B` has column-cut rank `k`, its row matrix is

```math
B_{\rm row}=\sum_{j=1}^{k}P_j\otimes Q_j
```

for endomorphisms `P_j,Q_j` of `V`. Choose any subspace `U` of
dimension `s`, where `0<=s<=d`, and set `L=\sum_j P_jU`.
Then `dim L<=ks`. Restrict the input to `U tensor V` and quotient
the first output factor by `L`. Every product operator vanishes.
Swap maps that input onto `V tensor U`, whose image in the quotient
has dimension `s(d-dim L)`. Restriction and quotient cannot increase
rank. Consequently,

```math
r(E_d-B)+s^2c(B)\ge sd.                               \tag{3}
```

This remains valid when `ks>d`; the resulting lower bound for the
rank before rearranging may be negative. The argument uses no
positivity, algebraic closure, or field-characteristic assumption.

For any subset `J` of pieces, apply (3) to `B=\sum_{i in J}A_i`.
Rank subadditivity gives

```math
\sum_{i\notin J}r_i+s^2\sum_{i\in J}c_i\ge sd.          \tag{4}
```

Choose `J` separately for each `s`, placing `i` in `J` exactly
when `s^2c_i<r_i`. This minimizes the left side term by term and
proves the first inequality of (1). Exchanging the indices `y,u`
fixes `E_d` and exchanges its two cuts, proving the second.

For `s=1` this recovers `sum min(r_i,c_i)>=d`. Both ordinary
rank budgets `sum r_i>=d^2` and `sum c_i>=d^2` are also retained.

## 2. A symmetric lower interpolant for the geometric mean

For an integer `k>=2`, define

```math
\begin{aligned}
a_k&=\frac1{2k-1},\qquad w_1=\frac16,\qquad
w_s=\frac{2}{(2s-1)(2s+1)}\quad(2\le s\le k-1),\\
Q_k(r,c)&=a_k(r+c)
  +\sum_{s=1}^{k-1}w_s
       \bigl(\min(r,s^2c)+\min(s^2r,c)\bigr).
\end{aligned}                                                   \tag{5}
```

If `r,c>0` and `k^(-2)<=r/c<=k^2`, then

```math
Q_k(r,c)\le\sqrt{rc}.                                 \tag{6}
```

To prove this, use homogeneity and symmetry to set `c=1` and
`1<=t=r<=k^2`. Telescoping gives `Q_k(1,1)=1` and, on the
interval `j^2<=t<=(j+1)^2`,

```math
\frac{d}{dt}Q_k(t,1)
=\frac1{2k-1}+\sum_{s=j+1}^{k-1}
       \frac2{(2s-1)(2s+1)}
=\frac1{2j+1}\qquad(1\le j<k).                       \tag{7}
```

Thus `Q_k(t,1)` is exactly the linear interpolation of `sqrt(t)`
at `1,4,9,...,k^2`. Concavity of the square root proves (6)
on this interval. Symmetry gives the reciprocal interval as well.
The unusual coefficient `w_1=1/6` makes the interpolation join
correctly at ratio one; substituting the formula for `s>=2` at
`s=1` would be incorrect.

All nonzero piece ranks lie in `[1,d^2]`, so (6) applies with
`k=d`. Sum (5), use (1), and use both ordinary rank budgets:

```math
S\ge\frac{2d^2}{2d-1}
       +2d\sum_{s=1}^{d-1}s w_s
 =d\left(2\sum_{s=1}^{d}\frac1{2s-1}-1\right)
 =d(2H_{2d}-H_d-1).                                    \tag{8}
```

This proves the all-dimension theorem (2). The asymptotic expression
in the opening section follows from the usual harmonic-sum estimate;
the finite certificate uses (8) exactly and needs no such estimate.

## 3. Truncate the interpolant at the largest actual imbalance

For a given decomposition, define

```math
R_* = \max_i\max\{r_i/c_i,c_i/r_i\},\qquad
k=\lceil\sqrt{R_*}\rceil.                              \tag{9}
```

Here `1<=k<=d`. If `k=1`, all `r_i=c_i`, so `S=sum r_i>=d^2`.
Suppose `k>=2`. Every piece meets the ratio restriction in (6),
so the same proof, with `k` in place of `d`, gives

```math
S\ge B(d,k):=
d(2H_{2k}-H_k-1)+\frac{2d(d-k)}{2k-1}.                 \tag{10}
```

If additionally `sum r_i+sum c_i>=2d^2+e`, then (5) gives
the stronger bound

```math
S\ge B(d,k)+\frac e{2k-1}.                            \tag{11}
```

There is an actual piece whose ranks, after exchanging the two cuts
if needed, are `(a,b)` with

```math
(k-1)^2b<a\le k^2b,\qquad a,b\le d^2.                 \tag{12}
```

The existence of this piece, rather than only a bound on the allowed
rank ratios, supplies `e` in the next section.

## 4. Charge the complement of that piece

Write `N=d^2` and define the nonnegative integer

```math
f_d(b)=\max_{0\le s\le d}s(d-bs).                      \tag{13}
```

Applying (3) to the selected piece gives
`r(E_d-A)>=f_d(b)`. Hence

```math
\sum_i r_i\ge a+f_d(b),\qquad
\sum_i r_i+\sum_i c_i
\ge2N+\max(0,a+f_d(b)-N).                              \tag{14}
```

For `b>=2`, (12) allows the lower replacement
`a_0=(k-1)^2b+1`. It also forces
`b<=floor((N-1)/(k-1)^2)`. Therefore a valid scalar bound is

```math
e_b(d,k)=\max\{0,(k-1)^2b+1+f_d(b)-N\},
\quad 2\le b\le\left\lfloor\frac{N-1}{(k-1)^2}\right\rfloor.
                                                               \tag{15}
```

The rank-one case has a stronger available bound. The all-field
[swap/product theorem](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/swap-product-rank-sum.md) proves that for
a column-rank-one piece of row rank `a`,

```math
r(E_d-A)\ge N-\frac{h(h+1)}2,
\qquad h\le\lfloor\sqrt a\rfloor.                    \tag{16}
```

Indeed the row matrix is `P tensor Q`; the stable invertible part
used in that theorem has dimension at most
`min(rank P,rank Q)<=sqrt(a)`.
For `(k-1)^2<a<k^2`, the integer square root is `k-1`, giving

```math
a-\frac{h(h+1)}2
\ge(k-1)^2+1-\frac{k(k-1)}2
=\frac{k^2-3k+4}{2}.
```

At `a=k^2`, (16) instead gives `k(k-1)/2`, which is no
smaller, since their difference is `k-2`. Thus put

```math
e_1(k)=\frac{k^2-3k+4}{2}.                             \tag{17}
```

This is a nonnegative integer. Define

```math
\begin{aligned}
e(d,k)&=\min\Bigl(\{e_1(k)\}\cup
  \{e_b(d,k):2\le b\le\lfloor(N-1)/(k-1)^2\rfloor\}\Bigr),\\
L_d&=\min\left(\{d^2\}\cup
  \left\{B(d,k)+\frac{e(d,k)}{2k-1}:2\le k\le d\right\}\right).
\end{aligned}                                                   \tag{18}
```

Every actual decomposition has one of these values of `k` and one
of the permitted extremal-piece cases. Equations (11), (15), and
(17) therefore prove `S>=L_d`. The minimization relaxes tensor
realizability, so extra scalar profiles can only weaken the result.

The whole certificate has only `O(d^2)` scalar cases per dimension:
`sum_{k=2}^d floor((d^2-1)/(k-1)^2)=O(d^2)`. Evaluate (13)
at the two integers surrounding `d/(2b)`; that quadratic is
concave and its maximum is attained at one of them. No rank oracle,
tensor search, floating-point optimization, or unverified numerical
minimum is involved.

## 5. Exact finite certificate and its stopping point

The accompanying standard-library checker computes (18) with exact
rational arithmetic. For every `2<=d<=46` it verifies

```math
L_d^{408}>d^{577},\qquad 577^2-2\cdot408^2=1.          \tag{19}
```

Thus `L_d>d^(577/408)>d^sqrt(2)`, excluding the strict target.
The final dimension has

```math
L_{46}=224.9098175102055\ldots,
\qquad46^{\sqrt2}=224.6439239425683\ldots.
```

These decimals are summaries of exact comparisons, not their evidence.
The minimum scalar case at `d=46` has `k=31`, `b=2`, and
zero extra rank budget. At `d=47` the certificate gives
`L_47=230.5574231437782...`, which is below `47^sqrt(2)`.
The checker certifies this failure using

```math
L_{47}^{985}<47^{1393},\qquad1393^2-2\cdot985^2=-1.    \tag{20}
```

So the first uncertified dimension is a limit of this scalar
certificate. It is not a passing tensor seed or evidence that a
forty-seven-label seed exists.

For controls, the checker also verifies the interpolant against
the geometric mean by exact squaring, tests the quotient inequality
on independently generated finite-field matrices, and checks all
weighted budgets on finite exact crossing decompositions. A rank
profile satisfying only the earlier unweighted budgets fails a
weighted budget, confirming that (1) adds a restriction.

```sh
python tools/check_weighted_cut_rank_budget.py \
  --output evidence/weighted-cut-rank-budget.json
```

The finite field examples check the transcription; Sections 1–4 prove
the all-field statements. No new Lean formalization accompanies this
result. The new weighted inequalities, interpolating dual, proof, and
checker used substantial OpenAI assistant support. The prior finite-seed
transfer and swap/product refinement remain credited to their pinned
repository notes above.
