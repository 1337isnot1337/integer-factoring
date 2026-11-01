# Geometric compression excludes exact rank-sum seeds through ninety-five labels

## Result

Every exact decomposition of the crossing over every field satisfies

```math
E_d=\sum_i A_i,
\qquad S:=\sum_i\sqrt{r(A_i)c(A_i)}>d^{\sqrt2}
\qquad(2\le d\le95).                                      \tag{1}
```

Thus a seed passing the strict rank-sum amplification criterion needs
at least **ninety-six labels**. In particular, a binary-label search must
start at dimension `128` or above. This excludes one sufficient seed
criterion; it does not prove a lower bound for factoring. The achieved
factoring exponent remains `1/2`.

The new ingredient is the
[universal geometric compression theorem](geometric-subspace-compression.md).
Explicit integer Jacobian minors give common input and output subspaces
for every tuple of matrices, over every field after scalar extension.
This improves the complement bounds used in the
[flagged rank-budget theorem](flagged-cut-rank-budget.md), which already
excluded dimensions through `82`. The present exact rational replay
checks every remaining scalar class for `83<=d<=95`.

## 1. The stronger complement table

For `m>=2`, define `F_geo(d,m)` as the maximum product `s*q` over positive
integers `s,q` and integers `2<=A<=2m-2`, `B=2m-A`, satisfying

```math
2s+Aq\le2d,\qquad Bs+2q\le2d.                         \tag{2}
```

Use zero when no positive pair is feasible. The geometric theorem proves
that every sum `B_tensor` of at most `m` product operators satisfies
`rank(E_d-B_tensor)>=F_geo(d,m)`. Even `A,B` in (2) come from disjoint
Jacobian blocks. Odd `A,B` come from a three-map system whose four basic
minors have determinant `1` or `-1`. Those integral certificates make the
proof independent of the characteristic.

For exact crossings retain the sharper one- and two-product theorem:

```math
F_d(m)=
\begin{cases}
d^2,&m=0,\\
d(d-1)/2,&m=1,2,\\
F_{\rm geo}(d,m),&m\ge3.
\end{cases}                                             \tag{3}
```

For every subset `J` of an exact crossing decomposition, concatenate
column-rank factorizations of its pieces. Their sum has a product
representation with at most `m=sum_(i in J)c_i` terms. Therefore

```math
\sum_{i\notin J}r_i\ge F_d(m).                           \tag{4}
```

The same statement holds with the cuts exchanged. No independence of the
product terms is required; padding a representation by zero terms is allowed.

The geometric table is zero for `m>=2d-1`, and can remain positive through
`m=2d-2`. Consequently the complete finite weighted minimum is

```math
G_d(s;b)=\min_{b\le m\le\max(b,2d-1)}
                  \bigl(F_d(m)+s^2m\bigr).               \tag{5}
```

The older cap `d+1` is insufficient for (3). The checker verifies the
correct closure by comparison with longer minima. For
`J={i:r_i>s^2c_i}`, equation (4) gives

```math
\sum_i\min(r_i,s^2c_i)\ge G_d(s;0),\qquad
\sum_i\min(s^2r_i,c_i)\ge G_d(s;0).                     \tag{6}
```

## 2. Conditioning on an actual piece

Discard zero pieces and choose one attaining the maximum of all rank
ratios and reciprocal rank ratios. Exchange the cuts so its ranks are
`(a,b)` with `a>=b`. If every ratio is one, `S>=d^2` settles (1).
Otherwise let `k` be the ceiling of the square root of the largest ratio.
Then

```math
a_0=(k-1)^2b+1\le a\le a_1=\min(k^2b,d^2),\qquad
2\le k\le d,\quad
1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor. \tag{7}
```

For `1<=s<k`, this piece belongs strictly to the minimizing subset in
the first inequality of (6). Its product count is at least `b`, so that
side improves to `G_d(s;b)`. The reverse cut keeps `G_d(s;0)`.

Write `n=d^2`. The residual pieces have rank totals at least

```math
R=\max(n-a_1,F_d(b)),\qquad C=n-b.                     \tag{8}
```

For `b=1`, additionally take the maximum of `R` with
`n-k(k+1)/2`. This is the exact stable-dimension swap/product bound
proved in the [residual-budget note](weighted-cut-residual-budget.md).
It also gives the following surplus in the sum of both ordinary rank
budgets:

```math
\sum_i(r_i+c_i)\ge2n+e,\qquad
e=
\begin{cases}
(k^2-3k+4)/2,&b=1,\\
\max(0,a_0+F_d(b)-n),&b\ge2.
\end{cases}                                             \tag{9}
```

All piece ratios lie in `[1/k^2,k^2]`. For each `1<=s<k`, put

```math
\mu_s=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_s=\frac{k^2-s^2}{k^4-1}.
```

The concave piecewise linear function `min(t,s^2)` lies above its
endpoint chord. Summing the resulting inequalities over the residual
pieces, then adding the selected piece, proves

```math
\begin{aligned}
\sum_i\min(r_i,s^2c_i)&\ge
 U_s:=\max\{G_d(s;b),s^2b+\mu_sR+\nu_sC\},\\
\sum_i\min(s^2r_i,c_i)&\ge
 V_s:=\max\{G_d(s;0),b+\nu_sR+\mu_sC\}.
\end{aligned}                                           \tag{10}
```

These statements refer to the same decomposition and chosen piece, so
their nonnegative combinations are valid simultaneously.

## 3. An exact scalar certificate

The [weighted cut interpolant](weighted-cut-rank-budget.md) gives

```math
\sqrt{rc}\ge\frac{r+c}{2k-1}
 +\sum_{s=1}^{k-1}w_s
   \bigl(\min(r,s^2c)+\min(s^2r,c)\bigr),
\quad
w_1=\frac16,\quad
w_s=\frac2{(2s-1)(2s+1)}\ (s\ge2),                    \tag{11}
```

whenever `1/k^2<=r/c<=k^2`. Combining (9)--(11) proves

```math
S\ge J(d,k,b):=\frac{2d^2+e}{2k-1}
                 +\sum_{s=1}^{k-1}w_s(U_s+V_s).        \tag{12}
```

Every term in (12) is rational. The replay enumerates all `167,984`
classes, including the equal-rank classes, for `83<=d<=95`. It checks

```math
J(d,k,b)^{408}>d^{577},\qquad577^2-2\cdot408^2=1.       \tag{13}
```

This proves the strict target comparison without floating-point input.
Together with the already verified exclusion through `82`, it proves (1).
At the final dimension the minimum is

```math
J(95,51,3)=\frac{203689540088152253}{324506874817125}
          =627.689444801257\ldots
          >95^{\sqrt2}=626.503563\ldots.              \tag{14}
```

At `d=96`, the class `(k,b)=(52,3)` gives a value below the target.
The checker proves that comparison using the lower Pell exponent
`1393/985`, with `1393^2-2*985^2=-1`. This records the limit of (12),
not a realizable seed or an impossibility of further refinement.

## 4. Replay and limits

Run

```sh
python3 tools/check_geometric_subspace_compression.py \
  --output evidence/geometric-subspace-compression.json
python3 tools/check_geometric_cut_rank_budget.py \
  --output evidence/geometric-cut-rank-budget.json
```

The first replay checks the integer geometric witnesses and finite-field
controls. The second checks all rational scalar classes, `9,690`
complement comparisons with the preceding flag bound, `38,000` weighted
budget comparisons, and `38,000` extended-minimum closure checks. Both use
only the Python standard library. Their source hashes are recorded in
the evidence files.

The [dyadic scalar-profile extension](geometric-cut-relaxation-limit.md)
shows that even these geometric subset budgets admit abstract rank
lists of cost `O(d log d)`. Realizing those lists by compatible crossing
pieces remains open. The separate
[formal-limit theorem](rank-one-formal-cutoff.md) reaches `71` labels;
the exact-only small-product bounds used here are not transferred to
Laurent-series degenerations.

This proof and the exact replay used substantial OpenAI assistant support.
