# A mixed-cut rank budget for every crossing decomposition

## Result and scope

The existing [geometric cut-budget theorem in PR #27](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md) already excludes exact finite seeds through 95 labels. This note gives a shorter independent argument for a smaller range; it does not extend that cutoff. A native finite-size exclusion does not by itself exclude asymptotic constructions obtained from tensor powers of a small polynomial source.

Every exact crossing decomposition on `q` labels, over every field,
satisfies

```math
\boxed{\quad
\sum_i\sqrt{r_i c_i}
\;\ge\;q-2+2\sqrt{q^2-q+1}.
\quad}
\tag{1}
```

Here `r_i` and `c_i` are the two ordinary cut ranks of its nonzero
pieces. In particular, **no seed on at most eleven labels** can
satisfy the strict finite-seed target

```math
\sum_i\sqrt{r_i c_i}<q^{\sqrt2}.
\tag{2}
```

The proof applies to arbitrary four-tensors and arbitrary coefficient
fields. It makes no assumption about nilpotence, permutation support,
polynomial paths, a coefficient span, or a chosen parameterization.
It strengthens the earlier
[universal finite-size gate](noncommutative-crossing-finite-gates.md),
which used the two cuts separately and excluded at most five labels.
The new ingredient is a budget that combines the two cuts.

This is an obstruction to the stated finite-seed sufficient criterion.
It is not a lower bound on factoring or on all possible clique evaluators.
It does not exclude larger seeds or an improvement through a different
algorithmic interface. The achieved factoring exponent remains `1/2`.

## 1. The two cut conventions

Let `V=F^q` and define

```math
E(x,y,u,v)=[x=v][y=u].
```

For a four-tensor `T`, use

```math
r(T)=\operatorname{rank}T_{(x,y),(u,v)},\qquad
c(T)=\operatorname{rank}T_{(x,u),(v,y)}.
\tag{3}
```

The row matrix of `E` is the flip `P` on `V\otimes V`,
`P(\xi\otimes\eta)=\eta\otimes\xi`. Both cuts of `E` have rank
`n=q^2`. Thus any decomposition `E=sum_i T_i` has

```math
\sum_i r_i\ge n,\qquad \sum_i c_i\ge n.
\tag{4}
```

Zero pieces may be removed, so all the remaining ranks are positive
integers at most `n`.

## 2. A mixed rank inequality

**Lemma.** If `E=A+B`, then

```math
r(A)+c(B)\ge q.
\tag{5}
```

Write `c=c(B)`. A column factorization of `B` gives matrices
`H_j,J_j` on `V` such that its row matrix is

```math
B_{\rm row}=\sum_{j=1}^c H_j\otimes J_j.
```

The reversal of the column index in (3) only transposes the second
factor in a rank factorization; it does not change this conclusion.
Fix any nonzero `xi in V`, and put

```math
W=\operatorname{span}\{H_1\xi,\ldots,H_c\xi\}.
```

Restrict the input to `xi tensor V`. Project the first output factor
onto `V/W`, retaining the second factor. This kills `B`. On the
same input, the flip becomes

```math
\eta\longmapsto[\eta]\otimes\xi.
```

This map has rank `q-dim(W)`, because `xi` is nonzero. Restriction
and projection cannot increase the rank of `A`. Therefore

```math
r(A)\ge q-\dim W\ge q-c(B),
```

which proves (5), including when its right side is nonpositive.
No extension of the coefficient field is needed.

For each piece of a decomposition, choose the cheaper cut. Put
the pieces with `r_i<=c_i` into `A`, and all other pieces into `B`.
Rank subadditivity and (5) give the new budget

```math
\boxed{\quad \sum_i\min(r_i,c_i)\ge q.\quad}
\tag{6}
```

The cut choice can vary from one piece to another. That freedom
is the reason (6) adds information to the two separate sums in (4).

## 3. Convert the three budgets to a rank-cost lower bound

Set

```math
x_i=\min(r_i,c_i),\quad d_i=|r_i-c_i|,\quad m=\sum_i x_i.
```

Equation (6) says `m>=q`. The cost of a piece is

```math
\sqrt{r_i c_i}
=x_i+f_{x_i}(d_i),\qquad
f_x(d)=\sqrt{x^2+xd}-x.
```

For fixed `d>=0`, the function `f_x(d)` increases with `x>0`.
For example, differentiating gives

```math
\frac{\partial f_x(d)}{\partial x}
=\frac{2x+d}{2\sqrt{x^2+xd}}-1\ge0,
```

because `(2x+d)^2-4(x^2+xd)=d^2>=0`. Since `x_i>=1`,

```math
f_{x_i}(d_i)\ge g(d_i),\qquad
g(d)=\sqrt{1+d}-1.
```

For nonnegative `a,b`, `g(a)+g(b)>=g(a+b)`. One direct proof
uses `g(0)=0` and the decreasing derivative
`g'(d)=1/(2 sqrt(1+d))`: the increment from `a` to `a+b` is at
most the increment from zero to `b`. Repeatedly applying this
inequality combines all surpluses on each side of the decomposition.

Let `R` be the sum of `d_i` over pieces with `r_i>=c_i`, and
let `C` be the corresponding sum over pieces with `c_i>r_i`.
The cut budgets (4) imply

```math
m+R\ge n,\qquad m+C\ge n.
```

If `q<=m<=n`, these inequalities give

```math
\begin{aligned}
S&=\sum_i\sqrt{r_i c_i}\\
 &\ge m+g(R)+g(C)\\
 &\ge m+2\sqrt{n-m+1}-2.
\end{aligned}
\tag{7}
```

The last expression increases with `m` on `[q,n]`, since its
derivative is `1-1/sqrt(n-m+1)>=0`. Its minimum occurs at `m=q`,
which is exactly (1).

If `m>=n`, simply use `S>=m>=n`. This also implies (1), since,
with `D=n-q+1`,

```math
n-\bigl(q-2+2\sqrt{n-q+1}\bigr)
=(\sqrt D-1)^2\ge0.
```

This completes the all-field proof.

### Sharpness for these numerical budgets

The abstract list of rank pairs

```math
(n-q+1,1),\quad(1,n-q+1),\quad
\underbrace{(1,1),\ldots,(1,1)}_{q-2\text{ pairs}}
\tag{8}
```

attains all three budgets and the lower bound (1). Thus no larger
bound follows solely by optimizing these three numerical constraints.
Equation (8) is not asserted to come from actual tensors whose sum
is the crossing. Tensor structure can impose additional restrictions.

For comparison, a simpler pointwise chord inequality gives

```math
\sqrt{rc}\ge
\frac{r+c+(q-1)\min(r,c)}{q+1},
```

and hence `S>=(3q^2-q)/(q+1)`. Formula (1) is the stronger version
used for the eleven-label cutoff.

## 4. Exact cutoff certificates

The rational inequalities

```math
\frac{140}{99}<\sqrt2<\frac{99}{70}
```

follow by squaring integers. The checker computes a rational lower
bound `L_q` for the right side of (1), using an integer square root,
and verifies

```math
L_q^{70}>q^{99},\qquad 2\le q\le11.
```

Thus (1) exceeds `q^sqrt(2)` throughout that range, without any
floating-point decision. At twelve labels the same lower bound is
strictly below `12^(140/99)`, also verified using integer powers.
This particular bound first leaves numerical room at twelve labels.
The [strict elementary refinement](strict-mixed-cut-gate.md) also excludes
twelve, while the prior geometric theorem excludes all sizes through 95.

For example, at eight labels (1) is `6+2 sqrt(57)`, greater than
21.09; the target is below 18.94. At eleven labels it is
`9+2 sqrt(111)`, greater than 30.07; the target is below 29.71.
The independently checked integer certificates are authoritative
for these comparisons.

## 5. Evidence and reproduction

Run

```sh
python3 tools/check_mixed_cut_rank_budget.py \
  --output evidence/mixed-cut-rank-budget.json
```

The checker independently reshapes all 65,536 binary four-tensors
on two labels and verifies (5) for every complementary pair.
It also checks (6), records an equality witness for the mixed
budget, verifies the abstract sharp profiles (8), and certifies
the cutoff with exact integer inequalities.

The negative control `E=I+(E+I)` in characteristic two has mixed
budget exactly two, despite both crossing cuts having rank four.
It prevents an unsupported replacement of the right side of (6)
by `q^2`. Finite checks support the stated controls; the proof above
establishes the uniform all-field result.
