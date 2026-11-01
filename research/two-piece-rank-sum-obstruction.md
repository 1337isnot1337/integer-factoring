# Two pieces cannot pass the finite rank-sum criterion

## All-field theorem

For any field, any `d>=2`, and any two-piece decomposition of the
crossing

\[
E_d(x,y,u,v)=[x=v][y=u]=A_1+A_2,
\]

let `r_i,c_i` be the ranks of `A_i` in the cuts
`(x,y)|(u,v)` and `(x,u)|(y,v)`. Then

\[
\boxed{\quad
\sqrt{r_1c_1}+\sqrt{r_2c_2}
\ge\sqrt{d^2-d+1}+\sqrt{(d^2-1)(d-1)}
>d^{3/2}>d^{\sqrt2}.
\quad}                                                        \tag{1}
\]

Thus every seed passing the strict
[finite rank-sum criterion](finite-seed-amplification.md)
has at least three nonzero pieces, regardless of its label dimension.
This complements the universal finite-dimensional cutoffs; it does
not exclude other amplification methods or give a factoring lower
bound.

## 1. Cross constraints from a quotient

If a tensor `B` has column rank `c`, its row matrix is a sum
`sum_(j=1)^c P_j tensor Q_j`. Fix a nonzero `u in F^d` and
quotient the first output factor by
`L=span(P_1u,...,P_cu)`. On the input `u tensor F^d`, all
product terms disappear. Swap still has rank at least `d-c`.
Consequently

\[
r(E_d-B)+c(B)\ge d.                                          \tag{2}
\]

This is also the `s=1` case of the
[higher-cut quotient bound](minimum-cut-rank-budget.md).
Apply (2) to both pieces. With `N=d^2`, rank subadditivity and
the two quotient inequalities give

\[
r_1+r_2\ge N,\quad c_1+c_2\ge N,\quad
r_1+c_2\ge d,\quad r_2+c_1\ge d.                              \tag{3}
\]

If one piece is zero, the other is the crossing and its cost is
`N`, which is at least the bound in (1). Hence assume all four
ranks are positive.

## 2. Reduce the scalar problem to exact totals

Holding the column ranks fixed, the row ranks can be decreased until
their sum is exactly `N` while preserving (3) and positivity.
Their lower bounds are
`max(1,d-c_2)` and `max(1,d-c_1)`, each at most `d-1`.
The sum of these lower bounds is at most `2d-2<N`, so the
decrease is feasible. It cannot increase the cost.

Now do the same to the column ranks, using the already adjusted row
ranks. This preserves both cross inequalities and leaves

\[
r_1+r_2=c_1+c_2=N.
\]

This is a scalar relaxation: no claim is made that the decreased
rank tuple comes from actual tensors.

Put `x=r_1` and `y=c_2`. Then `1<=x,y<=N-1`, and the
remaining constraints become

\[
d\le s:=x+y\le2N-d.
\]

The cost is

\[
f(x,y)=\sqrt{x(N-y)}+\sqrt{y(N-x)}.                            \tag{4}
\]

Replacing `(x,y)` by `(N-x,N-y)` preserves (4) and sends
`s` to `2N-s`. It is therefore enough to take `d<=s<=N`.

## 3. Two elementary concavity steps

At fixed `s<=N`, the admissible interval is `1<=x<=s-1` and

\[
f(x,s-x)=
\sqrt{x^2+(N-s)x}+
\sqrt{x^2-(N+s)x+Ns}.
\]

Each radical is concave on this interval. Its second derivative is
`-(N-s)^2/(4q(x)^(3/2))` for the corresponding positive
quadratic `q`. Thus the minimum is at `x=1` or `x=s-1`,
and both endpoints give

\[
g(s)=\sqrt{N-s+1}+\sqrt{(N-1)(s-1)}.                          \tag{5}
\]

Moreover `g` is nondecreasing for `s<=N`. Indeed, `g'(s)>=0`
is equivalent to

\[
(N-1)(N-s+1)\ge s-1,
\]

whose difference is `N(N-s)>=0`. Equations (4)--(5) are
therefore minimized at `s=d`. This proves the first inequality
in (1). The same calculation gives `g(d)<=g(N)=N`, as used
for a zero piece.

For the strict comparison, square the resulting expression:

\[
g(d)^2=d^3-2d+2+
2\sqrt{(d^2-d+1)(d^2-1)(d-1)}.
\]

The radical exceeds `d-1` because
`(d^2-d+1)(d^2-1)>d-1` for `d>=2`. Hence
`g(d)^2>d^3`. Finally `3/2>sqrt2` because `9>8`.
This completes the all-field proof.

## Reproduction and limits

The checker uses an exact integer comparison of sums of two square
roots. It verifies the scalar budget minimization for small dimensions
and the normalized scalar region through `d=12`. Independent
finite-field controls use arbitrary small matrices `A_1` and set
`A_2=E_d-A_1`; they check both cross constraints and (1).
These examples verify the transcription, while Sections 1--3 prove
the theorem for every dimension and every field.

```sh
python tools/check_two_piece_rank_sum.py \
  --output evidence/two-piece-rank-sum-obstruction.json
```
