# Residual cut budgets exclude seeds through fifty labels

## Result and scope

For every field and every exact decomposition of the `d`-label
crossing, let `r_i,c_i` be the two cut ranks and put
`S=sum_i sqrt(r_i c_i)`. Then

```math
S>d^{\sqrt2}\qquad(2\le d\le50).                       \tag{1}
```

A seed passing the strict sufficient rank-sum criterion therefore
needs at least **fifty-one labels**. This strengthens the
[forty-six-label certificate](weighted-cut-rank-budget.md).
The new step uses the ranks still needed after removing one actual
piece. It is a restriction on that finite seed criterion, not a
factoring lower bound, a new factoring exponent, or an assertion
that a fifty-one-label seed exists.

The proof retains all weighted cut budgets and the symmetric
square-root interpolant from the preceding note. Everything below
uses rational arithmetic, including the final finite certificate.

## 1. A chord lower bound for each remaining piece

Fix `k>=2`, and suppose `r,c>0` satisfy
`k^(-2)<=r/c<=k^2`. For `1<=s<k`, put

```math
\mu_{k,s}=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_{k,s}=\frac{k^2-s^2}{k^4-1}.                       \tag{2}
```

Both coefficients are nonnegative, and

```math
\begin{aligned}
\min(r,s^2c)&\ge\mu_{k,s}r+\nu_{k,s}c,\\
\min(s^2r,c)&\ge\nu_{k,s}r+\mu_{k,s}c.
\end{aligned}                                                   \tag{3}
```

To prove the first inequality, divide by `c` and set `t=r/c`.
The concave function `min(t,s^2)` lies above the chord connecting
its endpoint values at `t=1/k^2` and `t=k^2`. That chord has
exactly the slope and intercept (2). Exchange `r,c` for the
second inequality. There is no matrix-rank claim in this step;
it is an elementary inequality for two positive numbers.

## 2. An actual piece leaves nonzero rank budgets for the rest

Use the largest-imbalance classification from the preceding note.
The case `k=1` gives `S>=d^2` directly. For `k>=2`, choose an
actual piece with ranks `(a,b)`, exchanging the two cuts globally
if necessary, such that

```math
(k-1)^2b<a\le k^2b,\qquad a,b\le N:=d^2.
```

Set

```math
a_0=(k-1)^2b+1,\qquad a_1=\min(k^2b,N),
\qquad f_d(b)=\max_{0\le t\le d}t(d-bt).               \tag{4}
```

The integer parameter satisfies
`1<=b<=floor((N-1)/(k-1)^2)`. The sum of the other pieces is
`E_d-A`. Rank subadditivity and the quotient bound give

```math
\sum_{i\ne A}r_i\ge R_0:=\max(N-a_1,f_d(b)),\qquad
\sum_{i\ne A}c_i\ge C_0:=N-b.                         \tag{5}
```

When `b=1`, the all-field swap/product theorem also permits the
replacement

```math
R_0\longleftarrow\max\left(R_0,N-\frac{k(k+1)}2\right).
                                                               \tag{6}
```

Indeed `a<=k^2`, so the stable dimension in that theorem is at
most `k`. All rank totals in (5) concern the *other* pieces.
Every one of them still has ratio between `1/k^2` and `k^2`.
Summing (3) over those pieces and adding the selected piece gives

```math
\begin{aligned}
\sum_i\min(r_i,s^2c_i)
&\ge\min(a_0,s^2b)+\mu_{k,s}R_0+\nu_{k,s}C_0,\\
\sum_i\min(s^2r_i,c_i)
&\ge b+\nu_{k,s}R_0+\mu_{k,s}C_0.                    \tag{7}
\end{aligned}
```

For the second selected-piece term, `a>b` and `s>=1` imply
`min(s^2a,b)=b`. For the first, `a>=a_0` justifies the lower
replacement. Nonnegative coefficients in (2) justify substituting
the lower rank totals in (5).

Each left side in (7) also has the universal lower bound `sd`.
Taking the larger lower bound separately in each cut is valid.
The same ordinary rank surplus remains available:

```math
e_b(d,k)=
\begin{cases}
(k^2-3k+4)/2,&b=1,\\
\max(0,a_0+f_d(b)-N),&b\ge2,
\end{cases}
\qquad
\sum_i r_i+\sum_i c_i\ge2N+e_b(d,k).                  \tag{8}
```

The proof of (8), including the rank-one endpoint, is in Sections
3–4 of the weighted-cut note. Its use alongside (7) does not assume
that the different lower bounds are attained simultaneously. Each
is valid for the same decomposition, and all dual coefficients are
nonnegative.

## 3. The explicit rational certificate

Let `w_1=1/6` and
`w_s=2/((2s-1)(2s+1))` for `2<=s<k`. Define

```math
\begin{aligned}
U_s(d,k,b)&=\max\{sd,\min(a_0,s^2b)+\mu_{k,s}R_0+\nu_{k,s}C_0\},\\
V_s(d,k,b)&=\max\{sd,b+\nu_{k,s}R_0+\mu_{k,s}C_0\},\\
J(d,k,b)&=\frac{2d^2+e_b(d,k)}{2k-1}
        +\sum_{s=1}^{k-1}w_s\bigl(U_s(d,k,b)+V_s(d,k,b)\bigr).
\end{aligned}                                                   \tag{9}
```

The lower interpolant `Q_k(r,c)<=sqrt(rc)` from the preceding
note, summed over the pieces and combined with (7)–(8), proves
`S>=J(d,k,b)` for the actual extremal-piece class. Therefore

```math
\boxed{\quad
S\ge J_d:=\min\left(
\{d^2\}\cup
\left\{J(d,k,b):
2\le k\le d,\quad
1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor
\right\}\right).
\quad}                                                         \tag{10}
```

This formula is an all-field scalar relaxation. The actual piece
could satisfy stronger restrictions than its class permits; allowing
extra scalar cases only lowers `J_d` and preserves validity.
It has `O(d^2)` cases and `O(d)` summands per case.

The exact checker verifies, for `d=47,48,49,50`,

```math
J_d^{408}>d^{577},\qquad577^2-2\cdot408^2=1.            \tag{11}
```

The already proved certificate covers `2<=d<=46`. This proves (1).
For reference, the new lower bounds are approximately

| Labels | Certified lower bound |
|---:|---:|
| 47 | 234.7720984567 |
| 48 | 241.0500106747 |
| 49 | 247.0874738235 |
| 50 | 252.8726010042 |

The last bound exceeds `50^sqrt(2)=252.7588799437...`.
All comparisons use rational powers, not those decimal summaries.
For `d=51`, the same formula has value `259.2456221550...`,
below `51^sqrt(2)`. The checker proves that inequality using
the lower rational approximation `1393/985<sqrt(2)`.
This is the stopping point of (10), not a tensor construction.

## 4. Verification

Run

```sh
python tools/check_weighted_cut_residual_budget.py \
  --output evidence/weighted-cut-residual-budget.json
```

The checker verifies (3) by exact arithmetic, enumerates every scalar
class needed for (10), and checks the residual budgets against actual
small finite-field decompositions with low-rank pieces. The boundary
`k=1` is handled separately. It also records an explicit first
uncertified dimension, so a failed lower-bound certificate cannot be
mistaken for a successful seed.

The proof is ordinary mathematics. No new Lean formalization or
factoring-runtime improvement is asserted.
The refinement and its exact checker used substantial OpenAI assistant
support. The preceding note pins the source amplification criterion and
swap/product theorem used here.
