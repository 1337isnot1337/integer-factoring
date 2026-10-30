# Flag quotients exclude finite rank-sum seeds through eighty-two labels

## Result and scope

Every exact decomposition of the `d`-label crossing over every field satisfies

```math
S:=\sum_i\sqrt{r_i c_i}>d^{\sqrt2}\qquad(2\le d\le82).       \tag{1}
```

Thus a seed meeting the strict sufficient rank-sum criterion needs at least
**83 labels**. In particular, no 64-label tensor decomposition can meet that
criterion, regardless of its coefficient field, number of pieces, or source
construction. This is an exclusion of one sufficient amplification criterion,
not a lower bound on integer factoring. The achieved factoring exponent
remains `1/2`.

The new structural input is a stronger complement bound for **every** number
of product operators. Two of their first factors admit common input/output
flags, so their images cost one subspace dimension instead of two. Combining
this observation with the earlier sharp one- and two-product bounds improves
all weighted cut budgets. The unconditioned budgets alone prove (1) through
77 labels. Conditioning on an actual piece of largest imbalance extends the
exact certificate through 82.

The proof uses the [two-map flag lemma](two-product-complement.md), the
[weighted cut interpolant](weighted-cut-rank-budget.md), and the
[residual chord bounds](weighted-cut-residual-budget.md). It is ordinary
mathematics with exact integer and rational replay, without new Lean work or
a claim of literature priority. The construction and checker used substantial
OpenAI assistant support.

## 1. Saving one product in the quotient argument

Let `V=F^d`, let `E_d` swap the factors of `V tensor V`, and write

```math
B=\sum_{j=1}^m P_j\otimes Q_j.
```

For `m>=2`, extend scalars to an algebraic closure. The first two maps
`P_1,P_2` admit complete input and output flags `U_s,W_s` of dimension `s`
such that `P_1 U_s+P_2 U_s` is contained in `W_s`. This is simultaneous
equivalence of two maps; it does not assume simultaneous similarity or an
invertible matrix in their pencil. The inductive proof in the two-product
note also covers wholly singular pencils.

For a chosen `s`, put

```math
L_s=\sum_{j=1}^m P_jU_s.
```

The first two summands together have dimension at most `s`; each remaining
summand adds at most `s`. Therefore

```math
\dim L_s\le(m-1)s.
```

Restrict the input to `U_s tensor V` and quotient the first output factor
by `L_s`. Every summand of `B` vanishes. The swap sends this input to
`V tensor U_s`, whose image in the quotient has dimension
`s(d-dim L_s)`. No identification of the input and output flags is needed.
Restriction and quotient cannot increase rank. Hence

```math
\boxed{\quad
r(E_d-B)\ge\max_{0\le s\le d}s\bigl(d-(m-1)s\bigr),\qquad m\ge2.
\quad}                                                        \tag{2}
```

Ranks are unchanged by scalar extension, so (2) holds over the original
field. Zero summands and dependent product representations cause no problem.
The number `m` may be an upper bound on the actual column-cut rank.

## 2. A complement table valid for padded representations

Retain the sharper half-dimensional bound for one or two products. Define

```math
F_d(m)=
\begin{cases}
d^2,&m=0,\\
d(d-1)/2,&m=1\text{ or }m=2,\\
\displaystyle\max_{0\le t\le d}t\bigl(d-(m-1)t\bigr),&m\ge3.
\end{cases}                                                   \tag{3}
```

Every representation of `B` by at most `m` products satisfies

```math
r(E_d-B)\ge F_d(m).                                      \tag{4}
```

For `m=0` the statement is the rank of swap. For `m=1,2` it follows from
the earlier [two-product complement theorem](two-product-complement.md).
For `m>=3`, use (2). These bounds do not require `m` to equal the minimal
number of products.

For `m>=d+1`, the maximum in (3) is zero: `t=0` gives zero and every
integer `t>=1` gives a nonpositive value. Thus the table is eventually zero.
For computation, the integer maximum for `m>=3` is attained at one of
the two integers nearest `d/(2(m-1))`. The checker independently compares
that implementation with direct enumeration of all `t`.

## 3. Stronger weighted budgets for every decomposition

Suppose `E_d=sum_i A_i`, with positive cut ranks `(r_i,c_i)` after discarding
zero pieces. For a subset `J`, choose a column-rank factorization of every
`A_i` with `i in J`. Their concatenation represents `B=sum_(i in J)A_i`
by

```math
m=\sum_{i\in J}c_i
```

products. Equation (4) and rank subadditivity give

```math
\sum_{i\notin J}r_i\ge F_d(m).                         \tag{5}
```

For an integer `0<=s<=d` and a nonnegative integer `b`, define

```math
G_d(s;b)=\min_{m\ge b}\bigl(F_d(m)+s^2m\bigr),\qquad
G_d(s)=G_d(s;0).                                      \tag{6}
```

This is a finite calculation: the minimum may stop at `max(b,d+1)`.
Beyond that point `F_d(m)=0` and `s^2m` is nondecreasing, including at
`s=0`.

Choose `J={i:r_i>s^2c_i}`. This minimizes the cut cost term by term.
Combining (5) with (6) proves

```math
\sum_i\min(r_i,s^2c_i)\ge G_d(s),\qquad
\sum_i\min(s^2r_i,c_i)\ge G_d(s).                     \tag{7}
```

The second inequality follows by exchanging the two cuts, which fixes the
crossing. Both ordinary rank sums are still at least `d^2`.

A simple visible improvement is

```math
G_d(s)\ge sd+s^2,\qquad1\le s\le\lfloor(d-1)/2\rfloor.        \tag{8}
```

To see this directly, if the minimizing subset has `m>=2`, use (2) at
`t=s` to get `F_d(m)+s^2m>=sd+s^2`; the stronger cases in (3) remain
valid. For `m=1`, its cost is at least `d(d-1)/2+s^2>=sd+s^2`.
For `m=0`, it is `d^2>=sd+s^2`. The exact minimum (6), rather than
only (8), is used in the certificate.

## 4. An unconditioned bound through 77 labels

Use the nonnegative coefficients

```math
w_1=\frac16,\qquad
w_s=\frac2{(2s-1)(2s+1)}\quad(s\ge2).
```

The earlier interpolant proves that for positive ranks with ratio between
`1/k^2` and `k^2`,

```math
\frac{r+c}{2k-1}
 +\sum_{s=1}^{k-1}w_s\bigl(\min(r,s^2c)+\min(s^2r,c)\bigr)
 \le\sqrt{rc}.                                       \tag{9}
```

All `d`-label ranks lie in `[1,d^2]`. Take `k=d`, sum (9), and apply
(7) and the ordinary rank sums:

```math
S\ge H_d^{\rm flag}:=
\frac{2d^2}{2d-1}+2\sum_{s=1}^{d-1}w_s G_d(s).        \tag{10}
```

The exact checker verifies `(H_d^flag)^408>d^577` for every `2<=d<=77`.
Since `577^2-2*408^2=1`, this proves the strict comparison with `d^sqrt(2)`.
At `d=78`, (10) is below the target; that negative control is also proved
with rational powers. The next section supplies the additional argument.

## 5. Condition on the piece of largest imbalance

Let `k=ceil(sqrt(max_i max(r_i/c_i,c_i/r_i)))`. If `k=1`, every piece
has equal cut ranks and `S>=d^2`. Otherwise choose an actual piece `A`
of ranks `(a,b)`, exchanging the cuts globally if necessary, with

```math
a_0:=(k-1)^2b+1\le a\le a_1:=\min(k^2b,d^2),\qquad
1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor.        \tag{11}
```

For every `1<=s<k`, the inequality `a>s^2b` is strict. Thus `A` belongs
to the canonical minimizing subset in the first cut of (7). Its product
budget is at least `b`, so that particular inequality strengthens to

```math
\sum_i\min(r_i,s^2c_i)\ge G_d(s;b).                  \tag{12}
```

The other cut keeps the unconditioned budget `G_d(s)`.

The remaining pieces have row and column rank sums at least

```math
R_0=\max(d^2-a_1,F_d(b)),\qquad C_0=d^2-b.            \tag{13}
```

For `b=1`, also replace `R_0` by its maximum with
`d^2-k(k+1)/2`. This uses the stable-dimension swap/product theorem already
proved and applied in the [residual-budget note](weighted-cut-residual-budget.md).
The corresponding total-rank surplus is

```math
e(d,k,b)=
\begin{cases}
(k^2-3k+4)/2,&b=1,\\
\max(0,a_0+F_d(b)-d^2),&b\ge2,
\end{cases}
\qquad\sum_i(r_i+c_i)\ge2d^2+e(d,k,b).                \tag{14}
```

The rank-one formula in (14) retains the integer stable-dimension argument
from the earlier note; simply subtracting the upper-endpoint bound in (13)
would give a weaker formula.

Define the nonnegative chord coefficients

```math
\mu_{k,s}=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_{k,s}=\frac{k^2-s^2}{k^4-1}.
```

Every remaining piece has ratio in `[1/k^2,k^2]`, so
`min(r,s^2c)>=mu*r+nu*c` and its reversed version hold. Summing them,
adding the selected piece, and retaining (7) and (12) gives

```math
\begin{aligned}
U_s&=\max\{G_d(s;b),\min(a_0,s^2b)+\mu_{k,s}R_0+\nu_{k,s}C_0\},\\
V_s&=\max\{G_d(s),b+\nu_{k,s}R_0+\mu_{k,s}C_0\}.
\end{aligned}                                                   \tag{15}
```

All coefficients are nonnegative, and all inequalities hold for the same
decomposition. Summing (9) now gives the rational bound

```math
S\ge J^{\rm flag}(d,k,b):=
\frac{2d^2+e(d,k,b)}{2k-1}
 +\sum_{s=1}^{k-1}w_s(U_s+V_s).                       \tag{16}
```

The verifier enumerates every integer class (11) for `78<=d<=82` and
proves `(J^flag)^408>d^577` in every class. Together with (10), this
proves (1). Approximate summaries of the exact minima are:

| Labels | Certified lower bound | Minimizing class `(k,b)` |
|---:|---:|---:|
| 78 | 480.0474341923 | `(43,3)` |
| 79 | 487.5339500607 | `(43,3)` |
| 80 | 494.4818267325 | `(44,3)` |
| 81 | 501.9686331303 | `(44,3)` |
| 82 | 508.9302867823 | `(45,3)` |

The JSON contains the exact fractions and every per-`k` minimum; the
decimal summaries are not proof inputs. At dimension 83, the class
`(k,b)=(45,3)` gives a bound below `83^sqrt(2)`. This is a limit of this
certificate, not evidence of a realizable seed.

## 6. Replay and remaining scope

```sh
python3 tools/check_flagged_cut_rank_budget.py \
  --output evidence/flagged-cut-rank-budget.json
```

The checker uses only the Python standard library. It verifies the finite
minimum closure, the integer quadratic maxima, every required scalar class,
and the strict target comparisons. It also checks the subset budgets,
forced inclusion of the selected piece, and the full interpolant bound on
actual finite-field decompositions. These finite controls supplement the
all-field proof; they do not establish it by extrapolation.

The [dyadic scalar profiles](weighted-cut-relaxation-limit.md) still meet
the strengthened subset budgets with cost `O(d log d)`. Thus the new flag
constraints do not settle the criterion at every dimension. A compatible
tensor realization of those abstract profiles remains unknown.

The half-dimensional complement bounds in (3) concern exact crossings.
This note does not transfer the 82-label cutoff to formal degenerations
with Laurent poles; those require their own moving-lattice argument and
their own complement table.
