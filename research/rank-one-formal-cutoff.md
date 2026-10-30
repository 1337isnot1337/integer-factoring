# Factor kernels exclude formal crossing limits through seventy-one labels

## Result

For every field and every supplied formal decomposition

```math
T(t)=\sum_i A_i(t),\qquad T(t)\in\operatorname{Mat}_{d^2}(F[[t]]),
\qquad T(0)=E_d,
```

let `r_i,c_i` denote the two cut ranks of the Laurent-series pieces.
Then

```math
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}
\qquad(2\le d\le71).                               (1)
```

The [geometric formal-limit theorem](geometric-border-cutoff.md)
already proves this through sixty-four labels. The present extension
uses the kernels of the two factors of an opposite-cut-rank-one
piece. It applies with poles and does not import the stronger
exact-only swap/product or half-rank complement theorems.

A formal seed passing the strict rank-sum criterion therefore needs
at least **seventy-two labels**. This is not a claim that such a seed
exists, and it does not improve the factoring exponent above `1/2`.

## 1. A rank-one complement bound that survives poles

Suppose `T=A+B`, where `T(0)=E_d` and `B=P tensor Q` over
`F((t))`. If the first factor has rank `p`, its kernel has dimension
`d-p`. For every `s`, choose an `s`-dimensional subspace `U` with

```math
\dim P(U)=\max(0,s-d+p).                            (2)
```

Take as much of the kernel as possible, then extend using vectors
whose images are independent. Projecting out `P(U)` kills `B` on
`U tensor V`. The saturated input and output lattices from the
geometric theorem make the projection integral before reducing
`T`. Its reduction has rank

```math
s\bigl(d-\max(0,s-d+p)\bigr).
```

Maximizing over integer `s` gives

```math
r(T-B)\ge K_d(p):=\left\lfloor\frac{(2d-p)^2}{4}\right\rfloor.
                                                               (3)
```

Indeed, above `s=d-p` the expression is `s(2d-p-s)`; its vertex
lies between `d-p` and `d`. The part below `d-p` is no larger.
This argument is valid for moving subspaces and Laurent poles,
because the entire product term is annihilated before specialization.

Exchange the two tensor factors if necessary. If `a=r(B)` and
`c(B)=1`, then

```math
a=\operatorname{rank}P\,\operatorname{rank}Q,
\qquad\min(\operatorname{rank}P,\operatorname{rank}Q)
 \le\lfloor\sqrt a\rfloor.
```

Since `K_d(p)` decreases with `p`,

```math
\boxed{\quad r(T-B)\ge K_d(\lfloor\sqrt a\rfloor).\quad}          (4)
```

The same statement holds after exchanging the two cuts. Equation
(4) uses factor rank and kernel dimension only. It is distinct from
the exact crossing's stable-rank Fitting argument.

## 2. Classify an actual piece of maximum ratio

Discard zero pieces. Choose a piece that **attains the largest of
all ratios `r_i/c_i` and `c_i/r_i`**. Exchange the cuts globally so
its ranks are `(a,b)` with `a/b>=1`. This actual maximum, not an
arbitrary piece in the same rounded class, is used throughout.

If the maximum is one, all ranks are equal and the rank sum is at
least `d^2`. Otherwise choose the integer `k` with

```math
(k-1)^2<a/b\le k^2,
\qquad2\le k\le d,
\qquad1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor.
                                                               (5)
```

Every piece has rank ratio in `[1/k^2,k^2]`. In particular, any
piece whose smaller cut rank is one has the other cut rank at most
`k^2`; equation (4) applies with parameter `h=k`.

When `b=1`, split the class into two complete subcases:

| Subcase | Selected rank interval | Parameter `h` for all rank-one pieces |
|---|---|---:|
| Below the square | `(k-1)^2+1<=a<=k^2-1` | `k-1` |
| At the square | `a=k^2` | `k` |

In the first row, the selected ratio is the actual maximum and is
strictly below `k^2`. Thus every other piece of ranks `(z,1)` or
`(1,z)` has the integer rank `z<k^2`, so its smaller factor rank
is at most `k-1`. This is the reason the stronger value is valid
for every rank-one subset encountered below.

For `b>=2`, use the unsplit interval

```math
a_0=(k-1)^2b+1,\qquad a_1=\min(k^2b,d^2),\qquad h=k.              (6)
```

## 3. Class-dependent weighted budgets

Let `F_form` be the geometric formal complement function in the
preceding note. For the current class parameter `h`, put

```math
F_h(d,m)=
\begin{cases}
d^2,&m=0,\\
K_d(h),&m=1,\\
F_{\rm geo}(d,m),&m\ge2.
\end{cases}                                                    (7)
```

Why may the one-product entry depend on the class? A subset whose
sum of opposite-cut ranks is one contains exactly one nonzero
piece. Section 2 bounds the factor rank of that piece by `h`, so
(4) applies. For larger sums, use the universal geometric bound.

For `s<k` define

```math
G_h(d,s;b)=\min_{b\le m\le\max(b,2d-1)}
                 \bigl(F_h(d,m)+s^2m\bigr).          (8)
```

The finite range remains complete. The entries for `m>=2d-1`
are zero, and the linear term is nondecreasing there.

For the minimizing subset of `sum_i min(r_i,s^2c_i)`, the
selected piece belongs strictly to the subset, because
`a>(k-1)^2b>=s^2b`. Therefore its product count is at least `b`.
The opposite weighted cut has the ordinary minimum starting at
zero. Consequently

```math
\sum_i\min(r_i,s^2c_i)\ge G_h(d,s;b),\qquad
\sum_i\min(s^2r_i,c_i)\ge G_h(d,s;0).                (9)
```

Both sides use the class-dependent rank-one entry. Applying it
only to the selected piece's residual would discard a valid
strengthening.

## 4. Residual ranks and the rational class bound

For a class interval `[a_0,a_1]`, put `n=d^2` and

```math
f=\begin{cases}K_d(h),&b=1,\\F_{\rm geo}(d,b),&b\ge2,\end{cases}
\quad R=\max(n-a_1,f),\quad C=n-b,
\quad e=\max(0,a_0+f-n).                             (10)
```

The remaining row- and column-rank totals are at least `R,C`,
and the total row-rank sum is at least `n+e`. The usual rank
budgets give the other total at least `n`.

For `1<=s<k`, define

```math
\mu_s=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_s=\frac{k^2-s^2}{k^4-1}.                         (11)
```

Concavity of `min(t,s^2)` on `[1/k^2,k^2]` gives, for every
remaining piece,

```math
\min(r,s^2c)\ge\mu_s r+\nu_s c,\qquad
\min(s^2r,c)\ge\nu_s r+\mu_s c.                     (12)
```

The selected piece contributes exactly `s^2b` and `b` to these
two weighted cuts. Combining (9)--(12) with the positive rational
interpolant `Q_k` from the preceding note gives

```math
S\ge J(d,k,b,a_0,a_1,h)
:=\frac{2n+e}{2k-1}
 +\sum_{s=1}^{k-1}w_s\left(
   \max\{G_h(d,s;b),s^2b+\mu_sR+\nu_sC\}
  +\max\{G_h(d,s;0),b+\nu_sR+\mu_sC\}\right).         (13)
```

Every quantity in (13) is rational. The proof uses only generic
cut ranks, formal complement bounds, and the actual maximum ratio.

## 5. Complete finite replay and stopping point

The checker inherits and rechecks the rational bound through
sixty-four labels. It then enumerates every class in (5) for
`d=65,...,71`, including both rank-one subcases and the equal-rank
case. It verifies each value in (13) against the exact upper
exponent `577/408`, using `577^2-2*408^2=1`.

At seventy-two labels, the explicit class
`k=d=72,b=1,a_0=a_1=d^2,h=d` gives a value below
`72^sqrt(2)`. This comparison is also exact, using the lower
Pell exponent `1393/985`. It records the limit of this certificate,
not a realizable seed. No numerical optimizer or sampled rank
profile is part of the proof.

Run

```sh
python tools/check_rank_one_formal_cutoff.py \
  --output evidence/rank-one-formal-cutoff.json
```

The checker also verifies (3) by integer maximization, tests the
class partitions, and reuses moving-lattice controls with explicit
rank-deficient first factors. The stronger exact-decomposition
cutoff remains a separate result.

The construction, proof development, and replay used substantial OpenAI assistant support.
