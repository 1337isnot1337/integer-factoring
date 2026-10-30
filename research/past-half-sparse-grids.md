# Sparse grids and the square-root logarithm coefficient

## Question and result

A fully occupied square grid makes the parent have about the square of the
child's vertex count. Could leaving grid cells empty give more recursion
levels, while retaining the Appendix B crossing saving?

For regular supports and the same row/column crossing operation, it cannot.
Fewer occupied cells also mean fewer complete rectangles. The exact count
below gives a coefficient bound for **arbitrary changing densities**, rather
than only a fixed power schedule. It concerns the declared channels and a
scalar cost recurrence. It is not a lower bound on matrix ranks, the runtime
of a factoring machine, or all possible uses of sparse networks.

The crossing pieces and their declared factorizations are from Claude
(Anthropic), *A 2^o(n) Time Algorithm for 3-SAT*, at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The source's grid operation is in Section 2.5 and Appendix A of the
[paper](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
Its half-piece factors are equations (54)--(56) in Appendix B, also stated in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The counting and variable-density argument below are new ordinary proofs.
There is no new Lean formalization or improved factoring exponent here.

## 1. How many crossings survive missing cells?

Take a subset `S` of a rectangular grid, with one clique vertex in each
occupied cell. Write `N=|S|`, and let `R` be the number of complete
axis-aligned rectangles whose four corners lie in `S`. Suppose every row
has at most `d_r` occupied cells and every column at most `d_c`.
Then

```math
4R\le N(d_r-1)(d_c-1).
\tag{1}
```

For each occupied corner, choose the other corner in its row and the other
corner in its column. There are at most `(d_r-1)(d_c-1)` such choices.
A rectangle requires the fourth corner to be occupied too. Each actual
rectangle is counted at its four corners. This proves (1).

Two occupied cells in different rows and columns determine exactly one
potential rectangle. Their clique edge participates in the grid crossing
operation precisely when the other two corners are occupied. Each complete
rectangle pairs its two diagonal edges, so exactly `2R` clique edges are
processed by crossings. A missing fourth corner cannot be replaced by an
empty group: the diagonal endpoint tensor would then be absent.

For a `K` by `K` support with exactly `d` occupied cells in every row and
column, `N=Kd`, and (1) becomes

```math
R\le\frac{N(d-1)^2}{4}.
\tag{2}
```

If `d` divides `K`, a disjoint union of `K/d` complete `d` by `d` blocks
attains equality. Thus the rectangle bound cannot be improved for all
regular supports.

## 2. The exact declared volume of a balanced piece profile

Give every original clique edge `lambda` bits. Its directed channel volume
is `V=N(N-1)lambda`. Replace the diagonal pair in every complete rectangle
by a piece whose declared row and column widths are the same integer `c`,
where `0<=c<=2lambda`. Retain all other edges.

Let `D` count the edges joining different rows and columns whose opposite
corners are missing. There are `K*binom(d,2)` edges within rows and the same
number within columns. Therefore

```math
\binom N2=2K\binom d2+D+2R.
```

After merging rows, row-side edges disappear, column-side edges remain,
and each processed rectangle contributes one `c`-bit channel. The resulting
volume is

```math
\begin{aligned}
U
 &=2\lambda\left(D+K\binom d2\right)+2Rc\\
 &=N(N-1)\lambda-2\lambda K\binom d2-2R(2\lambda-c).
\end{aligned}
\tag{3}
```

Merging columns gives exactly the same number. In particular, choosing the
cheaper orientation does not improve this profile.

Put `r=c/(2lambda)`. By (2), equation (3) gives

```math
\frac UV\ge
1-\frac{d-1}{Kd-1}
 -(1-r)\frac{(d-1)^2}{Kd-1}.
\tag{4}
```

This statement is exact for the declared factorization and finite integers.
It does not assert that either declared channel count is the minimum rank
of the resulting tensor.

For the source's central piece `s=0`, both declared channel counts are
`R_0=sum_(j=0)^(2h)p(j)`. The source's choice is the least `h` with
`p(h)>=2^lambda`. Its partition estimates imply

```math
\frac{\min\{2\lambda,\lceil\log_2R_0\rceil\}}{2\lambda}
\longrightarrow\rho=\frac1{\sqrt2}
\quad\text{as }\lambda\longrightarrow\infty.
```

The optional trivial-width cap does not change the limit. Thus (4) gives
the following lower bound for the leading normalized coefficient of the
all-central profile:

```math
\eta(K,d)=
1-\frac{d-1}{Kd-1}
 -(1-\rho)\frac{(d-1)^2}{Kd-1}.
\tag{5}
```

Taking the leading coefficient means that the limit in the label width is
taken with this finite layout fixed. This is the coefficient relevant to
a bound uniform over every label width. No finite-width rank lower bound
or unproved error summation is being used.

## 3. Arbitrarily changing the density does not improve the coefficient rate

Consider an idealized recursion that pays no grouping or construction
costs. At step `j`, let its child have `K_j` vertices and let its parent
use a `d_j`-regular `K_j` by `K_j` support:

```math
K_{j+1}=K_jd_j,\qquad 2\le d_j\le K_j,\qquad K_0\ge4.
```

Suppose the scalar leading coefficient `e_j>0` is propagated by the
same declared-channel estimate, including the all-central profile. It must
therefore obey

```math
e_{j+1}\ge\eta(K_j,d_j)e_j.
\tag{6}
```

Then for every number of steps `ell`,

```math
e_\ell\sqrt{\log K_\ell}
\ge
\left(1-\frac{2\sqrt2}{K_0}\right)
 e_0\sqrt{\log K_0}>0.
\tag{7}
```

This rules out an improved power of the logarithm within (6), even when
`d_j` is chosen adaptively and approaches full density at an arbitrary rate.
It also rules out any unbounded extra logarithmic saving in this scalar
coefficient. It does not rule out a sharper analysis of the actual program.

### Proof of (7)

Fix `K>=4` and `2<=d<=K`, and put

```math
t=\frac{d-1}{K-1},\qquad
b=\frac{\log(Kd)}{\log K}=1+\frac{\log d}{\log K}.
```

The elementary rational bounds

```math
\frac{d-1}{Kd-1}\le\frac1K,\qquad
\frac{(d-1)^2}{Kd-1}\le\frac{d-1}{K-1}=t
```

give `eta(K,d)>=1-(1-rho)t-1/K`. Concavity of `log` on `[1,K]` gives
`log d>=t log K`, hence `b>=1+t`. The convex function
`f(t)=(1+t)^(-1/2)` lies below its chord from zero to one, so

```math
\eta(K,d)
\ge1-(1-\rho)t-\frac1K
\ge\frac1{\sqrt{1+t}}-\frac1K
\ge\frac1{\sqrt b}-\frac1K.
\tag{8}
```

As `b<=2`, multiplying (8) by `sqrt(b)` and applying (6) gives

```math
e_{j+1}\sqrt{\log K_{j+1}}
\ge\left(1-\frac{\sqrt2}{K_j}\right)
 e_j\sqrt{\log K_j}.
```

Also `K_j>=2^jK_0`. Every factor is positive, and the elementary finite
product inequality `prod(1-x_j)>=1-sum x_j`, for `0<=x_j<=1`, yields

```math
\prod_{j<\ell}\left(1-\frac{\sqrt2}{K_j}\right)
\ge1-\sum_{j<\ell}\frac{\sqrt2}{2^jK_0}
\ge1-\frac{2\sqrt2}{K_0}.
```

This proves (7). If a later grouping reduces the child below `K_j`, a
size-only bound of the proposed form `C/(log K)^alpha` has a larger child
coefficient, so that extra grouping does not repair the failed estimate.
Grouping overheads cannot improve this idealized calculation either.

## 4. A fixed sparse power loses even the geometric contraction

For intuition, take `d` of order `K^(b-1)` with a fixed `1<b<2`. Then the
parent has order `K^b` vertices, but (5) gives

```math
1-\eta(K,d)=O(K^{b-2}+K^{-1}).
```

Along a repeated fixed-power schedule, these losses are summable because
`K_j` grows doubly exponentially. Their product stays positive. Thus the
resulting coefficient need not even tend to zero under this estimate.
The larger number of recursive levels has been bought by making each
crossing layer asymptotically ineffective. The variable-density theorem
(7) gives the stronger conclusion for all the regular schedules at once.

## 5. Finite checks and scope of the remaining opening

Run:

```sh
python3 tools/check_sparse_grid_geometry.py
```

The checker enumerates every binary support through grid side four and
retains all 156 nonempty regular supports. It checks the exact rectangle
count and the channel-volume identities for 3,744 integer choices of
layout, old width, and balanced new width. It checks equality in (2) for
35 unions of complete blocks through side twelve. A three-corner example
checks that an incomplete rectangle contributes no crossing.

The output is recorded in
[sparse-grid-geometry.json](../evidence/sparse-grid-geometry.json).
These finite checks support the transcription of (1)--(4); the proof of
(7) is the ordinary argument above.

The statement is limited to regular grid supports, the original complete-
rectangle operation, the declared central channels, and a uniform scalar
leading coefficient that handles every piece profile. It leaves open a
new multi-orientation construction, a proof that relevant profiles vanish,
compression of the final pieces, or an input-dependent analysis using
properties that are lost in a scalar coefficient. It establishes no new
factoring bound.
