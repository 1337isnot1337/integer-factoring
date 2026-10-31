# Compatible compression and rational cuts exclude formal seeds through 73

## Result and scope

Let `F` be any field and suppose a supplied finite formal decomposition has

\[
T(t)=\sum_i A_i(t),\qquad T(t)\in\operatorname{Mat}_{d^2}(F[[t]]),
\qquad T(0)=E_d.
\]

The individual pieces may have Laurent-series poles. If `r_i,c_i` are their
two generic cut ranks over `F((t))`, then

\[
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}\qquad(2\le d\le73).
\]

The preceding [PR27 factor-kernel theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/rank-one-formal-cutoff.md)
proved this through 71. This note adds dimensions 72 and 73 using the
[compatible integer compression cores](compatible-compression-budget.md)
and [rational scalar cut meshes](rational-cut-mesh-budget.md).

This is a formal finite-seed rank-sum exclusion. It does not improve the
factoring exponent, which remains `1/2`, and it does not exclude larger
completed crossings obtained from tensor powers of a smaller polynomial
source. The separate exact finite-seed cutoff remains 97.

## 1. Which complement bounds remain valid with poles

### Compatible geometric bounds

For a sum of `m` product operators, the compatible-core construction
supplies common input and output subspaces satisfying
`P_i U subset W`, with `dim U=s` and `codim W=q`. The quotient compression
annihilates the entire product sum and has swap rank `sq`.

To use this over Laurent series, first choose the subspaces over an
algebraic closure of `F((t))`. Their finitely many coordinates lie in a
finite extension `L`. Let `O` be its valuation ring. The saturated
lattices `U_L intersect O^d` and `W_L intersect O^d` are direct summands
of the free module `O^d`. Integral inclusion and quotient matrices
therefore retain ranks `s` and `q` after reduction.

The product sum is killed before specialization, including all its
poles. The compressed complementary piece equals the compression of the
integral tensor `T(t)`. Its reduction is the corresponding restriction
and quotient of the swap, with rank `sq`. Thus the generic rank of the
complementary piece is at least `sq`.

This proves that all compatible-core geometric bounds used below are
valid for formal limits. It is the same saturated-lattice argument as in
the preceding geometric theorem, with stronger integer Jacobian cores.

### The class-dependent one-product bound

If a rank-one opposite-cut piece is `B=P tensor Q` and the smaller factor
rank is `p`, choose an input subspace `U` with

\[
\dim P(U)=\max(0,s-d+p).
\]

The same moving-subspace compression gives

\[
r(T-B)\ge K_d(p):=
\max_{0\le s\le d}s\bigl(d-\max(0,s-d+p)\bigr)
=\left\lfloor\frac{(2d-p)^2}{4}\right\rfloor.
\tag{1}
\]

The generic row rank of `B` is `rank(P)rank(Q)`. Hence a piece of ranks
`(a,1)` has smaller factor rank at most `floor(sqrt(a))`, and the
complement bound is at least `K_d(floor(sqrt(a)))`.

Only these kernel and geometric arguments are used. In particular, the
exact-only Fitting bound and the exact one/two-product half-rank bound
are not imported. At `m=2`, the formal table uses
`floor(d^2/4)`, the geometric two-map bound.

## 2. Complete classes and their formal budgets

Choose an actual piece attaining the largest rank ratio or reciprocal,
exchange cuts so its ranks are `(a,b)` with `a>=b`, and let
`k=ceil(sqrt(a/b))`. The equal-rank case has cost at least `d^2`.
The remaining classes are

\[
2\le k\le d,\qquad
1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor.
\]

For `b>=2`, set

\[
a_0=(k-1)^2b+1,\qquad a_1=\min(k^2b,d^2),\qquad h=k.
\]

For `b=1`, retain the two complete PR27 subcases:

| Subcase | Selected rank interval | Factor-rank parameter `h` |
| --- | --- | --- |
| Below the square | `(k-1)^2+1<=a<=k^2-1` | `k-1` |
| At the square | `a=k^2` | `k` |

In the first subcase, the selected ratio is the actual maximum. Every
other rank-one piece therefore also has its larger rank strictly below
`k^2`, giving the stronger factor-rank parameter `h=k-1`.

For the current class, define the formal complement table

\[
F_h(m)=\begin{cases}
d^2,&m=0,\\
K_d(h),&m=1,\\
F_{\rm compatible}(d,m),&m\ge2.
\end{cases}
\tag{2}
\]

Here the last entry combines the old geometric allocations with the
verified denominator-three, four, and five integer cores. A subset
whose opposite-cut rank sum is one contains exactly one nonzero piece,
so its factor rank is controlled by the current class. This justifies
using the improved one-product entry throughout the weighted subset
budgets, not only for the chosen piece.

Define

\[
G_h(s;b)=\min_{b\le m\le\max(b,2d-1)}\{F_h(m)+s^2m\}.
\]

The tail is complete because the geometric tables vanish from `2d-1`
onward. Let `f=F_h(b)`, and put

\[
R=\max(d^2-a_1,f),\qquad C=d^2-b,
\qquad e=\max(0,a_0+f-d^2).
\]

For a rational cut threshold `1<=s<k`, use

\[
\mu_s=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_s=\frac{k^2-s^2}{k^4-1}.
\]

The selected piece is guaranteed to belong to the first minimizing
subset only when `a_0>s^2b`. Set `b_s=b` in that case and `b_s=0`
otherwise. The weighted cut bounds are

\[
U_s=\max\{G_h(s;b_s),\min(a_0,s^2b)+\mu_sR+\nu_sC\},
\]

\[
V_s=\max\{G_h(s;0),b+\nu_sR+\mu_sC\}.
\tag{3}
\]

The positive rational chord interpolant of the square root combines
(3) with the ordinary rank budget `sum(r_i+c_i)>=2d^2+e`. This is the
same all-field scalar inequality used for exact seeds; the complement
table and residual parameters here are entirely formal-safe.

## 3. Exact replay at 72 and 73

The checker covers all subclasses at both dimensions. It first tries
the nonnegative base term, then an integer cut mesh, then a quarter
mesh, and finally a mesh of width `1/32` where needed.

| Dimension | Scalar subclasses | Base term | Integer mesh | Quarter mesh | Mesh `1/32` | Certified common floor |
| --- | --- | --- | --- | --- | --- | --- |
| 72 | 8,487 | 8,079 | 407 | 1 | 0 | 424 |
| 73 | 8,735 | 8,311 | 421 | 2 | 1 | `21583/50=431.66` |

The sole finest-mesh case is the full-square rank-one class
`d=k=73,b=1,a=d^2,h=d`. Its certified lower bound is

\[
\frac{431661183669}{10^9}=431.661183669.
\]

Each weighted summand is rational. The checker rounds each summand
**down** to a multiple of `10^-9` by exact integer division. Thus the
recorded sum remains a rigorous lower bound; no floating-point
comparison enters the certificate.

The target comparisons use

\[
3363^2-2\cdot2378^2=1,
\]

\[
424^{2378}>72^{3363},\qquad
21583^{2378}>73^{3363}\,50^{2378}.
\]

They prove that both common floors strictly exceed `d^sqrt(2)`.
Together with the preceding all-field theorem through 71, this gives
the formal cutoff through 73.

## 4. Reproduction and boundary

```sh
python tools/check_compatible_formal_cutoff.py \
  --cores evidence/compatible-compression-cores.json \
  --output evidence/compatible-formal-cutoff.json
```

The checker uses the preceding compatible-core and rational-mesh
modules. It recomputes all 50 integer cross-block certificates,
verifies 2,700 integer kernel maxima and 132,348 factor-rank controls,
and covers all 17,222 new scalar subclasses. Both rank-one subcases
are counted and their intervals are checked to form complete
partitions. The replay takes about eight seconds in the recorded
workspace.

Exploratory calculations at dimension 74 remain below the target even
with the dimension-feasible geometric relaxation and a fine cut mesh.
No statement about the existence or impossibility of a useful formal
seed at 74 is made. The result above is the proved finite exclusion,
with no stronger factoring score claimed.
