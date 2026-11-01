# A protected lower-space obstruction for coherent boundary retractions

## Result and scope

**Theorem.** In the eight-mode reflection construction, retain all 28
coordinate states of height below 16. On the entire five-dimensional
height-16 cell, allow an arbitrary homogeneous embedding and retraction
of any rank, over any field of characteristic two. The natural
common-scale two-cut Smith pressure, optimized over every common physical
type, has infimum at least

\[
\frac{17}{24}>\frac1{\sqrt2}.
\tag{1}
\]

Thus this whole boundary family cannot reach the strict profile-pressure
threshold needed by the present factoring transfer. The conclusion
includes the binary 32-letter optimum, every rank-three `F4` boundary
retraction, and arbitrary extensions of the coefficient field. It also
allows mixtures between the four two-particle height-16 states and the
one four-particle state.

The proof protects a fixed polynomial submatrix. It therefore avoids
enumerating all boundary projectors. It does not bound actual coefficient
ranks obtained by another method, additional physical gauges, unequal cut
scales, or spaces that remove or coherently change the lower 28 coordinates.
The useful next constructions must escape at least one of these hypotheses.
No stronger factoring exponent is established.

The finite input is the complete 784-order natural Smith profile of the
28-letter coordinate space. The original binary compound constructor and
the independently normalized field-valued constructor agree on both cuts,
checking 3,136 complete orders. Both use the existing native Smith engine;
this is not a claim of two independent Smith implementations. After that
finite profile, the all-scale certificate uses integer arithmetic only.

## 1. The fixed submatrix survives every boundary retraction

For an occupation word `a` on modes `0,...,7`, write

\[
h(a)=\sum_{i\in a}(2i+1),\qquad
U=\operatorname{span}\{e_a:h(a)<16\}.
\]

There are 28 such words. The full height-16 cell is

\[
W=\operatorname{span}(e_{15},e_{24},e_{36},e_{66},e_{129}).
\]

Choose any `P_W:K^r -> W` and `L_W:W -> K^r` with `L_W P_W=I_r`,
where `0<=r<=5`. Set

\[
P=I_U\oplus P_W,\qquad L=I_U\oplus L_W.
\tag{2}
\]

Every nonidentity bit deletion lowers height. Thus both map spaces are
deletion invariant, and the
[coherent compound reduction](coherent-hot-compound-retractions.md)
applies. Its algebraic reduction does not require occupation-homogeneous
labels; that additional grading is only used to speed up evaluation.
The two regular cut matrices are

\[
M_0=(L\otimes L)E_B(P\otimes P),\qquad
M_\infty=(L\otimes P^{\mathsf T})E_B(L^{\mathsf T}\otimes P).
\tag{3}
\]

Restricting every row and column label of either matrix to `U` yields
the same fixed 784-by-784 coordinate compound matrix `N`. Row and column
permutations used to express its crossing endpoint do not affect this
statement or its Smith orders. Hence `N` is a polynomial submatrix of
each matrix in (3), for every choice in (2).

All these matrices are generically invertible: their endpoint at `t=1`
is the crossing permutation because `LP=I`. The regular polynomial
matrices have nonnegative Smith orders. Opposite scalar shifts in the
original two tensor cuts cancel in their product counting function, so
the regular matrices are the appropriate objects for the natural pressure.

## 2. Smith counting functions increase under taking a larger matrix

**Lemma.** Let `N` be a submatrix of a polynomial matrix `M` over a field,
and let their finite Smith orders at zero be `n_i` and `m_i`. For every
`0<r<1`,

\[
\sum_i r^{m_i}\ \ge\ \sum_i r^{n_i}.
\tag{4}
\]

**Proof.** Multiplication by `M` modulo `t^k` is represented by a finite
block Toeplitz matrix. Multiplication by `N` modulo the same power is its
row-and-column submatrix. Their ranks therefore satisfy `R_k(M)>=R_k(N)`.
Invertible power-series row and column operations preserve these truncated
ranks, so the Smith normal form gives

\[
R_k(M)=\sum_i(k-m_i)_+.
\]

For a nonnegative integer `m`, the convergent geometric-series identity is

\[
\sum_{k\ge1}(k-m)_+r^{k-1}=\frac{r^m}{(1-r)^2}.
\]

Multiply each truncated-rank inequality by `r^(k-1)`, sum, and cancel
the positive factor `(1-r)^(-2)`. This proves (4). The proof also allows
rectangular matrices and zero Smith factors, which contribute zero. No
unproved componentwise interlacing of individual Smith orders is used.

Let

\[
Z(r)=\sum_{i=1}^{784}r^{n_i}
\]

for the protected matrix `N`. By (4), the two complete counting functions
of (3) obey

\[
Z_0(r)Z_\infty(r)\ge Z(r)^2.
\tag{5}
\]

The entries of `N` lie in `F_2[t]`. Its profile is unchanged under every
characteristic-two field extension: ranks of its finite Toeplitz matrices
are unchanged. Thus the same fixed `Z` proves (5) over all these fields.

## 3. Entropy comparison with the largest physical alphabet

The largest permitted physical height polynomial is

\[
\begin{aligned}
A_{33}(t)={}&1+t+t^3+t^4+t^5+t^6+t^7
  +2t^8+2t^9+2t^{10}+2t^{11}\\
&+3t^{12}+3t^{13}+3t^{14}+4t^{15}+5t^{16}.
\end{aligned}
\tag{6}
\]

Any choice of rank `r` in (2) has its own polynomial `A_q(t)<=A_33(t)`
for `0<t<=1`. Put `rho_0=17/24`. The exact certificate below proves

\[
Z(s^{17})^{12}>A_{33}(s^{24})^{17}
\qquad(0<s\le1).
\tag{7}
\]

Set `r=s^17`. Equations (5) and (7) give

\[
\log Z_0(r)+\log Z_\infty(r)
 >4\rho_0\log A_{33}(r^{1/\rho_0}).
\tag{8}
\]

For any common physical probability vector `p` with positive entropy,
let `H(p)=-sum p_j log p_j` and `E(p)=sum p_j h_j`. The elementary
Gibbs entropy inequality gives

\[
H(p)+\frac{\log r}{\rho_0}E(p)
\le\log A_q(r^{1/\rho_0})
\le\log A_{33}(r^{1/\rho_0}).
\]

Combining this with (8) proves the pointwise bound

\[
\frac{\log Z_0(r)+\log Z_\infty(r)-4\log r\,E(p)}{4H(p)}
>\rho_0.
\tag{9}
\]

The infimum over scales and types is consequently at least `rho_0`.
We assert a non-strict lower bound on that infimum; pointwise strictness
alone need not supply a uniform gap. Finally,

\[
2\cdot17^2-24^2=2>0,
\]

so this non-strict infimum bound already strictly exceeds the required
threshold. The entropy accounting retains every physical label, including
all five possible boundary directions.

## 4. Exact proof of the scalar inequality

The [checker](../tools/check_coherent_boundary_protected_floor.py) stores
the complete coefficient list of `Z`. It has degree 82, total coefficient
sum 784, and initial terms `1+r+O(r^3)`. The retained
[certificate](../evidence/coherent-boundary-protected-floor.json) records
the profile, both exact constructor comparisons, and every interval used
below.

For `0<s<=4/5`, put `t=s^24`. From (6),

\[
A_{33}(t)\le1+t+31t^3\le1+2t,
\]

because `31*(4/5)^48<1`. The coefficient inequalities
`binom(17,j)*2^j <= 34*32^(j-1)` imply

\[
(1+2s^{24})^{17}-1
\le \frac{34s^{24}}{1-32s^{24}}.
\]

On the other side, `Z(s^17)^12>=1+12s^17`. Exact integer comparisons give

\[
\begin{split}
32\cdot4^{24}&<5^{24},\\
12(5^{24}-32\cdot4^{24})-34\cdot4^7\cdot5^{17}
&=182169346247795596>0.
\end{split}
\]

They imply `12s^17 > 34s^24/(1-32s^24)` throughout this interval, proving
(7) there.

On `[4/5,1]`, both sides of (7) are increasing functions. The checker
constructs 142 rational intervals with maximum dyadic depth 11. On each
interval `[a,b]`, it compares the left side at `a` with the right side at
`b` using integers only. It also verifies that the intervals have no gaps
or overlaps and cover the whole range. The largest polynomial degrees
before denominator clearing are 16,728 and 6,528, so the denominator
power difference is exactly 10,200. All 142 comparisons are strict.

Run the independent profile reconstruction and the integer certificate:

```sh
python3 tools/check_coherent_boundary_protected_floor.py \
  --recompute-profiles --output evidence/coherent-boundary-protected-floor.json
```

Omitting `--recompute-profiles` replays only the standard-library scalar
certificate. The mathematical submatrix and entropy arguments are given
above, rather than replaced by finite sampling. The construction builds
on the existing coherent compound reduction and its native exact profile
engine. This work had substantial OpenAI assistant support and introduces
no Lean formalization.
