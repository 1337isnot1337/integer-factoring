# Partial transpose bounds the early ranks of nilpotent crossings

## Results and scope

Let `F` be any field. For an invertible operator `U` on `Mat_q(F)`, put
`n=q^2` and use the crossing reshuffling

\[
\Phi(U)_{(x,y),(u,v)}=U(E_{vy})_{xu}.
\]

Write `r_k=rank Phi(U^k)`. Then, for any integers `i,j`,

\[
\boxed{\ n\le r_i r_j r_{i+j}.\ }
\tag{1}
\]

In particular, if `r=rank Phi(U)`,

\[
\operatorname{rank}\Phi(U^2)\ge\left\lceil\frac{n}{r^2}\right\rceil.
\tag{2}
\]

For the nilpotent resolvent crossing with `S=I-U`, this gives the following
lower bound on its existing sufficient Smith-profile cost:

\[
\boxed{\ B(s)\ge nr+\operatorname{rank}(S)
 \max\left(0,\left\lceil\frac{n}{r^2}\right\rceil-2r\right),
 \qquad s\ge1.\ }
\tag{3}
\]

It applies to every nilpotent Jordan type. In the regular case,
`rank(S)=n-1`. Every regular nilpotent source with Choi rank at most two and
physical size `q>=5` therefore fails the Smith sufficient threshold
`B<q^(2 sqrt(2))`. This conclusion is about this local resolvent and its
stated sufficient bound. It is not a lower bound for arbitrary factoring
algorithms or arbitrary recombinations of crossing coefficients.

The underlying complex-matrix rank inequality is Theorem 9 of
[Song, Chen, Sun and Hu, arXiv:2105.07679](https://arxiv.org/abs/2105.07679).
The proof below gives its arbitrary-field version explicitly. It uses an
exterior-product argument for the block staircase and induction on ordinary
rank, so it does not require complex conjugation, positivity, or a
characteristic-zero specialization.

## 1. The arbitrary-field partial-transpose inequality

For a block matrix `M` with `a` by `b` blocks, each of size `c` by `d`, let
`K=OSR(M)` be the dimension of the span of its inner blocks. Equivalently,
`K` is the least number of terms in a decomposition
`M=sum A_t tensor B_t`. Let `Gamma` transpose every inner block. We prove

\[
\operatorname{rank}(M^\Gamma)\le K\operatorname{rank}(M).
\tag{4}
\]

Extending the field does not change any of these three ranks: all are ranks
of explicitly reshuffled matrices, detected by minors. We may consequently
work over an infinite extension of `F`. A rational-function extension is
enough. All coordinate changes below are invertible, and preserve both the
ordinary rank and the partial-transpose rank.

### 1.1 A block-space staircase

Let `W` be the span of all inner blocks. Choose a linear combination of the
outer block columns whose entries span a subspace of the largest possible
dimension `k_1`. Make it the first outer column. Outer row operations put a
basis of its block span in the first `k_1` positions and zero blocks below.
Call that span `W_1`.

Every block below outer row `k_1` now belongs to `W_1`. To see this, suppose
another outer column has a lower block `b` outside `W_1`. Mix that column into
the first with scalar `z`. If the original independent upper blocks are
`a_1,...,a_(k_1)`, the exterior product of their perturbed versions and the
new lower block has coefficient

\[
z\,a_1\wedge\cdots\wedge a_{k_1}\wedge b
\]

at its first nonzero degree. This coefficient is nonzero in every
characteristic. Choosing `z` outside the roots of this nonzero polynomial
increases the block-span dimension, contradicting maximality.

Repeat this construction in `W/W_1`, using the remaining outer columns and
row operations on only the first `k_1` outer rows. It produces `k_2<=k_1`
new independent blocks in outer column two. Continuing gives

\[
W_0=0\subset W_1\subset\cdots\subset W_p=W,
\qquad k_s=\dim(W_s/W_{s-1}),
\qquad k_1\ge\cdots\ge k_p>0,
\qquad\sum_s k_s=K.
\tag{5}
\]

The distinguished new blocks at stage `s` occupy the first `k_s` outer rows
of outer column `s`. In that column, every lower block lies in `W_(s-1)`.
In every later column, blocks below outer row `k_s` lie in `W_s`. Previous
stages remain valid because each subsequent row operation acts only within
their distinguished upper rows, and subsequent column operations leave
their selected columns unchanged. Set `k_0=a`.

### 1.2 Common input kernels and column strips

For each `s`, define

\[
D_s=\bigcap_{B\in W_s}\ker B\subseteq F^d,
\qquad D_0=F^d.
\]

Choose one basis adapted to the descending flag `D_0 superset D_1 superset
... superset D_p`. Write `I_s` for a basis of a complement of `D_s` inside
`D_(s-1)`, and put `h_s=|I_s|`. Every block in `W_(s-1)` vanishes on `I_s`.
The `k_s` distinguished new blocks, restricted to `I_s` and stacked
vertically, have column rank exactly `h_s`: a vector in their common kernel
would belong both to `D_s` and to its chosen complement.

Let `N` be the matrix after these outer and inner coordinate changes. Let
`N_s` retain only the inner input columns in `I_s`, setting the other inner
columns to zero in every block. Then

\[
N=\sum_s N_s,\qquad \operatorname{rank}N_s\le\operatorname{rank}N.
\tag{6}
\]

The final common kernel `D_p` contributes only zero columns. If `h_s=0`,
then `N_s=0` and the term can be omitted.

After deleting zero rows and columns, the support of a nonzero strip is

\[
N_s=\begin{pmatrix}A_s&C_s\\0&\omega_s\end{pmatrix}.
\tag{7}
\]

Here `A_s` is the vertical stack of the first `k_s` blocks in outer column
`s`, restricted to `I_s`, and has column rank `h_s`. The lower-right block
matrix `omega_s` occupies outer rows `k_s+1,...,k_(s-1)` and outer columns
`s+1,...,b`. Its inner blocks have size `c` by `h_s`. By the staircase, their
span is contained in the restriction of `W_s/W_(s-1)`, so

\[
\operatorname{OSR}(\omega_s)\le k_s,
\qquad h_s+\operatorname{rank}\omega_s\le\operatorname{rank}N_s.
\tag{8}
\]

After partial transpose, the upper part has only `k_s h_s` potentially
nonzero rows. Thus

\[
\operatorname{rank}(N_s^\Gamma)
 \le k_s h_s+\operatorname{rank}(\omega_s^\Gamma).
\tag{9}
\]

### 1.3 Induction on ordinary matrix rank

The zero-rank case is immediate. Suppose (4) holds for every block matrix
of ordinary rank less than `R`, with arbitrary block dimensions, and let
`rank M=R`. For every retained strip, `h_s>0`, and (8) gives
`rank omega_s<R`. Applying the induction hypothesis to (9), then (8), gives

\[
\begin{aligned}
\operatorname{rank}(N_s^\Gamma)
&\le k_s h_s+\operatorname{OSR}(\omega_s)\operatorname{rank}\omega_s\\
&\le k_s\bigl(h_s+\operatorname{rank}\omega_s\bigr)\\
&\le k_s\operatorname{rank}N_s\le k_sR.
\end{aligned}
\]

Sum over the strips and use `sum k_s=K` to obtain (4). This induction also
covers a staircase with only one stage; no separate full-column case is
needed. Since extension leaves the relevant ranks unchanged, (4) holds over
the original field.

## 2. Exact reshuffling and power bounds

Let `M(V)` denote the ordinary row-major matrix of a linear operator `V`,
and let `T` be the permutation that interchanges the two input indices.
The definitions give the entrywise identity

\[
\boxed{\ \Phi(V)=M(V)^{\Gamma_2}T.\ }
\tag{10}
\]

For example, at row `(x,y)` and column `(u,v)`, its right-hand side is
`M(V)_(x,u),(v,y)`, exactly the defining entry of `Phi(V)`.

A rank decomposition of `Phi(U)` is equivalent to

\[
U(X)=\sum_{a=1}^{r}P_a X^T Q_a.
\]

Composing two such decompositions gives

\[
U^2(X)=\sum_{a,b}P_aQ_b^T X P_b^TQ_a.
\tag{11}
\]

Consequently the ordinary matrix of `U^2` has operator Schmidt rank at most
`r^2`. Apply (4) to `M(U^2)^(Gamma_2)`, use the invariance of operator Schmidt
rank under partial transpose, and use invertibility:

\[
n=\operatorname{rank}M(U^2)
 \le r^2\operatorname{rank}\Phi(U^2).
\]

Composing decompositions of `U^i` and `U^j` proves (1) in the same way.

The square exponent on `r` cannot be removed in general. For `q=4`, rotate
the four binary matrix-unit label bits by one position. The resulting
invertible permutation operator satisfies `U^2=T` and has
`(rank Phi(U),rank Phi(U^2))=(4,1)`. Thus `r^2 r_2=n=16`, whereas the stronger
claim `r r_2>=n` would be false. This example is not a regular 16-dimensional
nilpotent source.

## 3. Consequences for every nilpotent Jordan type

Let `S=I-U` be nilpotent of index `m`, with Jordan block sizes `b_1,...,b_t`.
Use the [nilpotent resolvent construction](nilpotent-flag-crossings.md). If
`mu_i` are its row orders at infinity and `nu_j` its column orders at zero,
the sufficient profile cost is

\[
B(s)=\left(\sum_i s^{m-1-\mu_i}\right)
     \left(\sum_j s^{-\nu_j}\right),\qquad s\ge1.
\]
Its reversed row matrix is

\[
H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr)
 =\Phi(U)+w\Phi(U-U^2)+O(w^2).
\]

There are exactly `r` zero row-Smith orders. After their constant pivots are
removed, the first Schur coefficient is the compression of `-Phi(U^2)` from
`ker Phi(U)` to `coker Phi(U)`. Removing an `r`-dimensional domain quotient
and an `r`-dimensional target subspace lowers rank by at most `2r`. Hence
there are at least

\[
k=\max\left(0,\left\lceil n/r^2\right\rceil-2r\right)
\tag{12}
\]

row-Smith orders equal to one.

The column valuation sum is

\[
V(s)=\sum_{j=1}^t
 \left(s^{-(m-b_j)}+(b_j-1)s^{-m}\right).
\]

A row order zero contributes to `B(s)` at least

\[
s^{m-1}V(s)=\sum_j\left(s^{b_j-1}+(b_j-1)s^{-1}\right)\ge\sum_j b_j=n
\]

by arithmetic-geometric mean, including `b_j=1`. A row order one contributes
at least

\[
s^{m-2}V(s)=\sum_j\left(s^{b_j-2}+(b_j-1)s^{-2}\right)
 \ge\sum_j(b_j-1)=\operatorname{rank}S.
\]

For completeness, the scalar inequality used here is
`x^(b-2)+(b-1)x^(-2)>=b-1` for `x>=1` and `b>=2`.
The case `b=2` is immediate. For `a=b-2>0`, evaluating the unique stationary
point reduces the assertion to

\[
\left(1+\frac2a\right)^{(a+2)/2}\ge 2+\frac2a,
\]

which follows from Bernoulli's inequality. A block of size one contributes a
nonnegative term. All other row orders make positive contributions as well.
This proves (3).

### Regular rank-two sources of physical size at least five

For a regular block, `rank S=n-1`. If `r=2` and `n>=25`, (3) implies

\[
\begin{aligned}
B(s)&\ge2n+(n-1)\left(\frac n4-4\right)\\
&=\frac{n^2}4-\frac{9n}4+4
 \ge \frac{104}{625}n^2>n^{\sqrt2}.
\end{aligned}
\]

The penultimate inequality follows because
`1/4-9/(4n)+4/n^2` increases for `n>=25`. For the last inequality, it suffices
to check the endpoint `104>25^sqrt(2)` and use `2-sqrt(2)>0`. A purely integer
certificate is `104^7>25^10`, together with `sqrt(2)<10/7`.

Rank one is even more restrictive: (3) gives
`B>=n+(n-1)(n-2)`. Thus the conclusion applies to Choi rank at most two.
For `q=8,r=2`, the sharper integer floor is `B>=128+63*12=884`, compared with
the success threshold `64^sqrt(2)`, approximately 358.36.

The bound does not exclude regular rank-two sources at `q=3` or `q=4` by
itself. It also leaves larger Choi ranks open. At `q=8`, a regular source
passing this first test can have `r=3,4,5`; the respective maximum permitted
numbers of order-one pivots are two, one, and zero. Further row orders and
the complete construction still have to be checked.

## 4. Exact verification

```sh
python tools/check_partial_transpose_rank.py --output /tmp/partial-transpose-rank.json
```

The checker uses only the Python standard library. It verifies (4), (10), and
the square-rank bound on all 65,536 binary operators on `Mat_2`, including all
20,160 invertible maps and all 2,520 regular unipotent maps. It also checks
deterministic finite-field controls over characteristics two, three, and five;
the sharp four-bit square example; and the exact regular-rank-two threshold
certificates. These finite checks are separate from the all-field proof.
