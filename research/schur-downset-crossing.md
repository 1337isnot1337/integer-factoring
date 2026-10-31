# Sparse Boolean mode sets in the Schur crossing

## Result and scope

The complete Schur construction can use a Boolean down-set of fermion
indices and the complemented indices for the opposite species. The
resulting kernel has an explicit polynomial form and an exact crossing
endpoint on a fixed-energy label set.

The reduction differs from the
[full binary-interval Gaussian reduction](schur-gaussian-digit-reduction.md).
A proper down-set introduces a matrix-valued holonomy. Two matrices,
`H_S` and `F_S`, control its one-particle cut profiles. Their determinant
identities give asymptotic profile limits for fixed Hamming radius,
without a formula for every individual Smith order.

This note establishes the sparse family, its endpoint, effective finite
construction, and compact sufficient profiles. It does not establish an
improved factoring exponent. Its displayed degree gauge is a concrete
control; optimizing the profiles and the projected coefficient ranks
remains open.

The construction uses the complete **inhomogeneous** integral Schur dual
`T(s_x)=exp(L_1) omega(s_x)`. Its crossing multiplier is one at every
weight. The determinant formula is proved in
[the Schur query note](schur-fermion-queries.md), following
[PR #30 at `28105165`](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/integral-schur-crossing.md)
and the pinned Appendix B source cited there. The new sparse construction
and checks were developed with substantial OpenAI assistant support.

## 1. Labels from a down-set and its complement

Let `M=2^b`, put `T=M-1`, and let `S` be a nonempty down-set in the Boolean
lattice of masks `0,...,M-1`. Write `n=|S|` and `J(a)=T-a`.
Choose two `r`-subsets `A,B` of `S` with equal numerical sums. They specify
a Schur partition with arms `A` and legs `J(B)`. Its weight is

```math
\sum_{a\in A}a+\sum_{c\in B}(M-1-c)+r=Mr=:h.              \tag{1}
```

If `N_t` counts `r`-subsets of `S` with sum `t`, the number of external
labels is

```math
d=\sum_tN_t^2\ge\frac{\binom nr^2}{r(M-1)+1}.              \tag{2}
```

This is the same elementary Cauchy--Schwarz bound as for the full interval.
The complement in the leg indices makes equality of the two subset sums
the required energy condition.

## 2. Restricted matrices and the new holonomy

Work over `F_2(z)` and put `q=z+1`. All rows and columns below are indexed
by `S`. Define

```math
\begin{aligned}
P_{a,c}&=\mathbf1_{c\subseteq a},&
Q&=\mathrm{diag}(q^a),&
Q_J&=\mathrm{diag}(q^{T-a}),\\
G_{a,c}&=z^{a\mathbin{\&}c},&
H_{a,c}&=z^{T-(a\mathbin{|}c)},&
Z_{a,c}&=z^{a\mathbin{\&}\neg c}.
\end{aligned}                                                \tag{3}
```

Down-set closure gives the exact identities

```math
P^2=I,\qquad G=PQP^T,\qquad
H=z^T\mathrm{diag}(z^{-a})G\mathrm{diag}(z^{-a}). \tag{4}
```

For `G`, expand the product over the bits in `a intersection c`.
Every subset of this intersection lies in `S`, so the matrix product
in (4) retains every term. The identity for `P^2` similarly counts all
intermediate masks in a Boolean interval.

One must retain this closure hypothesis. For `S={0,3}` the claimed
factorization of `G` fails by the nonzero polynomial `z+z^2`.

The full interval also obeys `HG=q^(M-1)I`. A proper down-set need not.
For `M=4,S={0,1,2}`, the `(0,1)` entry of `HG` is `z`. Thus restricting
the full-interval scalar identities would lose a real matrix degree of
freedom.

More explicitly, a spanning-tree change of basis turns the rational
four-boundary kernel into the full-interval scalar pattern, except for
the `V,X` block. That block becomes

```math
q^{-1}G^{-1}H^{-1}.                                         \tag{5}
```

Equation (5) is the additional holonomy. No scalar replacement is made
in the construction or rank calculations below.

## 3. An explicit polynomial kernel

Use `z=w^2` and `s=w+1`, so `s^2=q`. Relabel every complemented index
`J(a)` by `a in S`. In boundary order `B,X,V,Y`, the polynomial
plus-to-minus kernel is

```math
\boxed{\quad
\mathcal A_S(w)=
\begin{pmatrix}
0&s z^M H&sP^TQ_J&I\\
s z^M G&0&z^M P&sG\\
sQP^T&z^M P&0&sQP^T\\
I&sH&sP^TQ_J&0
\end{pmatrix}.\quad}                                       \tag{6}
```

Every displayed entry is a polynomial in `w`. To obtain (6) from the
folded Schur kernel, multiply each `X` fermion by `z^M s`. Multiply a
`V` plus fermion of index `a` by `q^a s`, and a `V` minus fermion of
complemented index `J(c)` by `q^(T-c)s`.

For a selected label of rank `r` and weight `h=Mr`, these factors multiply
the determinant by

```math
z^{2Mr}q^{h+r}=z^{2h}q^{h+r},                               \tag{7}
```

exactly the prefactor in the complete Schur polynomial. Thus selecting
the occupied rows and columns of (6) gives `P_h(w^2)` directly. There is
no uncharged normalization or missing factorial.

The block occupations are

```math
\begin{array}{c|cc}
&\text{plus indices}&\text{minus indices}\\\hline
B&B_u&A_u\\
X&B_x&A_x\\
V&A_v&B_v\\
Y&A_y&B_y.
\end{array}                                                   \tag{8}
```

These conventions include the bra and conjugate-parameter species flips.

### Endpoint without a Gram-matrix assumption

At `w=1`, only `B,Y` and `X,V` survive. The first pairing is the identity.
The second is

```math
\det P_{B_x,B_v}\,\det P_{A_v,A_x}.                         \tag{9}
```

A nonzero minor of `P` requires a matching of column masks into row masks
by subset inclusion. Hence the first factor implies
`sum(B_v)<=sum(B_x)`, and the second implies
`sum(A_x)<=sum(A_v)`. The equal-sum condition of each label makes both
inequalities equalities. Each matched inclusion must then be equality,
so (9) is `delta_xv`.

Therefore the coefficient tensors of (6), restricted to the labels in
Section 1, sum exactly to `delta_xv delta_yu`. This proof applies to every
`r`; it does not require `Mr` to be a prime power.

If all pairs of `r`-subsets are retained without the energy condition,
the endpoint pairing in (9) remains invertible. Multiplying the `X` row
and column blocks by the constant matrix `P` normalizes it. That gives
a separate crossing on `binom(n,r)^2` occupation labels. Index-dependent
energy gauges below are used with the energy slice, where their shifts
are constant across labels.

## 4. Two small matrices control the cut profiles

Put

```math
L_0=I+qGH,\qquad L_\infty=z^M I+qGH,\qquad F=I+qZ^2.       \tag{10}
```

Write `A=sum_(a in S) a`. From (4),

```math
\det G=q^A,\qquad \det H=z^{n(M-1)-2A}q^A.                 \tag{11}
```

Both `G` and `q` are units at zero. Also

```math
z^{M-1}H^{-1}
=\mathrm{diag}(z^a)G^{-1}\mathrm{diag}(z^a)
\in\mathrm{Mat}_n(\mathbb F_2[[z]]).                  \tag{12}
```

If `nu_1,...,nu_n` are the zero Smith orders of `H`, equations
(11)--(12) prove

```math
0\le\nu_i\le M-1,\qquad
\sum_i\nu_i=n(M-1)-2A.                                     \tag{13}
```

The factorization

```math
L_\infty=qGH\left(I+z^M H^{-1}G^{-1}/q\right)               \tag{14}
```

has a final factor congruent to the identity modulo `z`. Thus the zero
Smith profile of `L_infinity` is exactly that of `H`.

For the Hamming ball

```math
S_{b,k}=\{a:\mathrm{popcount}(a)\le k\},\qquad 2k<b,  \tag{15}
```

every entry of `H` vanishes at zero: the union of two allowed masks cannot
contain every bit. Consequently `L_0` is a unit at zero.

### Reciprocal identity and the second profile

The identity

```math
HG=z^{M-1}Z(1/z)^2
```

gives

```math
\boxed{\quad z^M F(1/z)=L_\infty^T.\quad}                  \tag{16}
```

The transpose is part of the entrywise identity. It has no effect on
Smith orders or determinants.

Each summand in an entry of `Z^2` has exponent
`(a & ~c)+(c & ~b)`. The two masks are disjoint, so this exponent is at
most `M-1`. Hence `deg F<=M`. Also `F(1)=I`, so its determinant is nonzero.
Equation (16) and (13) imply

```math
\deg\det F=nM-\sum_i\nu_i=n+2A.                            \tag{17}
```

Let `phi_i` be the zero Smith orders of `F`. Since `Z(0)=P^T`, all
`phi_i` are positive. In particular,

```math
0\le\sum_i\phi_i=\mathrm{ord}_0\det F\le n+2A.       \tag{18}
```

These aggregate identities suffice for the limiting exterior pressures.
The [singleton profile note](sparse-holonomy-singleton-profile.md)
also proves an explicit individual profile at radius one, and a Gaussian
encoding for every radius.

## 5. From one-particle matrices to Gaussian profiles

For an alternating Gaussian kernel, split the fermions between the two
sides of a flattening. If its internal blocks are regular at the chosen
valuation point, their Gaussian factors act by invertible matrices over
the power-series ring. Removing them leaves the cross block `C`.

The remaining coefficient matrix is the direct sum of all exterior
powers of `C`. Therefore, if `a_1,...,a_m` are the one-particle Smith
orders, the full Gaussian tensor has the subset-sum profile

```math
\left\{\sum_{i\in I}a_i:I\subseteq\{1,\ldots,m\}\right\}.
\quad\text{Its profile sum is }\prod_i(1+t^{a_i}).           \tag{19}
```

This statement also works with Laurent cross blocks. The internal
Gaussian factors must be regular for this exact removal argument.

The plus/minus form of (6) makes its one-particle cross block the direct
sum, up to permutations, of two matrices of size `2n`. For the row cut,
they are

```math
\begin{pmatrix}z^MP&s z^MG\\sP^TQ_J&I\end{pmatrix},\qquad
\begin{pmatrix}z^MP&sQP^T\\s z^MH&I\end{pmatrix}.            \tag{20}
```

Their Schur complements are `z^(2M)P` and `z^M P L_0`.
For the column cut, the corresponding complements are `P` and
`P L_infinity`. Using Sections 2--4 gives the following exact
one-particle zero profiles under (15):

| Cut | Orders measured in `z` |
| --- | --- |
| Row | `0` repeated `2n` times, `M` repeated `n` times, `2M` repeated `n` times |
| Column | `0` repeated `3n` times, followed by the `n` orders of `H` |

The polynomial variable in (6) is `w`, so its ordinary integer orders
are twice those in the table. The full tensor profiles follow from (19).
The fixed-occupation and energy projections can reduce coefficient ranks;
the unprojected Smith profiles provide upper bounds after those projections.

## 6. A fully charged infinity gauge

The energy condition permits more than constant occupation gauges.
At a block with plus index `a` and minus index `c`, give the two fermions
heights `c_+ + alpha a` and `c_- - alpha c`. Their total on a selected
label is `r(c_++c_-)`, because the two index sums agree.

Here is one explicit gauge, with heights measured in powers of `w`:

| Boundary | Plus height | Minus height |
| --- | --- | --- |
| `B` | `0` | `-2M` |
| `X` | `-1-2a` | `-4M+1+2c` |
| `V` | `-1-2a` | `-2M+1+2c` |
| `Y` | `0` | `0` |

Every gauged kernel entry is regular at infinity. The total shift on
every selected external tuple is exactly `-8Mr`. Thus this gauge has
no label-dependent or uncharged coefficient shift.

Put `xi=1/z` and

```math
(P_\xi)_{a,c}=\mathbf1_{c\subseteq a}\xi^{a-c}.
```

Down-set closure gives `P_xi^2=I` and
`Z(xi)=P_xi diag((1+xi)^a)P^T`. Block elimination in the gauged kernel
then gives the following infinity profiles, expressed as orders in `xi`:

| Cut | Orders |
| --- | --- |
| Row | `4n` zeros |
| Column | `n` zeros, `2n` copies of `M`, and `M+phi_i` for each zero order `phi_i` of `F` |

The checker independently computes all four small cross-block Smith
profiles before comparing them with these formulas.

Combining the ungauged zero row profile with the gauged infinity column
profile gives the exact sufficient pressure

```math
\mathcal B=
2^{3n}t^{4r}(1+t^{-1})^{3n}(1+t^{-2})^n
\prod_{i=1}^n(1+t^{-1-\phi_i/M}),\qquad t>1.                \tag{21}
```

The coefficient-degree factors cancel because the gauge shift is charged.
Equation (21) is one certificate for this family, rather than an assertion
about its optimal rank cost.

## 7. The sparse asymptotic regime

For the set in (15),

```math
n=\sum_{j=0}^k\binom bj,\qquad
A=(M-1)\sum_{j=0}^{k-1}\binom{b-1}{j}.                      \tag{22}
```

Every bit occurs equally often. Since every mask has at most `k` bits,
`A/[n(M-1)]<=k/b`. Equations (13) and (18) therefore give

```math
\frac1n\sum_i\left(1-\frac{\nu_i}{M}\right)
\le\frac1M+\frac{2k}{b},\qquad
\frac1n\sum_i\frac{\phi_i}{M}
\le\frac1M+\frac{2k}{b}.                                   \tag{23}
```

For fixed `k`, the normalized `H` orders concentrate at one and the
normalized `F` orders concentrate at zero. Markov's inequality bounds
the exceptional fractions by the right side of (23) divided by any
fixed positive tolerance.

For fixed `t>1`, the function `log(1+t^(-x))` is Lipschitz with constant
at most `log t`. Thus (23) controls normalized log profile sums directly,
with error `O((k/b+1/M) log t)`. In particular, with `p=r/n`, the
per-mode limit of (21) is

```math
\mathcal B^{1/n}\longrightarrow
8t^{4p}(1+t^{-1})^4(1+t^{-2}).                               \tag{24}
```

Detailed radius-two or radius-three Smith recurrences are not needed
for this limit. Other gauges may lead to different sufficient pressures.
No strict post-half certificate is supplied here.

For fixed `k>=2`, `n` grows as a constant times `b^k`, while `log M=b=o(n)`.
With fixed occupation density away from zero and one, equation (2) gives

```math
\log_2d=2nH_2(p)-O(b+\log n).                               \tag{25}
```

The polynomial degree and the number of coefficient pieces are
`O(Mr)`, with logarithm `o(n)`. Each entry uses a determinant of order
at most `4r` with explicit polynomial entries. Coefficient interpolation
or direct polynomial arithmetic takes `2^O(b) poly(n)` operations.
Enumerating all tables and ordinary rank factorizations remains
`2^O(n)` bit operations.

The number of pieces here is subexponential in the label width, rather
than automatically polynomial in it. A strict rank-cost margin could
absorb this overhead and yield a fixed improving seed. The present note
does not skip that margin or claim the final uniform factoring theorem.

## 8. Exact checks

Run:

```sh
python3 tools/check_schur_downset_crossing.py \
  --output evidence/schur-downset-crossing.json
```

The checker verifies the down-set factorizations, the reciprocal identity
with its transpose, determinant-order bounds, and all eight small
one-particle cut profiles. It compares `2,064` selected determinant values
with the independently audited complete Schur evaluator. These checks
include weight `24`, where the separate homogeneous Gram construction
cannot be substituted for the complete dual.

It also verifies `27,957` endpoint pairing entries, including the
`167`-label example `b=4,k=2,r=2`, and two negative controls: an invalid
non-down-set and a proper down-set whose holonomy is not scalar.

Files: [checker](../tools/check_schur_downset_crossing.py) and
[recorded evidence](../evidence/schur-downset-crossing.json).
