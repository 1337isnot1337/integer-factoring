# A moving physical constraint with an eight-mode profile

## Result and scope

There is a uniform bounded-core calculation for a physical code whose
support grows with the reflection interval. For every power of two
`M >= 4`, forbid joint occupation of physical modes `0` and `M-1`,
and retain the full occupation alphabet on every middle mode. The
physical alphabet has `3*2^(M-2)` letters. Its complete two natural
Smith profiles are given by an eight-mode active matrix of order 256,
generic rank 144, and explicit Gaussian spectators.

This construction lies outside the fixed-support hypothesis of
[the surviving-band injectivity theorem](coherent-tail-injectivity.md).
Its constraint acts on the far mode that served as a free partner in
that proof. The construction is also outside the fixed bounded-height
supports in [the protected-height floor](coherent-protected-height25.md).
These facts do not imply a better pressure.

The exact profiles for `M=4,8,...,512` are supplied and replayed below.
At every one of these intervals, the numerical pressure is slightly
worse than the unmodified full Gaussian pressure. There is no new
factoring exponent or competition score. The useful result is an
explicit growing-support construction with complete finite profiles,
including its actual alphabet, height, degree, and regularity costs.

## 1. Physical tensor and two-cut reduction

Use the deflated reflection variable `t=w^2`. A physical occupation
word `A` has height

\[
h(A)=\sum_{a\in A}(2a+1).
\tag{1}
\]

Let

\[
\mathcal A_M=\{A\subseteq\{0,\ldots,M-1\}:\{0,M-1\}\nsubseteq A\}.
\tag{2}
\]

The embedding and retraction are the coordinate inclusion and
coordinate restriction for this alphabet. They satisfy `LP=I`,
are height homogeneous, and have deletion-closed images. Therefore
the internal Gaussian factors descend by the argument in
[the coherent compound reduction](coherent-hot-compound-retractions.md).
For example, if an output occupation lies in (2), every occupation
obtained by deleting some of its bits still lies in (2). The
restricted triangular Gaussian factor has diagonal one and a regular
triangular inverse.

The two regular matrices are consequently the corresponding
restrictions of the compound of

\[
B_M=\begin{pmatrix}D&J\\J&D\end{pmatrix},\qquad
D_{ab}=(1+t)t^{2(a\mathbin{|}b)+1},\qquad
J_{ab}=[a=b]t^{2a+1}.
\tag{3}
\]

Since the embedding is the transpose of the retraction, the zero and
infinity compound matrices agree up to fixed index permutations.
Their complete regular Smith profiles agree. This uses both map
orientations in the compound reduction; it is not an inference from
the zero cut alone.

The raw tensor is the constant restriction of the original polynomial
reflection tensor. Its conservative degree bound remains `8M^2` in
`w`, or `4M^2` in `t`. The natural scalar shifts are `+M^2` and
`-M^2` in `t` at the two cuts, and cancel in their product. No
additional physical normalization is imposed below.

The alphabet polynomial is exactly

\[
A_M(t)=(1+t+t^{2M-1})\prod_{j=1}^{M-2}(1+t^{2j+1}).
\tag{4}
\]

Thus `A_M(1)=3*2^(M-2)`. In particular, the high physical symbol
has height `2M-1`; it is not treated as a fixed low-energy letter.

## 2. Regular concentration

After swapping the two column mode groups, write

\[
B_{M,\mathrm{swapped}}=L\Delta R,
\qquad
L=\begin{pmatrix}I&0\\F&I\end{pmatrix},\quad
R=\begin{pmatrix}I&F^{\mathsf T}\\0&I\end{pmatrix},\quad
\Delta=\operatorname{diag}(J,t^{2M}J),
\tag{5}
\]

where `F_ab=(1+t)t^(2(a & ~b))`. This is the same polynomial
factorization as in [the fixed-tail calculation](coherent-full-tail-core.md).
Its proof uses `FJ=JF^T` and `(F^T)^2=(1+t^(2M))I`.

Put `x=t^2`, `T=t^(2M)`, `p=M-2`, and `e=M-1`. Use the
one-particle concentration `C` specified by

\[
\begin{aligned}
Ce_0&=e_0,\\
Ce_1&=\sum_{b=1}^{M-1}x^{b-1}e_b,\\
Ce_p&=\sum_{b=1}^{M-2}e_b,\\
Ce_e&=e_e.
\end{aligned}
\tag{6}
\]

For `2 <= j <= M-3`, its remaining columns are

\[
Ce_j=e_j+\alpha_j e_1+\beta_j e_p,
\quad
\alpha_j=\frac{1+x^{M-2-j}}{1+x^{M-3}},
\quad
\beta_j=\frac{x^{M-2-j}+x^{M-3}}{1+x^{M-3}}.
\tag{7}
\]

This is the `k=1` specialization of the concentration in
[the full physical tail calculation](coherent-physical-tail-low-band.md).
The denominator in (7) is a local unit at `t=0`. Both `C` and
`J^-1 C J` are regular and have invertible constant matrices.
For the weighted matrix, the factor `x^(M-2-j)` in `beta_j`
cancels exactly the negative weight of its row `p`, while all
other entries have nonnegative valuation. Modulo `t`, the columns
in (6)--(7) are triangular after elementary column additions;
the same is true after weighting. Their inverses are therefore
regular as well.

Apply `U=C^-1` on each mode group on the left and
`V=J^-1 C J` on each group on the right. Then

\[
UJV=J,\qquad U(TJ)V=TJ.
\tag{8}
\]

Their occupation actions are regular unimodular matrices, so (8)
introduces no Smith-order shift. Only the four mode indices
`0,1,p,e` in each of the two Gaussian groups remain active.
The concentration is a local calculation device; it can have a
pole at `t=1`. Generic rank is proved from the original polynomial
tensor, not by evaluating this rational change of basis there.

## 3. Complete transported operators

Write `c_z,a_z` for creation and annihilation on a Gaussian mode,
and abbreviate

\[
c=1+t,\qquad S_1=\sum_{j=0}^{M-2}x^j,\qquad
S_2=\sum_{j=0}^{M-3}x^j.
\tag{9}
\]

The eight active modes are ordered as
`x0,x1,xp,xe,y0,y1,yp,ye`. The row operators for the two selected
physical modes are

\[
\begin{array}{c|l|l}
&\text{creator}&\text{annihilator}\\ \hline
x_0&c_{x0}+c(c_{y0}+xc_{y1})&a_{x0}\\
x_e&c_{xe}+c(c_{y0}+c_{yp}+c_{ye})&(T/x^2)a_{x1}+a_{xe}\\
y_0&c_{y0}&a_{y0}+c(a_{x0}+S_1a_{x1}+a_{xe})\\
y_e&c_{ye}&(T/x^2)a_{y1}+a_{ye}
 +c((T/x)a_{x0}+(T/x^2)a_{x1}+xS_2a_{xp}+a_{xe}).
\end{array}
\tag{10}
\]

After the column swap in (5), the column operators are

\[
\begin{array}{c|l|l}
&\text{creator}&\text{annihilator}\\ \hline
v_0&c_{x0}&a_{x0}+c(a_{y0}+xS_1a_{y1}+(T/x)a_{ye})\\
v_e&c_{xe}&a_{x1}+a_{xe}+c(a_{y0}+a_{y1}+S_2a_{yp}+a_{ye})\\
u_0&c_{y0}+c(c_{x0}+c_{x1})&a_{y0}\\
u_e&c_{ye}+c((T/x)c_{x0}+xc_{xp}+c_{xe})&a_{y1}+a_{ye}.
\end{array}
\tag{11}
\]

Equations (10)--(11) follow by applying (6)--(7) to the columns
and rows of `F` at indices `0,e`. The two relevant columns are
geometric and constant, respectively; their weighted counterparts
use the same two finite sums in (9). One can check the formulas
without computing `C^-1`: multiply each asserted creator image by
`C`, or by `J^-1 C J` on the column side, and compare with the
unconcentrated creator. For annihilators, multiply the original
row by the appropriate concentration matrix. The checker performs
all these identities as exact binary polynomial identities after
clearing the common local unit.

All coefficients in (10)--(11) are polynomials for `M>=4`.
Their canonical pairings are
`a_i c_j+c_j a_i=[i=j]`. No physical mode or hidden tail coordinate
is omitted from these operators.

The projector onto (2), on either physical group, is

\[
\Pi=I+n_0n_e,\qquad n_a=c_a a_a.
\tag{12}
\]

Let the four copies of (12), using (10)--(11), be
`Pi_rx,Pi_ry,Pi_cv,Pi_cu`. The active matrix is

\[
E_M=\Pi_{rx}\Pi_{ry}\,\Gamma(\Delta_{\rm act})\,
     \Pi_{cv}\Pi_{cu},
\tag{13}
\]

with active diagonal orders

\[
(1,3,2M-3,2M-1,2M+1,2M+3,4M-3,4M-1).
\tag{14}
\]

Each omitted mode is a spectator with its original diagonal order:

\[
\{2j+1:2\le j\le M-3\}\ \cup\
\{2M+2j+1:2\le j\le M-3\}.
\tag{15}
\]

## 4. Rank, degree, and complete counting function

The original restricted polynomial tensor has a strict crossing
endpoint, so both cut matrices have generic rank

\[
q^2=9\,2^{2M-4}.
\tag{16}
\]

The spectator diagonal in (15) has generic rank `2^(2M-8)`.
Therefore (13) has generic rank exactly 144. All changes of basis
and projectors preserve total occupation. The original endpoint
rank generating polynomial is

\[
(1+2z)^2(1+z)^{2M-4}.
\]

Dividing out the spectator occupation polynomial gives the active
block ranks

\[
(1+2z)^2(1+z)^4
\quad\longleftrightarrow\quad
(1,8,26,44,41,20,4,0,0).
\tag{17}
\]

The largest ambient occupation block has order 70. Thus the full
calculation uses complete blocks, not a selected coefficient window
or a partial list of pivots. The implementation starts the native
series precision above the complete constructed coefficient degree,
and increases it if a promised pivot is missing.

A direct degree bound for (13) is `28M` in `t`: the diagonal
degree is `16M`, and summing the maximum degrees of the sixteen
operators in the nonidentity terms of (12) adds at most `12M`.
For the eight replayed values of `M`, the actual degree is `16M-4`.
The latter equality is reported as a finite observation, not used
as the uniform bound. These active degrees do not replace the raw
physical degree bound in Section 1.

If `nu_1,...,nu_144` are the nonzero Smith orders of (13), the
complete regular counting function at either cut is

\[
Z_M(t)=\left(\sum_{i=1}^{144}t^{\nu_i}\right)
       \prod_{d\in\mathcal S_M}(1+t^d),
\tag{18}
\]

where `S_M` is (15). The common-scale pressure equation, with all
physical heights charged, is

\[
Z_M(e^{-u})^2=A_M(e^{-u/\rho})^{4\rho}.
\tag{19}
\]

## 5. Replayed finite profiles and numerical pressure

The evidence contains every active order and every occupation block
for all eight intervals. The scalar values below come from a
logarithmic grid followed by local scalar refinements; they are not
certified global optima or exact rational pressure certificates.

| `M` | Moving-pair pressure | Full Gaussian pressure |
| ---: | ---: | ---: |
| 4 | 0.718893493435791 | 0.716653659728187 |
| 8 | 0.710993298405198 | 0.710604404520169 |
| 16 | 0.708381217463066 | 0.708315755938014 |
| 32 | 0.707517951673040 | 0.707507211314948 |
| 64 | 0.707236867559971 | 0.707235145177966 |
| 128 | 0.707147149530903 | 0.707146878995158 |
| 256 | 0.707119087874509 | 0.707119046170249 |
| 512 | 0.707110474250595 | 0.707110467928395 |

The change in physical entropy is explicit. Relative to the full
Gaussian physical alphabet,

\[
\frac{A_M(t)}{\prod_{j=0}^{M-1}(1+t^{2j+1})}
=1-\frac{t^{2M}}{(1+t)(1+t^{2M-1})}.
\tag{20}
\]

In all eight complete active profiles, the first nonzero coefficient
of the difference between the full active Gaussian counting
polynomial and `sum_i t^nu_i` occurs at degree `6M-2`, with
coefficient one. This is a finite profile observation. It explains
why merely moving the forbidden pair far away does not help at the
tested scales: the physical alphabet changes earlier than the
active crossing profile.

The independent checker verifies 21,840 weighted polynomial entries,
12 constant unit ranks, 8,064 full-mode operator entries, and 192
canonical pairings. At `M=4`, it checks all 20,736 literal endpoint
entries, computes 144 Smith orders using the independent Python
binary-series implementation, and compares both original compound
profiles. The evidence replay reconstructs all 1,152 active orders.

Reproduction:

```sh
python tools/coherent_moving_pair_core.py \
  --output evidence/coherent-moving-pair-core.json
python tools/check_coherent_moving_pair_core.py \
  --evidence evidence/coherent-moving-pair-core.json \
  --output evidence/coherent-moving-pair-controls.json
```

This result does not address arbitrary coherent maps on a moving
pair, maps on a growing number of modes, or parameter-dependent
physical maps with their complete degree and gauge costs. Those
remain different constructive directions.
