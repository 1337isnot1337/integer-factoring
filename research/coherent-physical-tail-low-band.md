# Arbitrary low physical maps and the surviving Gaussian band

## Results and limits

The full-tail reduction extends to arbitrary constant physical split
maps. It retains the internal Gaussian factors and the physical height
normalization. For a map on the first `k` modes, an explicit regular
concentration uses `8k` active Gaussian modes, independently of the
total number `M` of physical modes. Height homogeneity and deletion
invariance are not required for this larger reduction.

For the four occupation states on the lowest two physical modes, an
exact finite certificate proves a further result. Every constant
`F4` split code of positive dimension retains all sixteen channels
in the leading infinite-tail band, at both cuts. This includes every
rank-three and rank-two projector on those four states. The certificate
includes the ambient height corrections for unequal-height mixtures.

Consequently, a proper code of dimension `q_hot < 4` has a positive
first pressure correction on the near-critical tail branch. If
`L=log(4/q_hot)` and `rho_0=1/sqrt(2)`, then, uniformly for bounded `s`,

\[
\begin{aligned}
u&=\frac{\rho_0}{2M}(\log M+s),\\
\rho_M(u)
&=\rho_0+\frac{6\rho_0}{\pi^2M}
 \bigl[L(\log M+s)+e^{-s}\bigr]+o(M^{-1}).
\end{aligned}
\tag{1}
\]

The physical label gauges in this statement are fixed bounded
integers. The statement remains uniform if the maps and gauges are
chosen from any fixed finite set as `M` varies.

This is not an all-scale pressure floor. A finite interval or a
different scale regime is not excluded. The complete high-band Smith
profiles of the larger cores have not been computed. The finite
certificate concerns constant `F4` maps; it is not a certificate for
arbitrary characteristic-two fields or parameter-dependent maps.
No improvement to the factoring exponent is established.

## 1. The full regular Gaussian and its physical maps

Use the notation of [the full-tail compound reduction](coherent-full-tail-core.md):
`M` is the number of physical modes, the reflection interval is `2M`,
and the deflated variable is `t=w^2`. Let

\[
T=t^{2M},\qquad
h(A)=\sum_{a\in A}(2a+1),\qquad 0\le a<M.
\]

The physical height increases with the mode index. After swapping the
two column groups, the cross block is `L Delta R`, as in that note.
The within-row and within-column alternating kernels are both

\[
\mathcal A=
\begin{pmatrix}Q&H\\H&Q\end{pmatrix},\qquad
Q=(1+T)J_{\rm rev},\qquad
H_{ab}=(1+t)t^{2(a\mathbin{\&}b)}.
\tag{2}
\]

Here `J_rev` reverses the `M` indices. If `C_A` denotes multiplication
by the occupation Gaussian with kernel `A`, the complete normalized
coefficient matrix is

\[
C_{\mathcal A}\,\Gamma(L)\,\Gamma(\Delta)\,
\Gamma(R)\,C_{\mathcal A}^{\mathsf T}.
\tag{3}
\]

The factors `C_A` and their inverses are regular. In characteristic
two, a creator is unchanged by conjugation through `C_A`, while
an annihilator gains the creator with the corresponding column of
`A`. The transposed Gaussian gives the dual formula on the right.
These terms are required when a physical space is not invariant under
deletion.

The identities

\[
FJ_{\rm rev}=H,\qquad FH=(1+T)J_{\rm rev},
\tag{4}
\]

where `F_ab=(1+t)t^(2(a & ~b))`, simplify those conjugations. In
particular,

\[
L\mathcal A=\begin{pmatrix}Q&H\\TH&0\end{pmatrix},\qquad
\mathcal A R=\begin{pmatrix}Q&TH\\H&0\end{pmatrix}.
\tag{5}
\]

The elementary digit factorizations used for the cross block also
prove (4). The checker additionally compares every entry of (2)
with the original normalized reflection kernel.

### Physical height corrections at the two cuts

Let `P,L` be arbitrary constant maps on the first `k` modes, with
`LP=I_qhot`. Tensor them with the identity on the full tail. Put

\[
\Pi=PL,\qquad
D=\operatorname{diag}_{A\subseteq\{0,\ldots,k-1\}}t^{h(A)},
\qquad D_H=\operatorname{diag}_j t^{H_j},\qquad R_H=PD_HL.
\tag{6}
\]

The chosen physical label heights are `H_j` on the code, followed by
the usual heights on the tail. The natural ambient heights must be
removed on the hot coordinates before these new label heights are
inserted. In the group order `row x,row y,column v,column u` of the
regular zero kernel, the exact four operators are

| Cut | First row group | Second row group | First column group | Second column group |
| --- | --- | --- | --- | --- |
| Zero | `R_H D^-1` | `R_H D^-1` | `Pi` | `Pi` |
| Infinity | `Pi` | `R_H^T D^-1` | `D^-1 R_H` | `Pi^T` |

At infinity the groups in this table represent the original physical
roles `x,u,v,y`. Thus the ambient hot height corrections occur on
row `u` and column `v`. Assigning both corrections to the row groups
would give a different tensor when the physical maps mix heights.

These formulas follow by writing the compressed tensor before
pulling any height diagonal through a physical map. For example,
the zero row map is `D_H L D^-1`, while the infinity maps on the
second row and first column are `D_H P^T D^-1` and `D^-1 P D_H`.
Embedding the resulting compressed matrices with constant split
pairs gives the table. Constant complements show that the nonzero
Smith factors of these sandwiches are exactly the factors of the
compressed physical matrices.

The underlying physical restriction is polynomial and preserves the
crossing endpoint because `LP=I`. Its original degree is at most
`8M^2` in `w`, or `4M^2` in `t`. These are bounds for the actual
tensor, not for a locally normalized matrix. The ambient scalar
adds `+M^2` to the zero orders and `-M^2` to the infinity orders;
those two shifts cancel. Independent monomial clearings of local
Laurent matrices must instead be subtracted from their respective
computed profiles.

## 2. A concentration retaining both end blocks

Assume `M,k` are powers of two and `M >= 4k`. Write a mode as
`c+kB`, with `0 <= c < k`, `0 <= B < N=M/k`, and set

\[
x=t^{2k},\qquad
S'=\sum_{j=0}^{N-2}x^j=\frac{1+T/x}{1+x},\qquad
S''=\sum_{j=0}^{N-3}x^j=\frac{1+T/x^2}{1+x}.
\tag{7}
\]

For each residue keep `B=0` fixed. On its tail, define `C` by
the following columns:

\[
\begin{aligned}
C_1&=(1,x,\ldots,x^{N-2})^{\mathsf T},\\
C_{N-2}&=(1,\ldots,1,0)^{\mathsf T},\\
C_{N-1}&=e_{N-1},\\
C_j&=e_j+\alpha_j e_1+\beta_j e_{N-2}
       \quad(2\le j\le N-3),\\
\alpha_j&=\frac{1+x^{N-2-j}}{1+x^{N-3}},\qquad
\beta_j=\frac{x^{N-2-j}+x^{N-3}}{1+x^{N-3}}.
\end{aligned}
\tag{8}
\]

The last line is absent when `N=4`. Every denominator is a local
unit. Modulo `x`, `C` is upper triangular with diagonal one.
Moreover, each entry from column `j` to row `i` has valuation at
least `max(0,2k(i-j))`. For the only nontrivial lower entry, this
follows from the factor `x^(N-2-j)` in `beta_j`. Therefore `C`
and `J^-1 C J` are regular with unit determinant, and so are their
inverses.

Let `U=C^-1` and `V=Delta^-1 C Delta`, doubled over the two
Gaussian groups. Then `U Delta V=Delta`. These are regular
unimodular changes on the occupation spaces. No endpoint
regularity of (8) is needed or asserted.

The inactive columns in (8) lie in the common kernel of the tail
covectors

\[
(1,\ldots,1),\qquad (x^{-1},\ldots,x^{-(N-1)}),\qquad e_{N-1}^{\mathsf T}.
\]

The three relevant tail vectors lie in the active column span:
the geometric vector, the constant vector, and the last coordinate.
This simultaneously concentrates the additional creation and
annihilation terms in (5).

The retained physical indices are

\[
\{0,\ldots,2k-1\}\ \cup\ \{M-2k,\ldots,M-1\}.
\tag{9}
\]

Keeping these indices in both Gaussian groups gives `8k` active
modes. All other diagonal factors are unchanged spectators.

## 3. Explicit complete physical operator images

Write the four retained blocks of a group as `0,1,p,e`, standing
for `B=0,1,N-2,N-1`. Let `c_z,a_z` denote creation and
annihilation on an active Gaussian mode, and put

\[
f_{ca}=(1+t)t^{2(c\mathbin{\&}\neg a)},\qquad
h_{ca}=(1+t)t^{2(c\mathbin{\&}a)},\qquad \bar a=k-1-a.
\]

The row images are

\[
\begin{aligned}
C_{x_a}&=c_{x0a}+\sum_c f_{ca}(c_{y0c}+x c_{y1c}),\\
A_{x_a}&=a_{x0a}+(1+T)c_{xe\bar a}
 +T\sum_c h_{ca}(c_{y0c}+c_{ypc}+c_{yec}),\\
C_{y_a}&=c_{y0a},\\
A_{y_a}&=a_{y0a}
 +\sum_c f_{ac}(a_{x0c}+S'a_{x1c}+a_{xec})
 +\sum_c h_{ca}(c_{x0c}+c_{xpc}+c_{xec}).
\end{aligned}
\tag{10}
\]

The column images, after the column swap, are

\[
\begin{aligned}
C_{v_a}&=c_{x0a}+(1+T)(a_{x1\bar a}+a_{xe\bar a})
 +T\sum_c h_{ca}(a_{y0c}+a_{y1c}+S''a_{ypc}+a_{yec}),\\
A_{v_a}&=a_{x0a}
 +\sum_c f_{ca}(a_{y0c}+xS'a_{y1c}+(T/x)a_{yec}),\\
C_{u_a}&=c_{y0a}+\sum_c f_{ac}(c_{x0c}+c_{x1c})
 +\sum_c h_{ca}(a_{x0c}+a_{x1c}+S''a_{xpc}+a_{xec}),\\
A_{u_a}&=a_{y0a}.
\end{aligned}
\tag{11}
\]

All sums run over `0 <= c < k`. At finite `M`, (7) makes every
coefficient in (10)-(11) a complete polynomial. The denominators
used to construct (8) have canceled from these selected images.
The [constructor](../tools/coherent_physical_tail_core.py) implements
these formulas for both complete polynomials and bivariate field
specializations.

Every hot physical matrix is a linear combination of matrix units.
The matrix unit from occupation `B` to occupation `A` is
`c_A Omega a_B`, where `Omega=prod_i a_i c_i`. Substitution of
(10)-(11) therefore gives the complete active matrix for every
operator in the height-correction table. No deletion assumption
is used.

The active ambient order is `2^(8k)`. For equal hot dimensions
`q_hot` on the two paired alphabets, its generic rank is

\[
q_{\rm hot}^2 2^{6k}.
\tag{12}
\]

Indeed, the actual physical crossing has rank
`q_hot^2 2^(2M-2k)`, and there are `2M-8k` spectator modes.
The crossing endpoint proves that full generic rank before any
local change of basis, even if (8) has a pole at that endpoint.

There is also an explicit finite coefficient bound. The sum of
the active Gaussian orders is `16kM`. Each displayed operator
coefficient has degree at most `4M`. A physical matrix unit uses
at most `4k` creators or annihilators. Thus, when `0 <= H_j <= Hmax`,
clearing the two possible hot poles with `t^(2k^2)` gives the
conservative active coefficient-degree bound

\[
80kM+2H_{\max}+2k^2.
\tag{13}
\]

The clearing shift must be removed from the computed local
orders. For fixed `k` the matrix order is constant and these
degree bounds are linear in `M`, so complete exact polynomial
linear algebra is effective with a cost polynomial in `M`.
The constants are large: for `k=2`, the full ambient order is
65,536 and a rank-three code has generic active rank 36,864.
The present computation does not claim to have evaluated those
complete profiles.

## 4. The leading band and its exact `F4` certificate

In (7), now regard `t,T` as independent. All displayed
coefficients are polynomials in `T` with coefficients in
`F4(t)`. The active diagonal has `2k` modes of `T`-order zero;
the other active modes have positive `T`-order. Setting `T=0`
therefore leaves a central diagonal of rank `2^(2k)`.

For `k=2`, the input of that leading band has sixteen columns.
Equations (10)-(11) at `T=0` contain no operators on the two
terminal `y` blocks. Thus its left and transposed right images
use only twelve Gaussian modes, of ambient dimension 4,096.
This reduction is for the leading band only; it does not reduce
the full physical core to twelve modes.

On either side, the four selected physical creator/annihilator
pairs satisfy the canonical anticommutation relations and have
square zero. For the leading input matrix `V`, define

\[
W_A=\Omega\Bigl(\prod_{i\in A}a_i\Bigr)V,
\qquad \Omega=\prod_{i=0}^{3}a_i c_i.
\tag{14}
\]

The canonical relations give
`V=sum_A c_A W_A`. The matrices `W_A` lie in the common physical
vacuum space. Consequently, applying physical functionals
`alpha,beta` to the two two-bit groups is represented by

\[
\sum_{A=0}^{15}\alpha_{A\bmod4}\,
                    \beta_{\lfloor A/4\rfloor}W_A.
\tag{15}
\]

If a physical ambient-height correction is present on one group,
multiply its four coefficients by
`(1,t^-1,t^-3,t^-4)`. This includes every correction in the
table in section 1. New label gauges multiply entire physical
output rows or columns by monomials and do not change the rank
over `F4(t)`.

There are exactly 85 projective nonzero functionals in `F4^4`:
normalize the first nonzero coordinate to one. The exact
specialization uses the field of 256 elements with modulus

\[
z^8+z^4+z^3+z^2+1,
\]

sets `t=z`, and embeds `F4` by `omega=z^85`. The checker verifies
the field inverses and `omega^2=omega+1`. The denominators
`t` and `1+t^4` are nonzero at this specialization.

For each side and all four independent height-correction
patterns, every one of the `85^2` matrices in (15) has rank
sixteen. This gives 57,800 exact cases. A coordinate restriction
to 256 common-vacuum coordinates is used to find a nonsingular
sixteen-by-sixteen minor. A separate scalar C++ finite-field
engine checks every selected minor. The matrices are rebuilt
from the formulas; no unrecorded successful samples are used.

A nonsingular specialized minor proves that the corresponding
rational-function minor is not identically zero. Since the
leading band has only sixteen input columns, this proves exact
rank sixteen over `F4(t)` in every case.

### Consequence for every constant split code

Take any nonzero row of either physical retraction, or any
nonzero embedding column on the transposed right side. It is
one of the 85 functionals up to a nonzero scalar. Equation (15)
therefore proves injectivity of the leading input into that
single pair of physical outputs. Retaining more physical
labels preserves this injectivity. The constant embeddings
used to form the sandwiches are injective as well.

Both the left and right factors of the leading core consequently
have rank sixteen. Its central diagonal is invertible on those
sixteen channels. The complete leading core has rank sixteen
at both cuts, for every constant `F4` split pair of every
dimension `1 <= q_hot <= 4`. The two sides and all height
orientations were checked separately.

For each fixed code and bounded gauge, the remaining `T`-orders
are positive. After substituting `T=t^(2M)`, the Smith orders
are eventually affine functions of `2M`: every minor is a
fixed finite bivariate expression, so its smallest `T` power
dominates once `M` exceeds a fixed bound. Thus exactly sixteen
orders remain bounded as `M` grows. The argument does not
claim a numerical bound for that stabilization interval or
give the remaining affine orders.

## 5. The pressure correction on the near-critical branch

Let

\[
F(u)=\sum_{j\ge0}\log(1+e^{-(2j+1)u}),\qquad
a=\frac{\pi^2}{24},\qquad \rho_0=2^{-1/2}.
\]

The midpoint integral estimate, or Euler--Maclaurin with the
integrable derivative of `log(1+exp(-x))`, gives

\[
F(u)=\frac a u+O(u)\qquad(u\downarrow0).
\tag{16}
\]

The Gaussian active block also has sixteen channels at
`T=0`. The rank result in section 4 implies that, for a fixed
physical code, its active counting function divided by the
Gaussian active counting function has logarithm
`O(u)+O(exp(-2Mu))`. The same is true at the opposite cut.
Including the common spectators therefore gives

\[
\log Z_0(u)+\log Z_\infty(u)
=2F(u)+O(u)+O(e^{-2Mu})+O(e^{-4Mu}/u).
\tag{17}
\]

Here the profiles omit only the canceling ambient scalar
shifts. Bounded negative active orders from height gauges
are included in the `O(u)` term.

The complete physical partition, with `q=q_hot`, is

\[
A_M(v)=\Bigl(\sum_{j=1}^{q}e^{-vH_j}\Bigr)
        \prod_{b=2}^{M-1}(1+e^{-(2b+1)v}).
\]

For bounded physical heights and `Mv` large,

\[
\log A_M(v)
=F(v)-\log(4/q)
 -\frac{e^{-2Mv}}{2v}\bigl(1+O(v)+O(e^{-2Mv})\bigr)+O(v).
\tag{18}
\]

The constant loss `log(4/q)` counts the physical hot alphabet.
There is no corresponding constant loss in (17), since all
sixteen low channels survive.

The balanced common-type pressure equation is

\[
\log Z_0(u)+\log Z_\infty(u)
=4\rho\log A_M(u/\rho).
\tag{19}
\]

For example, whenever `u -> 0`, `M u^2 -> 0`, and
`exp(-2 sqrt(2) M u)=o(u)`, equations (16)-(19) give

\[
\rho_M(u)=\rho_0+
 \frac{12\log(4/q)}{\pi^2}\,u+o(u).
\tag{20}
\]

To retain the leading finite-tail correction, set
`u=rho_0(log M+s)/(2M)` with bounded `s`. Solving (19) first
gives `rho-rho_0=O(log M/M)`. This justifies replacing `rho`
by `rho_0` in the exponential, since the resulting exponent
error is `O((log M)^2/M)`. Equations (17)-(18) then give

\[
\rho_M(u)=\rho_0+
 \frac{12L}{\pi^2}u+
 \frac{6\rho_0}{\pi^2}e^{-2Mu/\rho_0}+o(M^{-1}),
\tag{21}
\]

which is (1). The discarded profile-tail term is
`O(M^(-sqrt(2)))` at this scale; the active high-band term
contributes `O(log(M) M^(-1-rho_0))` to the rate. Both are
`o(M^-1)`. The bounded-order and Euler--Maclaurin terms
contribute `O((log M)^2/M^2)`.

For `q<4`, minimizing the displayed correction over bounded
`s` gives `s=-log L`. On every fixed window containing that
point, the asymptotic minimum is

\[
\rho_0+\frac{6\rho_0 L}{\pi^2M}
          \bigl[\log(M/L)+1\bigr]+o(M^{-1}),
\qquad
u=\frac{\rho_0}{2M}\log(M/L).
\tag{22}
\]

This is an optimization within the specified drifting-scale
window, not a proof of the global minimum over all scales.

For the forbidden-pair code, `q=3` and the coefficient of `u`
in (20) is approximately `0.3497794571`. Formula (22) predicts
`u=0.0051681108` and `rho=0.7091560148` at `M=512`. The exact
profile calculation in the preceding note gives the numerical
values `u=0.0051590394` and `rho=0.7091235806`. The approximation
explains the slow approach from above; it is not a finite
strict-rate certificate.

The low-band certificate shows that choosing a different
constant rank-three or rank-two `F4` map on these same four
states cannot change this leading entropy loss. A search for
a lower surviving-band count must change the physical support,
the coefficient class, or the scale regime.

## 6. Reproduction and independent controls

The [full-kernel checker](../tools/check_coherent_physical_tail_core.py)
compares the transported operators with conjugation through
the complete original kernel and full concentration matrices.
It checks 2,704 weighted polynomial entries, fourteen constant
unit ranks, 20,224 full-mode operator entries, 11,776 bivariate
specialization entries, and 21,632 internal kernel entries.

A separate literal control uses a coherent, nonhomogeneous
one-bit code with `P=(1,omega)^T`, `L=(omega,omega)`. It constructs
the actual normalized four-leg tensor at `M=4`, applies the
physical maps and height corrections, and compares every
transformed entry with the active construction at two field
specializations and both cuts. All 262,144 comparisons pass.
This explicitly checks the unequal-height infinity orientation.
The counts are retained in
[the physical-core controls](../evidence/coherent-physical-tail-core-controls.json).

The [leading-band constructor](../tools/coherent_physical_tail_low_band.py)
retains all 85 projective functionals, the coordinate restrictions,
all eight rank histograms, and source hashes in
[its exact receipt](../evidence/coherent-physical-tail-low-band.json).
The [independent checker](../tools/check_coherent_physical_tail_low_band.py)
rebuilds every selected minor and checks it with
[the scalar field engine](../tools/coherent_physical_tail_low_band_check.cpp).
Its [controls](../evidence/coherent-physical-tail-low-band-controls.json)
record all 57,800 nonsingular minors and the digest of their
complete coefficient stream.

```sh
python3 tools/check_coherent_physical_tail_core.py \
  --output evidence/coherent-physical-tail-core-controls.json
python3 tools/coherent_physical_tail_low_band.py \
  --output evidence/coherent-physical-tail-low-band.json
python3 tools/check_coherent_physical_tail_low_band.py \
  --evidence evidence/coherent-physical-tail-low-band.json \
  --output evidence/coherent-physical-tail-low-band-controls.json
```
