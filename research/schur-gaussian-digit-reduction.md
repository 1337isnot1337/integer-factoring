# A fixed Gaussian seed inside the complete Schur crossing

## Result and scope

The characteristic-two Schur kernel on a binary interval of `M` fermion
indices has an exact reduction to `M` copies of one four-mode kernel.
The reduction uses explicit changes of basis at the four external
boundaries. It preserves the number of each fermion species at every
boundary.

This gives a concrete polynomial crossing on four labels. Its coherent
tensor powers can be restricted to fixed plus and minus occupations.
For the Schur family with Durfee rank `r` and weight `h=Mr`, the reduction
retains the energy projections and gives an explicit crossing decomposition
with only polynomial overhead in the rank cost.

No stronger factoring exponent follows yet. The bare Smith profiles of
the four-label seed have an all-parameter floor above the improvement
threshold. Physical degree gauges and fixed-occupation coefficient ranks
remain separate questions. A lower bound for the unprojected seed does
not become a lower bound for its energy projections.

The starting point is the complete **inhomogeneous** integral Schur dual
`T(s_x)=exp(L_1) omega(s_x)`, whose vacuum pairing is the identity at every
weight. This is separate from the homogeneous Gram-matrix construction;
no prime-power condition on `h=Mr` is needed here. We use the proved
[complete Schur determinant formula](schur-fermion-queries.md).
The construction and integral normalization come from
[PR #30 at `28105165`](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/integral-schur-crossing.md),
using the source identities in
[`anthropics/formal-math` at `16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The digit reduction, projection accounting, and implementation below were
developed with substantial OpenAI assistant support. The independent
Smith-profile audit also used the cancellation-module research lane.

## 1. The complete one-particle kernel

Work over `F_2(z)`, put `q=z+1`, and let `M=2^m`. All matrix indices lie
in `0,...,M-1`. Define

```math
\begin{aligned}
P_{k,a}&=\mathbf1_{a\subseteq k},&
D_{a,k}&=\mathbf1_{a\mathbin{\&}k=0},&
J_{a,k}&=\mathbf1_{a+k=M-1},\\
Z_{k,a}&=z^{k\mathbin{\&}\neg a},&
R_{a,b}&=\mathbf1_{a\subseteq b}\,q^{-b-1},&
S_{a,k}&=z^{-1-(a\mathbin{|}k)}.
\end{aligned}                                                \tag{1}
```

The subset relation is the bitwise subset relation. In the order
`B,X,V,Y` for bra, transformed first parameter, second parameter, and ket,
the one-particle contraction matrix is

```math
\mathcal K_M(z)=
\begin{pmatrix}
0&Z&P&I\\
Z^T&0&R&S\\
P^T&R^T&0&D\\
I&S&D&0
\end{pmatrix}.                                               \tag{2}
```

Each block row indexes plus fermions and each block column indexes minus
fermions. Selecting the occupied indices at all boundaries and taking
the corresponding determinant gives the Gaussian coefficient tensor.
The ordinary Schur entries select equal numbers of plus and minus
fermions at every boundary.

### Exact binary identities

All matrices in (1) factor over binary digits. In particular,

```math
P=\bigotimes_{j=0}^{m-1}\begin{pmatrix}1&0\\1&1\end{pmatrix},
\quad
D=\bigotimes_{j=0}^{m-1}\begin{pmatrix}1&1\\1&0\end{pmatrix},
\quad
Z=\bigotimes_{j=0}^{m-1}
\begin{pmatrix}1&1\\z^{2^j}&1\end{pmatrix}.                   \tag{3}
```

The two-by-two identities and `q^(2^j)=1+z^(2^j)` give

```math
\begin{gathered}
P^2=I,\qquad D^3=I,\qquad DP=J,\qquad JZ=Z^TJ,\\
Z^2=q^{M-1}I,\qquad S=z^{-M}JZ,\\
Z^TD=q^M R,\qquad SP=(q/z)^M R.                              \tag{4}
\end{gathered}
```

Thus `Z^-1=q^(1-M)Z`. No large determinant or partition enumeration
is required for these identities.

For example, the two crossing cuts of (2) have one-particle block
matrices `[[R,Z^T],[D,I]]` and `[[R,S],[P,I]]`. Their Schur complements
are `z^M R` and `z^(-M)R`. This identifies the exact cut behavior before
any Fock or energy projection.

## 2. Reduction to one polynomial four-mode kernel

Set `z=w^2`, `W=w^M`, and `s=W+1`. Then `s^2=q^M`. Apply independent
row and column changes at each physical boundary:

```math
\begin{array}{c|cc}
\text{boundary}&U&V\\\hline
B&I&z^{-M}J\\
X&s z^M Z^{-1}J&s Z^{-1}\\
V&s z^M D^{-1}&sP\\
Y&z^M J&I.
\end{array}                                                   \tag{5}
```

Equations (4) give the exact block identity

```math
\boxed{\quad
\mathrm{diag}(U_B,U_X,U_V,U_Y)\,
\mathcal K_M(w^2)\,
\mathrm{diag}(V_B,V_X,V_V,V_Y)
=K(W)\otimes I_M,\quad}                                    \tag{6}
```

where the fixed polynomial kernel is

```math
K(w)=
\begin{pmatrix}
0&w+1&w+1&1\\
w+1&0&w^2&w+1\\
w+1&w^2&0&w^2(w+1)\\
1&w+1&w^2(w+1)&0
\end{pmatrix}.                                               \tag{7}
```

For instance, the `X,V` block becomes `z^M I`, while the `V,Y` block
becomes `s z^M I`. These are the corresponding entries `W^2` and
`W^2(W+1)` of (7). The other ten off-diagonal blocks follow from (4).
All diagonal blocks remain zero.

### The four-label crossing

Use the four physical labels `0,1,2,3` for empty, plus, minus, and both
occupations of one mode. Let `flip` interchange plus and minus.
For external labels `(x,y,u,v)`, select occupations

```math
(B,X,V,Y)=(\mathrm{flip}(u),
\mathrm{flip}(x),v,y).                                 \tag{8}
```

The tensor entry `T(w;x,y,u,v)` is the selected determinant of (7),
or zero when the total selected row and column counts differ. Its degree
is at most eight. The full four-by-four determinant is
`1+w^4+w^8`, which attains this degree.

At `w=1`, only the `B,Y` and `X,V` pairing blocks survive. Therefore

```math
T(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.                    \tag{9}
```

Thus the nine coefficient tensors of `T` form an exact four-label crossing.
The construction uses no auxiliary dimension. The explicit seed is stored
in [the polynomial tensor table](../evidence/schur-gaussian-q4-seed.json).

## 3. Occupation sectors and Schur energy slices

Take the `M`-fold tensor power with the same coefficient variable at every
site, and retain exactly `r` plus occupations and `r` minus occupations
at each external boundary. Denote the resulting tensor polynomial by

```math
\Gamma_{M,r}(w)=\sum_{j=0}^{8M}\Gamma_{M,r,j}w^j.             \tag{10}
```

There are

```math
D_{M,r}=\binom Mr^2                                         \tag{11}
```

labels at one boundary. Equation (9) implies that the coefficient family
in (10) is already a crossing on these labels. This family can be studied
directly, without any symmetric-function conversion.

A Schur label of Durfee rank `r`, with all arms and legs less than `M`,
is a pair of `r`-subsets `A,B` of `0,...,M-1`. Its weight is
`sum(A)+sum(B)+r`. Put `h=Mr`, and let `N_t` count `r`-subsets with sum `t`.
Complementing the indices of `B` shows that the weight-`h` slice has size

```math
d_{M,r,h}=\sum_tN_t^2
\ge\frac{\binom Mr^2}{r(M-r)+1}.                            \tag{12}
```

The inequality is Cauchy--Schwarz. The denominator counts the possible
subset sums. Thus the fixed-energy condition costs only a polynomial
factor in label dimension.

| `M` | `r` | Weight `h` | Full occupation dimension | Fixed-energy dimension |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 1 | 2 | 4 | 2 |
| 4 | 2 | 8 | 36 | 8 |
| 8 | 4 | 32 | 4,900 | 390 |
| 16 | 8 | 128 | 165,636,900 | 4,836,452 |

These are exact dimension counts. The last two rows are not computations
of the corresponding full tensor ranks.

## 4. Keep the external maps and the energy projection

For matrices `A,B`, let `L(A,B)` be their action on the two exterior
powers of order `r`. Its entries are

```math
L(A,B)_{(S_+,S_-),(T_+,T_-)}
=\det A_{S_+,T_+}\,\det B_{T_-,S_-}.                        \tag{13}
```

This formula follows directly from Cauchy--Binet. It preserves the two
occupation numbers. Removing the scalar factors from the inverses of
(5) leaves the four boundary maps

```math
L_B=L(I,J),\quad L_X(z)=L(JZ,Z),\quad
L_V=L(D,P),\quad L_Y=L(J,I).                                 \tag{14}
```

The physical convention (8) conjugates the first two maps by fixed label
permutations. These permutations have no effect on either flattening rank.
Only `L_X` depends on `z`. Its polynomial degree is at most

```math
L=2r(M-1).                                                   \tag{15}
```

The product of the scalar factors removed from the four maps is
`z^(-2Mr)q^(-2Mr)`. In the complete Schur crossing polynomial, the
original prefactor is `z^(2h)q^(h+r)`. With `h=Mr`, their product is
`q^(-r(M-1))`.

Let `Pi_h` denote the coordinate projection onto the weight-`h` labels
at all four boundaries. For `M>=2`, the substitution `W=z^(M/2)` has
integer exponents. The exact tensor identity is

```math
P_h(z)=q^{-a}\,
\Pi_h\bigl(L_B\otimes L_X(z)\otimes L_V\otimes L_Y\bigr)
\Gamma_{M,r}(z^{M/2}),
\qquad a=r(M-1).                                            \tag{16}
```

Here `P_h(z)=z^(2h) sum_s A_s z^s` is the complete inhomogeneous Schur
polynomial with the integral dual `T(s_x)`, on the selected labels.
It has degree at most `4h`, and `P_h(1)` is the crossing tensor at every
weight. The local maps precede the energy projection in (16).
The derivation retains that order; it does not assume that a change of
one-particle basis preserves energy.

### An explicit finite crossing decomposition

Write `L_X(z)=sum_(i=0)^L L_i z^i`. Expanding the scalar at zero gives

```math
(1+z)^{-a}=\sum_{n\ge0}\binom{a+n-1}{n}z^n,                 \tag{17}
```

with the constant-series interpretation when `a=0`. Summing the
coefficients of `P_h` through its degree bound `4h` makes the scalar
multiplier of the pair `(i,j)` equal to

```math
b_{i,j}=\begin{cases}
\binom{a+N}{N}\bmod2,&N=4h-i-(M/2)j\ge0,\\
0,&N<0.
\end{cases}                                                 \tag{18}
```

The hockey-stick identity gives the binomial coefficient. Lucas's rule
reduces it to the test `a & N = 0`. Consequently the tensors

```math
B_{i,j}=b_{i,j}\,
\Pi_h\bigl(L_B\otimes L_i\otimes L_V\otimes L_Y\bigr)
\Gamma_{M,r,j}                                               \tag{19}
```

sum exactly to the crossing tensor. There are at most
`(2r(M-1)+1)(8M+1)` pieces. No infinite series is part of the output.

Every map in (19) acts on one physical leg, followed by coordinate
projections. Hence it cannot increase either flattening rank. If
`R_j,C_j` are the two ranks of `Gamma_(M,r,j)`, then (19) gives

```math
\sum_{i,j}\sqrt{\mathrm{rank}_R B_{i,j}
                       \mathrm{rank}_C B_{i,j}}
\le\bigl(2r(M-1)+1\bigr)\sum_j\sqrt{R_jC_j}.                \tag{20}
```

This is the promised polynomial overhead. It is a constructive upper-bound
transfer. Energy projection can still reduce actual ranks further.

All coefficient tables are effective in `2^O(M)` bit operations.
A coefficient of the seed tensor power is a product of `M` fixed-degree
polynomials. The exterior maps use matrices of dimension at most `M`
and polynomials of degree `O(M^2)`; division-free determinant algorithms
or subset dynamic programming give the required bound. The complete
tables and their ordinary binary rank factorizations have dimensions
`2^O(M)`. For fixed density `0<r/M<1`, the external label width is linear
in `M`, so this is the required exponential-in-width construction cost.

## 5. Bare profiles and the remaining balanced search

The seed's exact zero orders and infinity degrees are as follows.
Each list contains sixteen entries, one for each channel of its
sixteen-by-sixteen flattening.

| Cut | Zero orders | Infinity degrees |
| --- | --- | --- |
| `(x,y)` against `(v,u)` | `0` four times, `4` eight times, `8` four times | `8,7,7,6,5,5,4,4,4,4,3,3,2,1,1,0` |
| `(x,u)` against `(v,y)` | sixteen zeros | `8,6,6,4,2,2,0,0,0,0,-2,-2,-4,-6,-6,-8` |

For the useful ungauged orientation and an evaluation parameter `s>=1`,
the two profile sums are

```math
U(s)=(s^3+s^{-3})^2(s+s^{-1})^2\ge16,
\qquad V(s)=4(1+s^{-4})^2\ge4.                              \tag{21}
```

Their product is at least `64`, above `4^(2 sqrt(2))`. The other
orientation has product at least `256`. Thus the bare scalar Smith
criterion cannot beat the half-rate target for this seed. Tensor blocking
multiplies these scalar profiles and does not remove this particular
obstruction. These statements concern that certificate, not every
possible decomposition of its coefficient tensors.

Fixed occupation permits additional degree gauges. Give each of the eight
physical species an integer weight, and multiply a local tensor entry by
`w` to the sum of its occupied-species weights. On the balanced sector
the `M`-fold coefficient degree shifts by `r` times the sum of the eight
weights, independent of the external label.

For clarity, use a parameter `0<t<1`. Let `alpha` and `beta` be independent
eight-component gauges for the zero and infinity cuts, and let their
Smith orders and degrees be `sigma_R(alpha)` and `sigma_C(beta)`.
Writing `p=r/M`, the resulting pressure is

```math
\mathcal B(t,p;\alpha,\beta)=
t^{p(\sum\beta-\sum\alpha)}
\left(\sum_i t^{\sigma_{R,i}(\alpha)}\right)
\left(\sum_j t^{-\sigma_{C,j}(\beta)}\right).                \tag{22}
```

The coefficient-degree factors cancel between the two cuts. Since
`log_2 D_(M,r)=2M H_2(p)+O(log M)`, a strict certificate

```math
\mathcal B<2^{4\sqrt2 H_2(p)}                               \tag{23}
```

would meet the post-half rank-width target for the balanced Gaussian
family. The present note supplies no such certificate. More precise
coherent coefficient ranks can also lie below these Smith bounds.

## 6. Reproduction

The checker uses Python 3 and the standard library:

```sh
python3 tools/check_schur_gaussian_reduction.py \
  --output evidence/schur-gaussian-reduction.json \
  --seed-output evidence/schur-gaussian-q4-seed.json
```

It checks all `256` endpoint entries, `15,345` one-particle identities over
`F_128`, and `544` complete Schur entries after applying all exterior maps
and energy projections. The latter comparison uses the independently
audited Schur evaluator and includes `M=4,r=3,h=12`, a weight that is not
a prime power. It also verifies the dimension counts above.

For the four Smith profiles, `36` exact block-Toeplitz ranks determine
all orders through the largest one plus one. A separate implementation
in the cancellation-module lane independently obtained the same profiles.
The all-parameter lower bounds in (21) are algebraic consequences of those
exact lists, not extrapolations from sampled evaluation points.

See [the checker](../tools/check_schur_gaussian_reduction.py),
[verification output](../evidence/schur-gaussian-reduction.json), and
[complete seed coefficients](../evidence/schur-gaussian-q4-seed.json).
