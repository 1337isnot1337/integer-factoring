# Offset intervals give a two-scale Gaussian Schur kernel

## Result and scope

Using a short interval of arm indices and an offset interval of leg
indices gives a new two-parameter polynomial Gaussian crossing. The
one-particle kernel reduces exactly to copies of an explicit four-mode
matrix. The reduction retains the Schur energy projection and has
polynomial overhead when the two interval scales have polynomial ratio.

This is a constructive extension of the
[binary-interval reduction](schur-gaussian-digit-reduction.md) and a
special case of the [Boolean down-set construction](schur-downset-crossing.md).
It gives exact identities, effective coefficient tables, and a new family
for rank searches. It does not establish an improved factoring exponent.
No numerical pressure claim is needed for the reduction.

The construction uses the complete **inhomogeneous** Schur dual and the
[complete determinant formula](schur-fermion-queries.md), with the source
and normalization references recorded there. The offset reduction and
verification were developed with substantial OpenAI assistant support.

## 1. Separate the number of modes from their offset

Let `M>=N>=2` be powers of two. Put `q=1+z` and use the mode set
`S={0,...,N-1}`. A label consists of two `r`-subsets `A,B` of `S`
with `sum(A)=sum(B)`. Its Schur arms are `A`; its legs are
`{M-1-b:b in B}`. Every label has weight `h=Mr`.

Let `N_t` count the `r`-subsets of `S` with sum `t`. The full occupation
space and selected Schur space have dimensions

```math
D=\binom Nr^2,\qquad
d=\sum_tN_t^2\ge\frac{D}{r(N-r)+1}.                         \tag{1}
```

The denominator counts the possible subset sums. In particular, the
entropy loss depends on `N`, not on the offset `M`.

For indices `a,b` in `S`, write

```math
P_{a,b}=\mathbf1_{b\subseteq a},\qquad Q=\operatorname{diag}(q^a),
\qquad G_{a,b}=z^{a\mathbin{\&}b},\qquad
H_{M,a,b}=z^{M-1-(a\mathbin{|}b)}.                          \tag{2}
```

The ordinary full-interval identities give

```math
P^2=I,\quad G=PQP^T,\quad
H_M=z^{M-N}H_N,\quad GH_M=H_MG=z^{M-N}q^{N-1}I.           \tag{3}
```

Thus the sparse holonomy is scalar for this particular down-set, but
its scalar differs from the `M=N` identity. In the notation of the
down-set note, the three auxiliary matrices simplify to

```math
L_0=(1+z^{M-N}+z^M)I,\qquad
L_\infty=z^{M-N}I,\qquad F=z^NI.                           \tag{4}
```

Equation (4) also holds when `M=N`; then `L_0=z^M I`.

## 2. Exact reduction to a polynomial four-mode kernel

Relabel every complemented index `M-1-a` by `a`. In boundary order
`B,X,V,Y` for bra, first parameter, second parameter, and ket, the raw
plus-to-minus contraction matrix is

```math
\mathcal R_{M,N}(z)=
\begin{pmatrix}
0&H_M&P^T&I\\
G&0&q^{-M}PQ&z^{-M}G\\
P^T&q^{-1}Q^{-1}P&0&P^T\\
I&z^{-M}H_M&P^T&0
\end{pmatrix}.                                             \tag{5}
```

These are direct restrictions of the complete Schur contractions. The
complement relabelling makes this matrix nonsymmetric; the full fermionic
alternating matrix still has the blocks `[[0,R],[R^T,0]]`.

Apply these independent row and column changes at each boundary:

```math
\begin{array}{c|cc}
\text{boundary}&U&V\\\hline
B&I&z^{-M}I\\
X&q^{M/2}z^M G^{-1}&q^{N/2}H_M^{-1}\\
V&q^{N/2}z^M P^T&q^{M/2}P^T\\
Y&z^M I&I.
\end{array}                                                \tag{6}
```

All exponents are integral. Equations (3) prove

```math
\operatorname{diag}(U_B,U_X,U_V,U_Y)\,
\mathcal R_{M,N}(z)\,
\operatorname{diag}(V_B,V_X,V_V,V_Y)
=K_{M,N}(z)\otimes I_N,                                  \tag{7}
```

where

```math
\boxed{\quad
K_{M,N}(z)=
\begin{pmatrix}
0&q^{N/2}&q^{M/2}&1\\
q^{M/2}&0&z^M&q^{M/2}\\
q^{N/2}&z^N&0&z^M q^{N/2}\\
1&q^{N/2}&z^M q^{M/2}&0
\end{pmatrix}.\quad}                                     \tag{8}
```

For example, the `X,V` block becomes `z^M G^-1 P Q P^T=z^M I`.
The reversed block becomes
`z^M q^(N-1) G^-1 H_M^-1=z^N I`. Their different powers retain
the offset information. Every other block follows directly from (3).

### A crossing on four physical labels

Encode one site by `0=empty`, `1=plus`, `2=minus`, `3=both`.
For physical labels `(x,y,u,v)`, the boundary states are
`(flip(u),flip(x),v,y)`, where `flip` exchanges the two species.
Define the local tensor entry by the corresponding selected determinant
of (8), or zero when the total row and column counts differ.

At `z=1`, only the two pairings `B,Y` and `X,V` remain, both with
coefficient one. Therefore

```math
T_{M,N}(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u}.             \tag{9}
```

The tensor has degree exactly `3M+N`. Its full selected determinant is

```math
\det K_{M,N}
=1+z^{M+N}+z^{2M+N}+z^{3M}+z^{3M+N}.                       \tag{10}
```

When `M=N`, the two middle terms of equal degree cancel, and (8) is the
previous fixed kernel evaluated at `z^(N/2)`. When `M>N`, the two
crossing powers in (8) differ. Thus this is a two-scale extension of the
earlier seed, with no auxiliary dimension.

## 3. Occupation restriction and exact Schur projection

Take `N` copies of the four-label tensor with one common variable `z`.
Restrict to exactly `r` plus and `r` minus occupations at every physical
boundary. Write the resulting polynomial as

```math
\Gamma_{M,N,r}(z)=\sum_{j=0}^{N(3M+N)}\Gamma_jz^j.         \tag{11}
```

By (9), its coefficients already form a crossing on the `D` labels in
(1). This is a direct Gaussian construction, independent of its Schur
interpretation.

For matrices `A,B`, use the exterior action

```math
L(A,B)_{(S_+,S_-),(T_+,T_-)}
=\det A_{S_+,T_+}\det B_{T_-,S_-}.                         \tag{12}
```

Taking inverses of the changes in (6), and pulling
`H_M=z^(M-N)H_N` out of the `X` minus action, leaves

```math
L_X(z)=L(G,H_N),\qquad L_V=L(P^T,P^T),                     \tag{13}
```

with identity maps at `B,Y`. The polynomial degree of `L_X` is at most
`L=2r(N-1)`; `L_V` is constant. Fixed species-flip permutations convert
these block coordinates to the physical labels. They do not change
either flattening rank.

Before extracting the scalar from `H_M`, the inverse maps contribute
`z^(-2Mr)q^(-(M+N)r)`. The complete Schur polynomial has the prefactor
`z^(2Mr)q^(Mr+r)`. Including the extracted power of `z` proves

```math
\boxed{\quad
P_h(z)=z^\delta q^{-a}
\Pi_h\bigl(I\otimes L_X(z)\otimes L_V\otimes I\bigr)
\Gamma_{M,N,r}(z),\qquad
\delta=(M-N)r,\quad a=r(N-1),\quad h=Mr.
\quad}                                                    \tag{14}
```

Here `Pi_h` selects the equal-sum labels of Section 1 at every boundary.
The exterior maps precede the projection. No commutation with energy
projection is asserted. The left side is the complete inhomogeneous
Schur crossing, whose degree is at most `4h` and whose value at one is
the crossing tensor at every weight, including non-prime-power weights.

### Effective finite crossing pieces

Write `L_X(z)=sum_(i=0)^L L_i z^i`. The scalar `q^-a` has the ordinary
power-series expansion at zero. Sum the coefficients of (14) through
degree `4h`. For each `i,j`, set

```math
E_{i,j}=4h-\delta-i-j,\qquad
b_{i,j}=
\begin{cases}
\binom{a+E_{i,j}}{E_{i,j}}\bmod2,&E_{i,j}\ge0,\\
0,&E_{i,j}<0.
\end{cases}                                               \tag{15}
```

The hockey-stick identity gives this binomial coefficient. Lucas's rule
reduces it to the bit test `a & E_(i,j)=0`. Thus the explicitly finite
tensors

```math
B_{i,j}=b_{i,j}\,
\Pi_h\bigl(I\otimes L_i\otimes L_V\otimes I\bigr)\Gamma_j  \tag{16}
```

sum exactly to the selected Schur crossing. If `R_j,C_j` are the two
flattening ranks of `Gamma_j`, locality of every map gives

```math
\sum_{i,j}\sqrt{\operatorname{rank}_R B_{i,j}
                       \operatorname{rank}_C B_{i,j}}
\le \bigl(2r(N-1)+1\bigr)\sum_j\sqrt{R_jC_j}.              \tag{17}
```

The number of pieces is at most
`(2r(N-1)+1)(N(3M+N)+1)`. This is an upper-bound transfer. Actual
energy-projected ranks may be smaller.

## 4. A sparse limit with polynomial overhead

For example, take `M=N^2` and fixed density `r/N` strictly between zero
and one. Then every arm index divided by `M` is at most `1/N`, while
the leg indices lie in the last `N` positions before `M`. The ambient
energy scale is cubic in `N`, but the label width remains linear:

```math
\log_2 d=2N H_2(r/N)+O(\log N).                            \tag{18}
```

All polynomial degrees, coefficient counts, and the overhead in (17)
are polynomial in `N`. The complete coefficient tables and ordinary
binary rank factorizations are effective in `2^O(N)` bit operations:
the site tensor has constant dimension and polynomial degree, and the
boundary determinants use `N`-dimensional matrices of polynomial degree.

Consequently a strict asymptotic rank-cost margin for this Gaussian
family would survive the displayed maps and energy restriction. This
note supplies the construction and reduction, not such a margin.
Index-dependent gauges whose slopes grow with `M/N` are not controlled
by a bounded-slope profile approximation; no exclusion of those gauges
is claimed here.

## 5. Reproduction

Run the standard-library checker:

```sh
python3 tools/check_schur_offset_interval.py \
  --output evidence/schur-offset-interval.json \
  --seed-output evidence/schur-offset-q4-seed.json
```

It verifies the exact polynomial holonomies, all blocks of the
one-particle change of basis at three points of `F_128`, and every local
endpoint entry in eight parameter pairs. It also constructs full
occupation tensors, applies the exterior maps, and then compares the
projected entries with the independently implemented complete Schur
evaluator. These tests include `M=8,N=4,r=3,h=24`, which has
non-prime-power energy.

See the [checker](../tools/check_schur_offset_interval.py),
[verification output](../evidence/schur-offset-interval.json), and
[four-label seed for M=4,N=2](../evidence/schur-offset-q4-seed.json).
