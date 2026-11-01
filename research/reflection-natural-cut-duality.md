# Natural cut duality and downset compound blocks

## Result

The two natural weighted cut matrices of any paired physical restriction
of the uniform reflection tensor are identical up to a scalar monomial
and relabelling. If both physical codes contain the vacuum, their full
scalar content is a monomial. Their primitive natural Smith profiles
therefore have the form

```math
(a+\sigma_i)_i,\qquad(-a+\sigma_i)_i,
\qquad a\ge0,\quad \min_i\sigma_i=0.                  \tag{1}
```

If the two codes are closed under deleting occupied bits, the shared
profile is the profile of a restricted compound-minor matrix. That
matrix splits into blocks by total occupation number. This is an exact
computation shortcut, not an improved factoring exponent.

## 1. The natural kernels

Use the notation of
[stable-hot-reflection-codes.md](stable-hot-reflection-codes.md):
`M>=4` is a power of two, `m=M/2`, and hot index `a=m-1-i` has height
`h_a=4a+2` and physical energy `2h_a`. Put

```math
C_{ab}=(1+w^2)w^{4(a\mathbin{\&}b)},\qquad
D_{ab}=(1+w^2)w^{4(a\mathbin{|}b)+2},\qquad
J_{ab}=[a=b]w^{4a+2}.                                 \tag{2}
```

The regular zero-gauged kernel, in leg order `x,y,u,v`, has cross
blocks `xy=uv=C`, `xu=yv=D`, and `xv=yu=J`. The regular infinity
kernel interchanges `C` and `D`. Substitution in the source Laurent
exponents, followed by addition of the natural heights, proves these
identities directly.

The statement also holds outside the stabilized energy range. A
within-leg reflection edge has `a+b=m-1`, hence `h_a+h_b=2M`.
In either regular gauge its entry is `1+w^(2M)`, identically on all
four legs. Thus exchanging the entire `y` and `u` legs transforms the
zero kernel into the infinity kernel, including all within-leg blocks.

## 2. Exact matrix identity and scalar content

For physical states write `h(x)=sum_(a occupied in x) h_a`. The original
tensor has the common Gaussian scalar `w^(M^2/2)`. Its two regular
gauges consequently satisfy

```math
w^{h(u)+h(v)}T_M(w^{-1};x,y,u,v)
 =w^{-M^2}w^{h(x)+h(u)}T_M(w;x,u,y,v).                 \tag{3}
```

There is no sign ambiguity in characteristic two. Restrict `x,v` to
any code `C_X` and `y,u` to any code `C_Y`. The exchange preserves the
allowed state sets. The left-side natural infinity flattening and
right-side natural zero flattening are consequently the same matrix
with scalar ratio `w^(-M^2)`.

Suppose both codes contain state zero. The all-zero coefficient is
`w^(M^2/2)`. Every common scalar polynomial therefore divides that
monomial, and the monic gcd is `w^g`, where `0<=g<=M^2/2`.
For the primitive tensor define `a=M^2/2-g`. Its natural zero matrix
is `w^a` times the restricted Pfaffian matrix of the regular zero
kernel. Its infinity matrix is `w^(-a)` times exactly the same matrix.

The shared matrix is regular and has a vacuum entry one. Its minimum
Smith order is zero, proving (1). In particular,

```math
Z_0(t)Z_\infty(t)
 =\left(\sum_i e^{-t\sigma_i}\right)^2.                \tag{4}
```

Neither the monomial content nor the separate constant shifts are
needed to evaluate this natural pressure product.

## 3. Remove the two within-cut Gaussian factors

Partition the regular zero kernel into row variables `(x,y)` and
column variables `(u,v)`. Its cross block is

```math
B(w)=\begin{pmatrix}D(w)&J(w)\\J(w)&D(w)\end{pmatrix}. \tag{5}
```

For row and column occupation masks `R,S`, define

```math
E_B(R,S)=
\begin{cases}\det B[R,S],&|R|=|S|,\\0,&|R|\ne|S|.
\end{cases}                                           \tag{6}
```

Let `A` and `F` be the within-row and within-column kernel blocks.
The Pfaffian expansion gives a polynomial matrix factorization

```math
M_K=U E_B V,
```

where

```math
U(R,R')=[R'\subseteq R]\operatorname{Pf}A[R\setminus R'],
\quad
V(S',S)=[S'\subseteq S]\operatorname{Pf}F[S\setminus S'].
                                                               \tag{7}
```

Both matrices are triangular by inclusion, with diagonal one. Their
determinants are one, so they are invertible over the polynomial ring
and over the local power-series ring. Formula (7) is also obtained by
separating the within-row, cross, and within-column factors of the
Gaussian coefficient generating polynomial.

Assume `C_X,C_Y` are closed under deleting occupied bits. The selected
row masks are the product downset `C_X x C_Y`; column masks are
`C_Y x C_X`. Whenever a selected row `R` contributes through (7), its
intermediate mask `R'` is still selected. The same holds on the column
side. Thus the restriction of the full factorization is exactly

```math
M_K|_C=U_C(E_B|_C)V_C,                                \tag{8}
```

and the two restricted triangular factors still have determinant one.
The shared Smith profile in (1) is therefore precisely the profile
of the restricted matrix (6). This is a statement about one cut's
regular matrix factorization; it does not authorize ignoring a cost
at the opposite cut after changing the tensor itself.

## 4. Occupation blocks and a strict polynomial candidate

Equation (6) vanishes unless the row and column have the same total
occupation number. The compound matrix therefore splits into blocks.
If `a_j,b_j` count states of occupation `j` in the two codes, block `k`
has size

```math
N_k=\sum_j a_j b_{k-j}.                               \tag{9}
```

Row and column sizes agree because their codes are exchanged. At
`w=1`, equation (5) is the block swap: `D=0,J=I`. Every occupation
block is a permutation matrix there, proving its full generic rank.
The union of their Smith orders is exactly the shared profile in (1).

Only equal-occupation minors and smaller square profile problems need
be constructed. The full four-leg tensor and ambient alphabet are
unnecessary for this natural-profile computation, and (4) also removes
the second cut computation.

Viewed as a four-leg polynomial tensor, `E_B` itself has the literal
crossing endpoint and can be tested as a separate candidate. Its full
opposite-cut profile and actual degree must then be recomputed: (8)
preserves only the stated regular cut profile. Arbitrary additional
physical gauges can make the triangular factors nonregular, so the
shortcut requires a new integrality proof for such gauges.

## 5. Exact controls

```sh
python tools/check_reflection_natural_cut_duality.py \
  --output evidence/reflection-natural-cut-duality.json
```

The [checker](../tools/check_reflection_natural_cut_duality.py) passed
87,360 complete kernel identities at intervals 4 through 128, and
4,388 coefficient reciprocal identities. Independent polynomial
Pfaffian expansion and permutation determinants agree on 473
factorization entries, including unequal paired downsets. It also
checks 217 compound endpoints and 136 off-occupation zero entries.

A non-downset countercontrol fails the restricted factorization because
its omitted vacuum contributes through (7). The projection step is
therefore not asserted for arbitrary codes. The complete retained
sixteen-letter and twenty-four-letter profiles independently satisfy
the opposite-shift identity (1).

The [evidence](../evidence/reflection-natural-cut-duality.json) binds
all sources and inputs and records the exact non-downset difference.
