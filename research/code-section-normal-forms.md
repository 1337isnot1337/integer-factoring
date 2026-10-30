# Dual code sections: exact normal forms and complete finite screens

## Results and scope

This note gives an explicit family of non-coordinate dual embeddings of the
elementary binary crossing pencil. It reduces arbitrary ambient code lengths
to finite group-algebra data at the logical label width. The family permits
exact coefficient queries, exact factors in singly exponential time, and a
complete finite rank search. It does not improve the factoring exponent.

For coordinate sections, every ambient code length reduces, up to rank-preserving
coordinate permutations, to `2^(2^lambda-1-lambda)` parity masks. Complete
enumeration at logical widths one through five gives the following minimum
of the maximum coefficient cut-rank product:

| Logical width | Classes | Minimum maximum rank product |
| --- | ---: | ---: |
| 1 | 1 | 4 |
| 2 | 2 | 36 |
| 3 | 16 | 224 |
| 4 | 2,048 | 2,116 |
| 5 | 67,108,864 | 15,312 |

The ordinary elementary pencil attains every listed minimum. Through width
four, each individual coefficient rank, in each cut, is at least its ordinary
elementary value. At width five, the native checker proves the sharper
finite statement needed for the table: every class has degree-two rank
product at least 15,312. It does not compute all coefficient ranks there.

For arbitrary linear sections, the finite normal form has
`2^((lambda+1)(2^lambda-1-lambda))` gauge classes. Complete enumeration gives
minima 4, 36, and 224 at widths one, two, and three, respectively. The
width-three screen covers all 65,536 classes, of which 28 attain 224.

These screens are finite statements for this explicit elementary-pencil
family. The subsequent [frame-filtration theorem](frame-filtration-barrier.md)
proves an all-size obstruction for every affine-selector frame, including
arbitrary mixtures and recombinations of their coefficients. General crossing
decompositions and nonlinear selectors remain outside that theorem. No new
Lean claims are made, and the factoring score remains `1/2`.

The starting crossing is the elementary pencil already studied in the
repository, pinned to
[binary-coefficient-entropy.md at PR #23](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/binary-coefficient-entropy.md).
The subsequent embeddings and finite screens here are ordinary exact
mathematics and computation. No claim of literature priority is intended.

## 1. The pencil and a dual code section

All coefficients are in `F_2`. Set

```math
E(x,y,u,v)=[x=v][y=u],\qquad
R(x,y,u,v)=[x=y][u=v].
```

The one-bit pencil is `R+z(E+R)`. Its `t`-fold tensor power is `P_t(z)`;
`P_t(1)=E` on the label space `F_2^t`.

Let `Pi:F_2^t -> B=F_2^lambda` be a linear surjection, and choose a linear
section `J:B -> F_2^t`, so that `Pi J=I`. Put `A=ker(Pi)` and `k=t-lambda`.
Define two embeddings into the vector space on physical labels by

```math
U_b=\sum_{a\in A}e_{Jb+a},\qquad V_b=e_{Jb}.
```

Their pairing is `U^T V=I`. Applying `U` on the first two legs and `V` on
the last two therefore preserves the crossing at `z=1`.

Every nonzero local pencil entry satisfies additive conservation
`x+y=u+v`. Consequently the compressed pencil also has this conservation
law. In a row block let `S=x+y` and let the convolution coordinate be
`w=x+u`. The whole compressed tensor is specified by a group-algebra
kernel `K_S(z)` on `B`.

Write `g_p` for the group basis element at `p in B`, and

```math
N_p=1+g_p,\qquad
F_p(z)=g_p+z^{-1}N_p
      =1+\frac{1+z}{z}N_p.
```

The group algebra is the square-zero algebra

```math
\mathbf F_2[B]\cong
\mathbf F_2[N_1,\ldots,N_\lambda]/(N_i^2).
```

In particular `N_p^2=0` and `F_p^2=1`.

For physical coordinate `j`, let `p_j=Pi(e_j)` and let `a_j` be row `j`
of `J`. The exact normalized compressed kernel is

```math
Q_S(z)=z^\lambda g_S
       \prod_{j:\,a_j\cdot S=0} F_{p_j}(z).            (1)
```

The unnormalized compression equals `z^k Q_S(z)`.

To check (1), the one-bit group kernel is `z g_e` when the row difference
is one and `1+(1+z)g_e=z F_e` when it is zero. Multiply these coordinate
kernels, then push forward the group algebra through `Pi`. Summing over
the coset `A` is precisely this pushforward. The physical prefactor is
`z^t`, and the group element is `g_(Pi JS)=g_S`.

The apparent negative powers in (1) cancel. For a set of distinct nonzero
vectors `P`,

```math
\prod_{p\in P}N_p=
\begin{cases}
\displaystyle\sum_{w\in\langle P\rangle}g_w,
     &P\text{ linearly independent},\\
0,&P\text{ linearly dependent}.
\end{cases}                                            (2)
```

For an independent set this follows by expanding all subset sums. If the
set is dependent, one vector lies in the span of preceding independent
vectors; multiplying their norm by `N_p` gives zero. Thus products of
more than `lambda` distinct `N_p` vanish. Equation (1) is a polynomial
of degree at most `lambda`, and `Q_S(1)=g_S`, the required crossing.

Each norm in (2), for a subspace of dimension `j`, has convolution rank
`2^(lambda-j)`: its image is the space of functions constant on its cosets.
This rank fact alone does not determine the rank of a sum of norms.

## 2. Affine selector normal form

Because `F_p^2=1`, the physical coordinates with the same projection
column `p` matter only through their parity. Set

```math
\theta_p(S)=c_p+v_p\cdot S,
\qquad
c_p=\#\{j:p_j=p\}\pmod2,
\qquad
v_p=\sum_{j:p_j=p}a_j.
```

Then

```math
Q_S(z)=z^\lambda g_S\prod_{p\ne0}F_p(z)^{\theta_p(S)},
\qquad
M:=\sum_{p\ne0}p v_p^T=I.                             (3)
```

Conversely, every affine-selector descriptor satisfying `M=I` comes
from a dual linear section of finite length. For each `p`, add one
physical coordinate with row `e_i` for every nonzero bit `v_(p,i)`.
If needed, add one coordinate with row zero to obtain the requested
constant parity `c_p`. The resulting length is at most
`(lambda+1)(2^lambda-1)`, and `Pi J=M=I`. Hence (3) is an exact
description of this family, independent of ambient length.

There are `lambda(2^lambda-1-lambda)` free slope bits: choose `v_p`
arbitrarily for nonunit `p`, then solve uniquely for the unit-vector
slopes from `M=I`. Constants have `2^lambda-1` bits. A paired coordinate
translation removes `lambda` of these bits. Specifically, replacing
`c_p` by `c_p+v_p dot a` gives

```math
Q'_S=g_a Q_{S+a}.                                     (4)
```

On the four tensor legs this is the permutation `x -> x+a`,
`v -> v+a`, which preserves both cut ranks. The map from `a` to its
constant shift is injective because `M=I`. Choosing a linear complement
to its image gives exactly

```math
2^{(\lambda+1)(2^\lambda-1-\lambda)}                  (5)
```

representatives. The checker constructs this complement by exact binary
elimination, rather than assuming a particular set of pivot positions.

For a coordinate section of a systematic code, the slopes are fixed:
`v_(e_i)=e_i`, and every nonunit slope is zero. The original elementary
pencil has `c_(e_i)=1` and all other constants zero. Every code generator
adds one constant toggle `F_p`. Unit-vector toggles are absorbed by (4),
leaving exactly `2^(2^lambda-1-lambda)` possibilities.

This last reduction depends on the coordinate-section hypothesis. It
must not be applied to arbitrary sections by merely putting a code in
row-echelon form: that operation can change the chosen section.

## 3. Exact rank computation from small convolution blocks

Write `K_i(S,w)=[z^i]Q_S(z)` and `D=2^lambda`. In the row flattening,
the block at `S=x+y=u+v` is convolution by `K_i(S,.)` on `B`. In the
column flattening, the block at `w=x+u=y+v` is convolution by
`K_i(.,w)`. Therefore

```math
r_i=\sum_S\operatorname{rank}\operatorname{Conv}(K_i(S,\cdot)),
\qquad
c_i=\sum_w\operatorname{rank}\operatorname{Conv}(K_i(\cdot,w)).          (6)
```

These are exact identities. The implementation independently checks them
against both complete dense flattenings at small widths.

Multiplying a coordinate-section kernel by `F_p` acts on coefficients by

```math
K'_i=\operatorname{shift}_p(K_i)
     +(I+\operatorname{shift}_p)K_{i+1}.                (7)
```

Use `K_(lambda+1)=0`. Traversing the parity masks in Gray order changes
only one toggle at each step. A direct group-algebra routine separately
constructs (3), and the checker verifies agreement with (7).

The formulas are effective at every width. Constructing all kernels,
both cut ranks, and Gaussian factors takes `2^O(lambda)` bit time.
This supplies an explicit stable crossing family, but no improved
uniform rank profile.

## 4. Uniform augmentation and extreme-coefficient identities

The row-block augmentation of (3) is `z^lambda`, since the augmentation
of every `N_p` is zero. A second identity holds even if `M` is arbitrary:

```math
\sum_S Q_S(z)
 =\det(zI+(1+z)M)\sum_{w\in B}g_w.                  (8)
```

To prove it, expand in `u=(1+z)/z` and use (2). For an independent
`j`-element set `P` spanning `W`, the relevant term is

```math
N_W\sum_Sg_S\prod_{p\in P}\theta_p(S).
```

On any coset of `W`, summing this degree-`j` polynomial extracts its top
coefficient. In a basis given by `P`, that coefficient is the determinant
of the matrix `(v_p dot q)_(p,q in P)`. It is independent of the coset
and of the constants. Cauchy--Binet, summed over `P`, gives
`det(I+uM)`; multiplying by `z^lambda` proves (8).

Let `d=sum_p c_p p`. The top kernel is

```math
K_\lambda(S,w)=[w=(I+M)S+d].                          (9)
```

Thus its ranks are `r_lambda=D^2` and
`c_lambda=4^rank(I+M)`. For actual dual sections, `M=I`, so its ranks
are exactly `(D^2,1)`.

The bottom kernel is constant in its second argument:

```math
K_0(S,w)=b(S),\qquad
b(S)=\det\left(\sum_{p\ne0}\theta_p(S)pp^T\right).     (10)
```

This is Cauchy--Binet applied to the products of `lambda` selected
vectors in (2). If `M=I`, equation (8) gives `sum_S b(S)=1`. In the
local algebra `F_2[B]`, an element is invertible precisely when its
augmentation is one. Hence the bottom ranks are exactly

```math
(r_0,c_0)=(|\operatorname{supp}b|,D^2).                (11)
```

For a contraction `rho<1/sqrt(2)`, the bottom coefficient alone would
therefore require
`|supp(b)| <= D^(2sqrt(2)-2+o(1))`. This is a necessary condition for
the family, not an established obstruction. The determinant support can
be a singleton in the ordinary elementary pencil.

The symmetric affine matrix in (10) has diagonal `d+S`. The translation
gauge can normalize `d=0`. Its off-diagonal entries can be arbitrary
affine forms: pair vectors `e_i+e_j` prescribe these entries, and unit
vectors then repair the diagonal constraint. This offers a structured
way to investigate new sections without enumerating arbitrary codes.

## 5. The width-five convolution classifier

For `lambda=5`, a group kernel is a 32-bit word. The upper-subset
Möbius transform expresses it as a square-zero polynomial `f` in five
variables. The rank of multiplication by `f` is classified exactly:

| First nonzero part | Extra condition | Rank |
| --- | --- | ---: |
| Constant | — | 32 |
| Linear | — | 16 |
| Quadratic | alternating rank four | 12 |
| Quadratic | alternating rank two | 8 or 10 |
| Cubic | dual alternating rank two or four | 4 or 6 |
| Quartic | — | 2 |
| Quintic | — | 1 |
| Zero | — | 0 |

For a nonzero linear part, an algebra automorphism sends `f` to a
variable, proving rank 16. For a nonzero quadratic coefficient, split
two variables and write

```math
f=a+xb+yc+xyd,
```

where `d` is a unit and the other coefficients lie in the remaining
three-variable algebra. Because `f` has no linear part, `b,c` have
zero augmentation. Every unit in this characteristic-two square-zero
algebra satisfies `d^2=1`. The invertible triangular substitution

```math
x'=xd+c,\qquad y'=y+bd
```

gives

```math
f=x'y'+g,\qquad g=a+bcd.                              (12)
```

The residual has zero augmentation. If the remaining algebra has
dimension `m`, then `g^2=0`, and a block calculation gives

```math
\operatorname{rank}M_{xy+g}=m+2\operatorname{rank}M_g. (13)
```

Indeed, on the decomposition with basis `1,x,y,xy`, two blocks are
`M_g`; the other block is `[[M_g,0],[I,M_g]]`, of rank `m` because
`M_g^2=0`. If the original quadratic rank is four, the residual
quadratic on three variables has rank two and multiplication rank two,
giving 12. If the original rank is two, the residual is zero or the
top cubic, giving 8 or 10.

For a cubic leading part on five variables, degree zero contributes
rank one, degree one to degree four has the rank of the complementary
alternating two-form, and degree two surjects onto the top degree.
Higher-degree perturbations add no further rank. Quartic and quintic
leading parts follow by the same degree filtration.

The Python and C++ implementations use this classification, but neither
finite screen relies on it without independent checks. The Python
replay compares direct Gaussian ranks on all 65,536 embedded
four-variable kernels; on every higher-degree perturbation for
quadratic patterns 0, 1, 3, 11, and 129, including quadratic rank four;
and on 10,000 deterministic pseudorandom five-variable kernels.
The native program additionally verifies 100,000 direct ranks, 1,000
transposes, the ordinary coefficient profile, and all 31 involutions
before it processes any chunk.

## 6. Reproduction and finite certificate ledger

Python replay, including every coordinate-section class through width
four and every general-section class through width three:

```sh
python3 tools/check_code_sections.py --full-classifier \
  --output evidence/code-sections-audit.json
```

The width-five native replay is independent of the saved output:

```sh
g++ -O3 -march=native -std=c++20 -fopenmp \
  tools/check_code_kernel_width5.cpp -o /tmp/check-code-width5
/tmp/check-code-width5 0 64 4 /tmp/code-width5
python3 tools/check_code_sections.py --full-classifier \
  --native-prefix /tmp/code-width5 --output evidence/code-sections-audit.json
```

Each of the 64 native chunks covers a disjoint interval of `2^20`
Gray-ordered parity masks. The checkpoint records the interval,
coverage, violations, witness degrees, and a 64-bit additive checksum.
The Python audit checks that the intervals exhaust `2^26` masks and
that every case has a recorded witness. The source regeneration, rather
than the checksum alone, is the mathematical verification procedure.

The exact baseline at width five is

```math
(r_i,c_i)_{i=0}^5=
((1,1024),(12,496),(66,232),(232,66),(496,12),(1024,1)).
```

All 67,108,864 cases pass at degree two, with product at least
`66*232=15,312`; there are zero violations. The baseline attains
15,312, so this is the exact minimum maximum product in the stated
coordinate-section family. The audit JSON records exact coverage,
hashes of every implementation source, and the smaller full profiles.

The finite searches do not justify extrapolation to a uniform rank
barrier. Further progress requires either a new section whose rank
profile improves asymptotically, a different underlying pencil, or an
all-size argument explaining the observed coefficient growth.

The [combined audit](../evidence/code-sections-audit.json) includes the
Python checks and the complete native coverage. The 64 regenerated
chunk records are retained in
[`code-width5-chunks.json`](../evidence/code-width5-chunks.json).
Substantial OpenAI assistant support was used for the proofs,
implementations, and independent reviews.
