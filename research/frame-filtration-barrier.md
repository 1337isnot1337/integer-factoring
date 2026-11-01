# A common rank quotient for all affine-selector frames

## Result

Every dual linear-code section of the elementary binary crossing pencil has
a common algebraic obstruction. This includes arbitrary ambient code lengths,
arbitrary linear sections, all linear relabelings of the logical labels, and
arbitrary mixtures of coefficient tensors from different sections.

There is a linear map from this entire coefficient symmetry span onto the
ordinary elementary coefficient span which preserves the crossing and never
increases either cut rank. Thus every decomposition in this larger class,
with any number of pieces, has leading contraction at least

```math
H_2(1/4)=2-\frac34\log_2 3
        =0.811278124459\ldots>1/\sqrt2.                (1)
```

The elementary coefficient family attains this leading rate, so the rate is
sharp for the class. The finite width-five screen in
`code-section-normal-forms.md` is consistent with this theorem; the present
proof applies at every logical width and allows collective recombinations
across frames rather than only literal coefficients.

This is a construction-specific theorem, not a factoring lower bound. In
particular it does not cover nonlinear selectors, arbitrary nonlinear label
relabelings, arbitrary logical-label subsets, or different local pencils.
The factoring exponent remains `1/2`. The proof and replay are ordinary
mathematics and exact computation, without new Lean formalization or a
claim of literature priority.

## 1. Frames and the precise rank quotient

Work over a field `K` of characteristic two. Let `B=F_2^lambda` and
`D=2^lambda`. In the group algebra `K[B]`, write

```math
g_p g_q=g_{p+q},\qquad N_p=1+g_p.
```

Choose affine selectors

```math
\theta_p(S)=c_p+v_p\cdot S,
\qquad c_p\in\mathbf F_2,\quad v_p\in B,
\qquad \sum_{p\ne0}p v_p^T=I.
```

Put `d=sum_p c_p p`. The reciprocal form of the normalized dual-section
pencil is

```math
K_S(z)=g_d\prod_{p\ne0}(1+zN_p)^{\theta_p(S)}.         (2)
```

An exponent is a bit, so a factor can equivalently be written
`1+z theta_p(S)N_p`. Products of more than `lambda` of the `N_p` vanish.
Consequently (2) has degree at most `lambda`. At `z=1`,

```math
K_S(1)=g_d g_{\sum_p\theta_p(S)p}=g_S.                 (3)
```

As in the normal-form note, its four-leg tensor is

```math
T(z)(x,y,u,v)
 =[x+y=u+v]\,[g_{x+u}]K_{x+y}(z).                    (4)
```

Equation (3) gives `T(1)=E`, where `E=[x=v][y=u]`.

Let `W_lambda` be the `K`-linear span of coefficients of (2)--(4), for
**all** choices of selectors satisfying the displayed matrix identity,
together with their eight crossing-preserving leg permutations. No common
frame is required among pieces, and their number is unrestricted.

For the canonical frame, use only `p=e_j`, with `v_(e_j)=e_j` and all
constants zero. Its coefficients `B_i` have kernels

```math
B_i(S,N)=\Omega_i(S,N)
 :=\sum_{|A|=i}S_A N_A,\qquad 0\le i\le\lambda,       (5)
```

where `S_A=prod_(j in A)S_j` and `N_A=prod_(j in A)N_(e_j)`.
They sum to the crossing.

**Rank-quotient theorem.** There is a linear map

```math
\Pi:W_\lambda\longrightarrow\operatorname{span}_K\{B_0,\ldots,B_\lambda\}
```

which is the identity on the target span, satisfies `Pi(E)=E`, and obeys

```math
r(\Pi F)\le r(F),\qquad c(\Pi F)\le c(F).             (6)
```

The proof below constructs `Pi` from the actual tensor entries. It does
not depend on a possibly nonunique representation by frame coefficients.

## 2. Two algebra structures for the row blocks

Every tensor in (4) has row blocks indexed by `S=x+y=u+v`. The block
at `S` is convolution by its group kernel `K_S`. Therefore the row rank
is the rank of multiplication by the kernel in

```math
\mathcal A
 =K[S_1,\ldots,S_\lambda]/(S_j^2-S_j)
   \otimes K[N_1,\ldots,N_\lambda]/(N_j^2).           (7)
```

The first factor is the algebra of functions of `S`; evaluation splits
it into the `D` row blocks. The second is the group algebra in its
square-zero coordinates `N_j=1+g_(e_j)`. These are invertible changes of
basis and preserve rank.

For a basis monomial `S_B N_A`, give it weight

```math
\operatorname{wt}(S_B N_A)=|A|-|B|.                   (8)
```

Suppose a kernel polynomial `F` contains only monomials of nonnegative
weight. Multiplication by it does not decrease weight. Indeed, a term
`S_D N_C` multiplying `S_B N_A` is zero if the `N` variables overlap;
otherwise the weight increase is

```math
|C|-|D|+|D\cap B|\ge0.                               (9)
```

Order a basis by weight. The multiplication matrix is block triangular.
Its rank is at least the sum of the ranks of its diagonal weight blocks.
Equivalently, a filtered map has rank at least that of its associated
graded map. One can also prove this by an invertible diagonal scaling
over `K(t)` and specialization at `t=0`.

For weight zero, (9) forces `|C|=|D|` and `D` disjoint from `B`.
Thus Boolean overlaps disappear in the associated graded map. Its
diagonal blocks are exactly multiplication by the weight-zero part of
`F` in the square-zero algebra on **both** variable families:

```math
\mathcal E
 =K[S_1,\ldots,S_\lambda,N_1,\ldots,N_\lambda]/(S_j^2,N_j^2).          (10)
```

Consequently

```math
r(F)\ge\operatorname{rank}_{\mathcal E}M_{F_{(0)}}.   (11)
```

Here `F_(0)` means the part of weight zero, not its constant term.

## 3. Every affine frame has the same weight-zero coefficients

Let `F_i=[z^i]K_S(z)`. Every chosen factor contributes one `N_p`, so
every monomial of `F_i` has `N`-degree at least `i`. A product of `i`
affine selectors has Boolean `S`-degree at most `i`. The prefactor
`g_d=prod_j(1+N_j)^(d_j)` has constant term one and no `S` variables.
Thus all monomials have nonnegative weight, and a weight-zero monomial
must have both degrees exactly `i`.

To obtain it, take the linear `N` part `p dot N` of every `N_p` and
the linear `S` part `v_p dot S` of every selector. Repeated `S` variables
lower Boolean degree and do not contribute at weight zero. In the
associated graded square-zero algebra, the coefficient of `S_I N_J`,
with `|I|=|J|=i`, is

```math
\sum_{|P|=i}
 \det((v_p)_{p\in P,I})\det((p)_{J,p\in P})
 =\det\left(\sum_p p v_p^T\right)_{J,I}
 =[I=J].                                             (12)
```

This is Cauchy--Binet, valid in characteristic two. Constants in the
selectors and higher `N` terms have positive weight. Therefore

```math
(F_i)_{(0)}=\Omega_i.                                 (13)
```

For any linear combination from any number of frames, the entire
weight-zero part consequently has the unique form

```math
F_{(0)}=\sum_{i=0}^{\lambda}a_i\Omega_i.               (14)
```

Define `Pi(F)=sum_i a_i B_i`. The numbers `a_i` are intrinsic polynomial
coefficients of `F`, hence this is well-defined on actual tensors.

## 4. The column cut reverses the same diagonal coefficients

The column blocks are indexed by `w=x+u=y+v` and convolve in `S`.
Changing from the row kernel to the column kernel sends one monomial
`S_B N_A` to

```math
\left(\prod_{j\notin A}(1+w_j)\right)
 g'_B N'_{B^c}.                                      (15)
```

To verify this, the group coefficients of `N_A` are the indicator of
`w` supported in `A`, namely the first factor. Meanwhile

```math
\sum_S S_B g'_S=g'_B\prod_{j\notin B}(1+g'_j).
```

Every monomial of (15) has `N'`-degree at least `lambda-|B|` and
Boolean `w`-degree at most `lambda-|A|`. Nonnegative weight is preserved.
Equality is possible only when `|A|=|B|`; its extremal term is
`w_(A^c) N'_(B^c)`. In particular

```math
\Omega_i\longmapsto\Omega_{\lambda-i}
\quad\text{on weight-zero parts}.                    (16)
```

The column version of (11) therefore gives

```math
c(F)\ge\operatorname{rank}_{\mathcal E}
 M_{\sum_i a_i\Omega_{\lambda-i}}.                    (17)
```

Crossing-preserving leg permutations act on the conserved kernel by
either leaving it unchanged or swapping `S` and `w`. Formula (15)
shows that their images still have nonnegative weight and a diagonal
part of form (14). This proves that the construction of `Pi` covers
the full symmetry span stated in Section 1.

## 5. The exterior ranks are the exact elementary ranks

On one pair of variables, multiplication by `SN` in the Boolean-times-
square-zero algebra has rank one and square zero. In the basis

```math
1,\quad S+1,\quad N,\quad SN,
```

it sends only the first basis vector to the fourth. This is exactly
multiplication by `SN` in the two-variable square-zero algebra.
Tensoring these constant changes of basis simultaneously identifies
every elementary symmetric polynomial in the local operators.
Hence, for arbitrary `a_i`,

```math
r\left(\sum_i a_i B_i\right)
 =\operatorname{rank}_{\mathcal E}M_{\sum_i a_i\Omega_i}.             (18)
```

For the column cut, (15) applied to (5) gives

```math
g'_{\mathbf1}\sum_{|A|=\lambda-i}
                  \left(\prod_{j\in A}(1+w_j)\right)N'_A.
```

The common group element is a unit. Complementing all Boolean `w`
variables and using the same local change of basis gives

```math
c\left(\sum_i a_i B_i\right)
 =\operatorname{rank}_{\mathcal E}M_{\sum_i a_i\Omega_{\lambda-i}}.    (19)
```

Equations (11), (14), and (17)--(19) prove (6). On a canonical coefficient
`B_i`, the map is the identity. Since `E=sum_i B_i`, it also preserves
the crossing. This completes the rank-quotient theorem.

## 6. Exact elementary ranks and their entropy

The required rank input is pinned to
[`binary-coefficient-entropy.md` at PR #23](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/binary-coefficient-entropy.md).
That note derives the exact elementary ranks from Wilson's inclusion-matrix
formula and cites the primary proof in Peter Frankl,
[*Intersection Theorems and Mod p Rank of Inclusion Matrices*, Theorem 4.1](https://www.renyi.hu/~pfrankl/j62.pdf),
JCTA 54 (1990). The formulas needed here and their consequence are given
explicitly below; the companion checker compares them with direct ranks.

Let `u(k,i)` be the rank of multiplication by the elementary symmetric
polynomial `e_i` in the square-zero algebra on `k` variables. Since each
local operator `SN` has one Jordan block of size two and two zero blocks,

```math
R_{\lambda,i}:=\operatorname{rank}_{\mathcal E}M_{\Omega_i}
 =\sum_{k=i}^{\lambda}\binom\lambda k2^{\lambda-k}u(k,i).             (20)
```

Multiplication by `e_i` sends degree `j` to degree `j+i`; its matrix is
the subset inclusion matrix. With `b=min(j,k-i-j)`, Wilson's formula is

```math
\operatorname{rank}W_{j,j+i}(k)
 =\sum_{a=0}^b
  \left[\binom{b+i-a}{i}\text{ is odd}\right]
  \left(\binom ka-\binom k{a-1}\right).               (21)
```

Use `binom(k,-1)=0` and sum (21) over `0<=j<=k-i` to obtain `u(k,i)`.
The last summand `a=b` always survives. Consequently

```math
\frac1{k+1}\binom kb
 \le\operatorname{rank}W_{j,j+i}(k)\le\binom kb.       (22)
```

Put `B(lambda,t)=sum_(a=0)^t binom(2lambda,a)`. In (20), each possible
value of `b` occurs once or twice as `j` varies. Substituting (22) gives

```math
\frac{B(\lambda,\lambda-i)}{\lambda+1}
 \le R_{\lambda,i}\le2B(\lambda,\lambda-i).            (23)
```

For completeness, the counting identity behind (23), with `a=lambda-k`,
is

```math
\sum_{a+2b\le\lambda-i}
 \binom\lambda a2^a\binom{\lambda-a}b
 =\sum_{v\le\lambda-i}[z^v](1+2z+z^2)^\lambda
 =B(\lambda,\lambda-i).                              (24)
```

The binomial entropy estimate is uniform:

```math
\log_2 B(\lambda,t)
 =2\lambda H_2(t/(2\lambda))+O(\log(\lambda+2)),
 \qquad 0\le t\le\lambda.                            (25)
```

For example, the largest binomial summand is between its entropy power
divided by `2lambda+1` and that power; there are at most `lambda+1`
summands. Thus no asymptotic rank extrapolation is involved.

## 7. One common piece forces the sharp entropy rate

Let `E=sum_j F_j` be any decomposition from `W_lambda`, and apply `Pi`.
The canonical coefficients `B_i` are linearly independent, since their
polynomial parts have different equal bidegrees. In

```math
\sum_j\Pi(F_j)=E=\sum_{i=0}^{\lambda}B_i,
```

some one piece has a nonzero coefficient at
`i_0=floor(lambda/2)`. Write its image as `sum_i a_i B_i`, and let
`q=min{i:a_i!=0}`, `p=max{i:a_i!=0}`. Then `q<=i_0<=p`.

In the exterior algebra, scaling all variables and keeping the lowest
degree proves

```math
r(\Pi F_j)\ge R_{\lambda,q},\qquad
c(\Pi F_j)\ge R_{\lambda,\lambda-p}.                 (26)
```

The same specialization cannot increase rank. Apply (23) and monotonicity
of binomial tails, rather than assuming monotonicity of the exact ranks:

```math
r(F_j)\ge\frac{B(\lambda,\lceil\lambda/2\rceil)}{\lambda+1},
\qquad
c(F_j)\ge\frac{B(\lambda,\lfloor\lambda/2\rfloor)}{\lambda+1}.        (27)
```

Thus both logarithmic widths of this **same piece** are at least
`2lambda H_2(1/4)-O(log(lambda+2))`, proving (1). In exact product form,

```math
r(F_j)c(F_j)
 \ge\frac{B(\lambda,\lfloor\lambda/2\rfloor)
          B(\lambda,\lceil\lambda/2\rceil)}{(\lambda+1)^2}.          (28)
```

The canonical coefficient decomposition has ranks
`(R_(lambda,i),R_(lambda,lambda-i))`. Equations (23)--(25), concavity
and symmetry of binary entropy, and the arithmetic-geometric mean
inequality for the two widths give the matching upper leading rate
`H_2(1/4)`. The balanced central coefficient attains it. This proves
sharpness for the entire affine-frame span.

The argument is unchanged over a characteristic-two extension field:
the imported inclusion matrices have entries in the prime field, and
all filtered rank arguments are valid over the extension.

## 8. Replay and a constructive scope boundary

Run

```sh
python tools/check_frame_filtration.py \
  --output evidence/frame-filtration.json
```

The checker independently assembles canonical kernels, compares every
elementary coefficient combination through width five with exterior
multiplication, checks both polynomial-coordinate transforms, and tests
mixed-frame rank domination. It also verifies Wilson's formula against
direct inclusion matrices and the integer tail bounds at larger widths.

The nonlinear negative control is a valid crossing family outside the
quotient hypothesis. At width three, begin with selectors
`theta_(e_j)=S_j`. Toggle the selectors at `e_1,e_2,e_1+e_2` by
`h(S)=S_1S_3`. Their sum as vectors is zero, so the crossing endpoint
is unchanged. Their factor product is

```math
\prod_{p\in\{a,b,a+b\}}(1+zN_p)
   =1+z(1+z)N_aN_b.
```

At coefficient degree one the resulting polynomial has the nonzero
equal-degree term `S_1S_3 N_1N_2`. It is off the diagonal `A=B` and
therefore not a multiple of any `Omega_i`. The checker verifies this
explicitly while retaining the exact crossing endpoint. This marks a
genuine scope boundary; it does not claim a favorable rank profile for
the nonlinear family.

Substantial OpenAI assistant support was used in this proof and replay.
