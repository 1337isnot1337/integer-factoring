# Complementary-pair labels give exact subset-intersection blocks

## Result and scope

Overlapping reflected modes admit an explicit all-size reduction that
escapes the two-pairing support condition. For `K` complementary mode
pairs, use labels that occupy exactly `q` whole pairs. The label dimension
is `d=binom(K,q)`. Every complete crossing piece splits into
subset-intersection blocks, with the exact formula below.

This reduces rank computations from matrices of order `d^2` to matrices
of order

```math
\binom{2j}{j}\binom{K-2j}{q-j},
\qquad 0\le j\le\min(q,K-q),
\tag{1}
```

with known multiplicities. At `K=14,q=7`, the full flattening has order
11,778,624; the largest block has order 3,432.

The reduction is constructive and exact. The finite ranks obtained so
far do **not** establish an asymptotic improvement over `rho=1/sqrt(2)`.
This note retains that distinction: it provides a scalable family and
a substantially smaller exact rank problem, not a better factoring rate.

The operators and reflected duality are those of
[fermionic-piece-polynomial.md](fermionic-piece-polynomial.md), which
attributes the source construction to Anthropic's Appendix B at
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.

## 1. The label family

Choose a power of two `M>=2K`, and select `K` disjoint pairs

```math
\{m,2M-m\},\qquad m\text{ odd},\quad 0<m<M.
```

A label occupies both modes of exactly `q` selected pairs. Its weight is
`h=2Mq`, and its particle count is `2q`. Reflection permutes the two
members of every occupied pair, so the reflected unit dual equals the
label monomial itself. Thus `a_x=b_x=t^x` gives unit vacuum duality
on these labels. There is no requirement that `q` be a power of two.

The generalized locality parameter is `2h=4Mq`. The pieces satisfy the
exact crossing sum by the established operator identity. Taking the
least suitable `M` gives `h=O(Kq)`.

## 2. All surviving matches are equal or complementary

In the four-block Wick Pfaffian, use the original mode `x` to name an
`A` vertex whose operator mode is `2M-x`. Nonzero edges imply

```math
\begin{array}{c|c}
\text{edge}&\text{mode constraint}\\\hline
U_u B_b&u\ge b\\
A_x Y_y&x\ge y\\
A_x B_b&x\ge b\\
U_u Y_y&u=y\\
U_u A_x&u+x\ge2M\\
B_b Y_y&b+y\le2M.
\end{array}
\tag{2}
```

The last inequality follows from `(b-1) AND y=0` and the fact that
both modes are below `2M`. All four external states have `2q`
particles. Any perfect matching therefore contains the same number
of `U--A` and `B--Y` edges.

Subtract the total original weight of `B union Y` from that of
`A union U`. The difference is zero. For cross-class edges use the
nonnegative differences in the first four rows of (2). For internal
edges use `u+x-2M` and `2M-b-y`; the added constants cancel because
the two internal-edge counts agree. The sum of these nonnegative
quantities is zero. Hence every cross-class edge matches identical
modes, every `U--A` edge matches complementary modes, and every
`B--Y` edge matches complementary modes. Each surviving scalar
coefficient is one.

The matching graph consequently separates into independent
complementary mode pairs. This is stronger than a numerical rank
observation. It holds at every `K,q`.

## 3. A local Laurent tensor

Put `T=z^(2M)`. In the fraction field `F_2(z)`, divide each vertex by
the following factor:

```math
g_U(k)=z^k,\quad
g_B(k)=((1+z)/z)^k,\quad
g_A(x)=(1+z)^{2M-x},\quad
g_Y(k)=1.
```

The remaining edge values are

```math
UA=UB=UY=AY=1,\qquad BA=T/(1+T),\qquad BY=T,
\tag{3}
```

subject to their equal/complementary mode constraints. The product of
the vertex factors is `(1+T)^(2q)`, and the locality offset is
`T^(-2q)`. Assign the corresponding factor
`((1+T)/T)^(2x)` to a mode pair whose label bit is `x`.
The rational denominators cancel in each resulting local polynomial.

Write `Z=T^2=z^(4M)`. Each external label is now a weight-`q`
binary string of length `K`. The local tensor is zero unless
`x xor y=v xor u`. Its row flattening has these two two-by-two
blocks, indexed by `a=x xor y`:

```math
L_0(Z)=
\begin{pmatrix}1&1\\1+Z^{-1}&1+Z^{-1}+Z\end{pmatrix},
\qquad
L_1(Z)=
\begin{pmatrix}1&Z\\1+Z^{-1}&1\end{pmatrix}.
\tag{4}
```

Rows within a block are indexed by `x`; columns are indexed by `v`.
Equation (4) follows by listing the matchings on the two physical
modes. The finite checker independently compares it with the original
Wick polynomial and, at small sizes, both literal source operator orders.

## 4. Constant weight turns the local products into Johnson blocks

Fix an XOR support `a` of size `2j`. Since both row labels have weight
`q`, the first label has exactly `j` ones inside `a` and `k=q-j`
ones outside. The same is true of the column label. This proves the
block dimension (1); there are `binom(K,2j)` equal blocks.

Let `r_a` be half the distance between the two inside `j`-subsets,
and let `r_o` be half the distance between the outside `k`-subsets.
Multiplying the entries of (4) gives the complete block polynomial

```math
\mathcal B_{K,q,j}(Z)
=Z^{-k}(1+Z)^{r_a+r_o}(1+Z+Z^2)^{k-r_o}.
\tag{5}
```

Thus the entries depend only on the two intersection sizes. This is a
tensor product of two Johnson-scheme matrix algebras, with coefficient
extraction still performed on the **complete** polynomial.
Its support is contained in `-q,...,q`.

The column flattening has local blocks

```math
\begin{pmatrix}1&Z\\1+Z^{-1}&1+Z^{-1}+Z\end{pmatrix},
\qquad
\begin{pmatrix}1&1\\1+Z^{-1}&1\end{pmatrix}.
```

On the same constant-weight sets their product equals (5) with
`Z` replaced by `Z^(-1)`. Consequently the exact complete-piece ranks
are

```math
\begin{aligned}
r_s&=\sum_{j=0}^{\min(q,K-q)}
 \binom K{2j}\operatorname{rank}_{\mathbb F_2}
 [Z^s]\mathcal B_{K,q,j}(Z),\\
c_s&=r_{-s}.
\end{aligned}
\tag{6}
```

The physical source shift is `4Ms`. There is no uncomputed contribution
from any other shift.

The crossing sum can also be verified directly from (5). At `Z=1`,
the entry is one precisely when `r_a=r_o=0`; otherwise it is zero.
Every XOR block is therefore the identity. The sum of the coefficient
tensors is exactly `[x=v][y=u]`.

## 5. Construction costs and finite findings

For a fixed rational density `0<theta<1`, set `q=floor(theta K)` and
select the least `K` with `binom(K,q)>=2^lambda`. Standard binomial
estimates give `K=O(lambda)` and `h=O(lambda^2)`. Restricting the label
set to `2^lambda` labels preserves the exact crossing and cannot
increase either flattening rank.

There are at most `2q+1=O(lambda)` pieces. Formula (5) computes any
entry using polynomially many bit operations in `K`. Complete matrix
construction and Gaussian rank factorization over `F_2` take
`2^O(lambda)` bit time. Thus a future improved joint rank bound for
these matrices would supply effective crossing factors, including
uniform label enumeration, piece indexing, and factor queries.

The exact balanced-slice computations give:

| `K` | `q` | labels | maximum rank product | finite ratio `log2(product)/(4log2(labels))` |
|---:|---:|---:|---:|---:|
| 4 | 2 | 6 | 256 | 0.773706 |
| 6 | 3 | 20 | 9,996 | 0.768588 |
| 8 | 4 | 70 | 1,192,464 | 0.823323 |
| 10 | 5 | 252 | 75,189,960 | 0.819955 |
| 12 | 6 | 924 | 7,155,806,464 | 0.830727 |
| 14 | 7 | 3,432 | 578,697,595,476 | 0.831728 |

These are finite ratios, not asymptotic contraction constants.
In particular they provide no evidence of a rate below `1/sqrt(2)`
at balanced density. Other densities and constructions require their
own analysis.

For symbolic rank analysis, the simultaneous integral forms in
[Ducey et al., *Integer diagonal forms for subset intersection relations*](https://arxiv.org/html/2310.09227v2)
are relevant. Their Theorem 3.4 assumes `n>=3k-1`; it does not directly
apply to the inside middle-slice factor `n=2j,k=j`. The earlier Bier
triangular form applies, but off-diagonal nilpotent blocks cannot be
discarded in characteristic two. This note uses direct exact ranks
and does not claim an unsupported diagonalization.

## 6. Reproduction

Run the quick suite:

```sh
python3 tools/check_paired_mode_johnson.py \
  --output evidence/paired-mode-johnson.json
```

To reproduce the full table, include `--max-pairs 14`. The largest
case takes substantially longer than the quick suite. The implementation
is in [paired_mode_johnson.py](../tools/paired_mode_johnson.py).

The code and ordinary proofs used substantial OpenAI assistant support.
No new Lean formalization or improved factoring exponent is claimed.
