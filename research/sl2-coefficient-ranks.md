# Polynomial-size rank blocks and an exactly analyzed reciprocal pencil

## Result

Let `R,C,A,X` be the four width-two tensors in the
[relabeling block model](hamming-relabeling.md): two noncrossing matchings,
the all-distinct indicator, and the xor-zero indicator. Consider

```math
P_t(z)=\bigl(R+zA+(1+z+z^2)X+z^2C\bigr)^{\otimes t}.
```

At `z=1` it is the crossing on `[4]^t`, with binary width `lambda=2t`.
Let `Q_t^*` be the optimized maximum rank product over arbitrary exact
crossing decompositions in its entire coefficient symmetry span.
Over **any characteristic-two field**, even allowing arbitrary numbers
of pieces and extension-field recombinations,

```math
log2 Q_t^*=2t log2 L+O(log(t+2)),                         (1)
```

where

```math
L=\min_{x>1}\left(1+8/x+1/x^2+3x+3/x^3\right)
 =11.620139083019\ldots .
```

The minimizer is the unique positive root of
`3x^4-8x^2-2x-9=0`, lying between one and two. Thus this pencil's exact
optimal leading contraction is

```math
rho_*=(log2 L)/4=0.884638857886\ldots,
```

which is worse than the half threshold. Its ordinary transfer would give
only the weaker saving exponent `log2(1/rho_*)=0.17683948...`.
It does not improve the repo's overall exponent `1/2`.

The proof also gives a polynomial-time exact rank algorithm for symmetric
Boolean multiplication matrices of exponential dimension. Other
gradings of the relabeling seed and arbitrary domain restrictions are
not covered by (1). Substantial OpenAI assistant support was used. The
new applications await independent review and add no Lean formalization.

## 1. The modular representation facts used

Work over an algebraic closure of the characteristic-two field; extending
a field preserves matrix rank. Let `V` be the natural two-dimensional
`SL_2` module and `T(n)` its indecomposable tilting module of highest
weight `n`. Tensor powers of `V` are tilting; see
[Doty, Lemma 13](https://arxiv.org/pdf/math/0610591v3).

Write `n+1=2^j+u`, `0<=u<2^j`. In characteristic two, the tensor-product
formula in [Sutton--Tubbenhauer--Wedrich--Zhu, Proposition 3.4(d) and
Corollary 3.6](https://arxiv.org/pdf/2105.07724v2), together with its
simple-module formula, gives

```math
T(n)\cong\bigotimes_{i=0}^{j-1}V^{(i)}
       \otimes\bigotimes_{i:u_i=1}V^{(i)}.                (2)
```

Here `(i)` is Frobenius twist. Consequently `dim T(n)<= (n+1)^2`.
Characters of these modules have distinct highest weights and leading
coefficient one. The multiplicities in a known tilting module are
therefore recovered by descending character subtraction.

These are genuine direct sums of indecomposable tiltings. Replacing them
with a Weyl filtration would not preserve the rank computation: already
on `V tensor V`, the first unipotent coefficient has rank two, while
the corresponding ranks on its two Weyl factors sum only to one.
The checker retains this negative control.

## 2. Exact Boolean ranks in polynomial size

Let `B_m=F_2[N_1,...,N_m]/(N_i^2)` with its `2^m` monomial basis.
Multiplication by `sum_k alpha_k e_k(N)` is the coefficient-linear
combination of the unipotent polynomial

```math
\prod_i(I+sN_i)=\sum_k s^k e_k(N).
```

Decompose `V^(tensor m)` into copies of `T(n)`. On (2), the unipotent
acts as `prod_i(I+s^(w_i)N_i)`, where the weights `w_i` are the powers
of two listed there, each once or twice. Its matrix entry for a subset
of variables is the monomial whose exponent is their weight sum.
Thus the desired operator is computed on blocks of dimension at most
`(m+1)^2`, then ranked and weighted by the integer multiplicities.
There are at most `m+1` block types. Character subtraction, matrix
construction, Gaussian elimination, and multiplicity arithmetic all
take polynomially many bit operations in `m`.

This evaluates an exponential matrix without constructing it. The
checker compares every coefficient combination through `m=9` with
direct Boolean matrices, not just individual homogeneous coefficients.

## 3. Identify the reciprocal pencil

Normalize the local row matrix by right multiplication with `E`.
The constant block basis from the relabeling note gives two scalars
`z^2`, a block `z^2I+z(1+z)N`, and three blocks

```math
z^2I+z(1+z)U+(1+z)^2UV.
```

With `s=(1+z)/z`, the latter is
`z^2(I+sV)(I+s(U+V))`; the linear change from `U,V` to `V,U+V` is
invertible in the Boolean algebra. The entire local matrix is

```math
z^2\rho_W(u(s)),\qquad W=2\mathbf1\oplus V\oplus3(V\otimes V).
```

Its character is

```math
w(x)=8+x+x^{-1}+3x^2+3x^{-2}.
```

The normalized row polynomial is therefore
`z^(2t) rho_(W^tensor t)(u((1+z)/z))`. Its column polynomial is the
degree-`2t` reciprocal. The eight crossing symmetries exchange only
`R,C`, fixing `A,X`; hence the coefficient span is already closed under
them. Write `W^tensor t=direct_sum_n T(n)^mu_n` and
`c_n=[x^n]w(x)^t`. Then `mu_n<=c_n`, `n<=2t`.

The weights in (2) have every subset sum from zero to `n`. The row-entry
polynomial space of its block is exactly

```math
z^{2t-n}F[z]_{<=n};
```

the reciprocal column-entry space is `F[z]_{<=n}`.
The subset-sum claim follows because the first copy of every binary
weight gives the interval `0,...,2^j-1`; the extra weights have total
`u<2^j`, so their shifted intervals cover everything through `n`.

## 4. A common lower functional and a matching upper bound

For every `n>=t`, both spaces contain `z^t`. Coefficient pairing with
this monomial defines one actual matrix-entry functional in both cuts.
Its value on `P_t(1)=E` is one. A detected piece in any decomposition
therefore acts nontrivially on **every** copy of every `T(n)` with
`n>=t`, in both arrangements. Direct-sum additivity and the dimension
bound give

```math
r_row,r_col>=\sum_{n>=t}mu_n>=c_t/(2t+1)^2.               (3)
```

For the upper bound use the original `2t+1` coefficient pieces `G_j`.
The block `T(n)` contributes to `G_j` only when `n>=2t-j`. Hence, for
any fixed `x>1`,

```math
r_row(G_j)<= (2t+1)^2\sum_{n>=2t-j}c_n
 <=\frac{(2t+1)^2}{1-1/x}w(x)^t x^{-(2t-j)}.
```

Reciprocity gives the analogous column bound with exponent `-j`.
Multiplying yields, uniformly in `j`,

```math
r_row(G_j)r_col(G_j)
 <=\frac{(2t+1)^4}{(1-1/x)^2}\bigl(w(x)/x\bigr)^{2t}.    (4)
```

The minimum of `w(x)/x` is `L`. At its minimizer `x_0`, the positive
numbers `p_a=w_a x_0^a/w(x_0)` for `a=-2,-1,0,1,2` sum to one and
have mean one. Choose integer counts within a constant of `t p_a`
whose sum is `t` and weighted sum is `t`: round the counts at weights
`-2,-1,2`, then solve for weights `0,1`. All counts are nonnegative for
large `t`. The corresponding multinomial term in `c_t` has logarithm
`t log L+O(log t)` by the elementary factorial estimate
`log n!=n log n-n+O(log(n+2))`. Conversely `c_t<=L^t` by positivity.
Thus `log c_t=t log L+O(log t)`. Equations (3)--(4) prove (1).
Uniqueness of the positive saddle follows by rewriting its equation as
`3x^2=8+2/x+9/x^2`, whose sides strictly increase and decrease.

This uses one common functional and actual modular direct-sum blocks;
it does not assume characteristic-zero semisimplicity.

The standard odd-width padding also cannot evade the leading bound:
fixing its last local letter gives `z P_(t-1)`, whose full coefficient
span is unchanged by the scalar `z`. Arbitrary label subsets remain a
different question.

## 5. Replay

```sh
python3 tools/check_sl2_boolean_rank.py
python3 tools/check_sl2_boolean_rank.py --output evidence/sl2-coefficient-ranks.json
make verify
```

The source manifest pins both primary PDFs. The replay checks character
decompositions and dimensions, exhaustive small Boolean ranks, direct
global pencil matrices, exact finite full-span cutoffs, and rational
intervals for the saddle point and leading contraction. These support
the stated scope; they do not replace independent review of the
all-size argument.
