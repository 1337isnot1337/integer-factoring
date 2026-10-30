# Degree and quadratic-rank budgets for the full exterior model

## Result and scope

Use the exact characteristic-two crossing model from
[the bi-exterior construction](bi-exterior-crossing-algebra.md). Write

```math
\mathcal A_\lambda=F[s_1,n_1,\ldots,s_\lambda,n_\lambda]/(s_i^2,n_i^2),
\qquad d=2^\lambda,\qquad
g=\prod_{i=1}^\lambda(1+s_in_i).
```

For any finite decomposition `g=sum_i f_i`, its actual crossing pieces
have ranks `r_i=rank M_(f_i)`, `c_i=rank M_(Hf_i)`. Then

```math
\boxed{\quad\sum_i\sqrt{r_i c_i}\ge d(\lambda+1).\quad}       \tag{1}
```

This covers arbitrary coefficients, arbitrary mixtures of degrees, and
the entire `4^lambda`-dimensional algebra, without position symmetry.
The bound excludes the strict rank-sum seed criterion through `128`
labels inside this model. A second exact certificate gives

```math
\boxed{\quad S\ge\frac{51}{5}\,256=2611.2
                   >256^{\sqrt2}=2545.456152\ldots.\quad}    \tag{2}
```

Thus a binary seed in this model passing the criterion needs at least
**512 labels**. No such seed is constructed. The proof also covers
supplied formal limits whose Laurent-series pieces stay in this algebra.
It is not a restriction on arbitrary crossing tensors, nonlinear label
frames outside this algebra, or factoring algorithms. The proved
factoring exponent remains `1/2`.

## 1. A supported monomial gives a triangular minor

Let `f` be a nonzero polynomial in the square-zero algebra on `n=2lambda`
generators. Suppose the coefficient of `x_A` is nonzero. In the matrix
of multiplication by `f`, retain the columns `x_I` and rows `x_(A union J)`
for `I,J subset A^c`.

An entry can be nonzero only when `I subset J`: the input monomial must
be contained in the output monomial, and `I` is disjoint from `A`.
For `I=J`, the entry is precisely the nonzero coefficient of `x_A`.
Ordering the subsets by cardinality gives a triangular nonzero minor,
so

```math
\operatorname{rank}M_f\ge2^{n-|A|}.                       \tag{3}
```

Let `p` and `q` be the lowest and highest supported degrees of `f`.
The map `H` complements all variables and then permutes them in pairs.
Equation (3), first at degree `p` and then for `Hf` at degree `n-q`, gives

```math
r(f)\ge2^{n-p},\qquad c(f)\ge2^q,
\qquad\sqrt{r(f)c(f)}\ge d\,2^{(q-p)/2}.                  \tag{4}
```

If `f` has nonzero components in `k` distinct even degrees, then
`q-p>=2(k-1)`. Since `2^(k-1)>=k` for `k>=1`, its cost is at least
`d*k`. A piece with no even component contributes zero to this counting
bound. Every one of the `lambda+1` even components of `g` is nonzero,
and therefore occurs in some piece. Summing proves (1).

At `lambda=7`, equation (1) is `S>=1024`, whereas
`128^sqrt(2)=955.091363\ldots`. The numerical display is illustrative;
the replay verifies this and all smaller binary dimensions by integer
Pell-exponent comparisons.

## 2. Exact ranks of a quadratic form and its dual

Let a nonzero homogeneous quadratic form have alternating coefficient
matrix of rank `2a`. An invertible change of the generators brings it to

```math
q_a=x_1x_2+x_3x_4+\cdots+x_{2a-1}x_{2a}.
```

Both multiplication ranks are unchanged by this substitution, by the
duality identity in the construction note. Its ranks on `n` variables are

```math
R_n(a)=2^{n-1}-2^{n-a-1},\qquad C_n(a)=2a+2.              \tag{5}
```

For the first formula, tensor in the two-variable algebra one pair at a
time. The ranks on the active `2a` variables satisfy
`R_(2a)(a)=4^(a-1)+2R_(2a-2)(a-1)`, with zero base rank, giving
`2^(2a-1)-2^(a-1)`. The inactive variables multiply this by `2^(n-2a)`.
Equivalently, each active pair is a two-term square-zero differential;
the homology dimensions multiply, giving the same formula.

For the second formula, `Hq_a` has degree `n-2`. Its multiplication
operator has only three possibly nonzero graded blocks, from input
degrees zero, one, and two. Their ranks are `1`, `2a`, and `1`.
The middle block is the alternating coefficient matrix, up to paired
permutations. Their sum is `2a+2`.

If the quadratic component is the **lowest** component of an
inhomogeneous polynomial, its multiplication rank cannot exceed that
of the whole polynomial. One proof scales all generators by a formal
parameter, divides by the lowest power, and specializes the resulting
matrix; rank cannot increase under specialization. The same observation
applies to the lowest component after `H`. We use no such assertion for
an arbitrary component in the middle of the degree range.

## 3. Two additional rank budgets

For a piece `f_i`, let `2a_i` be the rank of its quadratic component's
alternating matrix, using zero if that component is absent. Similarly,
let `2b_i` be the alternating rank of the quadratic component of `Hf_i`.
The quadratic component of `g`, and also of `Hg=g`, is
`sum_j s_j n_j`, whose alternating rank is `2lambda`.
Rank subadditivity therefore gives

```math
\sum_i a_i\ge\lambda,\qquad\sum_i b_i\ge\lambda.         \tag{6}
```

These are ranks of coefficient matrices, used only as additive budgets.
They are not substituted for tensor-cut ranks.

## 4. A rational certificate at 256 labels

Set `lambda=8`, `n=16`, and `d=256`. Give weight one to the seven even
degrees

```math
\mathcal D=\{0,4,6,8,10,12,16\},                         \tag{7}
```

and weight zero to degrees two and fourteen. Let `k_i` count the degrees
in (7) where the piece `f_i` has a nonzero component. The exact scalar
certificate is

```math
\sqrt{r_i c_i}\ge d\left(k_i+\frac{a_i+b_i}{5}\right).    \tag{8}
```

Here is a finite, complete verification of (8). Let `[p,q]` be the
piece's supported degree endpoints. Enumerate `0<=p<=q<=16` and
`0<=a,b<=8`, with these consistency rules:

- `a=0` unless `p<=2<=q`; if `p=2` or `q=2`, require `a>=1`.
- `b=0` unless `p<=14<=q`; if `p=14` or `q=14`, require `b>=1`.

This permits additional classes that no polynomial realizes, which
only makes the verification stronger. Start with

```math
r_0=2^{16-p},\qquad c_0=2^q.                            \tag{9}
```

Use (5) only in the following endpoint cases, replacing a bound by the
maximum of itself and the listed value:

| Condition | Additional lower bound |
|---|---|
| `p=2` | `r_0>=R_16(a)` |
| `q=14` | `c_0>=R_16(b)` |
| `q=2` | `c_0>=2a+2` |
| `p=14` | `r_0>=2b+2` |

Let `k_0` count every degree in (7) lying between `p` and `q`. This
overestimates `k_i`, because an intervening component may be zero.
For all **1,368** allowed classes, the checker verifies the integer
inequality

```math
25r_0c_0\ge d^2(5k_0+a+b)^2.                           \tag{10}
```

All terms are nonnegative. Equations (9)--(10) therefore prove (8)
for every possible piece, including those with odd components.
Summing (8), using the seven required degrees in (7), and applying
(6) gives

```math
S\ge d\left(7+\frac{8+8}{5}\right)=\frac{51d}{5}.
```

The exact strict comparison in (2) is checked by
`(51d/5)^408>d^577`, using `577^2-2*408^2=1`.
No numerical optimization is part of the certificate.

At `lambda=9`, the same weights, omitting degrees two and sixteen,
remain valid on all `1,746` scalar classes, but give only
`S>=512*(58/5)=5939.2<512^sqrt(2)`. That exact negative control records
the limit of this certificate; it does not assert a realizable seed.

## 5. Supplied formal limits

Suppose `sum_i f_i(t)` is integral over `F[[t]]` and specializes to `g`,
while each individual polynomial has Laurent-series coefficients.
Apply the monomial and quadratic arguments over the field `F((t))`.
Every required even component is nonzero generically, since its sum has
a nonzero specialization. The summed quadratic coefficient matrices
have full generic rank `2lambda`, since they specialize to the
nondegenerate quadratic form. Thus the degree coverage and both budgets
in (6) still hold. The same proofs give (1)--(2) using generic cut ranks.

This extension keeps every piece in the fixed bi-exterior algebra. It
does not import the conclusion into arbitrary formal crossing families.

## 6. A limit of one-bit common-coefficient channel walks

The [channel-walk transfer in PR26, pinned at `ae5115b3`](https://github.com/1337isnot1337/integer-factoring/blob/ae5115b3d6b4175693eb0899288b71d7b9b4629f/research/twisted-channel-walk-transfer.md)
suggests lifting a one-bit scalar dictionary to matrix coefficients.
Consider the specific construction

```math
\sum_v\Phi(f_v)M_v=E_2 I_b,\qquad M_v\in\operatorname{Mat}_b(F),
                                                               \tag{11}
```

where `f_v` are arbitrary one-bit exterior polynomials. In both cut
factorizations, put the entire constant matrix `M_v` in each left
channel and keep the right channels scalar. Prune channel words only
by certified zero products of the `M_v`, including longer zero words
if desired. In this specific construction, **no such walk-count
certificate has asymptotic base below four**, for any matrix size `b`.

Let `e_v` count the nonzero coefficients of `1` and `sn` in `f_v`, so
`e_v` is zero, one, or two. Equation (4) implies

```math
w_v:=\sqrt{r_vc_v}\ge2e_v.                              \tag{12}
```

For every binary word of length `k` in the algebra letters `1,sn`, the
corresponding matrix coefficient of the ordered power of (11) is `I_b`.
At least one dictionary word contributing to it therefore has nonzero
full matrix product. A dictionary word can contribute to at most
`prod_j e_(v_j)` binary words. Let `W_k` be the dictionary words with
nonzero full matrix product. Union coverage and (12) give

```math
2^k\le\sum_{v\in W_k}\prod_j e_{v_j},\qquad
\sum_{v\in W_k}\prod_j w_{v_j}\ge4^k.                   \tag{13}
```

Every sound zero-product pruning rule must retain these words. For
arbitrary positive reciprocal weights `z_v`, its row and column walk
counts are at least

```math
P_R\ge\sum_{v\in W_k}\prod_j(r_{v_j}z_{v_j}),\qquad
P_C\ge\sum_{v\in W_k}\prod_j(c_{v_j}/z_{v_j}).
```

Cauchy--Schwarz therefore gives `sqrt(P_R P_C)>=4^k`. Any fixed
constant-prefactor exponential upper certificate for those counts
must have base at least four. For identical adjacency graphs this
also proves `spr(A diag(w))>=4`, and the symmetric choice
`z_v=sqrt(c_v/r_v)` attains that common weighting.

This obstruction concerns the stated common matrix coefficients and
their zero-product channel bounds. It does not exclude different
channel-dependent factorizations, extra rank cancellations after
type aggregation, or larger local exterior blocks. The fifteen
nonzero binary one-bit polynomials are included in the local exact
controls for (12); the coverage argument proves the all-field and
all-matrix-size statement.

## Replay and attribution

```sh
python3 tools/check_bi_exterior_degree_budget.py \
  --output evidence/bi-exterior-degree-budget.json
```

The standard-library checker verifies the integer scalar certificate,
the exact target comparisons, every nonzero polynomial on two and four
square-zero generators, and the quadratic normal-form ranks through ten
generators. Finite controls supplement the all-field ordinary proof.
The argument and exact replay used substantial OpenAI assistant support.
