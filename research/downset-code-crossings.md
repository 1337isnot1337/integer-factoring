# Nonlinear physical codes and an exact downset normal form

## Result and scope

Restrict the elementary characteristic-two crossing pencil to any
nonempty word code `C subset {0,1}^m`. Its value at one is still the
crossing on `C`. This changes the physical alphabet and is outside
the full-word hypothesis of the
[binary site-symmetry barrier](binary-site-symmetry-barrier.md).

When `C` is closed under deleting ones, there is an elementary exact
normal form for the **whole column polynomial matrix**:

> Its global Smith diagonal over `F[z]` consists of
> `z^(|x|+|u|)` for every ordered pair `x,u in C`, over every field of
> characteristic two. A separate **constant** basis change expresses
> arbitrary coefficient recombinations as a direct sum of weighted
> disjointness matrices on smaller downsets.

The second assertion concerns actual ranks, not just their Smith upper
bounds. It gives a structural way to test nonlinear codes. It does not
by itself give a favorable uniform rank upper bound. Full finite
optimization for no-adjacent-ones codes and several weight slices has
not produced a post-half rank-sum certificate. The factoring exponent
remains **1/2**; no Lean formalization or independent review is claimed.

## 1. The exact restricted pencil

At a single bit, let

```math
I(x,y,u,v)=[x=u][y=v],\qquad E(x,y,u,v)=[x=v][y=u],
\qquad P(z)=I+z(E+I).
```

For complete words in `C`, define `G_C(z)` by multiplying the `m`
local entries. Then `G_C(1)=E_C`. A nonzero polynomial entry must have
the same local multiset on `(x,y)` and `(u,v)` at every position.
Equivalently,

```math
x\cap y=u\cap v,\qquad x\mathbin\triangle y=u\mathbin\triangle v.
```

For such a tuple put `s=|x triangle y|`, `b=|x triangle u|`. Its entry is

```math
G_C(z;x,y,u,v)=z^b(1+z)^{s-b}.                      \tag{1}
```

This formula gives a direct sparse construction for arbitrary codes,
including codes that are not downsets. It also splits the row matrix
into blocks indexed by intersection and symmetric difference.

## 2. Column fibers of a downset

Now assume `C` is a downset. In the column cut fix

```math
D=x\mathbin\triangle u=y\mathbin\triangle v,\qquad
a=x\cap D,\qquad b=x\setminus D.
```

The outside labels form the downset

```math
C_a(D)=\{b\subseteq [m]\setminus D:
          a\cup b\in C,\ (D\setminus a)\cup b\in C\}.
```

It is unchanged on replacing `a` by `D\a`. The polynomial column
matrix has a block from the rows of orientation `a` to the columns of
orientation `D\a`, with entries

```math
z^{|D|} H_{C_a(D)}(z),\qquad
H_L(z)_{b,c}=(1+z)^{|b\mathbin\triangle c|}.          \tag{2}
```

All other orientation blocks vanish. Permuting the column blocks
makes this a direct sum. Empty fibers are omitted. Every ordered
pair `(x,u)` occurs in exactly one block and one outside label.

## 3. Two exact normal forms for the distance kernel

Let `L` be any finite nonempty downset, and write `w=1+z`.
The polynomial unitriangular matrix

```math
U(z)_{b,s}=[s\subseteq b]w^{|b|-|s|}
```

satisfies

```math
H_L(z)=U(z)\operatorname{diag}_{s\in L}(z^{2|s|})U(z)^T.
                                                               \tag{3}
```

To check an entry, factor out `w^(|b triangle c|)` and sum over all
subsets `s` of `b intersection c`. The remaining binomial sum is
`(w^2+z^2)^(|b intersection c|)=1` in characteristic two. All required
subsets belong to `L`. Since `det U=1`, (3) is a unimodular polynomial
equivalence to a pure-power diagonal, without a local-series truncation.

For actual coefficient ranks use a **constant** transform instead.
Let `Z_(b,s)=[s subset b]`. Downward closure and characteristic two
give `Z^2=I`. Direct subset summation yields

```math
(Z H_L(z) Z^T)_{b,c}
 =[b\cap c=\varnothing]z^{|b|+|c|}.                 \tag{4}
```

Indeed, a coordinate belonging to exactly one of `b,c` contributes
`1+w=z`; a coordinate in their intersection contributes
`1+w+w+1=0`. Thus (4) is a constant invertible change of basis and
commutes with every polynomial coefficient functional.

These two transformations have different uses: (3) proves the global
Smith factors, while (4) preserves the ranks of the actual recombined
coefficient matrices.

## 4. Consequences for the entire coefficient span

Combining (2) and (3), the diagonal orders of the full column matrix are

```math
|D|+2|b|=|x|+|u|\qquad ((x,u)\in C^2).
```

In particular, if `n=|C|`, then

```math
\det\operatorname{Col}(G_C(z))
 =z^{2n\sum_{x\in C}|x|}.                           \tag{5}
```

There are no other irreducible determinant factors. All identities
remain true over every characteristic-two extension field.

For arbitrary functional weights `alpha_j`, write
`A=sum_j alpha_j [z^j]G_C(z)`. Formula (4) gives the exact rank formula

```math
c(A)=\sum_{D,a}
 \operatorname{rank}\left(
 [b\cap c=\varnothing]\,
 \alpha_{|D|+|b|+|c|}
 \right)_{b,c\in C_a(D)}.                            \tag{6}
```

Thus a potentially large `n^2`-order matrix splits into smaller
disjointness matrices with occupation-indexed coefficients. The same
`alpha` is used in every block; independently optimizing each block
would not be a valid coefficient recombination.

Downward closure is essential. For example `C={1}` is not a downset.
Its polynomial kernel is the one-by-one matrix `(1)`, whereas the
invalid use of (5) would predict `z^2`. Constant-weight slices require
their original matrices unless a separate argument is supplied.

## 5. No-adjacent-ones codes

Take `C_m` to be binary words with no adjacent ones. In a nonempty row
block, a coordinate common to `x,y` cannot neighbor any occupied
coordinate. On each run of `S=x triangle y`, the two words must
alternate. Therefore a run of length `ell` has two orientations and
contributes the two-by-two polynomial matrix

```math
\begin{pmatrix}(1+z)^\ell&z^\ell\\
                z^\ell&(1+z)^\ell\end{pmatrix}.       \tag{7}
```

The full row block is the tensor product of (7) over the runs. This
description and (6) retain actual coefficients and their cancellations.
They can support more efficient searches than dense full-word matrices.

For `m>=1`, the extreme coefficients have simple all-length formulas:

```math
(r([z^0]G_{C_m}),c([z^0]G_{C_m}))=(|C_m|^2,1),
\qquad
(r([z^m]G_{C_m}),c([z^m]G_{C_m}))=(1,4).             \tag{8}
```

For the second formula, the local leading coefficient vanishes unless
`x,y` are complementary, and likewise `u,v`. A word and its complement
both avoid adjacent ones exactly when they alternate, leaving two
words. The row block is an all-ones matrix of order two; the column
matrix is a permutation of order four. A constant-cost leading piece
therefore coexists with the larger costs in the middle coefficients.

There is an effective candidate for every desired alphabet size.
The counts satisfy `|C_m|=|C_(m-1)|+|C_(m-2)|`, by splitting words
according to their initial zero or initial `10`. They grow exponentially
in `m`. Choose the least `m` with `|C_m|>=N` and retain the first `N`
words in increasing binary order. This prefix is still a downset,
because deleting ones decreases the integer value. For `N=2^lambda`,
`m=O(lambda)`. Exact minmax coefficient-span optimization and fixed-
precision rational rank-sum bounds are computable in `2^O(lambda)`
operations. This effectiveness statement is not a
uniform favorable-rank theorem or a factoring cost proof.

## 6. Exact evidence and remaining question

Run

```sh
make verify-downset-crossings
make verify-regular-code-crossings
```

The [normal-form checker](../tools/check_downset_crossing_normal_form.py)
examines all **167** nonempty downsets on four bits. It compares
1,322,030 literal entry polynomials, 3,866 actual coefficient-rank
queries, 24,604 local polynomial identities, and 12,302 predicted
Smith orders. These finite controls support the all-downset proof
above; they are not the proof of its arbitrary-size statement.

The [regular-code search](../tools/search_regular_code_crossings.py)
enumerates full `F_2` coefficient spaces for `C_m`, `1<=m<=7`, four
constant-weight slices at lengths eight and ten, and three exact-size
prefixes: **14** cases and **1,028** coefficient functionals in total.
It saves every rank pair, exact optimal maximum-rank-product
decompositions and separating covectors, and rational rank-sum bounds
with dual potentials. The normal form independently checks every
column-rank query in the downset cases.
Every saved rank-sum lower bound exceeds `|C|^(3/2)`, so the simple
repeated-block rank-sum envelope fails the half-rate threshold in all
14 finite cases. This comparison is checked using exact integer squares.

For scale, the full no-adjacent-ones codes at lengths four and seven
have sizes eight and 34 and optimal rank sums approximately 50.2492
and 772.2361. These are worse than the simple repeated-block target.
The [normal-form artifact](../evidence/downset-crossing-normal-form.json)
and [complete finite optimization](../evidence/regular-code-crossings.json)
state the exact tested scopes. No asymptotic exclusion of all nonlinear
codes, or all tensor lengths of this language, is asserted.

The open structural question is whether a different downset, a
different nonlinear language, or a new local polynomial can exploit
the disjointness-block formula to obtain a uniform rank **upper**
bound below the required threshold. A smaller Smith profile alone is
insufficient, by the [universal bare criterion floor](multipoint-smith-floor.md).

This construction and its proof used substantial Codex AI assistance.
