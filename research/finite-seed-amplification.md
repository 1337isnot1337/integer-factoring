# Amplifying a finite crossing certificate

## Result and scope

A finite exact crossing decomposition can give a uniform factoring bound
without extending its original operator construction to larger weights.
The sufficient quantity is a sum of square roots of rank products.

Let `d>=2`, let `E_d(x,y,u,v)=[x=v][y=u]`, and suppose

```math
E_d=\sum_{i=1}^k A_i
\quad\text{over a fixed finite field }\mathbf F_p.
```

Let `r_i,c_i` be the ranks of the `(x,y)|(u,v)` and `(x,u)|(y,v)`
flattenings. Omit zero tensors and put

```math
S=\sum_{i=1}^k\sqrt{r_i c_i}.
\tag{1}
```

If `S<d^2`, there is an effective scale-one crossing family at every
binary label length with

```math
\rho=\frac{\log S}{2\log d},\qquad
\text{width sum}\le4\rho\lambda+O(1).
\tag{2}
```

The family has polynomially many pieces. In particular, a verified seed
with `S<d^(sqrt(2))` would give a factoring exponent strictly above `1/2`
through the existing [fixed-field transfer](scale-one-modular-transfer.md).

There is an unconditional limitation on this particular certificate:

```math
S\ge2\sqrt{d^2-1}.
\tag{3}
```

Consequently **no seed on two, three, four, or five labels can certify a
strict half-exponent improvement by (1)**. This does not preclude a
better analysis of its tensor coefficients, a restriction of labels, or
another construction. A second limitation below concerns only a common
Smith-form estimate.

The later [six-label argument](six-label-rank-sum-obstruction.md)
extends the first exclusion to six labels. The
[streamed refinement](streamed-seed-transfer.md) constructs polynomial-time
factor queries for the amplified family and excludes the usual Smith
estimate on four labels as well. The Gaussian-elimination construction
below remains a valid, less space-efficient implementation.

The final section records useful structure of the newer three-label
certificate in [PR #24](https://github.com/1337isnot1337/integer-factoring/pull/24),
pinned here to commit `e96ee50560284301428d72421330d8773155d25b`.
Its finite decomposition is prior work from that branch. The amplification
theorems and the two common triangular flags are the contributions of
this note. None gives a new factoring exponent by itself.

## 1. Group tensor copies by their piece counts

Take `t` tensor copies of the seed. For each vector of nonnegative
integers `n=(n_1,...,n_k)` with `sum n_i=t`, define

```math
B_n=\sum_{\substack{w\in\{1,\ldots,k\}^t\\
                    |w^{-1}(i)|=n_i\ \forall i}}
       A_{w_1}\otimes\cdots\otimes A_{w_t}.
\tag{4}
```

External coordinates are the words in `[d]^t`. Each word of piece
indices occurs exactly once, so `sum_n B_n=E_(d^t)`. There is no division
by a factorial, and the identity is valid in every characteristic.
There are `binom(t+k-1,k-1)` possible count vectors.

The tensor-product rank rule and subadditivity give

```math
\begin{aligned}
\operatorname{rank}_{\rm row}B_n
&\le {t\choose n_1,\ldots,n_k}\prod_i r_i^{n_i},\\
\operatorname{rank}_{\rm col}B_n
&\le {t\choose n_1,\ldots,n_k}\prod_i c_i^{n_i}.
\end{aligned}
```

The multinomial theorem, applied to the positive real numbers
`sqrt(r_i c_i)`, therefore yields

```math
\operatorname{rank}_{\rm row}B_n\,
\operatorname{rank}_{\rm col}B_n
\le
\left({t\choose n_1,\ldots,n_k}
      \prod_i(r_i c_i)^{n_i/2}\right)^2
\le S^{2t}.
\tag{5}
```

The real numbers in this bound are used only in the proof. All tensor
entries and factors are computed over the fixed base field.

## 2. Binary labels, construction, and exact factors

For a positive binary label length `lambda`, choose the least `t` with
`d^t>=2^lambda`, and restrict all four external coordinates to the first
`2^lambda` words in a fixed order on `[d]^t`. The crossing sum remains
exact. Both ranks can only decrease. Minimality gives

```math
2^\lambda\le d^t<d\,2^\lambda,
\qquad t\log_2d=\lambda+O(1).
```

Compute each restricted piece matrix and its two actual ranks by exact
Gaussian elimination. A nonzero rank-`r` matrix has an exact factorization
through `r` coordinates; pad that factorization to `2^ceil(log2 r)`
coordinates. Zero pieces may be omitted. Thus each separate binary width
is at most `2lambda`, and (5) gives

```math
\lceil\log_2\operatorname{rank}_{\rm row}B_n\rceil+
\lceil\log_2\operatorname{rank}_{\rm col}B_n\rceil
\le2t\log_2S+2
=4\rho\lambda+O(1).
```

The seed and `k` are fixed constants. An entry of (4) is the coefficient
of `X_1^n_1...X_k^n_k` in a product of `t` linear forms, one per external
coordinate position. A coefficient dynamic program has at most
`(t+1)^k` states. Consequently the entire index set, all complete piece
tables, both ranks, and both factorizations can be constructed in
`poly(t)*(d^t)^O(1)=2^O(lambda)` bit operations. Indices need `O(log lambda)`
bits. No oracle for an uncomputed factorization is assumed.

All coefficients belong to `F_p` and the crossing multiplier is one.
The identities extend to every commutative `F_p`-algebra by scalar
extension. This verifies the exact effective-pieces interface used in
the fixed-field transfer, with (2).

An approximate numerical value of `S` is not itself a strict certificate.
For a positive claim one may, for example, prove `S<d^(2q)` for a rational
`q<1/sqrt(2)`. This is an inequality between computable algebraic numbers
after raising to the denominator of `q`.

## 3. A universal lower bound on the rank sum

Put `N=d^2`. The crossing has rank `N` in both flattenings, so

```math
\sum_i r_i\ge N,\qquad \sum_i c_i\ge N.
```

Every retained piece has `r_i,c_i>=1`. If `k>=N`, then `S>=k>=N`.
If `k=1`, both ranks equal `N` and `S=N`.

Suppose `2<=k<=N`. To lower-bound `S`, relax the ranks to positive real
numbers and decrease them, if necessary, so that both sums equal `N`
while every coordinate stays at least one. This can only decrease `S`.
The function `sum_i sqrt(r_i c_i)` is concave on the product of the two
resulting simplexes, so its minimum is attained at a vertex. At a vertex
one `r` coordinate equals `N-k+1` and the others equal one; the same is
true independently of the `c` coordinates.

If the two large coordinates coincide, the objective is `N`. Otherwise
it is

```math
2\sqrt{N-k+1}+k-2.
```

This is nondecreasing for integral `k` from `2` to `N`. Its minimum is
`2sqrt(N-1)`, attained in the relaxation at `k=2`. Also
`N>=2sqrt(N-1)`. This proves (3) for every exact crossing decomposition.
The argument does not assert that the relaxed minimum is realizable by
actual tensors.

For `d=2,3,4,5`, the checker verifies the exact integer inequalities

```math
\bigl(4(d^2-1)\bigr)^{70}>d^{198}.
```

Since `sqrt(2)<99/70`, these imply
`2sqrt(d^2-1)>d^(sqrt(2))`. At `d=6` this obstruction no longer excludes
the target; that is only a necessary-condition check.

## 4. A separate limitation of the usual Smith estimate

Consider any fixed polynomial tensor curve `P(z)` on `d` labels with
`P(1)=E_d`. Remove a common power of `z`, so its constant and leading
coefficients are both nonzero. Let its degree be `D>=1`, and let
`v_1,...,v_N` be the `z`-adic Smith valuations of the row matrix. Let
`w_1,...,w_N` be those of the reversed column matrix
`z^D P_col(1/z)`. Both matrices are invertible over `F_p(z)` because their
values at one are invertible. Each valuation list has a zero and has sum
at most `ND`, by the determinant degree bound.

The usual cumulative Smith estimate bounds the rank of coefficient `z^q`
of `P(z)^(tensor t)` by a polynomial factor times the number of Smith
state words with total valuation at most `q`. To see the polynomial
factor, write a Smith decomposition over `F_p[[z]]`; for each state word
of valuation at most `q`, its two outer series contribute at most `q+1`
rank-one coefficient products. The corresponding reverse-column
threshold is `Dt-q`.

Normalize the valuations by `D`. The exponential rate of the cumulative
count with mean-cost bound `u` is

```math
\Phi_v(u)=\inf_{0<z\le1}
 \left(\log\sum_j z^{v_j/D}-u\log z\right).
```

The worst combined row and column envelope is

```math
\max_{0\le u\le1}\bigl(\Phi_v(u)+\Phi_w(1-u)\bigr)
=\inf_{0<z\le1}\log
 \frac{(\sum_jz^{v_j/D})(\sum_jz^{w_j/D})}{z}.
\tag{6}
```

For completeness, an interior minimizing tilt has the sum of its two
expected costs equal to one; allocate `u` to the first expectation.
At the endpoint `z=1`, the two unrestricted mean costs sum to at most
one, so an allocation attaining both unrestricted entropies exists.
These observations prove (6), including the endpoint case.

For `0<z<=1`, convexity and the determinant bound imply

```math
\sum_jz^{v_j/D}\ge1+(N-1)z^{N/(N-1)},
```

and the same holds for `w`. Thus (6) is at least

```math
\inf_{0<z\le1}
 \log\frac{(1+(N-1)z^{N/(N-1)})^2}{z}
=2\log(2N)-\frac{N+1}{N}\log(N+1).
```

Accordingly this particular Smith envelope cannot certify a contraction
smaller than

```math
b(d)=\frac{2\log(2d^2)-((d^2+1)/d^2)\log(d^2+1)}
           {4\log d}.
\tag{7}
```

For three labels, `b(3)=0.7332695228...>1/sqrt(2)`. This is a bound on
the estimate, **not a lower bound on the actual coefficient ranks**.
It leaves open sharper cancellations, multigraded estimates, and label
restrictions. In particular it does not rule out a useful amplification
of the seed studied next.

## 5. Common flags of the twenty-piece three-label seed

The pinned PR #24 witness uses labels `(1,11),(3,9),(5,7)` and twenty
pieces. Its rank-pair multiset is five pairs with product three, two with
product two, and thirteen with product four. Hence its rank sum is

```math
S=26+5\sqrt3+2\sqrt2=37.48868116\ldots.
```

This is too large for (2), despite its small maximum rank product.
The new calculation reveals a different structure. Write its row
matrices as `A_i` in the order used by the pinned witness, numbered
`0,...,19`. Multiply every matrix on the right by the crossing permutation
`E_3`, so their sum becomes the identity. Do the same after reshuffling
to the column cut.

Both operator families are simultaneously upper triangular over `F_2`.
Here are explicit ordered bases, with each vector encoded by its nine
binary coordinates as an integer:

| Cut | Ordered basis |
|---|---|
| Row | `27,36,192,10,256,3,4,64,1` |
| Column | `27,192,36,10,3,256,64,4,1` |

The diagonal of each transformed operator contains only zeroes and ones.
Only four piece indices have any nonzero diagonal entries:

| Cut | Diagonal character multiplicities |
|---|---|
| Row | index `16`: four; `17`: two; `19`: two; `12`: one |
| Column | index `15`: four; `17`: two; `19`: two; `12`: one |

In particular the determinant of the normalized multivariate pencil is
`X_16^4 X_17^2 X_19^2 X_12` in the row cut and
`X_15^4 X_17^2 X_19^2 X_12` in the column cut.
The generated operator algebras have dimensions thirty and thirty-one,
respectively. These claims concern the selected twenty operators, not
every element of the larger source symmetry span.

The off-diagonal entries are essential. The scalar associated graded
alone gives neither an upper bound for the complete ranks nor a
factoring improvement. The two flags provide concrete data for a
subsequent weighted type or label-restriction analysis.

### An actual obstruction for unrestricted type pieces

The three common diagonal characters, indexed by `12,17,19`, have
multiplicities `1,2,2` in both cuts. For a type supported on these three
indices, each tensor-power matrix is upper triangular and has exactly

```math
{t\choose n_{12},n_{17},n_{19}}\,2^{n_{17}+n_{19}}
```

nonzero diagonal entries in each cut. A principal minor on these entries
is triangular with diagonal one, so this is a lower bound for both
actual ranks, regardless of off-diagonal cancellations.
The displayed quantities sum to `5^t` over at most `binom(t+2,2)` types.
Thus some literal unrestricted type piece has rank product at least

```math
\frac{25^t}{\binom{t+2}{2}^2}.
```

On the full `3^t` label domain this forces
`rho>=log(5)/(2log(3))=0.73248676...` in any uniform width bound for
these particular type pieces. This conclusion does not concern further
recombinations or restrictions of their labels.

### A label restriction that removes that obstruction

Let `U=span(e_0,e_1)` and `w=e_2`. The five common-character basis vectors
in either cut span

```math
(U\otimes w)\oplus(w\otimes U)\oplus\langle w\otimes w\rangle.
\tag{8}
```

Every vector in either stated flag is homogeneous for each of the two
individual indicators of a `w` coordinate. In particular the basis
changes preserve every prescribed number of `w` coordinates in each
external word.

Restrict all four labels to words with exactly `m=floor(t/3)` occurrences
of `w`. The number of labels is

```math
\binom tm\,2^{t-m}=\frac{3^t}{\operatorname{poly}(t)}.
```

For `t>=1` with `2m<t`, two such words cannot have a `w` in at least one
of their coordinates at every position. Hence their pair space contains
no basis word made entirely from (8). This removes all the nonzero
diagonal entries used in the preceding obstruction while losing only
a polynomial factor of label dimension. The crossing restriction stays
exact. It does **not** make the corresponding pieces zero.

The remaining ranks of this original-basis restriction are now analyzed
in [balanced-three-label-barrier.md](balanced-three-label-barrier.md).
An explicit degeneration gives a larger-rank type and rules out a
half-exponent improvement for this particular constant-weight family.
The [permutation-invariant defect theorem](fixed-defect-invariant-barrier.md)
also covers its sparse-defect endpoint, including arbitrary invariant
recombinations there. Neither theorem concerns general label codes made
from linear combinations of words, or arbitrary restrictions that break
the stated position symmetry.

## 6. Reproduce

```sh
python tools/check_finite_seed_amplification.py --output evidence/finite-seed-amplification.json
```

The checker includes the twenty explicit binary tensors, independently
checks their crossing sum and both ranks, verifies both stated flags and
algebra dimensions, checks the individual `U/w` homogeneity of each flag
vector, and constructs every two-copy type piece. The latter
check replays the crossing on nine labels and both rank bounds before
binary restriction. The small-seed exclusion uses the exact integer
inequalities above. The finite checks do not replace either uniform
proof.
