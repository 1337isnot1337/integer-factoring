# Amplify an exact finite crossing block to every label width

## Result

An exact finite block can certify a uniform crossing family. The cost of
grouping its tensor-power pieces is explicit. This supplies an alternative
to proving a continuous adaptive-capacity estimate at every entropy level.

Let `F` be a fixed explicitly represented finite field. Let `X,Y` be
finite label sets of sizes `d_X,d_Y>=2`, and suppose finite tensors satisfy

```math
\sum_{i=1}^{L} A_i(x,y,u,v)
   =\mathbf1_{x=v}\mathbf1_{y=u},
\qquad x,v\in X,\quad y,u\in Y.                         \tag{1}
```

Discard zero tensors. For each remaining piece let `r_i,c_i` be positive
integer upper bounds for its ranks across `(x,y)|(u,v)` and
`(x,u)|(y,v)`. Define

```math
d=d_Xd_Y,\qquad M=\max_i r_i c_i.
```

There is an effective all-width family with polynomially many pieces and
balanced parameter at most

```math
\rho_0=\frac{\log(M L^2)}{2\log d}.                     \tag{2}
```

The complete rank table gives the sharper parameter

```math
S=\sum_i\sqrt{r_i c_i},\qquad
\rho_* =\frac{\log S}{\log d}\le\rho_0.                \tag{2a}
```

Thus `S^2<d^(sqrt(2))` also certifies the strict threshold. It can
hold when the simpler worst-piece test below does not.

The rank bound has an additive constant at widths not equal to a power
of `d`. Consequently every strict upper bound on `rho_0` is a uniform
asymptotic crossing parameter. In particular, the exact finite condition

```math
M L^2<d^{\sqrt2}                                       \tag{3}
```

supplies a balanced crossing below `1/sqrt(2)`. The
[ordinary modular transfer](asymmetric-modular-transfer.md) then gives
some logarithmic saving exponent above one half, with the field and
finite block fixed before the factoring input is read.

No block satisfying (3) is asserted here. The result is a constructive
certification criterion, not a new factoring exponent.

## 1. Balance the external label sets

Use a second copy with its paired species exchanged. For `i,j` put

```math
B_{ij}\bigl((x,x'),(y,y'),(u,u'),(v,v')\bigr)
 =A_i(x,y,u,v)A_j(y',x',v',u').                          \tag{4}
```

Here `x',v'` belong to `Y`, and `y',u'` belong to `X`. Fix the obvious
coordinate-swap identifications between `X x Y` and `Y x X` so that
all four external alphabets use one labelled set of size `d`.
Equation (1) implies

```math
\sum_{i,j}B_{ij}=C_d,                                  \tag{5}
```

where `C_d` is the literal crossing tensor on that label set.

The first cut of (4) is a tensor product of the two first cuts, up to
row and column permutations. In the second cut the second factor is
transposed. Hence

```math
\operatorname{rank}_R(B_{ij})\le r_i r_j,
\qquad \operatorname{rank}_C(B_{ij})\le c_i c_j.          \tag{6}
```

There are `J=L^2` pieces. This explicit copy is essential when the two
original label sizes differ. It also gives the same bound (2) when
they are equal.

## 2. Group words by their composition

Take the `h`th tensor power of (5). For each vector
`v=(v_1,...,v_J)` of nonnegative integers summing to `h`, let

```math
Q_v=\sum_{\substack{j_1,\ldots,j_h\\
                    \#\{a:j_a=j\}=v_j\ \forall j}}
       B_{j_1}\otimes\cdots\otimes B_{j_h}.             \tag{7}
```

The sum over all compositions is exactly `C_(d^h)`. The number of
compositions is

```math
\binom{h+J-1}{J-1}\le(h+1)^J,                          \tag{8}
```

which is polynomial in `h` because the one finite block, and therefore
`J`, is fixed.

Write `R_j,C_j` for the bounds in (6), and
`N_v=h!/(v_1!...v_J!)` for the number of words in (7). Rank subadditivity
and tensor multiplicativity give

```math
\operatorname{rank}_R(Q_v)\le N_v\prod_jR_j^{v_j},
\qquad
\operatorname{rank}_C(Q_v)\le N_v\prod_jC_j^{v_j}.        \tag{9}
```

In particular, using `N_v<=J^h=L^(2h)`,

```math
\operatorname{rank}_R(Q_v)\operatorname{rank}_C(Q_v)
 \le N_v^2 M^{2h}
 \le(M^2L^4)^h.                                        \tag{10}
```

Every external label has `h log d` units of entropy. Dividing the
logarithm of (10) by `4h log d` proves (2). No balance between the two
individual cut ranks is assumed. Their summed logarithmic bound is
the balanced support-line condition used by the transfer theorem.

For the sharper bound, the multinomial theorem gives

```math
N_v\prod_j\sqrt{R_jC_j}^{\,v_j}
 \le\left(\sum_j\sqrt{R_jC_j}\right)^h
 =\left(\sum_i\sqrt{r_ic_i}\right)^{2h}=S^{2h}.          \tag{10a}
```

Squaring bounds the product of the two ranks by `S^(4h)` and proves
(2a). This retains each piece's actual rank charge before summation.

## 3. The decomposition and factors are effective

All data in (1) are finite. Exact field elimination supplies a rank
factorization of every cut of every `A_i` whenever only rank bounds
were initially recorded. These arrays are computed once and remain
fixed for the entire uniform family.

For a tensor word in (7), tensor the appropriate local factorizations.
To factor their sum, concatenate the first factor matrices and stack
the second ones. The intermediate index is a pair consisting of:

1. a word of the specified composition;
2. the tuple of its local intermediate rank indices.

The word can be encoded by its lexicographic number among the `N_v`
words. Multinomial counts decode that number in polynomial time in
`h`, with constants depending on `J`. The local rank tuple uses an
ordinary mixed-radix encoding. A factor entry is then a product of
`h` fixed local field entries. Padding the intermediate indices to
binary lengths changes each cut bound by at most one bit.

The coefficients of `Q_v` are also effective. For a requested external
tuple, extract the coefficient of `z^v` from

```math
\prod_{a=1}^{h}\left(\sum_{j=1}^{J}z_j B_j(e_a)\right).
```

A histogram dynamic program has at most `(h+1)^J` states. Thus all
piece and factor queries have polynomial bit cost in `h` for the fixed
finite certificate. Explicitly writing their exponentially many
possible entries, when required, costs `exp(O(h))`, with a fixed
block-dependent constant. This is within the effective crossing
interface; no table-query oracle is introduced.

For an arbitrary binary label width `lambda`, choose

```math
h=\left\lceil\frac{\lambda\log2}{\log d}\right\rceil
```

by comparing integer powers with `2^lambda`. Restrict all four external
alphabets to the first `2^lambda` labels. This keeps the literal
crossing endpoint and cannot increase either cut rank. The difference
between `h log d` and `lambda log2` is bounded by the fixed constant
`log d`. The number of pieces remains polynomial in `lambda`, and the
logarithmic rank bounds gain only an additive constant.

The transfer interface also requires each individual intermediate
width to be at most `2 lambda`. If the concatenated factorization is
larger, use the trivial factorization of that cut as its identity
matrix times the complete piece matrix. The histogram algorithm
already evaluates the latter matrix. Choose the smaller of these two
factorizations independently at each cut. This caps both widths at
`2 lambda` and can only improve their summed bound. The fixed additive
constant is bounded by `c log_2(2 lambda)` for a block-dependent `c`
at every positive width. All of the transfer's width and query
conditions therefore hold uniformly, including the small widths.

## 4. Exact arithmetic certificates

There is no need to certify a floating logarithm in (3). It suffices
to find positive integers `a,b` with

```math
a^2<2b^2,
\qquad (M L^2)^b<d^a.                                 \tag{11}
```

Both tests are exact integer inequalities. They give
`rho_0<a/(2b)<1/sqrt(2)` with a strict rational margin.
In the cited transfer theorem take its two rational support weights
to be `(1,1)` and its rational support parameter to be `a/b`.
The resulting intermediate widths satisfy
`r+c <= 2(a/b) lambda+c_0 log_2(2 lambda)`, with `a/b<sqrt(2)`.
Thus the transfer uses rational input data; no irrational weight is
silently supplied to its effective premise.

For the sharper bound, choose a positive rational upper bound
`U/V >= S`. It suffices to check

```math
a^2<2b^2,\qquad U^{2b}<V^{2b}d^a.                     \tag{11a}
```

For example, for any positive integer `V`, set
`U=sum_i ceil(sqrt(r_i c_i V^2))`. Integer square roots compute this
bound without floating arithmetic. Increasing `V` makes it converge
to `S`; every strict instance of the sharper criterion eventually
has a certificate of the form (11a).

There is an asymmetric version. Fix positive `a,b` with `a^2+b^2=1`,
and put

```math
W=\max_i r_i^a c_i^b.
```

The same construction gives normalized support at most

```math
R_0=\frac{\log W+(a+b)\log L}{\log d}.                 \tag{12}
```

The complete-table version is

```math
R_* =\frac{(a+b)\log\left(\sum_i
                 r_i^{a/(a+b)}c_i^{b/(a+b)}\right)}{\log d}
 \le R_0.                                             \tag{12a}
```

Indeed, apply the multinomial theorem to
`N_v prod_j (R_j^a C_j^b)^(v_j/(a+b))` and then raise to `a+b`.
The paired sum again factors as the square of the original sum.
For rational weights, rational upper bounds for the finitely many
positive algebraic roots in (12a) are obtained by exact integer-power
comparisons. A strict `R_*<1` therefore also has a finite exact
certificate and the same rational transfer interface.

A strict `R_0<1` meets the unit quarter-circle support criterion.
For rational Pythagorean weights `a=A/C,b=B/C`, an exact certificate is

```math
A^2+B^2=C^2,
\qquad r_i^A c_i^B L^{A+B}<d^C\quad\text{for every }i.  \tag{13}
```

The strict integer inequalities permit an effective choice of a
rational support parameter between `R_0` and one, using rational
bounds for the logarithms. The weights `A/C,B/C` are already rational.

The field, pieces, ranks, and integer certificate must all be fixed and
effective. The transfer's remaining uniform machine obligations are
discharged by Section 3 and the cited transfer theorem.

## 5. Application to a finite adaptive block

Suppose an exact typed tensor block of length `m` has label counts
`d_X,d_Y`, and its opposite-flag decomposition has at most `D m+1`
nonzero pieces. Compute valid finite rank bounds for every piece,
including endpoint pieces, zero jets, full-dimension caps, and every
positive constant-prefix cost. Then (11) or (13) is a complete
finite target for amplification.

When `log(d_Xd_Y)=mH+O(log m)` and `L=O(m)`, the additional balanced
cost in (2) is

```math
\frac{\log L}{\log(d_Xd_Y)}=O(\log m/m).                \tag{14}
```

The polynomial degree in (8) may be enormous. It is nevertheless a
constant after one finite block is certified. This is an asymptotic
existence route, and supplies no practical factoring speed claim.

A finite rank table alone is insufficient if its pieces do not
actually sum to the strict typed crossing, or if its local rank
estimates ignore a prefix or physical-type cost. Conversely, once
those exact finite obligations are proved, no continuous
capacity-curve certificate is needed for this amplification step.

## 6. Exact controls

The companion checker constructs real polynomial crossing coefficients,
groups their tensor-square words, and verifies both complete cut
factorizations and their multinomial bounds. It also checks the
rectangular pair swap on an explicit two-by-three restriction.

Run

```sh
python tools/check_finite_block_amplification.py \
  --output evidence/finite-block-crossing-amplification.json
```

The [checker](../tools/check_finite_block_amplification.py) passed all
30 complete grouped factorizations, covering 7,680 matrix entries;
3,840 independent histogram queries; 162 complete paired cut ranks;
104,976 paired tensor entries; and 1,552 crossing endpoint entries.
It checks both the worst-piece and rational square-root-sum bounds.
The [evidence](../evidence/finite-block-crossing-amplification.json)
binds the exact input and checker hashes and retains every grouped
composition, rank bound, and actual rank. These are construction
controls; none of these sample decompositions is claimed to satisfy
the strict factoring threshold.
