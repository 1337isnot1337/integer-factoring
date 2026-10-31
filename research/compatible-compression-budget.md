# Compatible integer compression cores extend the all-field cutoff to 97

## Result and scope

Every exact crossing decomposition over every field satisfies

\[
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}\qquad(2\le d\le97).
\]

The [rational-cut successor](rational-cut-mesh-budget.md) established the
all-field cutoff through 96 and the characteristic-two cutoff through 97.
This note removes the characteristic restriction at 97. Three small
families of integer Jacobians supply the missing all-field compression
bounds. Exact rational cut thresholds and ratio subintervals then cover
every finite rank class.

This remains a restriction on completed finite crossing seeds. It does
not improve the factoring exponent, which remains `1/2`, and it does not
exclude tensor-power constructions from smaller polynomial sources. The
separate numerical formal-limit cutoff is unchanged.

## 1. Compatible cores with rational allocations

For each `T` in `{3,4,5}`, the literal certificate file supplies three-map
families `A_t` and `D_t`, for every `1<=t<=T`, of sizes

\[
A_t:\quad t\longrightarrow\left\lceil\frac{(2T-1)t}{T}\right\rceil,
\qquad
D_t:\quad\left\lceil\frac{(T+1)t}{T}\right\rceil\longrightarrow t.
\tag{1}
\]

For **every** pair of types `(a,b)`, form the integer Jacobian

\[
\mathcal J_{a,b}(X,Y)
 =\bigl(D_{b,i}X-YA_{a,i}\bigr)_{i=1}^3.
\tag{2}
\]

The certificate gives maximal minors whose integer determinants have
greatest common divisor one. Consequently (2) has full row rank over
every field. This assertion covers all `3^2+4^2+5^2=50` cross blocks,
not only the equal-type pairs.

The largest square core sizes are 27, 48, and 75. Their determinants are
respectively `1,-1,1`. All cross blocks have a displayed unit minor
except the following three denominator-five blocks:

| `A` type | `D` type | Displayed determinants | Greatest common divisor |
| --- | --- | --- | --- |
| 4 | 4 | `-17,-19` | 1 |
| 4 | 5 | `25,-3` | 1 |
| 5 | 4 | `-29,9` | 1 |

Thus 53 exact integer determinant calculations verify the whole
certificate. A nonzero minor over one characteristic is not substituted
for this all-field argument.

### Why gcd one is sufficient

Let an integer `r` by `c` matrix `J` have selected `r` by `r` minors
`M_j` with determinants `d_j`. Embed each adjugate into the corresponding
input coordinates to obtain a `c` by `r` matrix `C_j`. Then
`J C_j=d_j I_r`. If `gcd(d_j)=1`, Bezout coefficients give an integer
linear combination of the `C_j` that is a right inverse of `J`.
Therefore `J` remains surjective after extension to every field.

### From finite cores to arbitrary dimensions

Write an arbitrary input size as a sum of blocks of size `T` and one
remainder in `{1,...,T-1}`. Taking the corresponding direct sum of the
`A_t` families gives total output size
`ceil((2T-1)s/T)`. Independently, direct sums of the `D_t` families with
total output size `q` have total input size `ceil((T+1)q/T)`.

Both `X` and `Y` split into rectangular blocks indexed by one `A` block
and one `D` block. Equation (2) acts separately on every such pair.
The complete Jacobian is therefore a direct sum of the certified
surjective cross blocks. This proves surjectivity for arbitrary `s,q`.
The checker independently constructs a nontrivial two-by-two block
array for each `T`, verifies each of its four cross blocks, and checks
that every off-block coordinate is zero.

Add `a` maps using disjoint `D`-only identity blocks and `b` maps using
disjoint `A`-only identity blocks. For `m=3+a+b`, the resulting all-field
compression conditions are

\[
s+a q+\left\lceil\frac{(T+1)q}{T}\right\rceil\le d,
\]

\[
q+b s+\left\lceil\frac{(2T-1)s}{T}\right\rceil\le d.
\tag{3}
\]

The chart and proper-projection argument in the
[interior-Jacobian note](interior-jacobian-budget.md) now proves that
**every** tuple of `m` endomorphisms admits subspaces `U,W` with
`dim U=s`, `codim W=q`, and `P_i U subset W`, after scalar extension.
The coefficient field is arbitrary because the Jacobians have integer
right inverses.

As before, restriction to `U tensor V` and projection of the first
output factor to `V/W` imply

\[
\operatorname{rank}\left(E_d-\sum_{i=1}^m P_i\otimes Q_i\right)\ge sq.
\tag{4}
\]

The construction also supplies geometric complement bounds for formal
degenerations, but this note computes only the exact finite-seed cutoff.

## 2. The strengthened all-field table at 97

Maximize `sq` over (3), over all three core families and all nonnegative
`a,b` with `a+b=m-3`. Combine these results with the existing all-field
geometric table and the exact one- and two-product bounds. Zero padding
allows the monotone envelope.

At `d=97`, the new strict improvements are:

| Product count `m` | Previous complement bound | New all-field bound |
| --- | --- | --- |
| 6 | 576 | 578 |
| 7 | 462 | 468 |
| 8 | 387 | 390 |
| 9 | 328 | 333 |
| 12 | 230 | 234 |
| 14 | 188 | 192 |
| 17 | 150 | 153 |
| 18 | 141 | 144 |
| 22 | 112 | 114 |

The evidence records an explicit `(T,a,b,s,q)` witness for every entry.
The weighted minima still stop at `max(b,2d-1)`: positive parameters in
(3) force `a,b<=d-3`, so these new families are already zero from
`m=2d-2` onward. The old geometric table is zero from `2d-1` onward.

## 3. Complete rational rank-class coverage

Use the same maximal-ratio classes `(k,b)` and exact subset budgets as in
the rational-cut note. The desired common lower bound is

\[
S\ge\frac{16131}{25}=645.24.
\tag{5}
\]

At dimension 97, there are 15,340 classes. The nonnegative base term
settles 14,825. The integer-threshold formula settles 511 more, and the
quarter-threshold formula settles another three. Only `(k,b)=(52,3)`
requires a finer calculation.

### Refining the last ratio class

If `theta` is the square root of the largest piece rank ratio, divide
`51<theta<=52` into the 100 rational intervals

\[
L_j=51+\frac j{100}<\theta\le
R_j=51+\frac{j+1}{100},\qquad0\le j<100.
\]

For the selected piece, the exact integer rank range becomes

\[
a_0=\lfloor L_j^2 b\rfloor+1,
\qquad a_1=\min(\lfloor R_j^2 b\rfloor,d^2).
\tag{6}
\]

Use `R_j` as the upper ratio endpoint and place square-root interpolation
knots at `1,1+1/16,...` below `R_j`, followed by `R_j` itself. For
consecutive positive rational knots `s_0=1<...<s_t=R_j`, the base
coefficient is `1/(s_(t-1)+s_t)`, and the hinge coefficients are

\[
w_0=\frac12-\frac1{1+s_1},\qquad
w_i=\frac1{s_{i-1}+s_i}-\frac1{s_i+s_{i+1}}.
\]

They are positive and define the exact chord interpolant of `sqrt(t)`
on the knot squares, together with the reciprocal half of the grid.
The proof is the same slope calculation as in the rational-cut note.

The residual parameters are recomputed from (6), and the chord
coefficients become

\[
\mu_s=\frac{s^2R_j^2-1}{R_j^4-1},\qquad
\nu_s=\frac{R_j^2-s^2}{R_j^4-1}.
\]

At each threshold, the selected piece is charged to the minimizing
subset only if `a_0>s^2b`. Otherwise the unconditioned subset budget
is used. This avoids assuming that the selected piece belongs to a
subset at every knot.

### Exact downward rounding

Every summand in the weighted scalar certificate is a nonnegative
rational number. To keep denominators small, the checker replaces each
summand `x` by

\[
\frac{\lfloor10^9x\rfloor}{10^9}\le x.
\]

This is exact integer division, not floating-point rounding. The sum
remains a rigorous lower bound. All 100 intervals pass (5). The
smallest recorded rounded bound is

\[
\frac{645242874434}{10^9}=645.242874434
\]

in interval `j=40`. Empty integer-rank intervals, if encountered, are
identified explicitly; no interval is dropped heuristically.

Finally,

\[
16131^{2378}>97^{3363}\,25^{2378},
\qquad3363^2-2\cdot2378^2=1.
\]

Thus `16131/25>97^sqrt(2)`. Together with the already established
all-field cases through 96, this proves the claimed cutoff through 97.
The sharper Pell fraction is necessary because the final margin is
small; the older upper approximation `577/408` is not used here.

## 4. Reproduction

```sh
python tools/check_compatible_compression_budget.py \
  --cores evidence/compatible-compression-cores.json \
  --output evidence/compatible-compression-budget.json
```

The checker imports the two preceding interior-Jacobian and rational-cut
checkers. It uses the Python standard library and takes about fourteen
seconds in the recorded workspace. Its proof inputs are the literal
integer core matrices and selected column minors. It recomputes every
determinant, checks all cross-block assemblies, rebuilds the complete
all-field complement table, covers every rank-ratio class and refined
interval, and proves the final comparison with integers.

The exploratory sparse searches and floating-point screens are not
needed for replay. No all-field result is inferred solely from a
finite-field rank experiment, and no improved factoring score is
claimed.
