# Minimum-channel jet factorizations with polylogarithmic workspace

## The algorithmic result

Let `R=F_q[t]/(t^H)`, where `q` is prime and `H>=1`. Suppose an
`m` by `n` matrix `A` is available through a coefficient oracle
returning the coefficient of `t^a` in any entry. With
`N=max(2,m,n)` and `L=HN`, one deterministic algorithm can:

- compute the minimum number of `R`-linear factor channels;
- query coefficients of an exact factorization `A=UV` with that
  minimum number of channels;
- test `Ax=b` and query coefficients of a canonical solution;
- recover the multiplicities of every nonzero Smith valuation.

The time and auxiliary workspace for one requested query are

```math
T=L^{O(\log L)}\operatorname{poly}(\log q)(1+T_{\rm entry}),
\qquad
S=S_{\rm entry}+O\bigl(\log L\,(\log L+\log q)\bigr).       (1)
```

Right-hand-side and other supplied coefficient-oracle costs are
included in the same manner. Factors and solutions are queried or
streamed. Retaining their full coefficient arrays would require
additional storage.

The algorithm combines the
[field rank and affine-query construction](polylog-space-rank-recovery.md)
with coefficient linearization and a generator selection in
`image(A)/t image(A)`. It never divides by a nonunit of `R` and
never retains a matrix-sized list of selected columns.

This applies to equal-characteristic truncated polynomial rings.
It does not treat the base-`p` digits of `Z/p^H Z` as a vector
space over `F_p`: addition in that ring has carries. It also does
not by itself improve the factoring exponent. In particular, if
`L=2^Theta(lambda)`, the time in (1) is generally
`2^O(lambda^2)`, not the `2^O(lambda)` entry budget in the simplest
crossing transfer. The matrix order, precision, and chosen time-space
tradeoff must all be charged by the caller.

## 1. Coefficient linearization

For `1<=s<=H`, reduce the entries modulo `t^s` and form an
`ms` by `ns` matrix over `F_q`:

```math
\mathcal T_s(A)_{(i,a),(j,b)}
 =\begin{cases}
   [t^{a-b}]A_{ij},&a\ge b,\\
   0,&a<b,
  \end{cases}
\quad 0\le a,b<s.                                      (2)
```

It is the matrix of multiplication by `A` on coefficient vectors.
A requested entry is zero or one call to the original oracle.
The indices and subtraction use `O(log L)` bits; no block matrix
is materialized.

Put

```math
R_s=\operatorname{rank}_{F_q}\mathcal T_s(A),\qquad R_0=0.
```

The field algorithm evaluates each `R_s` within (1). Field
extensions used inside that algorithm do not change the coefficient
ring or its solutions: they are only deterministic evaluation aids
for detecting the rank of the base-field matrix.

## 2. Smith valuations and the exact channel count

Elementary row and column operations over `R` put `A` in Smith
form with nonzero diagonal entries

```math
t^{e_1},\ldots,t^{e_r},\qquad 0\le e_1\le\cdots\le e_r<H.
```

Here is the needed existence argument. Select a nonzero entry with
minimum `t`-valuation, move it to the leading position, and remove
its unit factor. Every entry is divisible by its power of `t`, so
elementary operations clear its row and column. Repeat on the
remaining submatrix. These operations are invertible; they can be
used in this proof without being stored or executed by our algorithm.
The zero matrix has `r=0`.

After reduction modulo `t^s`, the diagonal map `t^e` has field
rank `max(s-e,0)`. Therefore

```math
R_s=\sum_{j=1}^r\max(s-e_j,0),\qquad
R_s-R_{s-1}=\#\{j:e_j<s\}.                            (3)
```

In particular,

```math
\boxed{r=R_H-R_{H-1}.}                                (4)
```

The number of valuations equal to `e` is

```math
\begin{cases}
R_1,&e=0,\\
(R_{e+1}-R_e)-(R_e-R_{e-1}),&1\le e<H.
\end{cases}                                           (5)
```

These counts can be queried or streamed using a constant number of
rank values and counters.

The integer `r` is exactly the minimum number of inner factor
channels, meaning the least `c` for which `A=UV` with `U` of size
`m` by `c` and `V` of size `c` by `n`. Smith form gives the upper
bound `c<=r`. For the lower bound, restrict a putative factorization
to the nonzero `r` by `r` Smith diagonal. If `c<r`, Smith reduction
of its right factor gives a kernel vector having a unit coordinate:
an unused free column of its Smith form supplies such a vector,
and an invertible change of coordinates preserves this property.
Multiplication by the nonzero diagonal `t^(e_j)` cannot annihilate
a vector with a unit coordinate. This contradicts the factorization.

Thus (4) computes the minimum channel count, not merely the number
of unit pivots. For example, the one-entry matrix `(t)` has one
channel over `F_q[t]/(t^2)`, although its reduction modulo `t` is zero.

## 3. Select generators modulo the full image's radical

Knowing the channel count does not justify choosing columns using
the channel counts of successive prefixes. For example, the columns
of `(t,1)` over `F_q[t]/(t^2)` have the same nonzero channel count
after either prefix, but the first column does not generate the full
image. The selection must use the quotient of the *full* image.

Let `M=image_R A` and let `a_j` be its original column generators.
The field coefficient space of `tM` is the image of
`mathcal T_H(tA)`, whose entries are obtained by shifting one
coefficient index in (2). Write

```math
Q=\mathcal T_H(tA),\qquad
C_j=(\operatorname{coeff}(a_1),\ldots,
                   \operatorname{coeff}(a_j)).
```

Each column of `C_j` is an `mH`-dimensional field vector, not all
of its shifts. Define

```math
s_j=\operatorname{rank}_{F_q}[Q\mid C_j]
                      -\operatorname{rank}_{F_q}Q,
\qquad s_0=0.                                          (6)
```

These are the ranks of the first `j` generator classes in the
field vector space `M/tM`. All classes are generated by the
original columns with constant coefficients: the higher powers
of `t` vanish in this quotient. Consequently `s_n=r`.
This also follows from Smith form, where each nonzero summand
`t^e R` contributes exactly one dimension modulo its own multiple
by `t`, including when `e=H-1`.

Choose the columns whose addition increases `s_j`. The index of
the `l`th chosen column is simply

```math
J(l)=\min\{j:s_j\ge l\}.                              (7)
```

All entries of the field matrices in (6) have direct coefficient
oracles for `A`. A linear scan or binary search computes (7) with
counters and recomputed ranks. No pivot mask is retained.

Let `U` consist of the selected original columns. Their classes
form a basis of `M/tM`, so

```math
M=\operatorname{image}_R U+tM.
```

Iterate this equality `H` times. Because `image U` is an
`R`-submodule and `t^H M=0`, it follows that

```math
M=\operatorname{image}_R U.                            (8)
```

This is the finite nilpotent version of the generator argument
usually called Nakayama's lemma. The iteration proves exactly
the case used here.

## 4. Query a minimum-channel factorization

A coefficient of `U_(i,l)` is obtained by computing `J(l)` from
(7) and querying `A_(i,J(l))`. For each original column `a_j`,
solve

```math
U v_j=a_j\quad\text{over }R.                           (9)
```

It is consistent by (8). After coefficient linearization, (9)
is an ordinary field system

```math
\mathcal T_H(U)\,\operatorname{coeff}(v_j)
                         =\operatorname{coeff}(a_j).  (10)
```

Use the canonical field solution that sets the nonpivot coefficient
coordinates to zero. The field affine-query algorithm returns any
requested coefficient of `v_j` without storing the whole vector.
Assemble these vectors as the columns of `V`. Equation (9) for
every column proves `A=UV`. There are exactly `r` channels, which
is minimal by Section 2.

This construction makes no assertion that a chosen solution map
on all possible right-hand sides is `R`-linear. It does not need
that property: the displayed matrix product is checked separately
on each original column, and the resulting entries of `V` belong
to `R`.

## 5. Affine solutions and the cost of nested queries

For a general right-hand side `b`, test membership of its coefficient
vector in the image of `mathcal T_H(A)`. Equivalently compare
the field ranks before and after adjoining that single vector.
When it is consistent, use the same field solver to query any
coefficient of a solution. Comparing only the ring channel counts
of `A` and `[A|b]` would be wrong: `(t)x=1` modulo `t^2` is
inconsistent even though both ring channel counts are one.

The field matrices in (2), (6), and (10) have dimensions `O(L)`.
A requested factor coefficient in (10) invokes the field solver.
Its virtual matrix entries invoke the generator selector (7).
That selector invokes ranks of (6), whose entries query the
original `A` directly. These dependencies have a fixed depth;
there is no recursive call to the same selector through itself.

Multiplying a fixed number of quasipolynomial time bounds and
adding their scalar recursion stacks preserves (1). The retained
items are a constant number of rank values, grade and matrix
indices, field elements, and the logarithmic recursion stacks.
No full coefficient vector, Smith basis, generator mask, or
factor matrix is retained.

The proof also works over any explicitly represented finite field
whose field arithmetic has polynomial bit time and linear space
in its element bit length. A tower extension in the rank detector
can be represented by polynomial coefficients over that field.

## Verification and use in the factoring project

The exact checker compares (3)--(5) with independent small Smith
reductions, constructs the generator quotient (6), and verifies
the resulting factor products and affine systems. Negative controls
include `(t,1)` and `(t)x=1`, so the field-rank test cannot be
silently replaced by a ring-channel-count test.

```sh
python3 tools/check_truncated_ring_factorization.py
```

The committed evidence in `evidence/truncated-ring-factorization.json`
covers 380 matrices, including all 256 two-by-two matrices over
`F_2[t]/(t^2)`. Independent enumeration of every one-channel product
checks minimality on that complete set. Bounded oracle examples also
execute the nested coefficient queries without retaining a pivot list.
The finite execution verifies the identities and negative controls;
the workspace bound is the scalar-oracle analysis above.

This provides a concrete algorithmic tool for the finite-precision
identities in [border seed extraction](border-and-nonunit-seed-transfer.md).
A smaller channel profile remains a separate mathematical input.
The procedure does not turn integer digit extraction into a linear
operation, and it does not establish a stronger factoring score.
