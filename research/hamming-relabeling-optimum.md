# Exact cutoff four in the relabeled width-two span

## Result

Let `S` be the 102-dimensional `F_2` span obtained by applying every common
permutation of the four labels to the old joint one/two-pair source span.
Its construction and the earlier six-piece product-16 witness are in
[the relabeling note](hamming-relabeling.md).

The exact minimum maximum rank product for a crossing decomposition
inside `S` is **four**. A retained **20-piece** witness attains it; every
piece has ranks `(1,4)` or `(4,1)`. Its maximum rounded width sum is
**two**, also optimal. This improves the old product-25 result in the
smaller 33-dimensional span. The six-piece product-16 witness remains a
different tradeoff with fewer summands.

The witness, source membership masks, and verification receipt are in
[the evidence](../evidence/hamming-relabeling-optimum.json). The upper
bound persists under characteristic-two field extension. The optimality
claim is **over `F_2` only**: the complete lower-bound classification
enumerates binary one-coordinate factors. No extension-field lower bound
is inferred from it.

This is a finite construction result, not a new factoring exponent.
The overall complete exponent remains `1/2`. The proof used substantial
OpenAI assistant support, adds no Lean formalization, and awaits
independent review.

## 1. A rank product below four factors off one coordinate

Write the two ranks as `r_row(F),r_col(F)`. A nonzero tensor with product
at most three has rank one in at least one arrangement. Suppose the row
rank is one, so

```math
F(x,y,u,v)=A(x,y)B(u,v).
```

Its column matrix is the Kronecker product of the two four-by-four
matrices `A,B`, up to row/column ordering. Thus

```math
r_col(F)=rank(A)rank(B).
```

If this is at most three, one of `A,B` has rank one. Consequently one
of the four individual coordinates factors off:

```math
F(x,y,u,v)=f(coordinate_i)G(other\ three\ coordinates).    (1)
```

The argument starting from column rank one is identical. This reduction
uses only elementary rank facts and holds over any field.

## 2. Complete slice intersections and a four-entry separator

Over `F_2`, the nonzero vector `f` in (1) has exactly 15 possible values.
For each of the four coordinate choices and each such vector, compute
the intersection of `S` with `f tensor F_2^64`. This requires only linear
equations in the 102 source-basis coefficients.

Explicitly choose any position `p` where `f_p=1`. For each assignment of
the other three coordinates, require

```math
F_i=f_i F_p\qquad(i\ne p).                              (2)
```

These 192 binary equations are necessary and sufficient for (1).
Exact elimination therefore gives the whole intersection, not a sample
of its rank-one elements. The 60 systems have dimensions

| Intersection dimension | Number of systems |
| ---: | ---: |
| 0 | 16 |
| 2 | 40 |
| 14 | 4 |

Their 136 kernel-basis vectors span a 60-dimensional tensor subspace.
Every one is annihilated by the same functional

```math
ell(F)=F(3,3,2,2)+F(3,3,2,3)
       +F(3,3,3,2)+F(3,3,3,3).                          (3)
```

The checker verifies this on the complete bases returned by (2).
Linearity covers every tensor in all 60 intersections. On the crossing,
only the last entry of (3) is one, so `ell(E)=1`.

Every tensor in `S` of rank product below four lies in one of these
intersections and hence has `ell=0`. Such tensors cannot sum to `E`.
This proves the lower bound, with 60 small linear systems and no SAT or
external solver dependency.

## 3. Upper bound and what remains open

The checker rebuilds the 196 distinct relabeled source images and their
102-dimensional span. For every retained witness piece it replays its
source-image coefficient mask, independently checks both dense ranks,
and verifies all 256 crossing entries of the sum. All 20 products are
four. Each width sum is `ceil(log2 1)+ceil(log2 4)=2`.
A width sum below two would force a product below four, already excluded.

```sh
python3 tools/check_hamming_relabeling_optimum.py
python3 tools/check_hamming_relabeling_optimum.py --output evidence/hamming-relabeling-optimum.json
make verify
```

No minimum number of pieces is claimed. A bounded integer-programming
probe for fewer than 20 pieces did not return a feasible solution within
its 120-second limit; it supplies no lower bound on the piece count.

As with the earlier finite seeds, tensoring this witness gives
exponentially many pieces. Collecting words into polynomial coefficients
can increase their ranks. The [universal rank-count restriction](post-half-directions.md#5-a-universal-restriction-within-polynomially-many-crossing-pieces)
already prevents extrapolating this finite normalized contraction of
`1/4` to a subexponential-piece family. In fact, if `M` nonzero pieces
have products at most `Q`, then

```math
2D^2<=\sum_i(r_row(F_i)+r_col(F_i))<=M(Q+1),              (4)
```

because positive integer ranks satisfy `r+c<=rc+1`. The required
unbounded-width theorem still has to account for both ranks and piece
count. The new finite optimum supplies a seed and a precise search
boundary; it does not supply that theorem.
