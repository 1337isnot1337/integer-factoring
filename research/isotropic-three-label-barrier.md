# A constant-composition obstruction after an isotropic basis change

## Result and scope

Start with the twenty-piece three-label crossing certificate from
[PR #24](https://github.com/1337isnot1337/integer-factoring/pull/24),
pinned at `e96ee50560284301428d72421330d8773155d25b`.
The matrices are stored in
[the finite-seed checker](../tools/check_finite_seed_amplification.py).
The label space is over `F_2`, with its standard pairing.

Change the ordered label basis to

```math
(e_0,e_2,e_0+e_1).
\tag{1}
```

The last vector is isotropic: its pairing with itself is zero.
Use the dual change of basis on the two matched external legs.
The crossing identity is unchanged. This construction is outside
the original-basis hypothesis of
[the earlier constant-weight obstruction](balanced-three-label-barrier.md).

Take `t=3m` tensor copies and retain words with exactly `m`
occurrences of each of the three new symbols. There are

```math
D_m=\frac{(3m)!}{(m!)^3}
\tag{2}
```

labels. Group tensor copies by their literal twenty piece counts.
The single type `16^m 17^m 19^m` has exact ranks

```math
r_m=D_m\sum_{k=0}^m\binom mk^4,
\qquad
c_m=D_m\binom{2m}{m}.
\tag{3}
```

Consequently this family requires

```math
\rho\ \ge\ \frac{\log 6}{2\log 3}=0.815464\ldots
\tag{4}
```

asymptotically. The same obstruction applies if only the count of
the last symbol is fixed to `m`: the domain then has
`binom(3m,m)2^(2m)=3^(3m)/poly(m)` labels, and contains the domain
in (2). Restriction cannot increase either rank.

The theorem concerns this basis and this density. It does not claim
an obstruction at every density, for arbitrary linear label codes,
or after arbitrary recombinations of the type tensors.

## 1. Two exact charge constraints

Let `S` have the three vectors in (1) as its columns. Its rows,
encoded in binary, are `(5,4,2)`. For a piece matrix with row index
`(x,y)` and column index `(u,v)`, the transformed matrix is

```math
A'_i=(S^{\mathsf T}\otimes S^{\mathsf T})
       A_i
       ((S^{-1})^{\mathsf T}\otimes(S^{-1})^{\mathsf T}).
\tag{5}
```

This is the same matched-leg change for all pieces, so their sum
remains `E_3(x,y,u,v)=[x=v][y=u]`.

For every nonzero entry of every transformed piece, the charge

```math
\Delta_2=[x=2]+[y=2]-[u=2]-[v=2]
\tag{6}
```

is nonpositive. In a tensor entry whose four external words have
the same count of symbol `2`, the sum of these charges is zero.
Every contributing local entry must therefore have `Delta_2=0`.

On those surviving entries, the second charge

```math
\Delta_1=[x=1]+[y=1]-[u=1]-[v=1]
\tag{7}
```

is also nonpositive. Equal counts of symbol `1` force it to vanish
at every position. These are exact restrictions on individual
terms; there is no limit, generic specialization, or cancellation
assumption.

After both restrictions, the only nonzero pieces have these supports:

| Piece | Allowed `(x,y,u,v)` |
| --- | --- |
| `12` | `(1,1,1,1)` |
| `15` | `x,y,u,v` in `{0,2}`, with `x!=y` and `u!=v` |
| `16` | `(a,b,a,b)`, for `a,b` in `{0,2}` |
| `17` | `(a,1,1,a)`, for `a` in `{0,2}` |
| `19` | `(1,a,a,1)`, for `a` in `{0,2}` |

Their sum is still the three-label crossing. The checker verifies
both charge inequalities directly on all `107` original nonzero
entries. There are `29` entries after the first constraint and
`13` after the second.

## 2. Exact rank in the row cut

For the type in (3), let `A,B,C` be the sets of positions carrying
pieces `17,19,16`, respectively. Each set has size `m`; together
they partition the `3m` positions. In the four external words,
the positions occupied by symbol `1` are

```math
x^{-1}(1)=v^{-1}(1)=B,
\qquad y^{-1}(1)=u^{-1}(1)=A.
\tag{8}
```

Thus a nonzero row `(x,y)` determines the partition `A,B,C`.
The local supports then determine `(u,v)` uniquely. The converse
also holds, so this flattening is a partial permutation matrix.

On `A`, write the common `x,v` word as `a`; on `B`, write the
common `y,u` word as `b`. On `C`, let the common `x,u` word be
`c`, and the common `y,v` word be `d`. These are binary words
on the alphabet `{0,2}`. The composition constraints say that,
for some `k`, the words `c,d` each contain `k` twos, while `a,b`
each contain `m-k` twos.

For a fixed partition and `k`, there are `binom(m,k)^4` nonzero
rows. There are `D_m` partitions. Counting nonzero rows proves
the formula for `r_m` in (3).

## 3. Exact rank in the column cut

Now use rows `(x,u)` and columns `(y,v)`. Equation (8) again
separates the partitions `A,B,C` into disjoint blocks.
Within one partition, fix `k` and the words `a,b`. The possible
rows are indexed by the `binom(m,k)` choices of `c`; the possible
columns are indexed by the same number of choices of `d`.
Every entry of this submatrix is one. Its rank is one over `F_2`.

Distinct choices of `k,a,b` use disjoint rows and columns.
The rank for one partition is therefore

```math
\sum_{k=0}^m\binom mk^2=\binom{2m}{m}.
```

Multiplying by `D_m` proves the formula for `c_m`.

## 4. The asymptotic consequence

The multinomial estimate and elementary binomial bounds give

```math
\log D_m=3m\log 3-O(\log(m+2)),
\quad
\log\sum_k\binom mk^4=m\log16+O(\log(m+2)),
\quad
\log\binom{2m}{m}=m\log4-O(\log(m+2)).
```

For example, the largest `binom(m,k)` is at least `2^m/(m+1)`;
these estimates need no unproved rank asymptotic. Hence

```math
\log(r_mc_m)=6m\log6+O(\log(m+2)).
```

A joint width bound has the form
`log2(r_m c_m)<=4 rho log2(D_m)+o(log D_m)`.
Substitution gives (4). It is strictly larger than `1/sqrt(2)`:
the exact integer comparisons `6^5>3^8` and `2*4^2>5^2`
show that it exceeds `4/5`, which exceeds `1/sqrt(2)`.

## 5. Verification

Run:

    python3 tools/check_isotropic_three_label_barrier.py \
      --output evidence/isotropic-three-label-barrier.json

The checker verifies the crossing identity after the dual basis
change and every local charge and support assertion. It constructs
the type matrices directly from the original transformed pieces,
without using the rank formulas, at `m=1,2`.

| `m` | Labels | Row rank | Column rank | Rank product |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 6 | 12 | 12 | 144 |
| 2 | 90 | 1620 | 540 | 874800 |

It also checks `256` exact instances of the binomial estimates and
the two integer comparisons above. The all-size result follows
from the charge constraints and the two explicit block decompositions.
