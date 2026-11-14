# Four-piece cross windows close q=29 through q=39

## Result and scope

For each odd `q` in `{29,31,33,35,37,39}`, no exact decomposition of the
full-conservation crossing `T_I=E_q` costs at most `q^sqrt(2)`.

The piece floor `sqrt(r_i c_i)>=q` leaves at most four pieces at these
labels. One and two pieces are excluded by `q^2` and the all-field
[two-piece bound](two-piece-rank-sum-obstruction.md); three pieces are
excluded by the [three-piece cross relaxation](full-conservation-two-piece-reduction.md#4-three-piece-cross-windows-close-q232527);
and `5q>q^sqrt(2)` excludes at least five pieces. The new content is the
four-piece case.

This is a finite cross-constraint closure. It does not cover even labels,
`q=22`, or `q>=41`, and it does not change the complete factoring score
`1/2`.

## 1. Four-piece reduction

Let `N=q^2`, `T=q^sqrt(2)`, and suppose
`E_q=A_1+A_2+A_3+A_4` has total cost at most `T`. The piece floor gives
`cost(A_i)>=q`, so every piece costs at most `T-3q`. Put

```math
m_i=\min(r_i,c_i),\qquad L_i=\max(r_i,c_i),\qquad
K=\lfloor T-3q\rfloor .
```

Then `m_i<=K`.

A rank-one cut piece cannot occur. If `m_i=1`, the floor gives
`L_i>=q^2=N`; the quotient bound applied to the small cut gives
`L_i<=N-f(1)<N`, a contradiction. Thus `m_i>=2`.

Write `m_i` as the row rank for the row-small pieces and as the column
rank for the column-small pieces. Apply the quotient bound
`r(E_q-B)>=f(c(B))` to the piece itself, with the cuts exchanged when the
small rank is a column rank:

```math
L_i\le N-f(m_i),\qquad f(m)=\max_{0\le s\le q}s(q-ms).
\tag{1}
```

The minimum-cut budget is

```math
m_1+m_2+m_3+m_4\ge q.
\tag{2}
```

If the row-small group is `R`, its small rows sum to `M_R=sum_{i in R}m_i`;
the other pieces have large rows summing to `N-M_R`. Thus the large
coordinates in the row-small group have total `N-M_R`, and symmetrically
the column-small group has large-coordinate total `N-M_C`. For fixed small
ranks and orientations, the cost

```math
\sum_i\sqrt{m_iL_i}
```

is a sum of concave functions of the large coordinates. Each group is a
box with one sum constraint, so its minimum is at a vertex: all but at
most one large coordinate are at a bound. The checker enumerates every
multiset of small ranks `2<=m_i<=K`, every orientation partition up to
transpose, and every bound vertex.

## 2. Checked minima

The smallest relaxed costs are:

```text
 q   target      four-piece minimum   small ranks        margin
29   116.988475       137.214539       (2,2,3,22)        20.2261
31   128.559429       148.768804       (2,2,3,24)        20.2094
33   140.443941       160.317723       (2,2,3,26)        19.8738
35   152.630700       172.350065       (2,2,3,28)        19.7194
37   165.109450       184.367341       (2,2,3,30)        19.2579
39   177.871046       196.375010       (2,2,3,32)        18.5040
```

Every minimum is above its target by more than `18`. The extremal
small-rank pattern is `(2,2,3,q-7)`, with the two `2` pieces row-small
and the `3` and `q-7` pieces column-small. The checker also verifies the
one-, two-, and three-piece exclusions for every listed label.

## 3. What is not claimed

- The four-piece relaxation is a necessary-condition lower bound. It
  does not construct a tensor or prove a rank lower bound for arbitrary
  non-conservation pieces.
- Odd `q=41` and larger are not checked. The same four-piece model gives
  `208.83` at `q=41`, above the target `190.91`, but that case is not in
  the checked receipt and is not claimed here.
- Even labels and `q=22` are outside the odd conservation floor.
- The complete bound remains `1/2`.

## 4. Verification

```sh
python3 tools/check_full_conservation_four_piece_cross.py
make verify-full-conservation-four-piece-cross
python3 tools/check_docs.py
```

The checker records every small-rank multiset, orientation partition,
bound vertex, the three-piece comparison, and the exact target
comparisons. It uses ordinary floating-point square roots with margins
above `18`.
