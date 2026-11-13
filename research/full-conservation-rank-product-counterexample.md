# A rank-two conservation kernel with sub-square-root rank cost

## Question and result

The characteristic-two rank-sum barrier proves
`sum_i sqrt(r_i c_i) >= q^(3/2)` for every exact conservation-kernel
decomposition on an elementary abelian two-group.  A tempting shortcut
to the same conclusion for every finite group and every field is the
piecewise inequality

```math
\sqrt{r(T_K)c(T_K)}\ge q\sqrt{\operatorname{rank}K}. \tag{1}
```

If (1) held, then any exact decomposition `I=sum_i K_i` would give

```math
\sum_i\sqrt{r_i c_i}
\ge q\sum_i\sqrt{\operatorname{rank}K_i}
\ge q\sqrt{\sum_i\operatorname{rank}K_i}
\ge q^{3/2},                                        \tag{2}
```

because `sum_i rank K_i >= rank I = q`.  This would close the whole
full-conservation family at every `q`, over every field, not only at
the characteristic-two two-groups where the degree proof applies.

Inequality (1) is **false**.  On `G=Z/3`, over `F_2`, take

```math
K=\begin{pmatrix}1&1&1\\1&1&1\\1&1&0\end{pmatrix}. \tag{3}
```

Then `rank K=2`, and the two conservation cuts have ranks

```math
r(T_K)=4,\qquad c(T_K)=4.                            \tag{4}
```

Thus `sqrt(r c)=4` while `q sqrt(rank K)=3 sqrt(2)=4.2426...`.
The checker also scans all `2^9` kernels on `Z/3` over `F_2`; the
example maximizes the gap `q^2 rank K-r c`.

## What this does and does not show

A weaker piecewise inequality also fails:

```math
\sqrt{r(T_K)c(T_K)}\ge\sqrt q\,\operatorname{rank}K. \tag{5}
```

For the full-rank kernel

```math
K=\begin{pmatrix}1&1&1\\1&0&1\\1&1&0\end{pmatrix}
```

on the same group, `rank K=3` and both cut ranks are `5`, so the cost
is `5` while `sqrt(3)*3=5.196...`.  Thus one cannot recover the
`q^(3/2)` decomposition bound from either per-piece inequality.

The failure of (1) and (5) is a **negative diagnostic for proposed
shortcuts**.  It does not contradict any proved statement:

- the all-field minimum-cut cutoff still closes every exact
  decomposition with `q<=21`;
- the characteristic-two p-group rank-sum theorem still proves the
  barrier for its stated groups and fields;
- the minimum-rank piece floor and the two-piece bounds are unchanged.

In particular, (3) and the full-rank kernel are single pieces, not
exact decompositions of the identity, so they are not seeds or
counterexamples to any crossing statement.  The `q>=22` full-conservation lane over a general field
remains open; the rank-product route would need a different lower
bound that survives this example.

## Reproduction

```sh
python3 tools/check_full_conservation_rank_product_counterexample.py
make verify-full-conservation-rank-product-counterexample
```

The checker uses exact Gaussian elimination over `F_2` and the
convolution-rank formula for the two conservation cuts.
