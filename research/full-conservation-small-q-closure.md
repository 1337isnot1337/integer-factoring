# Small odd conservation kernels cannot cross

## Result and scope

Combine two existing results for the full-conservation kernel

```math
T_K(x,y,u,v)=[x+y=u+v]\,K(x+y,x+u)
```

on `q` labels.  Worker 4's Discussion #51 checkpoint proves that every
nonzero piece has

```math
\sqrt{r(T_K)c(T_K)}\ge q                                      \tag{1}
```

over every coefficient field.  The all-field
[two-piece rank-sum theorem](two-piece-rank-sum-obstruction.md) proves
that any two-piece decomposition of the label crossing `E_d` costs at
least

```math
\sqrt{d^2-d+1}+\sqrt{(d^2-1)(d-1)}.                           \tag{2}
```

The full-conservation crossing is `T_I=E_q`, since

```math
T_I=[x+y=u+v][x+y=x+u]=[x+y=u+v][y=u]=[x=v][y=u].
```

At `q=5`, the floor (1) gives `2q>q^sqrt(2)`, and one piece costs
`q^2>q^sqrt(2)`.  At `q=7,11,13`, the floor gives
`2q<q^sqrt(2)<3q`: a crossing has at most two pieces, and (2) is already
larger than the target.  The numerical comparisons are

```text
 q   target q^sqrt2   floor cap   two-piece bound
 5       9.738518         1          14.380535
 7      15.672891         2          23.528001
11      29.699672         2          45.176670
13      37.614370         2          57.429853
```

Thus there is **no full-conservation-kernel crossing for q=5,7,11,13
over any field**.  At `q=17` the floor permits two or three pieces and
the two-piece bound no longer closes the lane by itself; the checker
labels that row ``not closed by this composition.''  The minimum-cut
cutoff below closes it (and `q=19`) independently.

This is a scope closure for the fixed kernel space, not a lower bound
on arbitrary crossing tensors.  In particular it does **not** close the
general non-conservation translation-invariant lane, and it does not
change the complete factoring score `1/2`.

## Relation to the all-decomposition cutoff

There is a stronger existing statement.  The
[minimum-cut rank budget](minimum-cut-rank-budget.md) proves

```math
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}
```

for **every** exact decomposition of `E_d`, with no conservation or
piece-count restriction, for all `2<=d<=21`.  Since `T_I=E_q`, this
already closes every conservation seed with `q<=21`, including
`q=5,7,11,13,17,19`, without the floor argument.  The composition above
is an independent, conservation-specific check of the same boundary,
not a new closure of q=13.  If a q=13 search has a different label
dimension or target, that scope should be stated explicitly; for the
full-conservation kernel it is already excluded.

## Relation to the parallel work

The floor (1) is worker 4's finite support lemma, posted at
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19505782).
The two-piece bound (2) is the repository's ordinary all-field theorem.
The composition is a conservation-specific check of the q=13
boundary, including dense kernels inside that kernel family: the floor
restricts such a crossing to two pieces before the two-piece theorem
applies.  The stronger minimum-cut theorem already implies the same
exclusion without using conservation.

## Reproduction

```sh
python3 tools/check_full_conservation_small_q_closure.py
make verify-full-conservation-small-q-closure
```

The checker evaluates (1), the target, and (2) for `q=5,7,11,13,17`
with ordinary floating-point margins far from zero; it labels `q=17` as
not closed by this composition, not as globally open.  It does not build
a tensor or a factor algorithm.
