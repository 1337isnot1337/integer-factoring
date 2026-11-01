# Ordered subgroup-return crossings

An auxiliary construction need not have fixed matrix dimension to be
effective. A group-word filter can preserve the exact crossing endpoint
and a rank-one starting cut at every length, while its scalar return
functional has infinite Hankel rank. Word reduction supplies the algorithm
without materializing an infinite auxiliary space.

This note constructs that family, gives a four-state permutation analogue,
and checks their complete binary coefficient spans at lengths one through
three. Neither finite test gives a half-rate improvement. A useful uniform
rank upper bound remains missing, and the factoring score stays `1/2`.

## 1. The subgroup-return construction

Fix a finite physical alphabet `Q`, a group `G`, elements `g_x` for
`x in Q`, and a subgroup `H` containing every commutator

```math
[g_x,g_y]=g_xg_yg_x^{-1}g_y^{-1}.
```

Assume a uniform algorithm decides membership in `H` for words of length
`L` in these fixed generators in time `2^O(L)`. Let `T(t)` be a fixed
scalar polynomial tensor over a fixed finite field, with

```math
T(0;x,y,u,v)=[x=u][y=v],\qquad
T(1;x,y,u,v)=[x=v][y=u].                                \tag{1}
```

Retain the site order and define

```math
W=\prod_{i=1}^m g_{x_i}g_{y_i}g_{u_i}^{-1}g_{v_i}^{-1},
\qquad
J_m(t;x,y,u,v)=[W\in H]\prod_{i=1}^m T(t;x_i,y_i,u_i,v_i).
                                                               \tag{2}
```

At `t=1`, every nonzero scalar entry has `u_i=y_i,v_i=x_i`.
Each local group word is then the identity. Consequently

```math
J_m(1;x,y,u,v)=[x=v][y=u].                               \tag{3}
```

At `t=0`, the support instead has `u_i=x_i,v_i=y_i`.
Every local word is a commutator in `H`, and their ordered product is
still in `H`. Thus

```math
J_m(0;x,y,u,v)=[x=u][y=v].                               \tag{4}
```

In the second cut `(x,u)|(y,v)`, (4) is the outer product of the two
diagonal indicators and has rank one for every length. The first cut
has rank `|Q|^(2m)`. No commutativity of the auxiliary generators is used.

More generally, replace the indicator in (2) by any effective function
`ell:G->F` whose value is one on `H`. Equations (3) and (4) still hold.
In a permutation action where each commutator fixes a distinguished
point `o`, any matrix coefficient `f^T P(W)e_o` with `f_o=1` has this
property. The left boundary can vary without changing either endpoint.

## 2. Complete cost, including the unbounded auxiliary

If the degree of `T` is `d`, (2) has degree at most `dm` and at most
`dm+1` coefficient pieces. A complete entry polynomial is computed by
constructing its word, testing subgroup membership, and multiplying `m`
fixed-degree scalar polynomials. This takes `2^O(m)` bit operations under
the stated membership hypothesis, and polynomial time in the examples
below.

There are `|Q|^(4m)=2^O(m)` complete entries. Writing every coefficient
tensor and Gaussian-eliminating both matrices of order `|Q|^(2m)` therefore
still takes `2^O(m)` bit operations. The resulting rank factors are
explicit. This argument charges the entire construction; it does not
assume that an auxiliary space can be enumerated for free.

For binary label width `lambda`, choose
`m=ceil(lambda/log2(|Q|))` and retain any explicitly indexed `2^lambda`
physical words on all four legs. Crossing identities survive restriction
and ranks cannot increase. Thus this is a uniform, polynomial-piece,
`2^O(lambda)` constructor. A sufficiently strong uniform bound on the
actual coefficient ranks would fit the existing
[scale-one finite-field transfer](scale-one-modular-transfer.md).
The construction alone supplies no such bound.

The connection with matrix auxiliaries is exact. Let `G` act by left
multiplication on its right cosets `G/H`. With vacuum `e_H`,

```math
e_H^T P(W)e_H=[W\in H].                                 \tag{5}
```

Every local cold commutator fixes that vacuum. For finite index this is a
finite permutation auxiliary. Infinite index is permitted because the
entry algorithm above does not enumerate the cosets.

## 3. A free-group example beyond fixed matrix dimension

Use the [three-letter Steiner cubic](steiner-cubic-crossings.md) over
`F2`, and take

```math
G=\langle a,b\rangle\text{ free},\qquad
(g_0,g_1,g_2)=(1,a,b),\qquad
c=aba^{-1}b^{-1},\quad H=\langle c\rangle.                \tag{6}
```

The nine pairwise commutators are `1,c,c^-1`, so the hypothesis of
section 1 holds. The commutator is cyclically reduced. Freely reduce
the input word with a stack; it is in `H` exactly when the result is
`c^k` or `(c^-1)^k` for some `k>=0`. This takes linear time in the word
length. A second recognizer is the folded four-cycle with successive
edge labels `a,b,a^-1,b^-1` and inverse reverse edges.

The subgroup is not normal: `aca^-1` is not a power of `c`. Its index
is infinite because the cosets `a^j H` are distinct. Indeed, a nonzero
power of `a` has nonzero `a` exponent sum, whereas all powers of `c`
have zero exponent sums. More precisely, the return functional has
arbitrarily large identity submatrices

```math
\bigl[ a^{-i}a^j\in H\bigr]_{0\le i,j<k}=I_k.           \tag{7}
```

Thus its Hankel rank is infinite over every field. It cannot be a
matrix coefficient of a fixed finite-dimensional representation of the
group. This is a statement about the group return functional; it does
not assert a lower bound on every possible representation of the
specialized crossing tensors.

The filter also breaks simultaneous site symmetry explicitly. Consider
the two local four-tuples

```text
(x,y,u,v) = (0,1,1,2) and (1,2,0,1).
```

Both scalar entries are `t^2+t^3`. Their group words are `b^-1` and
`aba^-1`. In the displayed order their product is not in `H`; in the
reverse order it is `c`. The complete two-position entry therefore
changes from zero to `t^4+t^6` when the positions are swapped.

## 4. A four-state noncommuting dilation

Let `a,b` act on `{0,1,2,3}` by the image tables

```text
a = (3,0,1,2),   b = (3,2,0,1).
```

They generate `S4`. Their commutator is `(0,2,3,1)`, which fixes zero
and cycles the other three states. With `H` the stabilizer of zero,
all nine required commutators lie in `H`, so the four-state vacuum
`e_0` gives another instance of (2). The generators themselves do not
fix zero and do not commute. The auxiliary effect is not identically
one. The same position-symmetry counterexample from section 3 applies.

This illustrates a difference from a scalar Heisenberg phase: local
cold words can act nontrivially away from the vacuum. They may alter
later transitions once a noncold factor has moved the auxiliary state.
Whether that can yield useful uniform coefficient-rank savings is open.

## 5. Complete finite coefficient-span tests

For each of the two stated filters and each `m=1,2,3`, write
`J_m(t)=sum_{j=0}^{3m}t^j A_j`. Enumerate all binary coefficient masks,
including dependencies, and compute both cut ranks of `A(w)=sum_j w_j A_j`.
The exact minimum possible largest rank product in an arbitrary
coefficient-span decomposition of the crossing is:

| Filter | Length 1 | Length 2 | Length 3 |
| --- | --- | --- | --- |
| Free commutator subgroup | 24 | 576 | 102,960 |
| Four-state return | 24 | 576 | 188,136 |

For both filters the lower detector is the coefficient coordinate `w_3`
at lengths one and two, and `w_5` at length three. Every detected mask
has at least the displayed rank product. In particular, a coefficient
relation with tensor value zero has detector zero; hence the detector
is well-defined on the actual tensor span. Its value on the crossing
is one, so every crossing decomposition contains a detected term.
This lower bound allows any number of pieces.

The evidence also stores a list of masks XORing to the endpoint, each
with rank product at most the table entry. This proves sharpness for
the stated finite minmax problem. It does not optimize total rank sum.

All six values are at least `24^m`. The exact integer comparisons
`24^7>3^20` and `2*5^2>7^2` imply normalized finite parameter
`rho>5/7>1/sqrt(2)`. They therefore fail even this permissive finite
half target. The length-three failures also show why the apparently
consistent first two normalized values cannot be treated as a uniform
rank formula.

The [constructor](../tools/subgroup_return_crossings.py),
[checker](../tools/check_subgroup_return_crossings.py), and
[evidence](../evidence/subgroup-return-crossings.json) retain all six
coefficient spans, lower detectors, upper mask lists, and source hashes.
The replay uses the opposite pivot order for **4,672 cut ranks** and
independently checks reduction and subgroup membership on **87,381 words**
of length at most eight. It also verifies the literal four-state group,
both cold cuts, all finite endpoints, and the position-symmetry witness.

```sh
make verify-subgroup-return
python3 tools/check_subgroup_return_crossings.py --generate
make verify
```

These checks are over `F2`; extension-field recombinations, other groups,
subgroups, boundary functionals, physical restrictions, and longer words
are outside the finite catalogue. The endpoint and cost proofs apply at
every length. A favorable uniform rank factorization remains the open
step. The construction and analysis used substantial OpenAI assistant
support and retain the Steiner seed's source attribution. No new Lean
work or stronger factoring theorem is claimed.
