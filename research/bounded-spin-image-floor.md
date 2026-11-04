# A rank floor for correction-module decompositions of the bounded-spin crossing

## Result and scope

The [correction-module lemma](bounded-spin-correction-modules.md) proposes
using the low-rank extreme tensor `T` to approximate the coherent
crossing. This note isolates an image obstruction for that proposal.
It does not close the full coefficient algebra, the central-idempotent
decomposition, or a changed local kernel.

Let `V` be a finite-dimensional vector space over a field, let `B` be a
unital subalgebra of `End(V)`, and let `T in End(V)` have rank at most
`a`. Put

```math
M=B(\operatorname{im}T),\qquad A=\operatorname{Alg}(B,T).
```

Write `N=dim V` and `m_k=dim M` in cut `k`. Then:

1. Every word in `A` containing `T` has image in `M`; words without `T`
   lie in `B`.
2. Consequently every `X in A` has a decomposition `X=Y+R` with
   `Y in B` and `im R subset M`.
3. If `X` is invertible, then
   ```math
   \operatorname{rank}_k(Y)\ge N-m_k.                 \tag{1}
   ```
   For a two-piece correction decomposition `X=Y+R`,
   ```math
   S=\sqrt{r(Y)c(Y)}+\sqrt{r(R)c(R)}
    \ge\sqrt{(N-m_1)(N-m_2)}.                          \tag{2}
   ```

The floor is exact and needs no tensor entry computation. It applies to
the single-`B`-part decomposition supplied by the correction-module
lemma, not to a sum of several `B`-pieces. The latter is the remaining
central-idempotent target.

## 1. The word-image argument

Expand `X` in words in the generators of `B` and `T`. A word with no
`T` is in `B`. In a word containing `T`, look at the last occurrence of
`T`. Its image is contained in `im T`; every later factor belongs to
`B`, so the image of the whole word is contained in `B(im T)=M`.
Linearity proves (1) for `R` and `Y`.

If `X` is invertible, rank subadditivity gives

```math
N=\operatorname{rank}(X)
 \le\operatorname{rank}(Y)+\operatorname{rank}(R)
 \le\operatorname{rank}(Y)+m_k.
```

This is (1). Taking the product of the two cut inequalities and using
`sqrt{r(Y)c(Y)} + sqrt{r(R)c(R)} >= sqrt{r(Y)c(Y)}` proves (2).

The obstruction is a statement about the **image**, not about the size
of the coefficient space. Adding more coefficient directions inside
`M` cannot remove it.

## 2. The published finite modules

The bounded-spin notes supply two exact modules:

| Length | `N` | `m` | target `d^(sqrt 2)` | floor (2) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 169 | 28 | `13^sqrt(2)=37.61` | 141 |
| 10 | 1156 | 170 | `34^sqrt(2)=146.50` | 986 |

The eight-site module is the invariant commutator module of
[bounded-spin-six-vertex.md](bounded-spin-six-vertex.md), section 4.
The ten-site module is the vertex module of
[bounded-spin-correction-modules.md](bounded-spin-correction-modules.md),
section 4. In both cases the floor already exceeds the finite target
by a wide margin. These are independent of the minmax certificates in
the rank-algebra notes; they exclude the single-module correction
decomposition even if its internal coefficient search is optimized.

At twelve sites the crossing has `N=89^2=7921` and target
`89^sqrt(2)=571.29`. The extreme tensor has rank `a_6=396`. If the
correction algebra has dimension `d_B`, then `m_k<=d_B a_6`. Thus any
`d_B<=18` gives `m_k<=7128` and floor at least `793>571`. The ten-site
original algebra has dimension eight. A polynomial dimension bound for
the twelve-site algebra would therefore close this route; a growing
`d_B` above eighteen leaves the image floor inconclusive.

## 3. Asymptotic consequence

The extreme tensor at length `n=2r` has rank

```math
a_r=\frac{(2+\sqrt2)^r+(2-\sqrt2)^r}{4}
    =O\bigl((2+\sqrt2)^{n/2}\bigr),
```

so its per-site base is `sqrt(2+sqrt2)=1.848...`. The cut length is
`L_n=Theta(phi^n)` with `phi=1.618...`, hence `N=L_n^2=Theta(phi^{2n})`.
If the correction algebra has dimension `d_B=exp(o(n))`, then

```math
m_k\le d_B a_{n/2}=o(N),
```

and the floor (2) is `N-o(N)`. The finite target is
`L_n^{sqrt2}=N^{1/sqrt2}`, and

```math
\frac{N}{L_n^{sqrt2}}
=L_n^{2-\sqrt2}
=\phi^{(2-\sqrt2)n+o(n)}\longrightarrow\infty.
```

Therefore no single low-rank-`T` correction decomposition can cross for
large `n`, even if its coefficient algebra has subexponential dimension.
The eight- and ten-site failures above are not finite-size accidents;
the single-module route is asymptotically closed. This does not close a
sum of several `B`-pieces, a different local kernel, or a non-low-rank
generator.

## 4. What remains

A decomposition of `X` into several `B`-pieces can have a much smaller
rank sum than one high-rank `B`-piece, because the pieces can cancel in
both cuts. The natural candidates are the central idempotents of the
original coefficient algebra. At ten sites the row primary dimensions
are

```text
34, 1, 225, 300, 35, 406, 155,
```

summing to `N=1156`; the row and column cuts pair them through the cut
involution. The finite minmax certificate already shows that every
decomposition from that coefficient span has square-root rank sum at
least `225`, above the `146.5` target. The open question is whether the
twelve-site primary dimensions and their two-cut pairing can give a
sum below `571.29`.

The [second-coefficient identity](bounded-spin-second-coefficient.md)
identifies the eight-site primary blocks with transposition class sums,
and the three-jet identities constrain the coefficient span at every
length. A representation-theoretic computation of the twelve-site
primary dimensions and cut pairing is the next concrete target. The
image floor above rules out the shorter single-module shortcut; it does
not rule out that central decomposition.

No factoring exponent changes here. The complete bound remains `1/2`.
