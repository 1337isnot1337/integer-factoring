# Target orders, pure relation rows, and common multipliers

**Status.** Ordinary proofs of the order-truncation lemma, the pure-row
determinant, the mixed-order-row determinant, and the common-multiplier
determinant identity; a method-scoped medium-order barrier for the unmixed
order/density hybrid; and one exact finite hybrid certificate. No uniform
multiplier selection, no mixed-density positivity theorem, and no new
complete factoring bound is proved. The complete saving exponent remains
**1/2**.

This continues the [smooth selector](smooth-leading-power-selection.md), the
[small-index density](small-index-smooth-density.md), and the
[pooled rank condition](clustered-hits-and-pooled-rank.md). The question here
is different from extending the even-product density: it is whether changing
the **target** class, rather than the density theorem, moves every needed
target out of the uncovered index range. The answer is a clean conditional
reduction plus a precise algebraic obligation.

## 1. The order row

Use the exact selector contract. Let `Delta<0` be a discriminant, `D=|Delta|`,
`H=floor(sqrt(D/3))`, and

\[
 B=\{p\le Y:\ p\nmid \Delta,\ (\Delta/p)=1\}
\]

a split-prime base, with prime form `g_p` and
`K_p=min{k>=1: p^k>H}`. For a target class `h` with an explicit base word
`w in Z^B`, a smooth hit at exponent `k` is a reduced form of class `h^k`
whose leading coefficient `A<=H` decodes exactly as

\[
 A=\prod_{q\in B}q^{|v_q|},\qquad h^k=\prod_{q\in B}g_q^{v_q}.
\]

The acquired relation is the integer row `r=k w-v`.

**Lemma 1.1 (identity hit).** Let `m=ord(h)`. At `k=m` the class is the
identity, the reduced identity form has leading coefficient `1`, and the
exact decoder returns `v=0`. Hence the relation

\[
 r=m w
\]

is acquired after at most `m` powers of `h`, with no smoothness search beyond
the identity test.

**Proof.** `h^m=1`. The identity form is reduced and has leading coefficient
`1<=H`; its oriented ideal factorization is empty. The decoder is exact, so
`v=0`. The relation is the group identity `m w=0`.

For a prime target `h=g_p` this gives the **pure row** `m e_p`. The row is
available at exponent `m` even when the index `[G:<g_p>]` is far above the
density-covered range.

**Corollary 1.2 (pure-row block).** Let `h_1,...,h_n` have base words `w_i`
and orders `m_i`. If the word matrix `W=(w_i)` is invertible over `Q`, the
order rows `m_i w_i` form the full-rank relation matrix

\[
 M=\operatorname{diag}(m_i)W,\qquad
 \det M=\Big(\prod_i m_i\Big)\det W\ne0.
\]

If the coordinate generators `g_q` generate `G`, then `|G|` divides
`|\det M|`, so `det M` annihilates `G`.

**Proof.** Let `phi:Z^n->G` be the coordinate surjection `e_i|->g_{p_i}`;
the base generating `G` makes it surjective. Every row `m_i w_i` lies in
`ker phi`. The row lattice `L` is a sublattice of `ker phi`, and its covolume
is `|det M|`. Hence `|G|=[Z^n:ker phi]` divides `[Z^n:L]=|det M|`.

For the unmixed prime targets, `W=I` and the determinant is `prod m_p`. This
is the diagonal relation matrix already used as a finite control in
[the class-group checker](../tools/check_arithmetic_class_groups.py). Its
index is generally much larger than the class number, but any nonzero
multiple of `|G|` is a valid annihilator for the existing extractor.

**Corollary 1.3 (mixed order rows).** Fix `z in Z_{\ge0}^B` and put
`h_i=p_i z`, with base word `w_i=e_i+z`. Then

\[
 W=I+\mathbf 1 z^T,\qquad \det W=1+\|z\|_1>0.
\]

If the mixed orders are `m_i`, the order rows `m_i(e_i+z)` form a full-rank
relation matrix with determinant

\[
 \Big(\prod_i m_i\Big)(1+\|z\|_1).
\]

Thus a common multiplier preserves the pure-row rank. The mixed target need
not generate `G`; the rows are relations in the original coordinate basis,
and the coordinate base is the generating set.

Corollary 1.3 is the only use of mixing that is unconditional here. It gives
a **mixed-order selector**: if the mixed targets' orders are all at most `W`,
the selector factors within `W` power tests per target. It does not say how
to choose `z`.

## 2. Smooth rows under a common multiplier

Now suppose a mixed target `h_i=p_i z` acquires a smooth hit at exponent
`k_i`, with nonnegative word `v_i` and row

\[
 r_i=k_i(e_i+z)-v_i.
\]

Set `D=diag(k_i)`, `V=(v_{ij})`, `A=I+\mathbf 1 z^T`, and `K=D`. The full
row matrix is

\[
 R=DA-V=(D-V)+D\mathbf 1 z^T.
\]

**Lemma 2.1 (common-multiplier determinant).** If `D-V` is invertible, then

\[
 \det R=\det(D-V)\Big(1+z^T(D-V)^{-1}D\mathbf 1\Big).
\]

**Proof.** Factor `R=(D-V)(I+(D-V)^{-1}D\mathbf 1 z^T)` and use
`det(I+uv^T)=1+v^Tu`.

**Lemma 2.2 (a sufficient positivity condition).** Suppose there is a
positive diagonal `W` such that `W(D-V)W^{-1}` is a nonsingular `M`-matrix.
Then `det R>0`. In particular, strict diagonal dominance
`k_i>sum_j v_{ij}` suffices, and more generally a positive diagonal scaling
with all row sums positive suffices.

**Proof.** Write

\[
 W R W^{-1}=M'+(WD\mathbf 1)(z^T W^{-1}),
 \qquad M'=W(D-V)W^{-1}.
\]

The rank-one determinant is
`det M' (1+z^T W^{-1}M'^{-1}WD\mathbf 1)`. Since `M'^{-1}>=0`, `WD\mathbf 1>0`,
and `z^T W^{-1}>=0`, the second factor is at least `1`.

The standard unmixed height bound gives `k_i log p_i>log A_i`, hence strict
diagonal dominance with weights `1/log p_i`; this is the existing selector
proof. For a mixed target the height bound is

\[
 k_i\log(p_i Z)>\log A_i,
\]

which is weaker than `k_i log p_i>log A_i` when `Z>1`. The extra
`k_i log Z` can break the unweighted or `p_i`-weighted row sum. Thus Lemma
2.2 is a real hypothesis, not an automatic consequence of the mixed height
bound. This is the main technical gap of the mixed-density route.

**Lemma 2.3 (hybrid block identity).** Partition the targets into `S` and
`T`. Use the unmixed order row `m_i e_i` for `i in S`. For `i in T`, fix a
common `z`, acquire rows `k_i(e_i+z)-v_i`, and let `D_T`, `V_{TT}`,
`K_T=diag(k_i)` be the restrictions to `T`. Then, after eliminating the `S`
columns with the pure rows,

\[
 \det M=\Big(\prod_{i\in S}m_i\Big)\det(D_T-V_{TT})
 \Big(1+z_T^T(D_T-V_{TT})^{-1}K_T\mathbf 1_T\Big).
\]

If `D_T-V_{TT}` is invertible and the last factor is nonzero, the hybrid row
matrix is full rank; if the base generates `G`, its determinant annihilates
`G`.

This is the algebraic form of the proposed hybrid: pure rows handle the
small-order targets, and a common multiplier is needed only on the remaining
targets. It reduces the independence question for the hard block to one
determinant.

## 3. The medium-order gap

The density theorem covers subgroups of index at most `D^kappa` for fixed
`kappa<1/4` (subgroups) and `kappa<1/64` (all cosets), with the even-product
defect

\[
 \gamma=(8t)^L,\qquad t=b^2Y^{-1/2},\qquad L\log Y=\log R-O(\log Y),
\]

so `log gamma=-(log D)/4+O(log b)` and the covered index is at most
`D^{1/4-o(1)}`. The order row covers targets
of order at most `W`. For a fixed score-one window

\[
 W=2^{c b/\log b},
\]

we have `log W=O(b/log b)` while `log D^{1/4}=Theta(b)`, so for large `b`

\[
 W<D^{1/4-o(1)}.
\]

**Proposition 3.1 (unmixed hybrid barrier).** Consider the selector that
uses only the unmixed order row and the unmixed density row, with window
`W=2^{c b/\log b}` for fixed `c`. A target whose order and index satisfy

\[
 W<\operatorname{ord}(g_p)<D^{1/4-o(1)},\qquad
 [G:\langle g_p\rangle]>D^{1/4-o(1)}
\]

is covered by neither route. The gap is nonempty for large `b`. A class
group with a cyclic factor having a divisor in the interval supplies such a
target; this is not asserted for every input. Hence the unmixed hybrid has no
uniform window bound from these two mechanisms alone.

The statement is method-scoped. It does not say that every input has such a
target, and it does not exclude a different density theorem, a different
target representation, or a direct medium-order algorithm.

The exact finite control below exhibits the logical gap at a deliberately
small window. For `D=-4*10403`, `h=90`, `D^{1/4}=14.282...`, and window
`W=2`:

| target | order | index | covered by order at `W=2`? | covered by density? |
| --- | ---: | ---: | :---: | :---: |
| 23 | 3 | 30 | no | no |
| 19 | 5 | 18 | no | no |
| 11 | 6 | 15 | no | no |
| 17 | 15 | 6 | no | yes |
| 3 | 30 | 3 | no | yes |

The window `W=2` is not an asymptotic window; the table only makes the two
intervals disjoint and nonempty on an actual group. At the certificate's
natural window `15`, target `17` is also order-covered.

## 4. Exact finite hybrid certificate

The following certificate is checked by
[`tools/check_target_order_mixed_rows.py`](../tools/check_target_order_mixed_rows.py).
Let `N=10403=101*103`, `D=-4N`, and use the prime forms of `11,17,19,23`.
The class group has order `90`. The target `19` has order `5` and index `18`;
multiplying by the base prime `23` gives a mixed target of order `15` and
index `6`. At exponent `k=2` the mixed class has the exact smooth hit

\[
 (19\cdot 23)^2 = 17
\]

in the class group, with word `v=e_17`, hence the row

\[
 r_{19}=2(e_{19}+e_{23})-e_{17}.
\]

Use the pure order rows for targets `23`, `17`, and `11`:

\[
 3e_{23},\qquad 15e_{17},\qquad 6e_{11}.
\]

In the coordinate order `(19,23,17,11)` the relation matrix is

\[
 M=
 \begin{pmatrix}
 2 & 2 & -1 & 0\\
 0 & 3 & 0 & 0\\
 0 & 0 & 15 & 0\\
 0 & 0 & 0 & 6
 \end{pmatrix},
 \qquad \det M=540=6|G|.
\]

The determinant is a nonzero multiple of `|G|`, and the checked two-part
extractor returns the form `(101,0,103)` and the factor `101`. The certificate
uses the natural window `15`; its mixed row is at exponent `2`, earlier than
the pure order of target `19`. It verifies the algebra and the extractor, not
a window improvement on an asymptotic family.

The finite table also separates the two mechanisms cleanly: target `17` is
density-covered at index `6`, target `3` is density-covered at index `3`, and
target `19` is a gap target at `W=2` whose mixed index is `6<=D^{1/4}`.

## 5. What would close the route

The remaining obligations are now separated:

1. **Multiplier selection.** Find a common `z` (or a small block of
   multipliers) such that every medium target either has mixed order at most
   `W` or mixed index at most `D^{1/4-o(1)}`. The counting obstruction is
   explicit: the small-order translates of a fixed target are a union of
   subgroups, and the union bound over many targets is not available. A
   positive result would need either a new equidistribution/avoidance lemma
   for prime classes or a structured block construction whose determinant
   factor is provably nonzero.
2. **Mixed-density positivity.** Prove the scaled `M`-matrix condition of
   Lemma 2.2 for the mixed smooth rows, or replace it by a direct
   determinant lower bound. The mixed height bound alone does not supply it.
3. **Index range.** Extend the even-product density beyond
   `D^{1/4-o(1)}` without losing a uniform transfer. The defect computation
   in [the density note](small-index-smooth-density.md) shows why the present
   product hierarchy stops there. A different representation, not a larger
   replay, is required.

None of these is supplied here. The useful content is the reduction: the
order-truncated rows are free, a common multiplier preserves their rank, and
the hard case is exactly the medium-order, large-index block.

## Reproduction

```sh
python3 tools/check_target_order_mixed_rows.py \
  --output /tmp/integer-factoring-target-order-mixed-rows.json
make verify-target-order-mixed-rows
```

The checker enumerates the reduced forms of the small discriminant, computes
orders and indices by exact cyclic enumeration, verifies the mixed identity
and the two-part extraction, and checks the determinant identity on synthetic
nonnegative matrices. It is a finite control; the enumeration cost is not a
factoring algorithm.
