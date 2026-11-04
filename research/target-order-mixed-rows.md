# Target orders, pure relation rows, and common multipliers

**Status.** Ordinary proofs of the order-truncation lemma, the pure-row
determinant, the mixed-order-row determinant, the common-multiplier
determinant identity, the order-index product law, and an order-only
annihilator criterion; a method-scoped medium-order barrier for the unmixed
order/density hybrid; and three exact finite certificates, one mixed and two
order-only (one of them cyclic). No uniform multiplier selection, no
mixed-density positivity
theorem, and no new complete factoring bound is proved. The complete saving
exponent remains **1/2**.

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
extractor returns the form `(101,0,103)` and the factor `101`. This mixed-base
certificate uses the natural window `15`; its mixed row is at exponent `2`,
earlier than the pure order of target `19`. It verifies the algebra and the
extractor, not a window improvement on an asymptotic family. Section 5 gives
a smaller three-target order-only certificate at window `6`, which does not
need the mixed row at all.

The finite table also separates the two mechanisms cleanly: target `17` is
density-covered at index `6`, target `3` is density-covered at index `3`, and
target `19` is a gap target at `W=2` whose mixed index is `6<=D^{1/4}`.

## 5. Order-only rows, the order-index tradeoff, and rank-only extraction

The finite example exposes a simpler mechanism than the mixed smooth row, and
a conservation law explains why target mixing alone cannot remove the middle
of the order/index range.

**Lemma 5.1 (order-index product).** For every `x in G`,

\[
 \operatorname{ord}(x)\,[G:\langle x\rangle]=|G|.
\]

**Proof.** Lagrange's theorem applied to the cyclic subgroup `\langle x\rangle`.

**Corollary 5.2 (coverage tradeoff).** Let the order row cover targets with
`ord(x)<=W`, and let the density row cover targets whose cyclic subgroup has
index at most `Q` (here `Q=D^{1/4-o(1)}`), in the sense that the density
theorem applies to that subgroup; a bound on the first hit exponent is a
separate obligation. A target with

\[
 \operatorname{ord}(x)>W,\qquad [G:\langle x\rangle]>Q
\]

is covered by neither mechanism, and Lemma 5.1 gives `|G|>WQ`. With
`W=2^{c b/\log b}` and `Q=D^{1/4-o(1)}`, this product is
`D^{1/4+o(1)}`, far below `|G|=D^{1/2+o(1)}`. In a cyclic group of order
`h`, a target of order `m|h` is uncovered exactly when

\[
 W<m<h/Q.
\]

For `D=-4*10403`, `h=90`, `W=2`, and `Q=D^{1/4}=14.282...`, the interval is
`(2,6.302...)`, whose divisors are exactly `3,5,6`: the three gap targets in
the table. At `W=6` the interval is empty, which is why the three-target
order-only certificate below succeeds at window `6`.

The product law also delimits target mixing. A common multiplier sends `x`
to another element of the same group, so the pair
`(\operatorname{ord}(xz),[G:\langle xz\rangle])` still lies on the hyperbola
`\operatorname{ord}\cdot[G:\langle\cdot\rangle]=|G|`. A non-generator
multiplier can trade order for index, as `19` to `19\cdot 23` trades
`(5,18)` for `(15,6)`, and can move an individual target from one window to
the other; it does not by itself give a uniform reason for all medium targets
to land in the covered ends. When `G` is cyclic, a **generator** multiplier
puts every image at index `1`, the density extreme, and removes the gap
wholesale. Finding a nonnegative generator word is therefore the real
multiplier obligation. For non-cyclic `G` there is no index-`1` element, and
the best available density extreme is the minimum index of a cyclic
subgroup.

**Lemma 5.3 (positive semigroup).** If `S` generates a finite group `G`, the
positive semigroup generated by `S` is all of `G`. In particular, if `G` is
cyclic, some nonnegative word in the base elements is a generator, so the
multiplier of Corollary 1.3 can be taken nonnegative and its determinant is
positive.

**Proof.** For `s in S`, `s^{-1}=s^{\operatorname{ord}(s)-1}` is a positive
power. The semigroup is closed under inverses, hence is the subgroup
generated by `S`. A subset-sum search is not enough: in `Z/30` the generating
pair `{2,3}` has no subset-sum unit, although `14\cdot 2+3\equiv 1`. This is
a search caution, not an impossibility.

**Lemma 5.4 (order-only annihilator).** Suppose `S={p_1,\dots,p_r}\subset B`
generates `G` and the orders `m_i=\operatorname{ord}(g_{p_i})` are known.
Then

\[
 A=\prod_{i=1}^r m_i
\]

is a multiple of `|G|`, and the diagonal rows `m_i e_{p_i}` are relations
with determinant `A`. Given the generators `g_{p_i}`, the annihilator `A`,
and an input-only bound `B=O(b)` on `\log_2|G|`, the existing two-part
extractor applies. Finding the orders costs at most `\sum_i m_i`
compositions, or `rW` if all `m_i<=W`; the remaining extraction has the
existing `2^{O(\sqrt b)}` meet-in-the-middle factor after small-factor
preprocessing, plus polynomial overhead in `b` and `\log A`. For a semiprime
`N` the ambiguous-form map is nontrivial, so a full two-part basis yields a
proper factor.

**Proof.** The map `\prod_i Z/m_i Z -> G`, `(a_i)|->\prod_i g_{p_i}^{a_i}`,
is surjective because `S` generates `G`. Its domain has order `A`, so `|G|`
divides `A`. Each row `m_i e_{p_i}` is in `\ker\phi`; the determinant is
`A`. Any multiple of `|G|` annihilates `G`, and the extractor accepts an
annihilator. The generation property is a hypothesis, not a consequence of
the row construction: if `S` generates only a proper subgroup, `A` still
annihilates that subgroup, and the extractor may return only the identity.
The finite checker verifies the generation property by enumerating the small
class group.

**Corollary 5.5 (cyclic order-only criterion).** If `G` is cyclic of order
`h`, the subgroup generated by all elements of order at most `W` is the
unique subgroup of order `\gcd(h,\operatorname{lcm}(1,\dots,W))`. It is all
of `G` exactly when every prime-power divisor of `h` is at most `W`. This is
the abstract criterion. A usable order-only certificate additionally needs
the small-order classes to have known base words; in the finite example they
are the target primes `23,19,11`. If the criterion holds but only non-target
classes are small-order, finding their base words is a discrete-log/search
obligation. For a non-cyclic group the analogous condition is that the
small-order elements generate `G`; this is a genuine group-structure
hypothesis, not a consequence of the product law.

**Exact order-only certificate.** For `D=-4*10403`, the class group has
order `90` and exponent `30`, so it is non-cyclic; the cyclic criterion of
Corollary 5.5 does not apply. The three targets `23,19,11` have orders
`3,5,6`, generate the class group, and give the diagonal relation matrix

\[
 \operatorname{diag}(3,5,6),\qquad \det=90=|G|.
\]

The checked two-part extractor returns the form `(101,0,103)` and the factor
`101`. This uses window `6`, no smooth hit, no mixed target, and no density
row. It is a strict finite improvement over the mixed certificate's window
`15`, but it is not a uniform bound: the existence of a small-order
generating subset is exactly the open hypothesis. The mixed hit at exponent
`2` remains a distinct row and is not needed for this order-only
certificate.

**Cyclic gap control.** The second checked group is
`N=3979=23\cdot 173`, `D=-4N`, with cyclic class group of order and exponent
`60`. Here `Q=D^{1/4}=11.232...`, `W=2`, and `h/Q=5.343...`; the divisors
of `h` in the gap interval `(2,5.343...)` are `3,4,5`. The target `29` has
order `5` and index `12>Q`, so it is a genuine gap target. The targets `7`
(order `12`, index `5`), `31` (order `10`, index `6`), and `41` (order
`30`, index `2`) are density-covered by index. The order-only subset
`{29,7}` has orders `5,12`, lcm `60`, determinant `60=|G|`, and generates
the group; the extractor returns factor `23` at window `12`. The abstract
cyclic criterion already holds at `W=5` because `60` divides
`lcm(1,...,5)`, but no target of order `3` or `4` was available in this
small prime window, so the target-only generating window is `12`. This is a
concrete instance of the base-word caveat in Corollary 5.5. It does not give
a uniform improvement: the partner `7` has order `12>Q` even though its
index is density-covered.

**Pooled-rank negative control.** One possible escape from the per-target
order window is to let smooth rows from other targets fill the gap column
collectively. The same group gives a bounded test. Take smooth bound `Y=20`,
which is strictly below the norm bound `H=72`, and the target set
`{5,7,11,13,17,19,29}`; target `29` is included even though it exceeds `Y`,
because the row label is a target and the smoothing primes are at most `Y`.
Enumerating all powers `g_p^k`, `k<=5`, with leading coefficient at most `H`
and a smooth word over primes at most `Y` gives verified rows
`k e_p-v`. The rank profile is

| window | rows | rank |
| ---: | ---: | ---: |
| 1 | 6 | 0 |
| 2 | 11 | 3 |
| 3 | 12 | 4 |
| 4 | 14 | 5 |
| 5 | 16 | 7 |

Full rank `7` first appears at window `5`, which equals `ord(29)=5`; at window
`4` the rank is `5<7`. The selected maximal minor is `900=15|G|`; using it as
an annihilator, the two-part extractor returns `23`. Thus in this example the
off-diagonal smooth rows do not fill the gap column before the gap target's
own order window. A regime with `Y>H` would make every small leading
coefficient automatically smooth; that is not evidence for the density
mechanism. A positive rank result needs a uniform lower bound on the rank of
the pooled smooth rows below the order window, not a single full-rank
instance at the window.

**Mixed-pair complete basis (positive finite control).** The negative control
above used only pure rows for powers of one target. Mixing targets pairwise
does much better on the same group. Take the same `Y=20` and
`H=72`, and all unordered pairs `p<=q` from `{5,7,11,13,17,19,29}`. For each
pair compute the class of `pq`; if its reduced leading coefficient is at most
`H` and is `20`-smooth, the smooth word `v` gives the relation
`e_p+e_q-v`. All `28` pairs were checked. Fifteen have an admissible smooth
word; five of those rows are zero (trivial exact identities such as
`5*5=5^2`), three duplicate an earlier row, and the remaining seven distinct
rows are

\[
\begin{pmatrix}
 1&-2& 0& 0& 0& 1& 0\\
 1& 1& 2& 0& 0& 0& 0\\
-1& 1& 0& 2& 0& 0& 0\\
 0& 0& 1& 0& 1& 0& 1\\
 1&-1& 0& 1& 0& 0& 1\\
 1& 0&-1& 0& 1& 1& 0\\
-1& 0&-1& 0& 0& 2& 0
\end{pmatrix},
\qquad \det=-60=-|G|.
\]

The seven target classes generate `G` (checked), so `[Z^7:\ker\phi]=60`.
The rows lie in the kernel, and their lattice has index `|\det|=60`; hence
the row lattice **equals** the full relation lattice in this example. The
two-part extractor with annihilator `60` returns `23`. This uses window `1`,
no order rows, and no smooth row for `29` itself: the gap target
`ord(29)=5`, index `12` is covered by pair rows. It is the smallest-window
positive finite certificate in this note and directly crosses the
order-index gap for this group.

The certificate is finite, not uniform. It does not show that the pair search
or the density estimate produces a complete relation basis for all targets.
It isolates the uniform obligation: prove that admissible mixed-pair rows at
bounded exponent generate `ker\phi`, or at least have index `|G|^{1+o(1)}`,
under the same smoothness hypotheses. The determinant equality also shows
that the hard part is not a sign or a multiple: it is the index of the
generated relation lattice.

**Second mixed-pair control (non-cyclic).** The same construction also works
on the non-cyclic group `D=-4*10403`, of order `90`. Take `Y=30<H=117` and
targets `{3,11,13,17,19,23}`. Allowing pair powers at exponents `1` and `2`
gives `11` distinct nonzero rows of rank `6`; the selected minor is

\[
\begin{pmatrix}
 2& 0& 2& 0& 1& 0\\
-1& 3& 2& 0& 0& 0\\
 6& 0& 0& 0& 2& 0\\
 1& 1& 0& 1& 0&-1\\
-1&-1& 0& 4& 0& 0\\
-1&-1& 0& 1& 1& 0
\end{pmatrix},
\qquad \det=180=2|G|.
\]

The two-part extractor with annihilator `180` returns `101` at window `2`,
below the order-only window `6`. Thus the pair mechanism is not an artifact
of cyclicity or of the determinant being exactly `|G|`: a full-rank minor
with index `2|G|` already gives an annihilator. The row lattice is not the
full kernel here, only a sublattice of index `2`; full rank is enough.

**Third mixed-pair control (larger gap).** A larger group shows the same
behavior. Take `N=1009*1013`, so `h=928`, `D^{1/4}=44.966...`, `H=1167`, and
`Y=31`. The target set is `{3,7,11,17,19,23,389}`. The small primes have
large orders (`58,116,464,464,464,464`), but small indices; the target `389`
has order `16` and index `58>D^{1/4}`, so it is a genuine gap target. Mixed
pair powers at exponents `1` and `2` give rank `2` at window `1` and full
rank `7` at window `2`. The selected minor has determinant

\[
1856=2|G|,
\]

and the two-part extractor returns `1009`. Window `2` is far below
`ord(389)=16`; the pair rows involving `389` need only smooth words over
primes at most `31`. This is the largest finite gap crossed here and shows
that the mixed-pair window need not track the target order.

The same pair construction has a checked 40-bit instance. For
`N=1000036000099`, cutoff `177`, the split primes through the cutoff lie in
the order-`198850` cycle generated by `5`. Two split primes outside the cutoff,
`1117` and `1289`, have orders `194` and `97` and indices `1025` and `2050`
above `D^{1/4}`. Pair rows at windows `2` and `3` give full rank `23`, selected
minors `4096h` and `1000h`, and extraction of `1000003`. See
[mixed-pair rows on a 40-bit cycle](mixed-pair-gap-cycle.md). This is a fixed
finite certificate, not a uniform rank theorem.

**Lemma 5.6 (maximal-minor annihilator).** Let `R` be an `r x n` integer
relation matrix of rank `n`, `r>=n`, with rows in `\ker\phi`. Every nonzero
`n x n` minor `\Delta` of `R` is a multiple of `|G|`, hence is an
annihilator of `G`.

**Proof.** Let `L` be the row lattice. Since `L` is a rank-`n` sublattice of
`\ker\phi`, `[Z^n:L]=[Z^n:\ker\phi][\ker\phi:L]`, so `|G|` divides
`[Z^n:L]`. The index `[Z^n:L]` is the gcd of the maximal minors of `R`, and
that gcd divides every maximal minor `\Delta`. Thus `|G|` divides `\Delta`.

Consequently a full-rank pooled relation matrix suffices for the extractor;
determinant positivity and the scaled `M`-matrix condition of Lemma 2.2 are
sufficient routes, not necessary ones. The hard problem is the rank of the
pooled matrix, not the sign of a particular mixed determinant.

## 6. What would close the route

The remaining obligations are now separated:

1. **Multiplier selection.** Find a common `z` (or a small block of
   multipliers) such that every medium target either has mixed order at most
   `W` or mixed index at most `D^{1/4-o(1)}`. By Corollary 5.2 a non-generator
   common multiplier stays on the order-index hyperbola; a generator
   multiplier is the clean escape. Lemma 5.3 supplies a nonnegative generator
   word in the cyclic case, but not an efficient search for it. The counting
   obstruction remains: the small-order translates of a fixed target are a
   union of subgroups, and the union bound over many targets is not
   available. A positive result needs a new equidistribution/avoidance lemma
   for prime classes or a structured block construction whose determinant
   factor is provably nonzero. The mixed-pair basis suggests weakening the
   target: do not select one common `z`; allow pair rows `p*q`. Then the
   obligation is a lower bound on the index of the pair-generated relation
   lattice, not a per-target order/multiplier selection.
2. **Mixed-density positivity.** Prove the scaled `M`-matrix condition of
   Lemma 2.2, or replace it by a direct determinant lower bound. By Lemma
   5.6 the route only needs a full-rank pooled matrix; positivity is a
   sufficient, not necessary, interface. The rank gap is the real
   obligation. The pooled-rank control above shows that a full-rank instance
   at the gap window is not by itself a window improvement; the rank must be
   full below that window, uniformly. The mixed-pair basis attains the
   strongest possible finite index (`|G|`) at window `1`, so the open part is
   a uniform relation-basis generation theorem, not a determinant-sign
   estimate. No such theorem is supplied here.
3. **Order-only generation.** Prove that the targets of order at most `W`
   generate the class group, or handle the quotient collectively. For a
   cyclic group Corollary 5.5 reduces this to the prime-power divisor
   condition; for a non-cyclic group even the finite example required a
   specific generating subset `{23,19,11}`, and no uniform criterion is
   supplied.
4. **Index range.** Extend the even-product density beyond
   `D^{1/4-o(1)}` without losing a uniform transfer. The defect computation
   in [the density note](small-index-smooth-density.md) shows why the present
   product hierarchy stops there. A different representation, not a larger
   replay, is required.

None of these is supplied here. The useful content is the reduction: the
order-truncated rows are free, a common multiplier preserves their rank, a
generating subset of small-order targets gives an order-only annihilator, and
the hard case is exactly the medium-order, large-index block. For the present
windows the product `WQ` is `D^{1/4+o(1)}` while `|G|=D^{1/2+o(1)}`, so the
two windows do not meet; covering the hyperbola's middle requires a
genuinely collective rank argument or a different arithmetic representation.
The mixed-pair certificate is finite evidence that the collective route can
attain a complete relation basis at window `1` in a gap example; it is not a
uniform bound and does not change the saving exponent.

## Reproduction

```sh
python3 tools/check_target_order_mixed_rows.py \
  --output /tmp/integer-factoring-target-order-mixed-rows.json
make verify-target-order-mixed-rows
```

The checker enumerates the reduced forms of both small discriminants,
computes orders and indices by exact cyclic enumeration, verifies the mixed
identity and the two-part extraction, verifies that the selected generating
sets generate the order-`90` and cyclic order-`60` groups, checks the
three-target and two-target order-only annihilators and their extraction,
and checks the determinant identity on synthetic nonnegative matrices. It is
a finite control; the enumeration cost is not a factoring algorithm.
