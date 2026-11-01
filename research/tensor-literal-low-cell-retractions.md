# Literal hot tensors for low-cell retractions without deletion closure

## Construction and scope

A homogeneous physical retraction can alter a low height cell while
retaining every higher cell. Such a map generally fails the deletion
conditions needed by the compressed-compound shortcut. The literal
regular Gaussian tensor still gives a valid finite crossing.

The new [constructor](../tools/tensor_literal_hot_pfaffians_f16.py)
evaluates that literal tensor over F16, including all within-cut Gaussian
terms. Its [native evaluator](../tools/tensor_literal_hot_pfaffians_f16.cpp)
uses a shared binary Pfaffian cache and exact finite jets. It accepts
arbitrary constant height-homogeneous maps satisfying `LP=I`; deletion
closure is not assumed.

For example, the four-dimensional height-15 cell has occupation words
`14,21,35,128`. Replacing this cell by a three-dimensional nondegenerate
subspace, while keeping all other words through height 20, gives **55
physical labels**. The coherent cubic-graph construction with the same
low modification has 51 labels because it imposes higher differential
relations as well. The literal route keeps those four additional labels
and charges their inherited heights.

This is a construction and validation result. No stronger factoring
exponent follows from its existence alone. Completed numerical profiles
must be assessed separately with their explicit capped-tail bounds.

## 1. The literal regular hot kernel

Choose a power-of-two reflection interval `M` with `M>=4m`, and take
the `m` hot modes at the end of each physical leg. In the deflated
variable `t=w^2`, define symmetric `m`-by-`m` matrices

\[
A_{ab}=(1+t)t^{2(a\mathbin{\&}b)},\qquad
D_{ab}=(1+t)t^{2(a\mathbin{|}b)+1},\qquad
J_{ab}=[a=b]t^{2a+1}.
\]

In physical leg order `x,y,u,v`, the regular zero-place alternating
kernel is

\[
\mathcal K_0=
\begin{pmatrix}
0&A&D&J\\
A&0&J&D\\
D&J&0&A\\
J&D&A&0
\end{pmatrix}.
\tag{1}
\]

This follows directly from the existing
[reflection kernel](../tools/reflection_dual_rail_gaussian.py). Put
`s=M/2` and replace each original logical index by `s-1-a`. For example,

\[
(s-1-b)\mathbin{\&}\neg(s-1-a)=a\mathbin{\&}\neg b,
\]

and the zero-place gauge on an active hot variable has height `2a+1`.
The resulting cross exponent is
`2(a&~b)+2b+1=2(a|b)+1`, which gives `D`. The other two nonzero block
types give `A` and `J`. The same-leg reflection partners lie outside
the chosen hot range because `M>=4m`, giving the zero diagonal blocks.

For physical occupation words `x,y,u,v`, the regular tensor coefficient
is the Pfaffian of the selected principal submatrix of (1). Its infinity
kernel is obtained by exchanging the complete `y` and `u` roles.
This statement retains the two `A` blocks; deleting them would silently
reintroduce the coherence assumption.

Every selected kernel entry has degree at most the sum of its two
vertex heights. Therefore every tensor coefficient has degree at most
`h(x)+h(y)+h(u)+h(v)`. It also satisfies the exact reciprocal identity

\[
T_{\infty}(t;x,y,u,v)
=t^{h(x)+h(y)+h(u)+h(v)}T_0(t^{-1};x,y,u,v).
\tag{2}
\]

For (1), this exchanges `A` and `D` and preserves the `J` blocks.
Homogeneous physical contractions preserve (2) term by term.

## 2. Physical compression and the endpoint

Let each embedding column `P_i` and corresponding retraction row `L_i`
have one common inherited height `h_i`, and assume `LP=I`. Apply
`L,L,P,P` to the physical roles `x,y,u,v` of the literal tensor.
The maps need not commute with any deletion.

At `t=1`, both `A` and `D` vanish and `J` is the identity. Hence the
only surviving perfect matching pairs `x` with `v` and `y` with `u`.
The ambient tensor endpoint is exactly

\[
[x=v][y=u].
\]

The actual sparse-map contraction then has the crossing endpoint by
`LP=I`. The native implementation checks this endpoint contraction
for every matrix entry, independently of the coefficient-jet precision.
It therefore knows the complete generic rank even when high coefficients
have been discarded.

The maps preserve the natural diagonal gauges. Multiplying by the
original common scalar monomial restores a polynomial crossing; the
two resulting scalar valuation shifts cancel between the two places.
Identity (2) gives the required infinity restoration. Thus the existing
type-pressure transfer applies with the inherited physical heights.

If `G=P^TP` is invertible and `L=G^{-1}P^T`, the zero and infinity
profiles agree even without deletion closure. To see this, form the
literal compression using `P^T` on every row role and `P` on every
column role. After the common column-label reversal, the actual two
cuts differ from this matrix only by constant factors involving `G^{-1}`
on their row or column tensor factors. These factors are invertible and
preserve height. This is the same constant-matrix argument used for
[orthogonal graph codes](tensor-cubic-graph-orthogonality.md), now applied
before any Gaussian factor is discarded.

## 3. A shared Pfaffian cache

For an even selected mask `S`, choose its first index `i`. The binary
Pfaffian recurrence is

\[
\mathrm{Pf}(S)=\sum_{j\in S\setminus\{i\}}
 \mathcal K_{0,ij}\,\mathrm{Pf}(S\setminus\{i,j\}).
\]

Each kernel entry is a monomial or a monomial times `1+t`. Every
recurrence term therefore requires only polynomial shifts and
exclusive-or operations. The cache is indexed by the complete selected
mask and is shared across matrix entries and parity blocks.

For a common physical height cutoff `H`, let `B_H` be the ambient
occupation downset and let `E_H,O_H` be its even and odd counts.
Every recursive submask still has each leg in `B_H`. Only masks of even
total occupation are cached. Excluding the empty mask, the cache has
the explicit bound

\[
\frac{(|B_H|)^4+(E_H-O_H)^4}{2}-1.
\tag{3}
\]

At height 20, `|B_H|=56`, `E_H=30`, and `O_H=26`, so (3) is
4,917,375 entries. The default budget of eight million entries covers
the entire 55-letter height-15 pilot. A larger requested search that
exceeds the budget fails explicitly; it does not return a partial
profile as a complete result.

For precision `T`, cache construction uses
`O(m |B_H|^4 ceil(T/64))` binary word operations and
`O(|B_H|^4 ceil(T/64))` words of cached polynomial storage, up to
dictionary overhead. Matrix construction and Smith reduction have
their separate, explicitly charged finite costs. This is a practical
constant-size tensor computation, not a new asymptotic factoring
algorithm.

## 4. Exact capped jets and independent controls

The Pfaffian recurrence can run directly modulo `t^T`. Minimum-valuation
elimination then returns each order strictly below `T` and the number
of orders at least `T`. The endpoint proof supplies the remaining
generic rank. The resulting upper counting polynomial replaces every
unresolved tail order by `T`; deleting that tail gives a lower counting
polynomial. The [capped-jet proof](tensor-f16-capped-search.md) applies
without a deletion-closure hypothesis.

Full polynomial exports use precision above the proved complete degree
bound. The independent
[checker](../tools/check_tensor_literal_hot_pfaffians_f16.py) extracts
the kernel from the original interval reflection formula and verifies:

- 3,608 original kernel entries across six hot-mode sizes;
- all 65,536 selected Pfaffians on four hot modes;
- 192 larger Pfaffians and their truncated jets;
- 22,064 complete F16 polynomial matrix entries;
- 272 independent binary Toeplitz-rank identities;
- 1,080 capped-order comparisons against complete profiles;
- 22,032 endpoint entries and the same number of exact infinity
  reciprocal identities.

The controls include a coherent six-letter map, whose complete orders
also agree with the established compound implementation, and a
deliberately noncoherent 12-letter map. The latter changes height 9
but retains height 10; its embedding and retraction each fail a
recorded deletion condition. The literal calculation remains valid.
The [evidence](../evidence/tensor-literal-hot-pfaffians-f16-controls.json)
retains the maps, complete control profiles, independent ranks,
and source hashes.

Reproduce the controls and a 55-letter literal pilot with:

```sh
python3 tools/check_tensor_literal_hot_pfaffians_f16.py
python3 tools/tensor_literal_boundary_profiles.py \
  --height 20 --boundary-height 15 --coefficients 1 1 0 \
  --output evidence/tensor-literal-height15-binary-cap20-profiles.json
```

The [boundary-map builder](../tools/tensor_hot_boundary_retractions.py)
reports deletion failures as diagnostics. It never treats them as
permission to use the compound shortcut. The
[literal search driver](../tools/tensor_literal_boundary_profiles.py)
records exact capped multiplicities, both counting-polynomial fits,
cache use, and whether a second cut was computed or derived from
the orthogonal role-swap identity. Its numerical search samples 81
logarithmic scales and refines detected local minima; no global
numerical optimum is asserted.

## 5. First completed noncoherent pilot

The [completed binary pilot](../evidence/tensor-literal-height15-binary-cap20-profiles.json)
retains every ambient word through height 20, except that its height-15
cell is replaced by

\[
e_{14}+e_{128},\qquad e_{21}+e_{128},\qquad e_{35}.
\]

The retraction is the actual orthogonal Gram inverse. Both map spaces
fail deletion closure, with five recorded failures each. The literal
calculation nevertheless has the exact crossing endpoint.

All 3,025 capped orders of the zero cut were computed. Exactly 2,540
orders are known below 64, and 485 remaining orders are certified to
be at least 64. The infinity profile follows from the orthogonal
role-swap identity. The cache reached precisely the bound of 4,917,375
nonempty masks in (3).

The numerical pressure fit to the upper counting polynomial is
`rho = 0.711115214706` at scale approximately `0.2939093743`.
Deleting the capped tail gives the numerical lower-polynomial fit
`0.711115005973`. These values are not certified global minima and
do not improve the existing pressure frontier. They establish that
the larger noncoherent search can be evaluated with all physical
heights and both places accounted for.

## 6. A complete pressure gate for the first 55-letter family

Before expanding the height-15 coefficient search, its entire natural
pressure family can be excluded from the current compact frontier.
Over any characteristic-two field, replace the height-15 cell by any
three-dimensional physical split `LP=I`, and retain every other ambient
word through height 20. The maps can be nonorthogonal and need not be
deletion closed. Then

\[
\rho\ge\frac{27}{38}=0.7105263157\ldots.
\tag{4}
\]

This is a pruning theorem for this specified family, not a general
lower bound on crossing constructions or factoring algorithms.

All such maps retain the 24 coordinate words of height below 15.
Their literal tensor has the same complete 576-order profile at both
places. Let `Z(r)` be its counting polynomial, and let

\[
A_{55}(r)=\sum_{h(w)\le20}r^{h(w)}-r^{15}.
\]

The protected-submatrix counting argument gives each complete candidate
channel sum at least `Z(r)`. The physical height partition is exactly
`A_55`. The [gate checker](../tools/check_tensor_literal_height15_pressure_gate.py)
proves the exact integer-polynomial comparison

\[
Z(s^{27})^{19}>A_{55}(s^{38})^{27},\qquad 0<s\le1.
\tag{5}
\]

For `0<s<=4/5`, write `A_55(s^38)<=1+s^38+53s^114<=1+2s^38`.
The binomial inequality
`(1+x)^27<=1+27x/(1-27x)` applies at `x=2s^38`.
Also `Z(s^27)^19>=1+19s^27`. The exact positive margin

\[
19\bigl(1-54(4/5)^{38}\bigr)-54(4/5)^{11}>0
\]

therefore establishes (5) on this first interval. The remaining interval
`[4/5,1]` is covered by 622 adjacent rational intervals. On each interval
`[a,b]`, the checker clears denominators and proves
`Z(a^27)^19>A_55(b^38)^27`. Positive coefficients extend each comparison
to the entire interval. The maximum dyadic depth is 15.

Put `r=s^27` in (5). It gives
`log Z(r)>2(27/38) log A_55(r^(38/27))`. Combining the two protected
channel sums with the Gibbs inequality proves (4), including independently
chosen paired physical types. The result uses the stated common natural
scale and inherited heights. Other gauges and coefficient-rank arguments
remain outside its scope.

The [receipt](../evidence/tensor-literal-height15-pressure-gate.json)
retains both core profiles, all intervals, and the source hashes.
Its recomputation checks 1,728 orders: both complete compound profiles
of the deletion-closed core and its separately constructed literal
zero profile. The literal profile has maximum order 88, strictly below
the control precision 128, so no unresolved tail enters (5).

```sh
python3 tools/check_tensor_literal_height15_pressure_gate.py --recompute-profiles
```

The [coordinate-hole control](../evidence/tensor-literal-height15-coordinate-cap20-profiles.json),
which deletes only word 128, has upper-polynomial fit `0.711115213048`.
The [full 56-letter Gaussian control](../evidence/tensor-full-height20-coordinate-profiles.json)
has upper-polynomial fit `0.710824337307`. Both are completed capped
calculations. The full control is better than the two 55-letter pilots;
none improves the existing compact frontier. The complete gate justified
moving the search to lower cells instead of expanding the excluded
height-15 coefficient grid.

## 7. Exhaustive F4 splits in earlier two-word cells

The new [rank-one cell driver](../tools/tensor_literal_rank_one_cell_search.py)
changes a two-word cell at height 8 or 9 and retains every other cell.
It includes all constant splits `LP=1`, including isotropic embeddings
whose Gram matrix cannot supply an orthogonal retraction.

After scaling the one-dimensional physical basis, every pair has one
of the forms

\[
P=(1,a)^T,\quad L=(1+ab,b),\qquad\text{or}\qquad
P=(0,1)^T,\quad L=(b,1).
\tag{6}
\]

Thus a field of order `q_field` has exactly
`q_field(q_field+1)` normalized pairs. Field Frobenius preserves both
profiles because the ambient Pfaffians have binary coefficients. There
are 13 representatives for all 20 F4 pairs, and 76 representatives for
all 272 F16 pairs. This enumeration does not assume that the retraction
is orthogonal.

The [completed F4 census](../evidence/tensor-literal-f4-rank-one-cells-height14-profiles.json)
retains all other words through height 14, giving 23 physical labels.
It computes every representative in both specified cells: 26 literal
profiles covering 40 normalized pairs. Orthogonal cases use the proved
role-swap identity for the second cut; all other second cuts are computed
directly. In total, 24,334 capped orders were directly computed.

| Altered cell | Words | Best retained P and L | Best numerical upper-polynomial fit |
|---|---|---|---:|
| Height 8 | `6,9` | `P=(1,0)^T`, `L=(1,1)` | 0.716170962458 |
| Height 9 | `7,16` | `P=L^T=(0,1)^T` | 0.714950105688 |

These finite fits do not improve the frontier. Across the entire F4
height-8 census, the first profile variation occurs at orders 38 and 40
in the two cuts. For height 9, it occurs at order 38 in each cut. The
coefficients therefore affect only relatively late orders in these
completed small controls. This observation is not an all-field theorem
or an exclusion of larger alphabets.

The independent [cell checker](../tools/check_tensor_literal_rank_one_cells.py)
enumerates normalized projective lines and all possible left inverses
directly, then compares the result with (6), over F2, F4, and F16.
Three nonorthogonal, noncoherent 12-letter controls are also checked
against the original interval kernel. The
[receipt](../evidence/tensor-literal-rank-one-cells-controls.json) records
62,240 complete polynomial entries, 484 independent Toeplitz-rank checks,
2,592 capped-order comparisons, and 62,208 endpoint identities with the
same number of exact infinity reciprocal identities.

```sh
python3 tools/check_tensor_literal_rank_one_cells.py
python3 tools/tensor_literal_rank_one_cell_search.py \
  --field-size 4 --height 14 --boundary-heights 8 9 \
  --output evidence/tensor-literal-f4-rank-one-cells-height14-profiles.json
```

No Lean formalization is added.
