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

No Lean formalization is added.
