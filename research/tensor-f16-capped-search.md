# Exact capped F16 profiles and constructive kernel searches

## Result and scope

The direct F16 constructor now has a separate accelerated implementation.
It computes the same literal polynomial contractions and exact Smith
orders below a chosen cap. For the searches here, a cap of 64 permits
one binary machine word per coefficient plane. The monomial height stays
separate from the coefficient series.

The completed three-pair pilot contains twelve profiles of 47-letter
codes over F16. Its best numerical fit is approximately
`rho = 0.71273154736`. This does not improve the existing pressure
frontier or the factoring exponent. Larger three-pair, four-pair, and
cubic kernels are separate constructive search inputs. A partial profile
file is not evidence that a full search has completed.

The [accelerator](../tools/tensor_hot_compounds_f16_fast.py) and its
[native implementation](../tools/tensor_hot_compounds_f16_fast.cpp) retain
the existing field arithmetic and minimum-valuation Smith algorithm.
They add a binary subset determinant, an exact capped-jet reduction,
and a lower-bound shortcut for repeated minimum pivot valuations.

## 1. Division-free binary minors

The underlying binary matrix has entries

\[
K_{ij}=
\begin{cases}
(1+t)t^{(i\bmod m)\mathbin{\mathrm{xor}}(j\bmod m)},
  &\lfloor i/m\rfloor=\lfloor j/m\rfloor,\\
[i\bmod m=j\bmod m],&\text{otherwise}.
\end{cases}
\]

For row indices `r_0,...,r_(d-1)` and column indices
`c_0,...,c_(d-1)`, set

\[
e=\sum_i(r_i\bmod m)+\sum_j(c_j\bmod m)+d.
\]

Every nonzero entry has degree at most its two mode indices plus one.
Consequently every complete or partial permutation product has degree
at most `e`. In particular, no intermediate product overflows a 64-bit
binary polynomial when `e<64`.

In characteristic two, determinant signs disappear. For a column subset
`S`, let `D[S]` be the sum of products assigning the first `|S|` rows
bijectively to `S`. Initialize `D[empty]=1` and apply

\[
D[S\cup\{j\}]\mathrel{+}=D[S]K_{r_{|S|},c_j},\qquad j\notin S.
\]

Each multiplication is one shift or the exclusive-or of two shifts.
The recurrence uses `O(d 2^d)` word operations and `O(2^d)` words of
temporary storage. The implementation uses this path for `d<=12`
and `e<64`; other minors use the established exact constructor.
The existing degree bound and endpoint checks remain active.

## 2. Why a 64-coefficient jet is sufficient

Each regular matrix entry has the form `t^e p(t)`, with
`deg p<=e`. The implementation stores `e` separately. Thus a large
regular degree does not require the same number of stored coefficient
bits: when `e<64`, the complete polynomial `p` fits in one word.
Its evaluation at one, used in the full-rank endpoint certificate,
is still exact.

Let `T` be the desired Smith-order cap. The values

\[
\min(n_1,T),\ldots,\min(n_N,T)
\]

depend only on the matrix modulo `t^T`, provided its generic rank is
known. Minimum-valuation elimination can therefore discard every term
of absolute order at least `T` from the start.

For completeness, suppose the selected pivot is `t^v u`, where `u`
is a unit. Every remaining entry has valuation at least `v`. Dividing
the pivot column by `t^v` in the elimination factors loses `v` digits,
but multiplying by a pivot-row entry of valuation at least `v`
restores that precision. All Schur complements are therefore exact
modulo `t^T`. If the residual matrix vanishes modulo `t^T`, every
remaining Smith order is at least `T`. Generic full rank, separately
certified at the endpoint, supplies the number of remaining channels.

The same argument shows that the next minimum valuation cannot be
smaller than the previous pivot. If the pivot scan encounters that
known lower bound, it can stop immediately. This changes only the
choice of an attaining pivot, not the Smith orders.

For `T<=64` and `e<64`, one word per binary plane supplies both the
complete unshifted minor and every required coefficient of the capped
jet. Full-order requests and literal coefficient exports keep their
larger precision rules. The implementation never labels a capped tail
as a list of exact high orders.

If `Z_<T(r)` counts the orders strictly below the cap and `n_tail`
is the number of remaining channels, then

\[
Z_{<T}(r)\le Z(r)\le Z_{<T}(r)+n_{\rm tail}r^T,
\qquad 0<r<1.
\]

The search reports floating pressure fits for both bounding
polynomials. The polynomial bounds are exact. The floating optimizer
does not certify their global minima.

## 3. Independent controls and the completed pilot

The [checker](../tools/check_tensor_hot_compounds_f16_fast.py) verifies:

- 252 binary minors against independent literal permutation sums,
  including 224 subset-path cases and 28 fallback cases;
- 6,044 literal polynomial entries;
- 938 independent binary Toeplitz-rank identities after restriction
  of scalars from F16 to F2;
- 1,664 complete or capped Smith-order comparisons with the previous
  native implementation.

The retained [control evidence](../evidence/tensor-hot-compounds-f16-fast-controls.json)
records every minor and the source hashes. A separate
[47-letter control profile](../evidence/tensor-f16-fast-three-pair-height20-control.json)
agrees with the previous implementation on all 2,209 computed capped
orders. The observed elapsed times were 20.48 seconds for the accelerated
run and 49.75 seconds for the earlier run. These measurements came from
a shared research session; they are not a controlled hardware benchmark.

The [three-pair pilot](../evidence/tensor-f16-three-pair-height20-profiles.json)
uses

\[
T_{12}=aD_0D_5+bD_1D_4+cD_2D_3.
\]

The six recorded coefficient triples are permutations of subsets of
the binary basis `1,2,4,8` of F16. They are linearly independent over
F2, so every nonempty subset has nonzero sum and the rank-one
all-subset Gram condition holds. Each triple is evaluated with the
coordinate retraction and the orthogonal retraction through height 20.
Every resulting physical alphabet has 47 letters.

| Retraction | Number of cases | Best numerical upper-polynomial fit |
|---|---:|---:|
| Coordinate | 6 | 0.71273154736 |
| Orthogonal | 6 | 0.71273157320 |

The twelve records contain 39,762 directly computed capped orders.
For orthogonal maps the infinity profile follows from the established
orthogonal identity; those six additional profiles were not recomputed.
Coefficient choices first change an orthogonal profile at order 57.
For coordinate maps, the first changes occur at orders 62 and 58 in
the two cuts. These differences do not improve the reported frontier.
This is a finite pilot, not an exhaustive F16 classification.

Reproduce the accelerator controls and the larger control profile with:

```sh
python3 tools/check_tensor_hot_compounds_f16_fast.py
python3 tools/tensor_f16_fast_kernel_profiles.py \
  --family three-pair --height 20 \
  --output evidence/tensor-f16-fast-three-pair-height20-control.json
```

The [constructive driver](../tools/tensor_f16_fast_kernel_profiles.py)
also accepts three-pair and four-pair kernels at higher height, and the
cubic kernel

\[
T_{15}=D_0D_1D_5+2D_0D_2D_4+4D_1D_2D_3.
\]

It checks deletion closure, homogeneous heights, `LP=I`, and the
actual Gram inverses before computing a profile. The supported larger
inputs include three-pair height 27 with 108 letters, four-pair height
27 with 115 letters, and cubic height 26 with 108 letters. Their largest
occupation blocks have sizes 3,886, 4,190, and 3,672, respectively.
The larger block limit is 4,608; callers must still budget memory.

## 4. A label-count gate before expensive profiling

The existing [protected height-25 certificate](coherent-protected-height25.md)
has a useful corollary that does not require every label to have height
at most 25.

Suppose each paired physical alphabet has at most 104 labels, its
maps preserve inherited height, and both maps act literally as the
identity on the 28 occupation words of height below 16. Then the same
certificate gives natural pressure at least `17/24`.

To see this, arrange an alphabet's inherited heights as
`h_1<=...<=h_q`. The rank identity `LP=I` bounds its multiplicity at
each height by the dimension of that ambient height cell. Thus its
ordered heights are componentwise no smaller than the first `q`
ambient occupation heights. When `q<=104`,

\[
\sum_{i=1}^q r^{h_i}\le A_{104}(r),\qquad 0<r\le1,
\]

where `A_104` is the existing partition polynomial of all words through
height 25. The protected 28-letter polynomial submatrix is unchanged.
Both inputs to the old exact integer certificate therefore remain
valid. The argument also permits independently chosen paired types.

This gate excludes the four-pair height-26 kernel with 103 labels,
despite its larger maximum height. Its expensive profile run was
stopped. The four-pair height-27 alphabet has 115 labels and escapes
this particular gate. The three-pair and cubic constructions already
change protected states below height 16, so this corollary does not
exclude them.

No new factoring exponent, global pressure optimum, or Lean
formalization is claimed here.
