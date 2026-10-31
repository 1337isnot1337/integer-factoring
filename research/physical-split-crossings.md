# Coherent physical projections of polynomial crossings

A constant split projection on a physical label space preserves an exact
crossing endpoint. It can introduce coefficient cancellation that a coordinate
alphabet restriction cannot produce. This note supplies the precise endpoint
convention, a reusable exact implementation, and finite controls over `F4`.
The recorded searches do **not** improve the factoring exponent.

## 1. The split projection identity

Let a polynomial crossing `K(z;X,Y,U,V)` over a field have alphabets
`A,B,B,A` and endpoint

\[
 K(1;X,Y,U,V)=[X=V][Y=U].
\]

Choose matrices

\[
 L_X\in F^{a\times A},\quad R_X\in F^{A\times a},\quad
 L_Y\in F^{b\times B},\quad R_Y\in F^{B\times b}
\]

with `L_X R_X=I_a` and `L_Y R_Y=I_b`. Define

\[
 K'(z;x,y,u,v)=
 \sum_{X,Y,U,V}
 (L_X)_{xX}(L_Y)_{yY}
 K(z;X,Y,U,V)(R_Y)_{Uu}(R_X)_{Vv}.
\tag{1}
\]

Then

\[
 K'(1;x,y,u,v)
 =(L_XR_X)_{xv}(L_YR_Y)_{yu}
 =[x=v][y=u].
\tag{2}
\]

The projection cannot increase the polynomial degree. For fixed ambient
dimensions, constructing all its entries costs a fixed number of field
operations per entry. Thus it preserves the polynomial-family and effective
construction hypotheses of the existing crossing transfer.

Equation (1) explains the transpose convention when every leg is stored as
an output-by-input map: the four maps are
`L_X, L_Y, R_Y^T, R_X^T`. Applying `L` to every leg is valid only under an
additional identity such as `L L^T=I`; that identity is not assumed here.

The two complete cost flattenings have order `ab`. Their generic ranks are
exactly `ab`, because evaluation at one is a permutation matrix. A finite
series profile routine may use this proved rank as its pivot promise.

## 2. Gauge accounting after a projection

When the maps mix only symbols of the same inherited physical energy, they
commute with that diagonal height operator. Its restriction therefore gives
the same gauge accounting on the projected symbols.

An arbitrary constant split is also allowed. In that case, form the actual
new polynomial tensor first, choose heights for its new physical symbols,
and recompute both complete flattening profiles. An old energy cannot simply
be assigned to a mixture of unequal-energy symbols without this recomputation.

The implementation permits separate energy lists `E_X,E_Y` and uses

\[
 h_{0,x}=E_X/2,\quad h_{0,y}=E_Y/2,\qquad
 h_{\infty,u}=E_Y/2,\quad h_{\infty,v}=E_X/2,
\tag{3}
\]

with all other heights zero. Thus the total paired physical costs are exactly
`G_X=E_X` and `G_Y=E_Y`. Let `nu,mu` be the restored zero and infinity
orders after (3). For `t=log(s)`, the independently optimized paired types
are described by the scalar equation

\[
 \log\sum_i e^{-t\nu_i}+\log\sum_i e^{-t\mu_i}
 =R\left(\log\sum_xe^{-tE_X(x)/R}
         +\log\sum_ye^{-tE_Y(y)/R}\right),\qquad \rho=R/2.
\tag{4}
\]

This is the existing
[paired-type formula](paired-type-crossing-pressure.md). Different alphabet
sizes are converted to a square crossing by pairing with the swapped copy,
with separate coefficient indices. The native search minimizes (4) over
`10^-5 <= t <= 2`. Its scale optimization is numerical; an improvement would
still need a strict exact rational certificate.

## 3. Two complete finite projection families

The base tensor is the strict
[reflection crossing](reflection-dual-rail-compression.md).

### Rank-three projections of the four-letter space

For interval `M=4`, enumerate every nonzero functional `ell` on `F4^4`,
normalized so its first nonzero coefficient is one. There are

\[
 (4^4-1)/(4-1)=85
\]

such functionals. For each, enumerate the 64 vectors `r` with `ell(r)=1`.
The projection `I-r ell` has image `ker(ell)` and kernel spanned by `r`.
This gives all 5,440 rank-three projections, with a specified basis of each
image. If `j` is the normalized pivot, use basis vectors
`e_i+ell_i e_j` for `i!=j`, and let the left map extract the nonpivot
coordinates of `(I-r ell)v`. Their product is the identity on the
three-dimensional new alphabet.

For the displayed basis, the search assigns the original natural energies
`(0,12,4,16)` of its three nonpivot indices, then recomputes both profiles.
This is one valid height choice for each seed, not optimization over all
possible bases or all heights. The 5,440 cases give 294 distinct
profile-and-energy classes. Their smallest sampled optimized value is

\[
 \rho\approx0.7416128814255847,
\]

attained by the coordinate alphabet `{0,1,2}`. The original four-letter
crossing has the smaller value `0.7310809359736005` at its natural gauge.

### The degenerate physical energy block at interval eight

For `M=8`, the bit energies are `(28,20,12,4)`. The only nonsingleton
energy eigenspace of the sixteen-letter physical alphabet is spanned by
symbols `6` and `9`, both of energy 32.

There are five one-dimensional subspaces of `F4^2`, and four linear
functionals that take value one on a chosen normalized generator. Thus there
are twenty rank-one split projections of this two-dimensional block.
Each produces a fifteen-letter physical crossing.

Applying the same projection to both crossing pairs gives eight distinct
profile classes. The best numerical value is

\[
 \rho\approx0.7186154898658788.
\]

For example, the combined state may use right vector `(1,0)` and left vector
`(1,1)`. This strictly improves the corresponding coordinate-deletion search
value, approximately `0.7186166115706545`, but remains above the full
sixteen-letter value `0.7166536597281878`.

Applying the projection only to the `X,V` pair gives a `15`-by-`16`
rectangular physical alphabet. The twenty projections give three profile
classes. Their best value is approximately `0.717704215927076`, again above
the full alphabet baseline.

These statements concern the finite projection families and stated height
choices. They are not obstructions for larger coherent codes, other physical
height choices, or projections on longer intervals. The earlier
[four-mode resonance theorem](binary-four-mode-resonance.md) changes an
invertible virtual occupation block while retaining every physical state;
it does not directly cover these rank-reducing physical maps.

## 4. Independent exact controls

The constructor checks `LR=I` and every endpoint entry. The three recorded
searches check respectively 440,640, 1,012,500, and 1,152,000 complete
endpoint entries.

The independent checker also reconstructs each selected tensor with a literal
four-leg contraction, multiplying individual `F4` coefficients rather than
using the constructor's successive-leg bit-plane transforms. Its eighteen
controls include unequal physical alphabets, unequal heights, and 506
polynomial entries with nonbinary coefficients.

For its profile check, restrict scalars from `F4` to `F2`. Multiplication by
`a+b omega`, where `omega^2=omega+1`, becomes

\[
 \begin{pmatrix}a&b\\b&a+b\end{pmatrix}.
\tag{5}
\]

An `F4[[z]]` Smith summand `z^k` becomes two `F2[[z]]` summands of order `k`.
Consequently each complete native Smith order must occur twice in the binary
expanded matrix. The checker proves this equality by independent Python
elimination and then checks the block-Toeplitz rank formula at every distinct
order threshold.

All controls pass: 1,323 complete endpoint entries, 36 independently computed
Smith profiles, and 284 Toeplitz rank checks. A deliberately invalid split
with `LR=0` is rejected.

## 5. Reproduction

Run from the repository root:

```sh
python3 tools/check_physical_split_crossings.py \
  --output evidence/physical-split-crossings.json
python3 tools/search_physical_split_crossings.py --suite hyperplanes \
  --output evidence/physical-split-hyperplanes.json
python3 tools/search_physical_split_crossings.py --suite energy \
  --output evidence/physical-split-energy.json
python3 tools/search_physical_split_crossings.py --suite asymmetric-energy \
  --output evidence/physical-split-asymmetric-energy.json
```

The implementation is
[`physical_split_crossings.py`](../tools/physical_split_crossings.py).
The searches require NumPy, SciPy, and the existing optional C++17 exact
polynomial engine. The checker uses the existing Python binary polynomial
engine as an independent profile calculation.

The `pairs` search suite additionally enumerates all 120 physical-symbol
pairs at `M=8`, three specified split types per pair, and both original
energies and their midpoint as the replacement height. Its scope remains
finite; the three searches summarized above are the completed checkpoints
used in this note.
