# Odd cyclic crossings from matrix polynomials

## Result and scope

For every odd alphabet size `q`, an arbitrary `q`-square matrix polynomial
`A(z)` with `A(1)=I` gives a strict crossing tensor. One cost flattening
is a direct sum of `q` copies of `A`; the other is diagonalized by a
constant Fourier transform. This reduces both complete local profiles to
one matrix Smith calculation and `q^2` scalar polynomials.

There is also an all-degree obstruction: the **unweighted two-place
profile certificate for this family has `rho >= 3/4`**. This is above the
required `1/sqrt(2)`. Independent physical weights and physical types are
outside that obstruction. A finite search of 133 five-letter weighted
examples did not cross the threshold either; that search is not an
all-family lower bound.

The construction is related to the odd-dimensional phase-point transform
used in discrete Wigner functions; see Gross,
[Hudson's theorem for finite-dimensional quantum systems](https://arxiv.org/abs/quant-ph/0602001).
The crossing normal form, profile formulas, and the elementary obstruction
below are stated directly, so no quantum-information result is needed in
the proof. The factoring score supported by this work remains `1/2`.

## 1. A literal matrix normal form

Let `q >= 3` be odd. Work over a field `K` whose characteristic does not
divide `q`, with physical labels in `Z/qZ`. Division by two in a label
means multiplication by the inverse of two modulo `q`.

For any matrix `A`, define

\[
 T_A(x,y,u,v)
 =[x+y=u+v]\,
 A_{(x+u)/2,\ v-(x-u)/2}.
 \tag{1}
\]

All label arithmetic in (1) is modulo `q`. Apply this map coefficient by
coefficient to a polynomial or Laurent-polynomial matrix `A(z)`.

**Endpoint.** If `A(1)=I`, then

\[
 T_A(1;x,y,u,v)=[x=v][y=u].
 \tag{2}
\]

Indeed, equality of the two matrix indices in (1) implies `x=v`.
Together with the charge constraint, this implies `y=u`. The converse is
immediate. The endpoint is entrywise, before any tensor power or type
projection.

### The column cut

In the cut `(x,u)|(v,y)`, nonzero entries have
`c=x-u=v-y`. For each `c`, put `u=x-c` and `y=v-c`. Its block is

\[
 A_{x-c/2,\,v-c/2}.
 \tag{3}
\]

Independent constant row and column permutations identify this block
with `A`. Thus the entire flattening is `q` identical blocks. In
particular, each Smith order of `A` occurs exactly `q` times, and its
ordinary rank is `q rank(A)`.

### The row cut

In the cut `(x,y)|(v,u)`, nonzero entries have
`a=x+y=v+u`. With row coordinate `x` and column coordinate `v`, the block
entry depends only on `r=x-v`:

\[
 k_a(r)=A_{(a+r)/2,\,(a-r)/2}.
 \tag{4}
\]

Pass, if necessary, to a fixed field extension containing a primitive
`q`th root `zeta`. The Fourier vectors diagonalize each circulant block.
Its eigenpolynomials are

\[
 s_{a,b}(z)=\sum_{x\in\mathbb Z/q\mathbb Z}
 \zeta^{-b(2x-a)} A(z)_{x,a-x},
 \qquad a,b\in\mathbb Z/q\mathbb Z.
 \tag{5}
\]

The Fourier matrices are constant invertible matrices, so they preserve
both local Smith profiles. Field extension also preserves ranks and
orders of nonzero minors. Consequently the complete zero profile of
this cut is the multiset of the orders of the `q^2` polynomials (5),
and its infinity degrees are their polynomial degrees. Every symbol
satisfies `s_(a,b)(1)=1`.

These facts hold for composite odd `q` as well: no field structure on
the label group was used. The characteristic restriction is essential
to Fourier invertibility.

## 2. Complete two-place pressure

First remove a common Laurent monomial from all entries of `A`. Let `D`
be the actual degree of the resulting polynomial matrix. Write

\[
 \mu_i=\operatorname{SmithOrd}_{0,i}(A),\qquad
 \delta_i=D-\operatorname{SmithOrd}_{0,i}(z^D A(z^{-1})),
\]

\[
 \nu_{a,b}=\operatorname{ord}_0 s_{a,b},\qquad
 d_{a,b}=\deg s_{a,b}.
 \tag{6}
\]

The matrix infinity degrees `delta_i` can be negative. They must not be
replaced by entry degrees or truncated at zero. The checker includes
an explicit delayed-pivot example with a Smith order larger than every
entry degree.

For scale `s >= 1`, the two possible orientations have products

\[
 B_0(s)=q\left(\sum_i s^{\delta_i}\right)
                 \left(\sum_{a,b}s^{-\nu_{a,b}}\right),
\]

\[
 B_1(s)=q\left(\sum_{a,b}s^{d_{a,b}}\right)
                 \left(\sum_i s^{-\mu_i}\right).
 \tag{7}
\]

The ordinary unweighted graded-power criterion uses
`rho=log(B)/(4 log q)`, minimized over the scale. A common monomial
shifts the two profiles in opposite directions, leaving each product
unchanged. This justifies the initial normalization.

## 3. An all-degree bare-profile floor

**Theorem.** For every matrix polynomial with endpoint `I`, every scale
`s >= 1`, and both orientations in (7),

\[
 B_0(s)\ge q^3,\qquad B_1(s)\ge q^3.
 \tag{8}
\]

**Proof.** Normalize as above. Some matrix entry has nonzero constant
term, so the least zero Smith order is zero. Some entry has degree
exactly `D`, so the reversed polynomial matrix has an entry with
nonzero constant term. Its least Smith order is also zero. Hence one
matrix infinity degree is exactly `D`.

Each of the `q^2` symbols in (5) is a nonzero polynomial of degree at
most `D`. Thus `0 <= nu_(a,b) <= D`, and

\[
 q\sum_i s^{\delta_i}\ge q s^D,
 \qquad \sum_{a,b}s^{-\nu_{a,b}}\ge q^2s^{-D}.
\]

This proves the first inequality. All symbol degrees are nonnegative,
so their positive powers sum to at least `q^2`. The zero matrix profile
contains an order zero, so its negative powers sum to at least one.
This proves the second inequality. Therefore `rho >= 3/4`.

The proof uses neither a degree bound uniform in `q` nor a restriction
to a prime coefficient field. It allows every regular row and column
frame because complete Smith profiles are invariant under those frames.
It concerns the unweighted profile certificate only. It does not rule
out sharper bounds on the actual coefficient ranks, physical height
gauges, additional reconstruction places, or other crossing families.

## 4. Exact finite controls and a weighted search

The implementation is
[`odd_cyclic_crossing.py`](../tools/odd_cyclic_crossing.py). Its polynomial
matrix elimination starts at precision `q D+1`, exceeding the determinant
degree bound. Every common valuation division reduces the retained
precision. Thus the algorithm does not assume that entry degree bounds
individual pivot orders.

Run the independent controls with:

```sh
python3 tools/check_odd_cyclic_crossing.py
```

The [recorded evidence](../evidence/odd-cyclic-matrix-crossings.json)
checks 9,668 endpoint entries, 9,792 entries in each block identity,
64 literal rank comparisons, 65 independent block-Toeplitz profile
comparisons, and 256 exact rational pressure inequalities. The control
alphabets include `3,5,7,9`. These finite calculations check the
implementation; the all-degree theorem is proved in Section 3.

The weighted driver is
[`search_odd_cyclic_gauges.py`](../tools/search_odd_cyclic_gauges.py).
It starts in characteristic two with the reflection matrix `R` and
the endpoint pencil `R+z(I+R)`. Corrections of the form
`(z^a+z^b)Q` preserve the endpoint. The retained finite catalogue includes
single-entry corrections over the four-element field and random
low-rank corrections. Each seed is searched in both orientations, with
eight independent physical height lists and two separately optimized
paired physical types.

The [five-symbol run](../evidence/odd-cyclic-five-symbol-search.json)
completed 133 seeds and 201,800 native profile queries. Its best floating
ratio was approximately `0.8112781512`, at a distribution close to a
two-letter boundary. It supplies no post-half witness. The driver records
25,059 rejected evaluations: 25,044 unresolved finite-precision profiles,
nine very small type entropies, and six unbracketed positive Gibbs roots.
They are recorded as failures, never as feasible pressure values. The
native precision can be exhausted by large normalized height differences;
these failures are another reason the run is not an exclusion theorem.
The driver retains the complete best matrix coefficients, optimizer outcomes, physical
parameters, distributions, both profiles, and elapsed budget. A numerical
candidate, if found, would still require an exact rational replay and the
existing [effective transfer](graded-even-power-transfer.md).

The finite catalogue is not exhaustive over matrix polynomials, field
coefficients, degrees, physical heights, or stationary points. No
weighted impossibility theorem is inferred from this run.
