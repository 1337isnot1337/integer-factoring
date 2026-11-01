# F16 coefficient constructions and exact capped profile computation

## Purpose and status

The coherent-compound implementation now supports the field

```math
\mathbb F_{16}=\mathbb F_4[z]/(z^2+z+\omega),\qquad
\omega^2=\omega+1.
```

This permits orthogonal rank-one coefficient systems with three or four
nonzero columns, and an orthogonal scaling of the six-column F4 hyperoval
after extension to F16. Both are unavailable over F4 in the stated
coefficient interfaces. The implementation computes their actual
homogeneous kernel maps and crossing profiles, with every field
operation exact.

The search targets natural crossing pressure. No new factoring exponent
is asserted. Floating pressure fits are discovery data; a strict rate
claim requires an exact scalar certificate with the charged physical
height budget.

## 1. Explicit field and code constructions

An element is encoded as `a+b*z`, where a and b are F4 elements encoded
by two bits. The four binary basis elements have integer encodings
`1,2,4,8`, representing `1,omega,z,omega*z`. The polynomial defining z
is irreducible because omega has absolute trace one over F2. This is
the second level of the
[explicit quadratic tower](tensor-orthogonal-scaling.md).

For paired differential operators, use coefficient rows

```math
(1,\omega,z),\qquad(1,\omega,z,\omega z).              \tag{1}
```

Their entries are independent over F2. Every nonempty subset therefore
has a nonzero coefficient sum, and hence a nonzero sum of squared
coefficients. All subset Gram conditions hold. The three pairs
`(0,5),(1,4),(2,3)` have common height twelve; the four pairs
`(0,7),(1,6),(2,5),(3,4)` have common height sixteen.

The associated differential kernels and all homogeneous height
truncations admit the orthogonal retraction
`L=(P^T P)^(-1)P^T`. The driver also allows coordinate retractions,
computing their two cuts separately.

A separate explicit extension-field hyperoval construction is

```math
C=\begin{pmatrix}
1&1&1&\omega&0&0\\
0&1&\omega&1&z&0\\
0&1&1+\omega&1+\omega&0&1+z
\end{pmatrix}.                                      \tag{2}
```

It is obtained from the usual F4 hyperoval by column scales
`(1,1,1,omega,z,1+z)`. Every three columns are independent, and all
64 subset Gram ranks equal the selected-column rank. The exact checker
verifies both statements. Thus the
[F4-only hyperoval obstruction](tensor-f4-hyperoval-orthogonal-obstruction.md)
is escaped by a quadratic extension. This is an explicit code
construction, not a pressure certificate for that code.

## 2. Exact finite-field compound evaluation

The new C++ implementation reuses the established binary hot minors
and their polynomial degree bounds. A weighted contraction accumulates
each binary polynomial into four coefficient planes instead of two.
The field multiplication is a quadratic-tower Karatsuba multiplication.
All input scalar coefficients, endpoint identities and polynomial degree
bounds are checked.

The unweighted Smith elimination uses integer valuations, exact field
arithmetic and unit inverses in a finite power-series window. It increases
precision when a complete requested profile is unresolved. Every block
is separately checked to be the identity at the crossing endpoint, so
its full generic rank is established independently of a series window.

The Python layer constructs sparse kernel bases and coordinate or
orthogonal retractions. It checks `LP=I`, homogeneity in both height
and occupation, and deletion closure of both map spaces. This interface
does not yet allow physical labels mixing different occupation numbers.

## 3. Capped orders give a conservative pressure bound

For a positive integer T, the implementation can return

```math
\widehat\nu_i=\min(\nu_i,T),                           \tag{3}
```

instead of resolving every order above T. It records these as capped
orders. An order equal to T is a lower bound for that channel's actual
order, rather than a claim that the actual order equals T.

During exact minimum-valuation elimination, every retained pivot has
the least valuation in the residual matrix. The recorded row operations
are units over the power-series ring. Clearing the pivot row by column
operations leaves the same residual block. If every residual entry is
divisible by `t^T`, then every remaining Smith order is at least T.
The endpoint identity supplies the number of remaining channels.

The truncated arithmetic preserves all residual coefficients below
the absolute precision cap. The routine stops only when that cap is
at least T. A residual zero in the window consequently certifies the
same lower bound T; it is not treated as a smaller generic rank.
These facts prove (3).

For `0<r<1`, let Z be the actual profile polynomial, let Z_low retain
only the certified orders below T, and let n_tail count the remaining
channels. Then

```math
Z_{\rm low}(r)\le Z(r)
\le Z_{\rm low}(r)+n_{\rm tail}r^T.                    \tag{4}
```

The upper bound is the counting polynomial of (3). It can be used in
any sufficient pressure certificate: replacing Z by a larger counting
polynomial makes the claimed rate weaker. The lower bound can screen
uncompetitive candidates. The driver reports numerical fits to both
bounds and keeps the underlying integer multiplicities. It does not
label either fit as an exact global optimization.

For orthogonal retractions, the established equality of the two natural
profiles permits computation of one cut. The second profile is explicitly
marked as derived. For coordinate retractions both cuts are computed.

## 4. Independent exact controls

Run

```sh
python3 tools/check_tensor_hot_compounds_f16.py
```

The [checker](../tools/check_tensor_hot_compounds_f16.py) constructs
literal compound entries using binary integer polynomial arithmetic
and permutation determinants. It then checks complete Smith orders
through block-Toeplitz ranks after restriction of scalars from F16 to F2.
Those rank calculations use ordinary binary Gaussian elimination and
do not call either native Smith implementation.

The [receipt](../evidence/tensor-hot-compounds-f16-controls.json) records:

- All 256 scalar multiplication comparisons with the quadratic tower.
- 6044 independently constructed literal polynomial entries.
- 938 independent block-Toeplitz rank checks, covering 416 complete orders.
- 1248 capped-order comparisons against those complete profiles.
- 144 complete-order comparisons with the embedded F4 implementation.
- All 64 subset conditions for the explicit hyperoval matrix (2).

These controls include genuine F16 coefficients in six- and eight-letter
codes, along with F4 subfield controls. They validate the new arithmetic
and profile path; they do not establish a stronger factoring bound.

Implementation files:

- [F16 Python construction layer](../tools/tensor_hot_compounds_f16.py)
- [F16 native compound and Smith engine](../tools/tensor_hot_compounds_f16.cpp)
- [Three-/four-pair search driver](../tools/tensor_f16_rank_one_pair_profiles.py)

## 5. A pruned coefficient family and continuing searches

The full height-eight/height-twelve coupled family has 3375 nonzero F16
coefficient triples. Frobenius conjugation reduces these to 851 exact
profile classes. The
[driver](../tools/tensor_f16_coupled_quadratic_profiles.py) starts with
classes represented by coefficients from the binary basis and records
every completed profile. Its source records distinguish an incomplete
bounded pilot from an exhausted finite catalogue.

The first 24 classes did not improve the binary coefficient reference.
An exact [protected-core floor](tensor-coupled-quadratic-field-floor.md)
then proved `rho>=5/7` for this whole coefficient family under the
specified retractions, over every characteristic-two extension. The
remaining classes were not computed. Search effort moved to the
three-/four-pair constructions in (1), including cutoffs above the
height-25 protected-boundary theorem.
