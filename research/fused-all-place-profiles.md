# Complete singular-place profiles for finite fused crossings

## Result and scope

The [literal fused six-vertex constructor](finite-fused-six-vertex-crossings.md)
has checked polynomial crossings on three, four, five, and six physical
letters. This catalogue computes both complete Smith profiles at every
singular place, including points outside the original coefficient field.
The profiles can be used in the
[mixed-place reconstruction theorem](mixed-place-hermite-crossings.md).

This is exact finite spectral data and an effective seed construction.
It does not supply a crossing parameter below `1/sqrt(2)` or a stronger
factoring exponent. Numerical allocation searches must state their own
scope and reproduce their complete degree and physical costs.

## 1. Literal seeds

Use the existing field
`F16=F4[u]/(u^2+u+omega)`, where `omega^2=omega+1`, with the same
integer coefficient encoding as the original constructor. Each instance
is constructed in the full space on `2m` bits. Its selected subspace is
reconstructed exactly after every complete fusion word, the common
polynomial factor is removed exactly, and the remaining endpoint at
one is checked against a nonzero scalar identity. Scalar normalization
and the original physical relabeling give the literal crossing.

| `m` | `theta` encoding | Alphabet `q=m+1` | Actual degree `D` | Degree of either homogeneous determinant `q^2D` |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 2 | 3 | 2 | 18 |
| 3 | 4 | 4 | 3 | 48 |
| 4 | 4 | 5 | 4 | 100 |
| 5 | 4 | 6 | 5 | 180 |

The last two rows are additional successful finite instances of the
same constructor. An unchecked all-spin fusion assertion is unnecessary
for these four seeds. The source's rejected endpoint cases remain
rejected.

## 2. Exact block minors at rational places

The support identity `x+y=u+v` splits both cost cuts into square
blocks. Their largest order is `q`. The checker constructs the literal
bipartite support components, so the decomposition uses row and column
permutations and loses no channels.

For each block, it computes every nonzero minor by complete polynomial
determinants over `F16`. At a finite place `a`, the first `k` Smith
orders sum to the minimum order of these `k`-by-`k` minors at `a`.
At infinity, a minor of the degree-`D` homogeneous section has order

```math
kD-\deg(\Delta).
```

Successive determinantal minima therefore recover the complete
ordered profiles. Independent ranks of local lower block-Toeplitz
matrices check every recovered order.

The rational-place catalogue checks all sixteen finite elements of
`F16` and infinity. It contains 136 complete cut profiles, 4,404
nonzero minors, 1,434 independent block-Toeplitz ranks, and 2,258
literal endpoint entries. Place one has only zero orders in both cuts,
as required by the strict crossing endpoint.

## 3. The rational-place catalogue misses useful column capacity

Summing determinant orders over rational places gives the following
degree accounting:

| Alphabet | First-cut degree covered over `F16` | Second-cut degree covered | Full degree | Missing second-cut degree |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 18 | 18 | 18 | 0 |
| 4 | 48 | 36 | 48 | 12 |
| 5 | 100 | 88 | 100 | 12 |
| 6 | 180 | 124 | 180 | 56 |

The missing degree consists of genuine irreducible finite-place factors.
It is not an infinity normalization error. In particular, searching
only the original field omits available column places in three seeds.

The complete checker factors every support-block determinant over
`F16` by exact Berlekamp linear algebra, with recursive handling of
repeated factors and inseparable polynomial powers. Every proposed
factor is checked independently for irreducibility using Frobenius
gcd tests. Multiplying the factors, with their full multiplicities,
recovers the original determinant including its leading scalar.

For every irreducible factor `f`, the checker again takes the minimum
`f`-adic valuation of every nonzero minor. This recovers both complete
profiles at that closed point. If `e=deg f`, independent multiplication
matrices modulo `f^L` satisfy

```math
rank_{F16}(M\bmod f^L)=e\sum_i(L-\nu_i)_+ .                \tag{1}
```

The implementation checks (1) at enough precisions to recover every
reported order. This verification uses the original polynomial matrix;
it does not require numerical roots or a chosen splitting-field basis.

## 4. Passing to one fixed splitting field

Finite fields are perfect. An irreducible factor of degree `e` splits
into `e` distinct roots over a degree-`e` extension, and each root has
the same minor valuations. Thus one closed-point profile becomes `e`
distinct rational-place profiles in a splitting field. The two cuts
must use the same root labels when allocating shared singular places.

| Alphabet | Distinct finite irreducible factors | Split singular places, including infinity | Sufficient extension degree over `F16` |
| ---: | ---: | ---: | ---: |
| 3 | 3 | 4 | 1 |
| 4 | 10 | 15 | 4 |
| 5 | 12 | 17 | 4 |
| 6 | 21 | 37 | 12 |

After expansion, both determinant ledgers equal `q^2D` exactly in
every row. Thus the catalogue includes every singular place over an
algebraic closure; all omitted places have entirely zero profiles.
Every newly discovered nonlinear place is nonsingular in the first
cut and singular in the second.

The extension degree is fixed once for each fixed seed. It does not
grow with the tensor-power length or factoring input. Constructing a
fixed field and its finitely many roots can be done deterministically
once, including by finite enumeration. The resulting constant cost
does not change the `2^O(lambda)` effective-construction bound. The
catalogue records closed polynomials and abstract conjugate labels;
it does not claim that the splitting fields were explicitly enumerated
during the finite verification.

## 5. Reproduction

```sh
python3 tools/check_fused_multipoint_profiles.py \
  --output evidence/fused-multipoint-profiles.json
python3 tools/check_fused_closed_profiles.py \
  --output evidence/fused-closed-profiles.json
```

The programs retain exact constructor hashes, complete profiles, and
the determinant-degree ledgers. Their imports use the existing fused
constructor and its installed numerical dependencies, but all profile,
factorization, endpoint, and rank checks here use exact finite-field
arithmetic.

The fusion construction is prior repository work. The all-place
enumeration, closed-point completion, and verification code used
OpenAI assistant support. No Lean changes accompany this catalogue.
