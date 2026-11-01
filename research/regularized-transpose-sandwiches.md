# Weighted transpose words around a rank-one bottleneck

**Constructive interface and registered pilot, 2026-11-01.** This supplies
literal F4 polynomial crossings with controlled leading ranks at both
places. Those ranks alone do not give favorable complete profiles. No
improving pressure, factoring exponent or new Lean formalization is claimed.
The complete score remains **1/2**.

## Regularized transpose factors

For an invertible `q`-by-`q` matrix `A` over a characteristic-two field,
let `S_A(M)=A M^T A^(-T)`. It is an involution. With `N_A=I+S_A`,
`N_A^2=0`, and `N_A` is conjugate by left multiplication with `A` to
`X -> X+X^T`. In particular its rank is `q(q-1)/2`.

For a positive integer `k` define

```
F_A,k(z) = I+(1+z^k) S_A = z^k I+(1+z^k) N_A.
```

Each factor has endpoint `F_A,k(1)=I`. Its determinant is
`z^(k*q^2)` and its inverse is `z^(-2k) F_A,k`. Products are polynomial
and have the literal crossing endpoint in operator order `(X,Y)|(V,U)`.
The opposite matrix must still be reshuffled explicitly.

The product has degree `D=sum k`: its leading operator coefficient is
the product of invertible full transpose channels. If the **number of
factors** is odd, that product is another full transpose channel. Its
opposite realignment therefore has rank one, with both associated
`q`-by-`q` matrices invertible. Degree parity is irrelevant when powers vary.

## A three-factor bottleneck for every alphabet size

Over F4 with `omega^2=omega+1`, take diagonal matrices

```
D0 = diag(1, 1, 1, ..., 1),
D1 = diag(1, omega+1, 1, ..., 1),
D2 = diag(1, omega, omega+1, ..., omega+1).
```

Then `(I+S_D0)(I+S_D1)(I+S_D2)` has rank exactly one, for every `q>=2`.
To see this, diagonal coordinates vanish under each `I+S_D`. Each
off-diagonal unordered pair is a two-dimensional block with rank-one
factor `(d_i,d_j)^T (d_i^(-1),d_j^(-1))`. For pairs `(0,j>=2)`, the
first intervening scalar is zero. For pairs `(1,j>=2)`, the second is
zero. Pairs among labels at least two also vanish. Both scalars are
nonzero on the single pair `(0,1)`. Repeated other-label classes cause
no problem; their codewords need not be distinct.

Consequently `B(z)=F_D0,k0 F_D1,k1 F_D2,k2` has rank one at zero for
any three positive powers. Multiply this bottleneck on either side by
further regularized factors. Its zero coefficient retains rank at most
one; the complete constructor reports whether it vanishes. With an even
number of additional factors, its leading opposite coefficient remains
rank one. Non-diagonal `A` can change the associated vector matrices and
leave the [coordinate support class](rectangular-coordinate-support.md).
They do not automatically improve higher jets or the pressure.

Both complete generic ranks are `q^2`, by the endpoint. The zero-cut
determinant has order `D*q^2`, and every unweighted zero order is at most
`2D`, since the inverse has pole at most `2D`. These are useful checks,
not an opposite-profile formula or a rate certificate.

## Literal implementation and independent controls

[The constructor](../tools/regularized_transpose_words.py) keeps both F4
polynomial bit planes and all four indices. It builds an invertible `A`
as a product of elementary shears, transports both input or output legs,
and applies `T -> T+(1+z^k) S_A T` or its right analogue. NumPy XOR/shift
operations compute the literal coefficient tensor. Total degree at most
62 avoids uint64 overflow. No Gaussian or occupation-closure formula is used.

[Independent controls](../tools/check_regularized_transpose_words.py) build
the shear matrices and inverses explicitly and propagate packed F4
polynomials with dense matrix multiplications. They compare all entries
at three and four letters and 576 selected entries at twelve letters,
**913 complete polynomial entries** in all. They check 21,073 endpoints,
50 independent untruncated small-cut orders, and 338 complete native orders
including all 288 orders of a twelve-letter degree-23 example.

Separate binary restriction-of-scalars ranks verify all six leading
matrices and their vector-matrix ranks, plus all 256 two-by-two F4 matrices
at each of the three sizes. The [receipt](../evidence/regularized-transpose-word-controls.json)
retains the complete word parameters, native profiles, source hashes and
explicit replay scope. The twelve-letter profile does not have a second
full untruncated replay; its independent coefficient controls are selected.

## Registered pressure pilot

The [manifest](../evidence/regularized-transpose-word-pilot-manifest.json)
freezes eight twelve-letter inputs before their pressure fits:

- Two diagonal controls, with base powers `(1,1,1)` and `(6,8,8)`.
- One nearby append pair with base `(6,6,6)` and powers two, degree 22.
- Two-sided nearby pairs at declared powers, each of degree 24.
- Generic non-diagonal four-factor words on each side, degree 14.
- Generic eight-factor words on each side, degree 22.
- The prevalidated full degree-23 control fixture, whose literal digest
  is already known and must match.

A nearby pair uses the same 48-shear matrix and that matrix followed by
one extra shear. Generic factors use 48 recorded random shears each.
All OS-random seeds are frozen in the manifest. No seed is discarded for
poor leading-vector rank or pressure; all eight outcomes are retained.

For each input compute both complete unit profiles, a balanced bare scalar
fit, and three fixed-profile normalized-gap angle fits. Then attempt one
fixed-balanced and one free-angle whole-height fit, with all 88 heights
and both physical costs. Starts have a recorded small random jitter to
leave large tie faces. Each free trial has at most 100 SLSQP iterations
and 15 seconds; the full pilot has 600 seconds. Heights lie in
`[-128,128]`, log scale in `[-9,log(4)]` and log tangent in `[-2,2]`.
Balanced trials fix log tangent to zero.

Every unresolved native query is rejected with coordinates and reason.
Retain all complete unit profiles, initial/best/last free parameters,
best/last native details, query counts, partial statuses and source words.
Intermediate successful free queries are counted rather than stored.
The best literal tensor is saved. An improving floating point requires
fresh independent exact profiles and outward pressure before any claim.
This pilot is not an exhaustive word search or a global gauge exclusion.

```
make verify-regularized-transpose-words
python3 tools/screen_regularized_transpose_words.py --manifest evidence/regularized-transpose-word-pilot-manifest.json --output /path/word-pilot.json
```

Credit: the existing transpose-channel/cofactor language, F4 arithmetic,
complete rectangular Smith backend, Gibbs/asymmetric pressure and
normalized-gap heuristic are reused. The contribution is this literal
weighted sandwich interface and its controls and pilot. Directed by Cinna
Davis / 1337isnot1337 with substantial OpenAI assistant support; no claim
of external novelty or independent mathematical acceptance is made.
