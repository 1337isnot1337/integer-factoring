# Complete weighted local profiles at arbitrary closed points over F4

The local profile search now accepts any explicitly supplied irreducible
polynomial over `F4`, including points outside `F16`. The construction
keeps every Smith channel, all physical heights, and the original paired
types. It supports unequal physical alphabets as well as square seeds.

This is an exact reduction of the local linear algebra. Its floating
pressure interface is a search tool; a proposed strict crossing bound
still needs an independent exact certificate.

## 1. Expand at a root and restrict coefficients

Let `f` be a monic irreducible polynomial of degree `e>=1` over `F4`,
and write

```math
E=F_4[\alpha]/(f(\alpha)).
```

For a polynomial cut matrix `A(w)` of size `n`, form `A(alpha+t)` over
`E[t]`. Every irreducible polynomial over a finite field is separable,
so `f(alpha+t)` is `t` times a unit. Thus the orders at `t=0` are
exactly the local orders at the closed point `f`.

Use the basis `1,alpha,...,alpha^(e-1)` of `E/F4`. Replace each
coefficient by its `e` by `e` multiplication matrix over `F4`. This
gives an `ne` by `ne` polynomial matrix `Res(A)`.

Suppose the complete local Smith orders of `A` are
`nu_1,...,nu_n`. The local invertible row and column operations over
`E[[t]]` restrict to invertible operations over `F4[[t]]`. Multiplication
by `t^nu_i` on one `E` coordinate restricts to `t^nu_i I_e`. Therefore
the complete profile of `Res(A)` is

```math
\underbrace{\nu_1,\ldots,\nu_1}_{e},\ldots,
\underbrace{\nu_n,\ldots,\nu_n}_{e}.                    \tag{1}
```

This argument remains valid with rational row and column weights.
After a common ramification, multiply each original row or column by
its uniformizer power, and repeat that weight on its `e` restricted
coordinates. The scalar restriction commutes with those diagonal
operations, so (1) still holds for the weighted orders.

At infinity the implementation uses the stated homogeneous degree `D`
and expands `t^D A(1/t)`. No coefficient-field enlargement is required
there. The exact degree correction is retained in the local profile.

## 2. Charge channels and conjugates separately

Let `sigma_i` denote the weighted orders, including the selected
intrinsic scale. The native restricted matrix returns

```math
\log\sum_{j=1}^{ne}e^{-\widetilde\sigma_j}
   =\log e+\log\sum_{i=1}^n e^{-\sigma_i}.               \tag{2}
```

The adapter subtracts `log e` and returns one copy of each order.
Repeating the physical row and column heights supplies the corresponding
gradient, and the paired type mean is restored on the original physical
alphabet. No extra physical letters are introduced by the coefficient
basis.

This normalization does not make the Hermite degree budget free.
The factor `f` has `e` conjugate roots. If a place group uses `j` of
those roots and has allocation `a`, its local Chernoff scale is
`outer_scale*j/a`. The remaining conjugate roots may be allocated to
the other cut. A pressure implementation must retain their distinct
labels and enforce that the two sets of interpolation roots are
disjoint. The same local order computation applies to all conjugates
because the coefficient matrix and the physical weights are defined
over the base field.

Consequently `e` appears in two different parts of the calculation:
one divides out the repeated channels in (2), while the number of
used conjugate roots is still charged in the allocation scale.

## 3. Explicit computation and bit costs

For each monomial, the implementation uses

```math
[(\alpha+t)^k]_{t^j}
  =\binom{k}{j}\alpha^{k-j}.                            \tag{3}
```

The binomial coefficient is zero or one according to its parity.
All powers of `alpha` are reduced using the supplied monic polynomial.
No splitting field, root enumeration, or floating algebraic number is
needed. Direct Frobenius and greatest-common-divisor checks verify
irreducibility before the matrix is used.

Constructing the expanded matrix takes polynomial work in the cut
size, the entry degree, and `e`. Its size is explicitly `ne`; its
entry degree is at most the original degree. Exact polynomial
elimination and rational physical-height comparisons give a
deterministic polynomial-time local-order algorithm in these input
parameters and their coefficient bit lengths. The implementation's
bounded native precision raises an error if it cannot resolve a full
profile. Such an error is not treated as a rank saving.

The independent exact path computes valuations of the original
polynomial matrix at `f` directly. It does not build the restricted
matrix or use the native pivot sequence. This makes it suitable for
replaying a numerical candidate with rationalized heights.

## 4. Implementation and verification

The [module](../tools/f4_closed_place_oracle.py) exports
`F4ClosedPlaceOracle`. Its paired dimensions are
`(q_x,q_y,q_y,q_x)`, and all physical-height vectors follow that order.
The finite point argument is the coefficient tuple of a monic
irreducible polynomial over `F4`; `None` denotes homogeneous infinity.

The [checker](../tools/check_f4_closed_place_oracle.py) uses the complete
[restricted-core catalogue](../evidence/finite-fiber-restriction-catalogue.json)
to select factors of every degree from one through seven. Its controls
include both equal and unequal paired alphabets.

It reconstructs every selected extension-field coefficient by independent
Horner multiplication by `alpha+t`, compares complete native profiles
with direct rational-height valuations, and checks all physical-height
and intrinsic-scale derivatives. It also agrees with the prior
rational-place adapter and rejects reducible or invalid place inputs.
The [evidence](../evidence/f4-closed-place-oracle-controls.json) records
all exact orders, coefficient controls, and source bindings.

```sh
python3 tools/check_f4_closed_place_oracle.py \
  --output evidence/f4-closed-place-oracle-controls.json
```

The oracle supplies local data. It makes no claim that a particular
place partition, gauge, or common type improves the global crossing
bound.
