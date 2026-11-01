# Exact local Smith profiles over binary extension fields

## Construction

The [new wrapper](../tools/binary_extension_smith.py) computes unweighted
local Smith orders of polynomial matrices over

\[
K=\mathbb F_2[z]/(f(z)),
\]

where f is a supplied irreducible binary polynomial of degree d. It uses
the existing binary-compatible native F4 matrix engine after restricting
scalars. The d-fold matrix enlargement and all coefficient conversions
are explicit. This provides a checked F8 and F16 interface without
silently treating field arithmetic as free.

The current implementation accepts degrees through 16 and requires
`d*n <= 1024`. It works at zero and infinity. Recentring at other places
must be supplied separately. No new crossing or factoring bound follows
from this utility alone.

## 1. Restriction of scalars

Use the basis `1,z,...,z^(d-1)`. For `a in K`, let R(a) be the binary
matrix of multiplication by a. Replace every coefficient in an n-by-n
polynomial matrix A(t) by R(a), obtaining a `dn`-by-`dn` binary polynomial
matrix R(A).

**Lemma.** If the finite Smith orders of A at zero are
`nu_1,...,nu_r`, then those of R(A) are each repeated d times. The same
statement holds at infinity, including negative Laurent orders.

**Proof.** A local Smith factorization over K is

\[
A=U\,\mathrm{diag}(t^{\nu_1},\ldots,t^{\nu_r},0)\,V,
\]

where U and V are invertible over the local power-series ring. Restriction
of scalars preserves products and inverses. It replaces each diagonal
`t^nu` by `t^nu I_d`. The restricted U and V are invertible over the binary
power-series ring, proving the claim. At infinity use the local parameter
`t inverse`. Extending the resulting binary matrix to F4 does not change
any determinantal valuation, so the existing backend returns the same
orders.

The wrapper checks that each returned order has multiplicity divisible
by d before dividing the multiplicities. It does not infer a generic rank
from an insufficient series window. The generic rank is an explicit
caller-supplied promise; precision exhaustion raises an error. When the
caller supplies a full-rank endpoint, that endpoint proves the promise.

## 2. Exact implementation boundaries

The modulus is checked for irreducibility by binary polynomial Frobenius
and gcd tests. Coefficients are binary integers in the chosen basis.
Multiplication matrices are cached by coefficient value. No factorization
of the modulus or primitive-element search is used.

The dense matrix has `d^2 n^2` scalar polynomial entries. Building it
processes every supplied coefficient and every multiplication-matrix bit.
The native elimination acts on dimension `dn`, not n. Precision starts
large enough to hold the input polynomials and can increase through 8192.
These bounds are implementation limits, not asymptotic speedup claims.

This wrapper intentionally exposes only unweighted integer local orders.
Exploratory floating-weight extensions require their own controls and an
exact replay before any pressure claim.

## 3. Independent finite controls

The [checker](../tools/check_binary_extension_smith.py) tests 64 full-rank
polynomial matrices of orders one through four over F2, F4, F8 and F16.
Each test has identity endpoint at one. Four additional diagonal examples
check an explicitly known generic rank smaller than the matrix order.

It compares against:

- 340 field products computed by a separate carryless-product reduction;
- 640 direct finite-field block-Toeplitz ranks;
- 128 direct small determinant valuation/degree identities;
- eight rejected malformed inputs or reducible moduli.

For a regularized local matrix with orders nu, the Toeplitz test uses

\[
\mathrm{rank}\,\mathcal T_h(A)=\sum_i\max(0,h-\nu_i).
\]

These reference ranks use ordinary finite-field Gaussian elimination,
not the native series-profile algorithm. The
[receipt](../evidence/binary-extension-smith-controls.json) records every
case and all relevant source hashes.
