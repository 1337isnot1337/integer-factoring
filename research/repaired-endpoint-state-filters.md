# Occupation filters with a diagonal endpoint repair

## Result and scope

The occupation-filter search can allow any nonzero endpoint value in
each occupation state. A constant diagonal map on the physical `Y`
alphabet restores the original endpoint and preserves both Smith
profiles. This enlarges the polynomial synthesis problem: the endpoint
image only needs a vector with no zero coordinate, rather than the
specified all-ones vector.

For one prime-state direction, there is an exact compressed synthesis
algorithm over the original prime field. No field extension is needed.
It returns explicit polynomials, charges their complete infinity
degrees, and retains the exact zero-profile formula. An explicit
qutrit example satisfies degree/order bounds that are impossible with
an unrepaired endpoint-one filter.

These are construction and search results. No strict critical-pressure
witness or improved factoring exponent is claimed.

## 1. The endpoint repair and its type interpretation

Use any of the selections proved in
[the coupled identity reduction](coupled-identity-filter-reduction.md),
Section 7. Let the selected identity directions be `Q_j`, and let
`d` index their joint occupation states. Replace their Gaussian
occupation coefficients by polynomials or Laurent polynomials
`alpha_d(x)`, where `x=-v` and the crossing endpoint is `x=1`. Assume

\[
\lambda_d=\alpha_d(1)\ne0
\quad\hbox{for every state }d.
\tag{1}
\]

At the endpoint, the selected left form in `Q_j` is exactly `Y_j`.
Removing the corresponding Gaussian factors leaves no dependence on
those selected `Y` variables. Consequently, the raw filtered endpoint
coefficient is the original endpoint coefficient multiplied by
`lambda_d`, where `d` is its selected physical `Y` occupation.

Apply the constant physical-`Y` diagonal map with entry `lambda_d^-1`.
The repaired tensor has exactly the original endpoint coefficient by
coefficient. This diagonal map commutes with every diagonal physical
height gauge. It and its inverse are regular at both cuts, so neither
Smith profile changes. It also leaves ordinary polynomial degree
unchanged. This argument permits any finite coefficient field.

For a tensor power restricted to a fixed full physical `Y` type, the
unrepaired scalar is explicitly

\[
\prod_d\lambda_d^{N_d},
\tag{2}
\]

where `N_d` counts positions with selected occupation `d`. It is
nonzero. Per-mode occupation totals alone need not determine (2),
because the dependence on `d` may be nonlinear. The local diagonal
repair avoids that restriction: it restores the original endpoint
pointwise, so the original fixed-mode-total normalization and
triangular-term cancellations apply without change.

The implementation exposes this distinction. The tensor constructors
accept `strict_endpoint=False` to permit (1).
`polynomial_seed(..., strict_endpoint=False, repair_endpoint=True)`
also applies the physical diagonal repair explicitly.

## 2. Exact scalar profiles remain available

For one direction over `F_p`, put `m=nP`, `u=x^m`, and

\[
\beta_k(x)=\sum_{j=0}^{k}(-1)^{k-j}\binom{k}{j}
 (1-x^m)^j\alpha_{p-1-k+j}(x),\qquad 0\le k<p.
\tag{3}
\]

At `x=1`, this gives
`beta_k(1)=(-1)^k alpha_(p-1-k)(1)`, which is nonzero by (1).
Write `d_j=deg_x alpha_j` and `b_k=ord_0 beta_k`. For the initial
identity direction, let `a=P(1-J)` and `kappa=P(1-2n)`. The exact
active profiles are

\[
\{a(r-k)+2b_k:0\le r,k<p\},
\qquad
\{\kappa j-2d_j:0\le j<p\}.
\tag{4}
\]

For the terminal direction at a prime-power interval length, replace
`a` by `P(2n-1-J)` and `kappa` by `-P`. The complementary zero and
infinity modes remain Gaussian. These formulas follow from
[the prime identity diagonalization](prime-identity-filter-diagonalization.md)
and the terminal concentration theorem; the endpoint repair preserves
them.

Signed monomials are now admissible in every prime characteristic:

\[
\alpha_j(x)=c_jx^{-e_j},\qquad c_j\in\mathbb F_p^*.
\tag{5}
\]

Their infinity orders are `kappa*j+2e_j`. Their zero orders come from
the exact coefficients of (3), including cancellations between the
nonzero `c_j`. A common monomial clearing cancels from the two-cut
partition product. A nonmonomial denominator still incurs the full
degree-minus-zero-order clearing cost described in
[polynomial synthesis](polynomial-state-filter-synthesis.md).

## 3. The repaired synthesis criterion

Supply nonnegative degree upper bounds `D_j` and zero-order lower
bounds `B_k`. Let `W` be the vector space of polynomial tuples that
satisfy

\[
\deg\alpha_j\le D_j,\qquad
\mathrm{ord}_0\beta_k\ge B_k.
\tag{6}
\]

The endpoint image is the linear subspace

\[
E=\{(\alpha_0(1),\ldots,\alpha_{p-1}(1)):\alpha\in W\}
 \subseteq\mathbb F_p^p.
\tag{7}
\]

**Theorem.** There is a polynomial tuple satisfying (1) and (6) if
and only if every coordinate projection of `E` is nonzero. If such a
tuple exists over any extension field, one already exists over `F_p`.

Necessity is immediate. For sufficiency, keep a vector `v` in `E`
whose nonzero-coordinate support is growing. If a basis vector `b`
is supported on every coordinate, use it and finish. Otherwise choose
a basis vector that has a nonzero coordinate outside the support of
`v`, and replace `v` by `v+t b`.

Require each coordinate in the union of their supports to stay
nonzero. A coordinate where `b` is zero excludes no value of `t`.
Each coordinate where `b` is nonzero excludes at most one value.
Because `b` is not fully supported, fewer than `p` field values are
excluded. An allowed `t` therefore exists in `F_p`, and the support
strictly increases. Repeating gives the required full-support vector.
Keeping the same combinations of coefficient lifts constructs the
polynomials as well.

Field extension cannot turn an identically zero coordinate projection
into a nonzero one, since all homogeneous equations lie in `F_p`.
The criterion and the construction therefore establish the last claim.

This theorem uses the fact that the number of scalar occupation states
is `p`. For `p^r` joint states with `r>1`, the same coordinate criterion
is sufficient over a fixed sufficiently large extension field, but
need not be sufficient over `F_p`. If `E` has a basis
`b_0,...,b_(h-1)`, the vector `sum_l theta^l b_l` in `F_(p^h)` has
full support whenever every coordinate projection is nonzero: each
coordinate is a nonzero polynomial in `theta` of degree below `h`.
This observation gives a finite-field construction, with its fixed
extension overhead charged in any resulting algorithm. The scalar
implementation does not need that extension.

### Residue compression

The transform (3) depends on `x^m`. Split the polynomials into
residues modulo `m`, exactly as in the endpoint-one synthesis. Only
the cut positions `D_j mod m+1` and `B_k mod m` can change a local
coefficient system. Thus at most `2p+1` residue patterns suffice.

The new algorithm computes a basis of the total endpoint image (7),
together with coefficient lifts. It first tests coordinate support,
then performs the support-growing construction above. For fixed `p`
and bounded ratios `D_j/m`, the field-operation count remains
independent of `m`. Integer arithmetic uses the bit lengths of the
degree and order bounds. Every returned witness is replayed against
the literal transform (3).

## 4. A strictly broader feasible example

Take `p=3`, any `m>1`, and prescribe

\[
(D_0,D_1,D_2)=(0,0,1),\qquad
(B_0,B_1,B_2)=(1,0,1).
\tag{8}
\]

The repaired construction returns

\[
(\alpha_0,\alpha_1,\alpha_2)=(1,-1,x).
\tag{9}
\]

Its endpoint values are `(1,-1,1)`, so its repair is diagonal and
invertible. With `C=1-x^m`, its transformed polynomials are

\[
\beta_0=x,\qquad
\beta_1=1+Cx,\qquad
\beta_2=1-C+C^2x.
\tag{10}
\]

Their zero orders are exactly `(1,0,1)`, and the degrees in (9) are
exactly `(0,0,1)`.

An endpoint-one solution with (8) cannot exist. Its first two
constant polynomials must both equal one. The condition
`ord beta_0>=1` makes `alpha_2` a scalar multiple of `x`, and its
endpoint makes that scalar one. The constant coefficient of
`beta_2=1+C+C^2x` is then `2`, contradicting `ord beta_2>=1`.

This proves that free diagonal repair enlarges the exact search
region. It does not by itself prove a pressure improvement: all
degrees and all zero orders in (10) must enter the full cost.

## 5. Exact verification and reproducibility

The implementations are
[signed scalar profiles](../tools/repaired_endpoint_state_filters.py)
and
[repaired polynomial synthesis](../tools/repaired_endpoint_polynomial_synthesis.py).
Their independent checkers and records are:

- [Endpoint and signed-profile checker](../tools/check_repaired_endpoint_state_filters.py)
  and [its evidence](../evidence/repaired-endpoint-state-filters.json).
  It compares 60 scalar predictions with literal active tensors,
  replays 14 complete raw/repaired physical cut pairs, checks 94,076
  endpoint coefficient identities, and checks 157,427 tensor-power
  entries after fixed-mode-total selection. It also exhibits the
  difference between a full type and mode totals for the raw scalar.
- [Compressed-synthesis checker](../tools/check_repaired_endpoint_polynomial_synthesis.py)
  and [its evidence](../evidence/repaired-endpoint-polynomial-synthesis.json).
  It compares 400 bound queries with 1,700 independent uncompressed
  affine-coordinate systems in characteristics 2, 3, 5, and 7. It
  enumerates 17,496 small nonzero-endpoint polynomial vectors and
  compares 5,103 bound queries with that complete catalog.
- The synthesis checker verifies (8)--(10) at scales including
  `m=2^201+1`. It also constructs two nonconstant repaired polynomial
  tensors, compares both complete physical Smith profiles, and
  checks all 6,642 endpoint coefficients against the original
  Gaussian tensor. The larger full matrix is `81` by `81`.

Run the checks from the repository root:

```bash
python tools/check_repaired_endpoint_state_filters.py \
  --output evidence/repaired-endpoint-state-filters.json
python tools/check_repaired_endpoint_polynomial_synthesis.py \
  --output evidence/repaired-endpoint-polynomial-synthesis.json
```

A single synthesis example is:

```bash
python tools/repaired_endpoint_polynomial_synthesis.py \
  --prime 3 --m 16 --degrees 0,0,1 --orders 1,0,1
```

## 6. Extension-field coefficients and finite searches

The binary coupled construction is also implemented over
`F4=F2[omega]/(omega^2+omega+1)`. The geometric tensors remain over
`F2`; arbitrary nonzero `F4` occupation coefficients enter by linear
combination. The endpoint repair from Section 1 applies over this
field without any profile or degree shift. Since the extension has
fixed degree two, it has a fixed arithmetic overhead.

The implementation is
[the repaired F4 wrapper](../tools/f4_repaired_identity_filters.py).
Its [independent physical checker](../tools/check_f4_repaired_identity_filters.py)
and [exact record](../evidence/f4-repaired-identity-filters.json) include
12 complete `F2` profile regressions, six full original physical tensor
cases through `256` by `256`, and 279,552 raw/repaired endpoint
coefficient equalities. Both physical cuts are checked before and
after repair. Conjugating every field coefficient preserves the
profiles in all six cases. Two explicit native-engine controls require
second-bit-plane cancellation through degree 600; a complete source
beyond the supported precision is rejected.

The [signed scalar search driver](../tools/search_repaired_endpoint_state_filters.py)
enumerates every normalized nonzero coefficient vector in its stated
finite exponent family. At `p=5,n=25,P=1`, each of the initial and
terminal geometries includes 115 exponent vectors and all 256
normalized coefficient vectors: 29,440 coefficient cases per
geometry, with 569 distinct degree/order profiles. The
[initial record](../evidence/repaired-endpoint-p5-initial-n25-screen.json)
and [terminal record](../evidence/repaired-endpoint-p5-terminal-n25-screen.json)
find no same-gauge gain or negative critical pressure on their recorded
grids, at tolerance `10^-9`. The smallest sampled critical pressure is
approximately `0.04557627` in each record. These are finite floating
screens of exact profiles, not global nonexistence results.

The [F4 coupled driver](../tools/search_f4_repaired_identity_filters.py)
similarly records its complete coefficient and exponent catalogs,
exact zero charge blocks, pressure grid, and bounded local refinements.
It keeps the characteristic-two physical alphabet while varying the
nonzero extension-field coefficients.

At `n=32` in the mixed two-mode geometry, the
[complete record](../evidence/f4-repaired-mixed-n32-screen.json)
covers 19,683 coefficient/exponent cases and 1,034 distinct exact
profiles. Its best refined critical pressure is the Gaussian value,
approximately `0.0330679625084`; no same-gauge gain below `10^-9`
was found. This remains a finite floating pressure screen.

The [compiled-profile checker](../tools/check_f4_weyl_filter_compiler.py)
and [its evidence](../evidence/f4-weyl-filter-compiler.json) independently
compare 30 colored reduced profiles with the old complete active
matrices, replay a coefficient at degree 600, and compare 48 pressure
values with expanded full profiles and physical partitions. The largest
recorded pressure discrepancy is below `6e-15`.

```bash
python tools/check_f4_repaired_identity_filters.py \
  --output evidence/f4-repaired-identity-filters.json
python tools/search_repaired_endpoint_state_filters.py \
  --prime 5 --n 25 --geometry initial \
  --output evidence/repaired-endpoint-p5-initial-n25-screen.json
python tools/search_repaired_endpoint_state_filters.py \
  --prime 5 --n 25 --geometry terminal \
  --output evidence/repaired-endpoint-p5-terminal-n25-screen.json
python tools/search_f4_repaired_identity_filters.py \
  --n 32 --geometry mixed \
  --output evidence/f4-repaired-mixed-n32-screen.json
```
