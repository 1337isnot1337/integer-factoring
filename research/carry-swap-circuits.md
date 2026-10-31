# Carry-based row-unimodular swap circuits

These constructions use addition modulo a power of two, so actual carries
enter the label permutations. They give polynomial crossings outside the
bitwise three-CNOT path. Direct circuits, their polynomial row inverses,
and a compressed binary-power variant are row-unimodular at every size.

The checked instances at four, eight, and sixteen labels do not improve
the factoring exponent. Exact low-order column data already exclude their
local Smith-profile criterion. For the inverse paths at eight and sixteen
labels, two explicit minors also exclude every nonzero finite column
basepoint from ordinary coefficient tensor-power constructions.

## 1. An all-size carry circuit

Let `q=2^r` and take labels in `Z/qZ`. Define permutation matrices by

\[
A(u,v)=(u+v,v),\qquad
B(u,v)=(u,v-u),\qquad
C(u,v)=(u,-v).
\]

Applying the rightmost operation first gives

\[
CABA(u,v)=(v,u).
\]

Over a field of characteristic two, put

\[
F_P(z)=I+z(I+P),\qquad
K(z)=F_C(z)F_A(z)F_B(z)F_A(z).                            \tag{1}
\]

The orders of `A,B` divide `q`, and the order of `C` divides two.
For a permutation of order dividing a power of two `m`,

\[
(I+P)^m=0,
\qquad
F_P(z)^{-1}=\sum_{j=0}^{m-1}z^j(I+P)^j.                  \tag{2}
\]

Thus every factor is row-unimodular. Formula (1) has degree at most four,
has `K(0)=I`, and has `K(1)=E`, the crossing. Its polynomial row inverse
also has endpoint `E`, because `E^2=I`. The checker constructs inverses
by (2) and verifies each inverse-factor identity directly.

The all-size degree bound for the inverse is `3(q-1)+1`. The smaller
degrees in the table below are exact finite computations; no unproved
all-size cancellation formula is assumed.

## 2. Compressing the binary powers

For a permutation of order dividing `q=2^r`, define

\[
G_P(z)=\prod_{j=0}^{r-1}\bigl(I+z(I+P^{2^j})\bigr).
\]

Each factor is row-unimodular, and

\[
G_P(1)=P^{1+2+\cdots+2^{r-1}}=P^{q-1}=P^{-1}.
\]

Replace each occurrence of `A,B,A` in (1) by
`G_(A^-1),G_(B^-1),G_(A^-1)`, respectively, and retain `F_C`.
The resulting crossing has degree at most `3r+1` and remains
row-unimodular. It differs from the inverse geometric series because
every binary power receives the same scalar parameter, rather than
parameters of degrees `1,2,4,...`.

The highest coefficient of each individual compressed gate is
`(I+P)^(q-1)`, which can be much thinner than its low-order coefficients.
This explains why several compressed crossings pass the highest
coefficient test even though their subsequent local profiles fail.

## 3. Exact finite profiles

The entries below are ordinary integer ranks over `F2`. The same ranks
hold after any characteristic-two field extension.

| q | Path | Degree | Top row rank | Top column rank | Column orders 0 and 1 |
|---:|---|---:|---:|---:|---|
| 4 | Direct | 4 | 3 | 3 | `0^1,1^3` |
| 4 | Row inverse | 7 | 1 | 1 | `0^1,1^3` |
| 4 | Binary compressed | 4 | 3 | 5 | `0^1,1^3` |
| 8 | Direct | 4 | 19 | 20 | `0^1,1^7` |
| 8 | Row inverse | 15 | 1 | 1 | `0^1,1^7` |
| 8 | Binary compressed | 8 | 1 | 1 | `0^1,1^7` |
| 16 | Direct | 4 | 98 | 98 | `0^1,1^15` |
| 16 | Row inverse | 31 | 1 | 1 | `0^1,1^15` |
| 16 | Binary compressed | 11 | 1 | 1 | `0^1,1^15` |

The checker obtains the order-zero, order-one, and order-two
multiplicities independently from block-Toeplitz ranks. If `T_j`
denotes the rank at precision `j`, those multiplicities are
`T1`, `T2-2T1`, and `T3-2T2+T1`.

At eight and sixteen labels, the direct paths already fail through
their actual highest coefficients:

\[
19\cdot20=380>8^{99/35},\qquad
98^2=9604>16^{99/35}.
\]

The highest coefficient of every ordinary tensor power is exactly
the tensor power of this coefficient, so these two exclusions cover
every such tensor power.

## 4. A cheap rigorous local-profile gate

Let a row-unimodular polynomial matrix have side `n=q^2`, degree `L`,
and top-coefficient row rank `r0`. Write `d_i` for its row Cartan
degrees at infinity. Its constant determinant gives `sum_i d_i=0`,
and `r0` of its degrees equal `L`. Hence, for `s>=1`,

\[
U(s)=\sum_i s^{d_i}\ge\max(n,r_0s^L).                    \tag{3}
\]

If `mu_j` column Smith orders equal `j`, retain any known lower
counts with `j<L`. The formal profile product satisfies

\[
B(s)\ge\max(n,r_0s^L)\sum_j\mu_j s^{-j}.                 \tag{4}
\]

The right side decreases up to
`s0=(n/r0)^(1/L)` and increases afterward. If a rational `a/b`
is at least `s0`, then

\[
\min_{s\ge1}B(s)\ge n\sum_j\mu_j(b/a)^j.                \tag{5}
\]

All these comparisons are verified with integers. The following
small subsets of the column orders suffice:

| q | Path | Counts retained `(mu0,mu1,mu2)` | Rational upper bound for s0 | Lower bound for B |
|---:|---|---|---|---|
| 4 | Direct or binary compressed | `(1,3,1)` | `8/5` | `209/4` |
| 4 | Row inverse | `(1,3,1)` | `3/2` | `496/9` |
| 8 | Row inverse | `(1,7,0)` | `4/3` | `400` |
| 8 | Binary compressed | `(1,7,2)` | `7/4` | `17728/49` |
| 16 | Row inverse | `(1,15,0)` | `6/5` | `3456` |
| 16 | Binary compressed | `(1,15,0)` | `5/3` | `2560` |

Each displayed lower bound exceeds `q^(99/35)`, which is strictly
larger than `q^(2 sqrt(2))`. For a bound `a/b`, the checker proves
the exact comparison `a^35 > q^99 b^35`.

This excludes the sufficient Smith-profile construction at the stated
basepoint. Except for the protected highest coefficients noted above,
it does not assert that every possible recombination of actual
coefficient tensors is excluded.

## 5. Every other finite basepoint of the inverse paths

The evidence file gives two exact column minors for each inverse path.
At `q=8` they have size six and their polynomial gcd is `z^12`.
At `q=16` they have size ten and their gcd is `z^37`.
The checker reconstructs the polynomial crossing, evaluates each
minor determinant, and checks both gcds exactly.

Therefore, over an algebraic closure and at every finite `a!=0`,

\[
c(K_8(a)^{-1})\ge6,\qquad
c(K_{16}(a)^{-1})\ge10.
\]

Here `K_q^-1` means the polynomial row inverse of (1). Its row
matrix is invertible at every finite point, so the corresponding
rank products are at least 384 and 2560. Both exceed `q^(99/35)`.
For any affine change of parameter sending such an `a` to zero,
the constant coefficient of every ordinary tensor power alone
exceeds the strict cost target. This is an all-field finite-basepoint
statement for these two inverse paths.

## Reproduction

Run:

```sh
python3 tools/check_carry_swap_circuits.py \
  --witnesses evidence/carry-swap-finite-minors.json \
  --output /tmp/carry-swap-circuit-audit.json
```

The checker needs only the standard library and the adjacent
`check_nilpotent_flag_crossings.py` helper. It avoids full Smith
calculations and checks the necessary finite jets and minors directly.
