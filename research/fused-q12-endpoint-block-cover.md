# Endpoint normalization leaves a global bare-pressure floor on fused q12

The twelve-letter compressed six-vertex family has bare two-place Smith
pressure rate **`rho>=3/4`**, for every nonzero F16 parameter, every regular
base-field endpoint normalization, every choice of two algebraic spectral
places and every thermal scale. The parameters and endpoints remain in the
specified base-field panel; the spectral places may lie over any extension.

The mechanism tested was joint endpoint normalization. Its inverse is a
unit on one cut and can change the opposite cut, so the existing natural
endpoint profiles did not settle this question. Polynomial block covers
now force opposite rank at least twelve everywhere for this panel. The
proved obstruction ends this bare endpoint sweep; larger numerical sweeps
under these same assumptions are unwarranted.

This is a floor for a pressure bound, not for actual coefficient ranks,
all decompositions, or integer factoring. Physical heights, type selections,
new joint maps and endpoints outside the base-field panel remain outside
the claim. The complete factoring exponent stays **1/2**.

The earlier [box-minor theorem](bare-seven-to-twelve-letters.md) already
excludes a post-half bare certificate for every twelve-letter strict crossing.
This source-specific rank cover proves the stronger floor `3/4` for its
panel; it does not make other bare twelve-letter endpoint normalizations
open post-half targets. Charged physical types and different coefficient
rank estimates require their own analysis.

## 1. Exact source and finite endpoint coverage

Reuse the [compressed fusion constructor](compressed-fusion-root-resonance.md)
with `m=11`, physical alphabet `q=12`, and every `theta in F16*`. The
field encoding, source row order, checked symmetry identities and exact
common-factor removal are unchanged. This contribution does not assume a
new all-spin fusion identity.

Let `F_theta(z)` be the resulting primitive polynomial matrix, of actual
degree `D`. For each `e in P1(F16)` with invertible homogeneous fibre
`E=F_theta(e)`, define `C_e(z)=F_theta(z)E^(-1)`. The inverse and product
are computed in the exact integer-charge blocks. Relabel

```
T_xyuv(z) = (C_e(z))_(12*y+x,12*u+v).
```

At `e` this is the strict crossing. A three-point projective coordinate
change sends `e` to one and the two selected distinct spectral places to
zero and infinity. Its complete homogeneous degree remains `D`.

| Source degree | Parameters | Regular endpoints per parameter | Instances |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 17 | 17 |
| 2 | 2 | 15 | 30 |
| 4 | 4 | 13 | 52 |
| 11 | 8 | 3 | 24 |
| Total | 15 | | 123 |

All 255 parameter/fibre possibilities are checked. The 132 singular
fibres are excluded by the stated inverse-normalization definition. These
123 instances need not be distinct polynomial tensors.

## 2. Disjoint polynomial blocks cover every finite place

Integer charge conservation gives `x+y=u+v`. The opposite matrix therefore
splits into 23 square blocks, indexed by

```
d=u-x=y-v,       -11<=d<=11.
```

For each block the packet retains an explicit polynomial identity

```
sum_j a_j(z) T_(x_j,y_j,u_j,v_j)(z) = g_d(z),           (1)
```

with all indicated entries in that block. Each nonzero `g_d` is completely
factored over F16. The checker verifies (1), factor reconstruction and
Rabin irreducibility certificates, without rerunning the gcd or factor
discovery algorithms.

At a finite point `a`, the block can vanish only if `g_d(a)=0`. Distinct
monic irreducible base-field polynomials remain coprime in every extension.
If `h` is the largest number of block polynomials containing any one factor,
at least `23-h` blocks are nonzero at every finite point. Each nonzero
block contributes at least one to rank, and their rows and columns are
disjoint. Every retained normalization has `23-h>=12`.

At infinity a separate collection of nonzero minors in the leading
coefficient blocks certifies rank at least twelve. Their determinants are
checked by division-free subset expansion. The finite Bézout identities
are not treated as identities of constant degree at infinity.

Consequently the opposite matrix has rank at least twelve at **every**
algebraic spectral place. Constant paired physical basis changes preserve
that rank and both unweighted local profiles.

## 3. The all-scale consequence uses only the determinant budget

For any row-place profile `r` over an extension field, the nonzero row
determinant has homogeneous degree at most `q^2 D`. Hence
`sum_i r_i<=q^2 D`. Jensen's inequality gives, for `s>1`,

```
s^D sum_i s^(-r_i) >= q^2.
```

At any opposite place, at least twelve orders are zero, so
`sum_i s^(-c_i)>=q`. The restored bare product therefore satisfies

```
B(s)=s^D (sum_i s^(-r_i))(sum_i s^(-c_i)) >= q^3,
rho_bare=log(B)/(4 log q) >= 3/4 > 1/sqrt(2).
```

The constant-degree-zero case satisfies the same bound immediately. This
uses the usual bare determinant and entropy framework; the new content is
the exact global block cover after joint endpoint normalization.

## 4. Reproduction, credit and stopping condition

Run `make verify-fused-q12-endpoints`. The
[generator](../tools/fused_q12_endpoint_block_cover.py) constructs the cover;
the [checker](../tools/check_fused_q12_endpoint_block_cover.py) verifies
5,865 endpoint block determinants independently, 2,550,528 strict endpoint
entries, 2,829 polynomial Bézout identities and 607 independent infinity
minors. It also checks the full coefficient identity `C_e(z)E=F_theta(z)`.
No complete Smith-profile window or numerical optimizer is used.

A fresh full replay at exact `acc5a4e36305b03fdff6f7d1e12a0eb7fdad6737`
passed the default gate and nineteen focused targets: 370 Python command
blocks, 66 unit tests, five crossing tests and links in 1,102 Markdown files,
in 924.17 seconds. The frozen checkout stayed clean. The
[full receipt](../evidence/fused-endpoint-full-verification.json) pins the
command and log hash. Later source additions are outside that full-run pin.

The [compressed proof packet](../evidence/fused-q12-endpoint-block-cover.json.xz)
retains all witnesses, and the [receipt](../evidence/fused-q12-endpoint-block-cover-controls.json)
pins the source files and coverage. To regenerate both:

```
python3 tools/check_fused_q12_endpoint_block_cover.py --generate \
  --output evidence/fused-q12-endpoint-block-cover-controls.json
```

The compressed fusion constructor, field representation and original
invariance proofs are credited to their existing contributions. The new
checker reconstructs those seeds and checks their normalized polynomial
identities; it does not repeat an ambient `2^22`-state fusion expansion.

Further work must change an assumption with a specified mechanism. This
result does not justify more bare base-field endpoint searches on these
seeds. Confluent degenerations, different tensor representations or fully
charged physical types require their own construction and cost arguments.
No stronger complete bound, new Lean proof or independent mathematical
review is claimed. Substantial OpenAI assistant support was used.
