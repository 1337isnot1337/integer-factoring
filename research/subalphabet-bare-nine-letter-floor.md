# Rectangular cofactors close the bare interface through nine letters

**Later refinement:** the
[first possible rectangular cofactor jet](rectangular-adjugate-first-jet.md)
improves the nine-letter floor to `709/1000`, with a separate short exact
entropy witness. Every certificate below remains valid.

For scalar polynomial strict crossings over every field, the bare local
Smith pressure satisfies the following **all-degree, all-scale** bounds:

| Physical alphabet | Certified bare rate |
| ---: | ---: |
| 7 | `rho > 719/1000` |
| 8 | `rho > 713/1000` |
| 9 | `rho > 17/24` |

All three rates exceed `1/sqrt(2)`. Together with the preceding
[floor through six letters](shared-minor-smith-floor.md), these exclude
the bare interface through **nine** physical letters. They improve the
earlier [seven-letter entropy witness](adjugate-seven-letter-pressure-floor.md)
and extend its alphabet coverage by two letters.

This is a limitation on one construction interface, not a lower bound on
integer factoring. Independent physical gauges, typed endpoint projections,
and other algorithmic mechanisms remain outside the scalar statement.
The complete factoring score remains **1/2**. The usual constant auxiliary
projection transfer, described below, also applies.

## 1. Homogeneous profile constraints

Use the notation of the seven-letter note: after removing a common
monomial, let the actual degree be `M>0`, let `mu,nu` be the increasing
complete opposite-cut orders at infinity and zero, and set

```math
z_i=(\mu_i+\nu_i)/(2M),\qquad t=M\log s>0.
```

Here `z_1=0`. Averaging valid opposite orientations preserves every
linear constraint. The shared-prefix theorem gives

```math
P_q(z):=\sum_{i\le q}z_i\le q/2.                       \tag{1}
```

Write `n=q^2`. The whole-alphabet adjugate supplies

```math
D_q(z):=(n-1)z_2+n\sum_{i<n}z_i\le n(n-1),             \tag{2}
```

```math
A_q(z):=(q-1)z_2+q\sum_{i\le n}z_i-Top_q(z)
                                            \le q(n-1). \tag{3}
```

The new input is the rectangular
[subalphabet cofactor constraint](subalphabet-adjugate-caps.md).
For `p<q`, put `N=pq`, `g=N-1-p^2>0`. Its shared line of length `q`
gives

```math
R_{p,q}(z):=qg z_2+q\sum_{i\le N}z_i
                    -\sum_{i=N-q+1}^N z_i
                                         \le q(N-1).   \tag{4}
```

This concerns the first `N` full orders. The proof uses individual Smith
interlacing and nonnegative absolute-order coefficients, so it does not
assume that restricted profiles have the same minimum as the full ones.

## 2. Small integer witnesses

The following positive integer weights combine (1)--(4). A blank cell
means zero weight.

| Constraint | `q=7` | `q=8` | `q=9` |
| --- | ---: | ---: | ---: |
| `P_q` | 53987 | 49183 | 45117 |
| `D_q` | 29 | 14 | 7 |
| `A_q` | 1197 | 753 | 504 |
| `R_(3,q)` | 213 | | |
| `R_(4,q)` | 995 | 549 | 236 |
| `R_(5,q)` | 513 | 422 | 309 |
| `R_(6,q)` | | 266 | 226 |
| `R_(7,q)` | | | 148 |

Let `V dot z <= W` be the resulting combined inequality. Divide all
coefficients other than the first by `2W` and set

```math
\pi_i=V_i/(2W)\ (i\ge2),\qquad
\pi_1=1-\sum_{i\ge2}\pi_i.
```

The [exact witness packet](../evidence/subalphabet-bare-entropy-witnesses.json)
and checker verify that all these rational probabilities are strictly
positive and sum to one. Since `z_1=0`, they satisfy

```math
\sum_i\pi_i z_i\le1/2.                                \tag{5}
```

One complete probability vector per alphabet, its repeated runs and every
exact rational normalization are retained in the
[receipt](../evidence/subalphabet-bare-entropy-controls.json).
The weights are fixed constants, independent of the field, kernel,
polynomial degree and thermal scale.

## 3. Uniform entropy proof

Let `B(s)=s^M (sum s^(-mu_i))(sum s^(-nu_i))` be the restored bare
pressure product. Cauchy--Schwarz and weighted arithmetic--geometric mean
give, using (5),

```math
\sqrt{B(s)}\ge\sum_i e^{t/2-tz_i}
 \ge\prod_i(e^{t/2-tz_i}/\pi_i)^{\pi_i}
 =\exp\{H(\pi)+t(1/2-\pi\mathbin\cdot z)\}
 \ge e^{H(\pi)}.                                      \tag{6}
```

The checker certifies `H(pi)>2*a*log(q)` for the stated rational rate
`a`, by 100-bit outward dyadic arithmetic. The lower margins are:

| Alphabet | Certified entropy margin above `2*a*log(q)` |
| ---: | ---: |
| 7 | `>0.00073403` |
| 8 | `>0.00231373` |
| 9 | `>0.00124207` |

It follows from (6) that `B(s)>q^(4a)` at every scale, hence
`rho_bare>a`. In particular the nine-letter comparison to the half-rate
threshold is exact: `2*17^2=578>576=24^2`. For degree zero the product is
`q^4` and the conclusion is immediate. Removing a common scalar monomial
preserves the restored product.

## 4. Auxiliary transfer and exact replay

The [seven-letter proof](adjugate-seven-letter-pressure-floor.md#6-transfer-to-fixed-auxiliary-bare-profiles)
shows that a constant endpoint contraction cannot increase either scalar
Smith pressure relative to its auxiliary cut. This follows by expressing
pressure as a positive weighted sum of truncated multiplication ranks.
Its argument is independent of alphabet size. Consequently the same bare
floors hold for any fixed auxiliary dimension when its auxiliary indices
stay together on the designated side of each cut and a constant contraction
produces the scalar strict crossing. This is a pressure comparison, not a
claim that all scalar cofactor inequalities hold directly for rectangular
auxiliary matrices.

Run `make verify-subalphabet-bare-entropy`. The
[checker](../tools/check_subalphabet_bare_entropy.py) independently constructs
the displayed bare rows, compares them with the general graded implementation,
checks the rational probability identities, and verifies the outward entropy
comparisons and exact half-rate comparisons. It does not invoke a numerical
optimizer. The ordinary all-field proof rests on the preceding cell/shared,
graded cofactor, rectangular rank-gap and interlacing arguments.

The optimization used to discover the weights is exploratory only. A
ten-letter point of the **linear** relaxation falls below the critical
threshold, but fails already known nonlinear second-jet constraints; it is
not a surviving tensor target. The ten-letter bare interface remains open
here, as do independently gauged constructions on smaller alphabets.
Substantial OpenAI assistant support was used; no new Lean formalization
or improved factoring algorithm is included.
