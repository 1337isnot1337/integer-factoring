# Multichannel cofactors close the bare interface through eleven letters

For scalar polynomial strict crossings over every field, the bare Smith
pressure satisfies the following **all-degree, all-scale** bounds:

| Physical alphabet | Certified bare rate |
| ---: | ---: |
| 10 | `rho > 71/100` |
| 11 | `rho > 177/250` |

Both exceed `1/sqrt(2)`. Together with the credited
[preceding bounds through nine letters](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/subalphabet-bare-nine-letter-floor.md),
this closes the scalar **bare** interface through eleven letters. The new
input is the [multichannel cofactor inequality](multichannel-rectangular-cofactors.md).
The entropy transfer and earlier prefix/adjugate/first-jet rows are retained
from the preceding work.

The earlier [box-minor budgets](bare-seven-to-twelve-letters.md) provide
stronger uniform bare coverage through twelve letters. The cofactor
certificate here is an alternative argument through eleven, rather than
the current frontier for all available necessary constraints.

This restricts one crossing construction interface. Independently gauged
energies, other endpoints and other factoring mechanisms remain open.
The complete factoring score remains **1/2**.

## 1. Homogeneous rows

Remove a common scalar monomial, let the remaining degree be `M>0`, and
let `mu,nu` be the complete increasing opposite-cut orders in the bare
convention of the preceding proof. Set `z=(mu+nu)/(2M)`, so `z1=0`.
Write `S_k` for prefixes of `z`. We use four kinds of necessary rows:

```
P_q = S_q <= q/2,
A_q = (q-1)S_(q*q) + S_(q*q-q) + (q-1)z2 <= q(q*q-1),
J_(p,q) = n S_(n-1) + (n-p*p)(n-1)z2 <= n(n-1),  n=pq,
K_(p,q) = (q-1)S_n + S_(n-q) + q L_(p,q)(z) <= q(n-1).
```

Here `J` is the preceding rectangular first-jet determinant row, and `K`
is the new multichannel shared-line row. The `L` term selects the `n-1`
smallest members of the multiset repeating every order `p^2` times.
Opposite orientations can be averaged because every row is linear.
Restrictions to smaller physical alphabets are already justified by the
nonnegative absolute-order coefficients and individual Smith interlacing.

## 2. Eight short integer weights per alphabet

Combine the following rows with the displayed integer weights:

| Row | `q=10` | `q=11` |
| --- | ---: | ---: |
| `P_q` | 2148 | 1872 |
| `A_q` | 33 | 23 |
| `J_(6,q)` | 2 | |
| `J_(7,q)` | 2 | 1 |
| `J_(8,q)` | 1 | 1 |
| `J_(9,q)` | | 1 |
| `K_(2,q)` | 64 | 52 |
| `K_(3,q)` | 32 | 27 |
| `K_(4,q)` | 18 | 17 |

Let the resulting inequality be `V dot z<=W`. Define
`pi_i=V_i/(2W)` for `i>=2`, and complete `pi_1=1-sum_(i>=2) pi_i`.
The [exact witnesses](../evidence/multichannel-bare-entropy-witnesses.json)
and [replay](../tools/check_multichannel_bare_entropy.py) verify that every
probability is strictly positive, their sum is one, and
`pi dot z<=1/2`. The complete probability vectors are retained in the
[receipt](../evidence/multichannel-bare-entropy-controls.json).

The weights were discovered numerically. Their validity, rational
normalization and entropy comparison are verified without an optimizer.
The checker independently constructs each displayed row and compares it
with the general cap compilers.

## 3. Uniform pressure consequence

As in the preceding entropy proof, put `t=M log s>0` and write the restored
bare product as `B(s)=s^M (sum s^(-mu_i))(sum s^(-nu_i))`.
Cauchy--Schwarz followed by weighted arithmetic--geometric mean gives

```
sqrt(B(s)) >= sum_i exp(t/2-t*z_i)
           >= exp(H(pi)+t*(1/2-pi dot z)) >= exp(H(pi)).
```

The 100-bit outward replay certifies `H(pi)>2*a*log(q)` for the stated
rational floor `a`. Therefore `B(s)>q^(4a)` and `rho_bare>a` at every
thermal scale and every positive degree. Degree zero has product `q^4`.
The half-rate comparisons are exact: `2*71^2>100^2` and
`2*177^2>250^2`.

The existing fixed auxiliary projection transfer also applies when a
constant contraction produces the scalar strict crossing and its auxiliary
indices stay together on the designated side of each cut. This uses the
previous pressure comparison; it does not assert scalar cofactor identities
for arbitrary matrix-valued entries or rule out independently gauged designs.

Run `make verify-multichannel-bare-entropy`. The receipt pins all compiler,
arithmetic, witness and checker sources. The later
[clean full replay](../evidence/higher-words-full-verification.json) at
`269efb8` includes this entropy package and passed defaults plus fifteen
joint/cofactor targets.
These witnesses only cover through eleven; the linked box theorem already
settles the desired bare exclusion through twelve.
No new Lean formalization, independent mathematical review or improved
factoring algorithm is claimed. Substantial OpenAI assistant support was used.
