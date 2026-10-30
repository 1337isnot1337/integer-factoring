# Arithmetic initial forms of rational crossing paths

## Result and scope

A rational polynomial path whose constant term is the crossing tensor can
be reduced in any prime characteristic after an exact valuation rescaling.
Crossing-preserving diagonal weights and a rational rescaling exponent are
allowed. The construction uses rational arithmetic and integer divisibility;
it does not require a numerical approximation in a valued field.

This gives, in particular, characteristic-two crossing paths from normalized
higher-spin polynomials even when the usual normalizing factorial is zero
modulo two. The normalization is performed over the rationals, all coefficient
valuations are charged, and only integral rescaled coefficients are reduced.
There is no division by a zero field element.

The checked examples do not improve the factoring exponent. Their local
Smith certificates are substantially above the required threshold. The
construction below is an additional exact source of finite candidate tensors,
not a claim that those candidates have sufficiently small ranks.

The local rank criterion used for the experiments is the
[Smith-profile transfer at 28105165](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md).
For a local alphabet of size `q`, its required inequality is
`B(s) < q^(2 sqrt(2))`, where `s >= 1` and `B` is the product of the two
specified positive Laurent sums. A failure of this sufficient certificate
is not a lower bound on the actual coefficient ranks.

## 1. A general exact construction

Let `q >= 2` and let

```math
P(u;x,y,a,b)=\sum_{d=0}^L c_{d;x,y,a,b}u^d
\quad\text{with }c_{d;x,y,a,b}\in\mathbb Q
```

be a completely specified polynomial four-tensor. Assume

```math
P(0;x,y,a,b)=[x=b][y=a].                         \tag{1}
```

Fix a prime `p` and two rational height vectors `h,g` on the labels. Put

```math
w(x,y,a,b)=h_x+g_y-g_a-h_b.                      \tag{2}
```

This weight is zero on every nonzero entry of the crossing. For a nonzero
rational number `c`, let `v_p(c)` denote its integer prime valuation.
If `P` has a nonconstant coefficient, define the rational number

```math
\sigma=\max_{d>0,\ c_{d;x,y,a,b}\ne0}
\frac{-v_p(c_{d;x,y,a,b})-w(x,y,a,b)}d.          \tag{3}
```

For a constant path use `sigma=0`. Every nonzero coefficient then satisfies

```math
v_p(c_{d;x,y,a,b})+w(x,y,a,b)+d\sigma\ge0.       \tag{4}
```

For each coefficient write `c=p^v r`, where `r` has numerator and
denominator coprime to `p`. Define a polynomial over `F_p` by retaining
the residue of `r` exactly when equality holds in (4), and setting the
coefficient to zero when the inequality is strict. Call the result `C(z)`.

**Theorem.** The polynomial `C` satisfies `C(0)=E_q`. Consequently the
coefficients of `K(z)=C(1-z)` sum to `E_q`.

**Proof.** At degree zero, (1) leaves only coefficient one on a crossing
entry. Both its valuation and its weight are zero, so it is retained as
one. Every other constant coefficient remains zero. Therefore `C(0)=E_q`.
Evaluation of `C(1-z)` at `z=1` is `C(0)`, which proves the coefficient
sum identity. This argument proves the endpoint directly. The following
valuation interpretation justifies the term "initial form."

Choose a positive integer `e` clearing the denominators of all heights and
of `sigma`, and use a uniformizer `pi` with `pi^e=p` in the finite extension
of `Q_p` defined by the Eisenstein polynomial `X^e-p`. Multiply tensor entry
`(x,y,a,b)` by `pi^(e w(x,y,a,b))` and substitute `u=pi^(e sigma) z`.
Its coefficient has valuation `e(v+w+d sigma)`. By (4) it is integral.
Reduction modulo `pi` gives exactly the rule above, with residue field
`F_p`. The degree-zero crossing is unchanged by the diagonal weights.
This interpretation works when individual heights or `sigma` are negative;
integrality is required for the complete coefficients in (4). No valued
extension needs to be constructed by the implementation.

All operations in (2)--(4) are finite and exact. In particular, a fractional
critical slope is not rounded to an integer. Rounding upward can erase the
nonconstant coefficients that the construction is intended to retain.

## 2. Higher-spin rational inputs

For concrete inputs let `q=k+1` and work on
`Sym^k(Q^2) tensor Sym^k(Q^2)`. Let `P_j` be its total-spin projector for
spin `j`, with rank `2j+1`. Define

```math
R_k(u)=\sum_{j=0}^k
\left(\prod_{l=1}^j(u+l)\right)
\left(\prod_{l=j+1}^k(u-l)\right)P_j.            \tag{5}
```

This spectral pattern comes from the rational higher-spin `sl_2` matrix;
compare equations (13) and (16) of
[A. G. Bytsko, On higher spin Uq(sl2)-invariant R-matrices](https://arxiv.org/pdf/math/0412482).
That paper uses the checked normalization. Multiplication by the permutation,
a change of sign of the parameter, and a scalar denominator clearing give
(5). No Yang--Baxter identity is needed for the present construction.

The permutation acts on spin `j` by `(-1)^(k-j)`. Evaluating (5) at zero
therefore gives

```math
R_k(0)=k!\,E_q.                                \tag{6}
```

Set `P(u)=R_k(u)/k!` over `Q`, then apply Section 1. This remains legal
when `p <= k`: the negative valuations contributed by `1/k!` occur in
(3), before any reduction.

The implementation constructs the rational projectors within each fixed
total-weight block. On the monomial basis `e_i=x^(k-i)y^i`, the operators
are `E e_i=i e_(i-1)`, `F e_i=(k-i)e_(i+1)`, and
`H e_i=(k-2i)e_i`. The total Casimir is

```math
\mathcal C=(H_{\rm tot}/2)(H_{\rm tot}/2+1)+F_{\rm tot}E_{\rm tot},
```

with eigenvalues `j(j+1)`. Lagrange interpolation gives the exact rational
projectors. A fixed-weight block has size at most `q`, so the implementation
does not invert a dense `q^2` by `q^2` rational matrix.

## 3. Independent finite controls

The default audit checks four sources and their complete reduced tensors:

| Alphabet | Prime | Critical slope | Reduced degree | Best tested Smith `rho`, approximately |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 2 | 0 | 1 | 0.811278 |
| 3 | 3 | 0 | 2 | 0.919776 |
| 5 | 2 | 2 | 1 | 0.969209 |
| 6 | 2 | 2 | 1 | 0.924496 |

The last row uses `h=g=(0,0,0,0,1,1)`; the other rows use zero heights.
The needed value is below `1/sqrt(2)`. The displayed decimals only describe
the profiles; they do not certify a minimum.

The rational source construction is checked against a separately constructed
highest-weight basis. For spin `j`, begin with
`(x_1 y_2-y_1 x_2)^(k-j) x_1^j x_2^j`, then repeatedly apply the total
lowering operator. The resulting basis is independent modulo `1000003`.
Every polynomial coefficient is checked on every basis vector against its
explicit eigenvalue in (5). This is independent of the Casimir projector
algorithm that constructs the source coefficients.

The four examples pass 11,176 rational spectral entry checks and 2,018
literal coefficient-sum endpoint checks. Ordinary block-Toeplitz ranks
independently verify the reported local Smith spectra in 207 checks.
These Toeplitz matrices use scalar modular elimination, not the series
elimination routine that supplies the candidate Smith orders.

An additional synthetic input `E_2+(u^2/3) I_4` has critical slope `1/2`.
Its reduced polynomial is `E_2+z^2 I_4` over `F_3`, providing a control for
ramified slopes. Malformed constant endpoints and composite moduli are
rejected.

For the four recorded families, both cut orientations and every distinct
allowed pair of projective parameter points over the stated prime field
are excluded by exact integer inequalities for this Smith certificate.
Here "allowed" excludes the fixed endpoint `u=0`; the two selected points
must be distinct so the parameter change is invertible.

To certify this without numerical minimization, expand `B(s)=sum_a b_a s^a`.
For `a>0`, pair exponents `a` and `-a`. If `b_a >= b_(-a)`, the pair is at
least `b_a+b_(-a)` for `s>=1`; otherwise it is at least
`floor(2 sqrt(b_a b_(-a)))`. Adding these integers and `b_0` gives an exact
integer lower bound `M`. The audit verifies `M^35 > q^99`. Since
`sqrt(2)<99/70`, this implies `B(s)>q^(2 sqrt(2))` for all `s>=1`.

The optional weighted search at `q=6,p=2` records 813 distinct initial forms
from 2,000 deterministic height choices. Its best recorded profile is the
last row above. This is a finite search, not a classification of rational
height vectors or arithmetic initial forms.

## Reproduction

From the repository root:

```sh
python3 tools/check_arithmetic_crossing.py --output evidence/arithmetic-crossing-audit.json
python3 tools/fused_rmatrix_probe.py --initial --max-q 9 --output evidence/fused-rmatrix-initial-small.json
python3 tools/search_arithmetic_gauges.py --q 6 --p 2 --count 2000 --output evidence/arithmetic-gauges-q6p2.json
```

The scripts use Python's standard library. The general construction is the
proof in Section 1; the finite checks establish only their stated instances.
No Lean proof or stronger numerical factoring exponent accompanies this work.

The construction, implementation, and checks used substantial OpenAI
assistant support. The higher-spin spectral input is attributed above.
