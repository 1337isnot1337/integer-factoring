# Complete natural profiles from occupation blocks

## Result and scope

The natural profiles of a bit-deletion downset can be computed directly
from small selected compound matrices. The ambient reflection tensor
and its complete four-leg polynomial array need not be constructed.
The new evaluator also computes the two natural cuts only once.

This is an exact finite-seed computation improvement. It does not change
the factoring exponent or assert that a natural gauge is globally
optimal. The reduction applies to the natural gauges proved in
[reflection-natural-cut-duality.md](reflection-natural-cut-duality.md).
Additional physical gauges require a separate integrality argument.

## 1. The selected matrices

Represent a physical letter by the set of its occupied hot modes.
Let `C_X,C_Y` be finite nonempty sets of such letters, each closed under
deleting an occupied bit. A sufficiently large reflection interval
realizes these codes as strict paired restrictions; the
[stabilization theorem](stable-hot-reflection-codes.md) makes their
primitive tensor independent of further interval growth.

For `s` hot modes, work in the deflated variable `t=w^2` and form

```math
B(t)=\begin{pmatrix}D(t)&J(t)\\J(t)&D(t)\end{pmatrix},
\qquad
D_{ab}(t)=(1+t)t^{2(a\mathbin{|}b)+1},
\qquad
J_{ab}(t)=[a=b]t^{2a+1}.                         \tag{1}
```

For every occupation number `k`, use row masks
`x | (y << s)` and column masks `y | (x << s)`, where
`x in C_X`, `y in C_Y`, and `|x|+|y|=k`. The block entry is the
complete polynomial determinant of the corresponding `k`-minor of
`B`. The column order makes the block endpoint exactly the identity:
at `t=1`, `B` exchanges its two mode sets.

If `a_j,b_j` are the occupation distributions of the two physical
codes, the block sizes are

```math
N_k=\sum_j a_jb_{k-j},\qquad \sum_kN_k=|C_X||C_Y|. \tag{2}
```

The downset factorization removes the regular within-row and
within-column Gaussian factors by invertible triangular polynomial
matrices. Therefore the union of these block Smith orders, multiplied
by two to restore the original `w` valuation, is exactly the common
profile `sigma` in the natural-cut duality theorem.

For the primitive original tensor its two natural profiles are
`a+sigma` and `-a+sigma` for one scalar shift `a`. In particular,

```math
Z_0(u)Z_\infty(u)
 =\left(\sum_i e^{-u\sigma_i}\right)^2.           \tag{3}
```

The shift cancels by an exact matrix identity. No polynomial degree
or scalar-content saving is inferred from changing a flattening.

## 2. Complete coefficients and complete profiles

[`hot_code_charge_profiles.py`](../tools/hot_code_charge_profiles.py)
provides the following interfaces:

| Interface | Result |
|---|---|
| `one_particle(s)` | The complete binary polynomial matrix (1) |
| `compound_block(B, rows, columns)` | Every ordinary Smith order of one selected block |
| `compound_block(..., coefficients=True)` | The same profile and every complete polynomial block entry |
| `natural_profiles(C_X, C_Y)` | All blocks and the common profile in the original `w` valuation |

The default interface retains one occupation block at a time in the
native evaluator. It returns orders without writing the full
four-leg tensor into Python. All coefficient arrays remain available
when independently required.

The [native binding](../tools/hot_code_charge_profiles.cpp) includes the
existing selected determinant and finite-field DVR kernels unchanged.
Before computing a selected determinant, the determinant kernel bounds
its polynomial degree by the sum of its row degrees. That bound is
strictly below the series precision. Formal-series inversion therefore
recovers the entire polynomial determinant; no unknown tail is dropped.
The Python wrapper also enforces the conservative bound
`k * max_degree(B)`.

Every selected block checks every endpoint coefficient and verifies
that each row has exactly its expected endpoint match. This proves
full generic rank independently of the finite series window. The
Smith evaluator must resolve that many pivots. If a required pivot is
not resolved, it returns an error and the wrapper retries at a larger
precision. It never reports the resolved prefix as a complete profile.
All weights here are zero and all pivot valuations are integers, so
the numerical pivot-comparison ambiguity of general real gauges is
absent.

The current implementation supports up to sixteen occupied hot-mode
positions and block order at most 4,096. The precision budget is at
most 8,192 deflated coefficients. These are implementation limits,
not restrictions in the mathematical reduction.

## 3. Independent exact controls

Run

```sh
OPENBLAS_NUM_THREADS=1 python tools/check_hot_code_charge_profiles.py \
  --output evidence/hot-code-charge-profile-controls.json
```

The [checker](../tools/check_hot_code_charge_profiles.py) and
[evidence](../evidence/hot-code-charge-profile-controls.json) record:

| Control | Count |
|---|---:|
| Complete selected polynomials checked by independent Laplace expansion | 94 |
| Block endpoint entries | 94 |
| Independent binary quotient-ring rank checks | 48 |
| Retained natural-profile channels, across eleven previous code cases | 5,920 |
| One-particle channels checked through the separate matrix interface | 272 |
| Invalid downset, charge, or endpoint inputs rejected | 5 |

The small controls include unequal paired alphabets in both orientations
and omitted hot-mode positions. Quotient-ring ranks verify
`rank(M mod t^ell)=sum_i max(ell-nu_i,0)` at the positive profile orders
and one level beyond the last order. Complete retained profiles agree
after subtracting their independently recorded scalar shifts.

The separate one-particle controls observe orders
`1,3,...,4s-1` in `t` for each `1<=s<=16`. This finite observation is
not used as a general Smith formula in the evaluator.

For the sixteen-letter stabilized code, block sizes are
`1,12,50,88,73,28,4`. Thus the largest profile matrix has order 88,
while the original cut has order 256. A local replay after native
compilation took approximately 0.04 seconds. Timings are operational
measurements; the exact coefficient and profile comparisons are the
verification gate.

## 4. Subsequent search

The evaluator supports arbitrary downsets, including those that are
not prefixes of natural energy. The bounded worker
[`search_hot_downsets.py`](../tools/search_hot_downsets.py) exchanges
a removable maximal letter with an admissible new letter. Every
proposal passes the downset validation, and every reported cost uses
the complete occupation-block profile.

Its outer scale minimization is numerical. Source hashes, completed
trials, unresolved calls, and time limits are recorded separately.
A finite search minimum does not prove a global lower bound. A proposed
improvement of the factoring threshold still requires an exact
pressure certificate and the full crossing transfer.
