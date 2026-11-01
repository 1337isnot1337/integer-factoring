# Sharp gauge profiles for two uniform ternary bases

## Result and scope

For each of two specified physical bases of the ternary crossing, the
known gauges simultaneously minimize **every local Chernoff profile sum**
over all rational physical-letter gauges, when the paired types are
uniform. The statement applies separately at the row places `2,infinity`
and the column place `0`. Thus allowing a different gauge for each jet
does not improve these local exponential bounds.

This is a sharp statement about the specified bases, types and local
Smith method. It is not an optimum over other physical bases, nonuniform
types, actual typed jet ranks, polynomial physical frames, or arbitrary
crossing decompositions. The factoring exponent remains `1/2`.

The result arose from the
[weighted gauge-catalogue construction](weighted-adaptive-gauge-catalogues.md).
A numerical probe found no improvement, and small exact minor frames
then proved the stronger uniform statement below. The proof uses those
finite frames, not extrapolation from the numerical search.

## 1. The fixed tensors and reference profiles

Start with

\[
 K(z;x,y,u,v)=[x+y=u+v]\bigl(1-(x-u-z(v-u))^2\bigr)
\]

over `F3`. Apply `A,A,A^-T,A^-T` to the four physical legs, using either

\[
 A_{\rm ev}=\begin{pmatrix}1&0&0\\0&2&1\\2&2&2\end{pmatrix},
 \qquad
 A_{\rm Pa}=\begin{pmatrix}1&0&0\\2&1&0\\1&1&1\end{pmatrix}.
\]

Both maps preserve the exact crossing endpoint. The paired types put
probability `1/3` on each symbol of each leg. In particular, their entropy
sum is `H=2 log3`.

Write

\[
 \gamma=(0,\tfrac12,\tfrac12,1,1,1,\tfrac32,\tfrac32,2).
 \tag{1}
\]

For both bases, the sharp **effective** profiles, after subtracting the
uniform physical mean, are:

| Place | Cut | Effective profile |
| --- | --- | --- |
| `2` | `(x,y)\mid(u,v)` | `gamma` |
| infinity | `(x,y)\mid(u,v)` | `gamma` |
| `0` | `(x,u)\mid(y,v)` | `2 gamma` |

Infinity uses the degree-two homogeneous local matrix
`t^2 K(t^-1)`. It uses the same degree budget as the opposite-flag
construction, rather than the raw-infinity convention of the two-place
supporting-line note.

The attaining slopes on `x,y,u,v` are, respectively,
`(0,0,1/2,1/2)`, `(-1/2,-1/2,0,0)`, and `(-1,-1,1,1)`.
The height of symbol `i` is its slope times the ordinary integer `i`.
The exact checker obtains all six complete profiles from independent
ramified block-Toeplitz ranks. These are the profiles of the
[certified ternary witness](ternary-basis-mixed-crossing.md).

## 2. A certificate that holds for every physical gauge

Set the height of symbol zero on each leg to zero. This loses no
generality: subtracting a constant height from a leg shifts every Smith
order and the physical mean by the same amount.
Let `h` be the remaining eight heights and put

\[
 A(h)=\frac13\sum_{j=1}^{8}h_j.
 \tag{2}
\]

For a nonzero `k`-minor `M`, let `v_M` be its ungauged local valuation.
Let `c_M` be its eight physical-letter counts, from its selected rows
and columns. Its gauged valuation is exactly

\[
 v_M+c_M\mathbin{\cdot}h.
 \tag{3}
\]

If `nu_1(h)<=...<=nu_9(h)` are the complete Smith orders, the
determinantal characterization gives

\[
 \sum_{i=1}^k\nu_i(h)=\min_M\bigl(v_M+c_M\mathbin{\cdot}h\bigr).
 \tag{4}
\]

Rational heights can be cleared by one ramification; dividing back
gives the same formula. No sign or valuation estimate for an individual
matrix entry substitutes for the complete minor minimum.

For each `k=1,...,9`, the certificate supplies nonnegative rational
weights `lambda_M` on a small list of nonzero minors such that

\[
 \sum_M\lambda_M=1,\qquad
 \sum_M\lambda_M c_M=(k/3,\ldots,k/3),\qquad
 \sum_M\lambda_M v_M=C_k.                              \tag{5}
\]

Here the sharp row prefix bounds are

\[
 (C_1,\ldots,C_9)=(0,\tfrac12,1,2,3,4,\tfrac{11}2,7,9),
 \tag{6}
\]

and the column bounds are twice (6). Averaging (3) with (5), and using
that a minimum is at most this average, proves

\[
 \sum_{i=1}^k\bigl(\nu_i(h)-A(h)\bigr)\le C_k
 \qquad\text{for every rational }h.                    \tag{7}
\]

For `k=9` equality holds. Direct determinants are `1+z^9` in the first
cut and `2z^18` in the second, for both bases. The full determinant's
physical count is three of every nonzero letter, so its gauge shift is
exactly `9A(h)`. The attaining profiles in Section 1 reach every bound
in (7) simultaneously.

The certificate contains 54 frames across the two bases and three
locations. The verifier independently recomputes all 234 listed minor
determinants by fraction-free elimination over `F3[z]`, verifies their
local valuations, and checks all rational equalities (5). The optional
LP proposer is not trusted by this verification.

## 3. Majorization gives every local pressure bound

Let `sigma_i=nu_i(h)-A(h)` and let `gamma_i` be the appropriate
reference profile. Equations (7) and the equal total sum imply

\[
 \sum_i e^{-s\sigma_i}\ge\sum_i e^{-s\gamma_i}
 \qquad(s\ge0).                                         \tag{8}
\]

This is the convex-majorization inequality. It can also be seen from
the elementary hinge formula

\[
 \sum_i(t-\sigma_i)_+
   =\max_{0\le k\le9}\left(kt-\sum_{i=1}^k\sigma_i\right).
\]

Every hinge sum for `sigma` is at least the corresponding sum for
`gamma`. Integrating these inequalities against the nonnegative second
derivative of `exp(-s t)`, with the equal sums handling the affine
term, proves (8).

At jet density `d`, every local gauge therefore has Chernoff bound

\[
 \inf_{s\ge0}\left\{\log\sum_i e^{-s\sigma_i}+sd\right\}
 \ge
 \inf_{s\ge0}\left\{\log\sum_i e^{-s\gamma_i}+sd\right\}.
 \tag{9}
\]

The clean gauge attains the right side. The type-dimension cap `H` is
common to every gauge and preserves the comparison. Taking the minimum
over any finite catalogue also preserves it. The result covers arbitrary
nonlinear assignments of heights to the three letters, not only the
linear slopes displayed in Section 1.

The complete determinant formulas also show that other places are
nonsingular in the relevant cut. At such a place the effective orders
sum to zero. Jensen's inequality gives a profile sum at least nine for
every `s>=0`, so those places supply no cheaper local exponential rate
at a nonnegative jet density under these uniform types.

## 4. Consequence for the inherited adaptive rate calculation

Put `Z(s)=1+exp(-s)+exp(-2s)`. The optimal local rate functions are

\[
 \Phi_R(d)=\inf_{s\ge0}\{2\log Z(s/2)+sd\},\qquad
 \Phi_C(d)=\inf_{s\ge0}\{2\log Z(s)+sd\}.
 \tag{10}
\]

The two row places have the same nondecreasing function. For a total
row-jet budget `x`, splitting it equally minimizes their larger rate.
Thus the continuous adaptive row and column rates are `h(x)` and
`h(2-x)`, where

\[
 h(x)=\inf_{s\ge0}\{2\log Z(s)+sx\}.
\]

Concavity gives `h(x)+h(2-x)<=2h(1)`, with equality at `x=1`.
Consequently the balanced rate calculation remains

\[
 \rho=\frac{h(1)}{2\log3}
      =\frac{\log\bigl(s_*^2(1+s_*^{-1}+s_*^{-2})^4\bigr)}{4\log3},
 \qquad s_*=(1+\sqrt{13})/2,
 \tag{11}
\]

approximately `0.8203391769`. This sharpness concerns the local
Smith/Chernoff rate calculation. Finite jet costs can still improve
their polynomial prefactors, as the gauge-catalogue checker demonstrates.
An improvement using a different common type or basis remains outside
this certificate.

For asymmetric supporting weights `a,b`, the same sharp local curves
reduce the inherited calculation to

\[
 \inf_{t\ge0}\frac{t+a\log Z(t/a)+b\log Z(t/b)}{\log3}.
 \tag{12}
\]

No global minimum over `a,b` is asserted here. Additional gauges cannot
lower (12) within the stated fixed-basis, uniform-type method.

## 5. Evidence and reproduction

```sh
python3 tools/check_ternary_gauge_majorization.py \
  --output /tmp/ternary-gauge-majorization-check.json
```

The default verifier needs only the standard library and the existing
exact polynomial helpers. `--propose` regenerates the small rational
frames using NumPy/SciPy and immediately verifies them exactly.

Files:

- [Balanced-minor certificate](../evidence/ternary-uniform-gauge-majorization.json).
- [Exact verification record](../evidence/ternary-uniform-gauge-majorization-check.json).
- [Verifier and optional proposer](../tools/check_ternary_gauge_majorization.py).

The preceding numerical catalogue probe completed 3,031 optimizations
in 205.181871588 recorded seconds. It was stopped after these exact
certificates made further trials in the same scope unnecessary.
The [complete probe record](../evidence/ternary-uniform-catalogue-probe.json)
retains its trials and stopping reason. The result does not depend on
that numerical record.
