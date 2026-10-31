# Arithmetic initial forms of quantum-spin crossing paths

## Result and scope

This note extends the rational crossing generator to a quantum-spin
parameter. It gives an exact polynomial path in every dimension `q>=2`,
then applies the arithmetic initial-form theorem to every prime,
including primes dividing the deformation parameter or a normalization
denominator. The construction and its endpoint are ordinary algebraic
statements. No Yang--Baxter identity is needed for either statement.

The finite screen contains 422 parameter-and-gauge cases. All 34,944
recorded choices of two distinct parameter points and cut orientation
miss the strict Smith-certificate target. These results do not improve
the factoring exponent beyond `1/2`. They do not exclude other gauges,
extension-field parameter points, or coefficient recombinations.

The spectral construction is motivated by the quantum-spin formulas in
Bytsko, *On higher spin Uq(sl_2)-invariant R-matrices*,
[arXiv:math/0412482](https://arxiv.org/abs/math/0412482), especially the
projector description and equation (19). The formulas below specify the
normalization and parameter convention used by the implementation.

## 1. A rational family with the crossing endpoint

Fix an integer `k>=1`, put `q=k+1`, and choose an integer `t>=2`.
Write `Q=t^2` and

\[
 [a]_Q=\frac{Q^a-Q^{-a}}{Q-Q^{-1}}.
\]

On the basis `e_0,...,e_k`, define

\[
 E e_i=[i]_Qe_{i-1},\qquad
 F e_i=[k-i]_Qe_{i+1},\qquad
 H e_i=(k/2-i)e_i,
\]

where out-of-range basis vectors are zero. On the tensor square use

\[
 E_\Delta=E\otimes Q^{-H}+Q^H\otimes E,
 \qquad
 F_\Delta=F\otimes Q^{-H}+Q^H\otimes F.
\]

Every entry is rational: `Q^(k/2-i)=t^(k-2i)`. If `w=i+j`, the total
weight is `k-w`. Define the total Casimir by

\[
 C=F_\Delta E_\Delta+[H_\Delta]_Q[H_\Delta+1]_Q.
 \tag{1}
\]

Its eigenspaces are indexed by `j=0,...,k`, with eigenvalue
`lambda_j=[j]_Q[j+1]_Q` and dimension `2j+1`. The eigenvalues are
distinct because `Q>1`. Let `P_j` be the rational spectral projector.

For completeness, these assertions can be obtained without an assumed
diagonalization algorithm. In total weight `j`, set `w=k-j` and form

\[
 v_j=\sum_{i=0}^{w}c_i e_i\otimes e_{w-i},\qquad
 c_0=1,\qquad
 c_{i+1}=-c_i\frac{[w-i]_Q}{[i+1]_Q}t^{2(j+1)}.
 \tag{2}
\]

Substitution gives `E_Delta v_j=0`. The defining matrices give the usual
commutator and Casimir identities by direct multiplication. In
particular,

\[
 E_\Delta F_\Delta^r v_j
 =[r]_Q[2j-r+1]_QF_\Delta^{r-1}v_j.
 \tag{3}
\]

For `0<=r<=2j`, these vectors are nonzero and are Casimir eigenvectors
of eigenvalue `lambda_j`. Within one chain they have different weights;
different chains have different eigenvalues. Thus they are independent.
Their total number is `sum_j(2j+1)=q^2`, so they form a basis.
This also proves that each chain ends after `2j+1` vectors.

Now set

\[
 D=\prod_{\ell=1}^k(Q^{2\ell}-1),
 \qquad
 f_j(u)=
 \prod_{\ell=1}^j(Q^{2\ell}-1-u)
 \prod_{\ell=j+1}^k(Q^{2\ell}-1+Q^{2\ell}u).
 \tag{4}
\]

Let `E_cross` denote the flip matrix on the tensor square. Define

\[
 R_{k,t}(u)=E_{\rm cross}\sum_{j=0}^k\frac{f_j(u)}D P_j.
 \tag{5}
\]

This is a rational matrix polynomial of degree at most `k`. Since
`f_j(0)=D` and `sum_j P_j=I`,

\[
 R_{k,t}(0)=E_{\rm cross}.
 \tag{6}
\]

Both the Casimir and the flip preserve total weight. Hence its entries
vanish unless `x+y=u+v`; this reduces both cut computations to small
blocks. This integer-weight conservation is not an identification with
the translation-invariant binary conservation spaces studied elsewhere
in the repository.

## 2. Reduction at arbitrary primes

Apply the theorem in
[arithmetic-crossing-initial-forms.md](arithmetic-crossing-initial-forms.md)
to (5). If `c_{xyuv,d}` is a nonzero coefficient and `h,g` are rational
diagonal weights, use

\[
 \sigma=\max_{d>0,c_{xyuv,d}\ne0}
 \frac{-v_p(c_{xyuv,d})-h_x-g_y+g_u+h_v}{d}.
 \tag{7}
\]

All adjusted coefficient valuations are nonnegative. Retain the residue
of the unit part exactly when the adjusted valuation is zero. The
resulting polynomial `C(u)` over `F_p` still has `C(0)=E_cross`;
therefore `K(z)=C(1-z)` is a crossing path with `K(1)=E_cross`.

The construction takes place over the rationals before reduction.
There is no requirement that `p` avoid `t`, `D`, or the denominators of
the projectors. Their valuations contribute to (7). A common denominator
of the rational weights and slope gives the formal ramification degree,
as in the general theorem. The implementation does not construct a
ramified number field.

The gauges in the recorded screen use zero weights, ordinary factorial
valuations, quantum factorial valuations

\[
 a_i=\sum_{\ell=1}^i v_p([\ell]_Q),
\]

and quantum binomial valuations `a_k-a_i-a_(k-i)`. Both signs and equal
or opposite weights on the two label legs are included. Duplicate
weight pairs are removed within each `(q,t,p)` case. The count of 422
refers to cases, not to 422 distinct resulting tensors.

## 3. Exact controls and finite results

The source checker constructs the vectors (2), applies the lowering
operator directly, and checks every spectral coefficient on every
resulting weight block. It does not use the Casimir projectors to build
these test vectors. A nonzero determinant modulo `1000003` independently
certifies each rational basis.

The source audit covers `q=2,...,9` at `t=2,3`. Five arithmetic cases
then check the crossing endpoint, both cut spectra at one and infinity,
and the spectra through independent block-Toeplitz ranks. These controls
include `p|t`, `p` dividing a normalization factor, and characteristic
two at an odd deformation parameter. Invalid source parameters are
rejected explicitly.

The larger screen is:

| Deformation | Label sizes | Primes | Gauge cases |
| --- | --- | --- | ---: |
| `t=2` | `3,...,9` | `3,5,7,11,13` | 227 |
| `t=2` | `3,...,9` | `2` | 91 |
| `t=3` | `6,...,9` | `2,5,7` | 104 |

For each case, the screen considers both cut orientations and all
ordered pairs of distinct points in `P^1(F_p)` other than the fixed
crossing endpoint. The cost is the same Smith-spectrum expression used
by the rational arithmetic screen. A floating-point minimizer ranks
candidates. Every reported exclusion instead uses the integer
AM--GM lower certificate, compared through
`lower_bound^35 > q^99`; since `sqrt(2)<99/70`, this implies failure of
the strict target `B<q^(2sqrt(2))` for every real scale `s>=1`.

All 34,944 parameter-pair cases have that exact exclusion certificate.
The smallest displayed contraction is approximately `0.919776`, at
`q=3,t=2,p=3`. Among `q>=6`, the smallest is approximately `0.924496`,
at `q=6,t=3,p=2` with equal factorial weights. Both are above the
required `1/sqrt(2)`.

## 4. Reproduction

Run the independent audit with:

```sh
python3 tools/check_quantum_crossing.py --output /tmp/quantum-crossing-audit.json
```

The three screen commands are:

```sh
python3 tools/quantum_crossing_initial_forms.py --q 3,4,5,6,7,8,9 --t 2 --primes 3,5,7,11,13 --output /tmp/quantum-t2-unit.json
python3 tools/quantum_crossing_initial_forms.py --q 3,4,5,6,7,8,9 --t 2 --primes 2 --output /tmp/quantum-t2-nonunit.json
python3 tools/quantum_crossing_initial_forms.py --q 6,7,8,9 --t 3 --primes 2,5,7 --output /tmp/quantum-t3.json
```

The combined records are in
[`quantum-crossing-screen.json`](../evidence/quantum-crossing-screen.json).
The independent controls are in
[`quantum-crossing-audit.json`](../evidence/quantum-crossing-audit.json).
All coefficient arithmetic is exact rational or finite-field arithmetic.
