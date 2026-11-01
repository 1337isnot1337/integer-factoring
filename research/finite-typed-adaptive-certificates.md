# Exact finite certificates for weighted adaptive crossings

## Result and scope

One finite typed block can now be checked through an exact integer gate.
The [weighted adaptive catalogue construction](weighted-adaptive-gauge-catalogues.md)
supplies a rank bound for each reconstructed piece. The
[finite-block amplification theorem](finite-block-crossing-amplification.md)
then turns a sufficiently strong finite rank table into an effective
crossing family at every external label width.

The implementation combines these two steps. It accepts unequal paired
physical alphabets, recomputes every type mean from integral ramified
heights, includes every jet and endpoint piece, caps ranks by the exact
type-class dimension, and tests the amplification inequality using
integers. Its stronger rank-table test uses the sum of square roots of
the individual rank products.

No finite block passing the factoring threshold is supplied here. The
real tensor controls below all fail that threshold. A separate artificial
rank table tests the arithmetic distinction between the coarse and
stronger gates; it is not asserted to be a crossing decomposition.

## 1. Exact local input

Fix a strict polynomial crossing of degree at most `D`, with leg
alphabets `(X,Y,Y,X)` of sizes `q_X,q_Y`. Fix one pair of physical type
counts `n_X,n_Y`, both summing to the positive tensor length `m`. Its
external type-class sizes and cut dimension are

```math
d_X=\binom{m}{n_X},\qquad d_Y=\binom{m}{n_Y},
\qquad d=d_Xd_Y.                                      \tag{1}
```

The implementation requires `d_X,d_Y>=2` when applying amplification.
The local rank estimate itself does not need positive entropy.

All geometric places are rational over one fixed coefficient extension.
The row and column place sets are disjoint and avoid the strict crossing
endpoint. A closed point of residue degree `e` supplies `e` conjugate
geometric places. Selected conjugates appear as separate labelled places
with identical profiles. A profile has exactly `q_X q_Y` Smith channels;
an unnormalized restriction-of-scalars profile with every channel
repeated `e` times is not the required input.

For each gauge, choose an integer ramification `L>=1` and integral
ramified physical heights `h_x,h_y,h_u,h_v`. The complete ramified Smith
orders are `nu_1,...,nu_(q_X q_Y)`, including any negative orders. Compute

```math
A=\frac1m\left(
 \sum_a n_X(a)(h_x(a)+h_v(a))
 +\sum_b n_Y(b)(h_y(b)+h_u(b))\right).                \tag{2}
```

Thus `mA` is an integer without an additional rounding convention. The
original order-`j` Hasse jet is the coefficient of order `Lj+mA` after
this local gauge. The global degree budget remains `Dm`.

The tensor, physical bases, types, geometric-place distinctness, and
complete local profiles must be certified independently. The wrapper
does not establish them from numerical optimizer records. Its input
parser rejects repeated place labels, but labels alone are not an
algebraic proof that two places are distinct. The command-line output
retains the supplied local provenance and input digest for this reason.

## 2. Every finite jet and reconstructed piece

For a channel word `w` of length `m`, write
`nu(w)=sum_t nu_(w_t)`. The exact local upper bound is

```math
c_g(j)=\sum_w\max\{Lj+mA-\nu(w)+1,0\}.                \tag{3}
```

The implementation evaluates (3) by the exact integer convolution from
the weighted catalogue module. Negative orders are translated together
with the mean. At each geometric place it takes the minimum of (3)
over all supplied gauges. This minimum is nondecreasing in `j`.

The greedy merge selects the cheapest `Dm+1` available initial jets on
each side. Let `R_k` be the sum of the first `k+1` row costs, capped by
`d`, and let `C_k` be the analogous column bound. The actual opposite-flag
pieces have bounds

```math
\operatorname{rank}_R A_k\le R_k,\qquad
\operatorname{rank}_C A_k\le C_{Dm-k},
\qquad 0\le k\le Dm.                                 \tag{4}
```

All indices, including both endpoints and the first positive jet after
any zero prefix, are retained. If either bound in (4) is zero, the piece
is zero and may be omitted from the amplification sum. If a nonzero
bound happens to correspond to an actually zero piece, retaining it is
conservative.

This calculation avoids a numerical convention for the value of an
inverse capacity curve at a jump. It does not replace the local
Smith-profile proofs or the opposite-flag reconstruction theorem.

## 3. The stronger rank-table gate

Let `r_i,c_i` be the remaining positive bounds in (4). Define

```math
S=\sum_i\sqrt{r_i c_i}.                                \tag{5}
```

The refined finite-block amplification theorem gives

```math
\rho\le\frac{\log S}{\log d}.                          \tag{6}
```

The key refinement is the multinomial theorem. For a composition `v`
of the tensor-power pieces, its multinomial factor and product of local
square-root rank bounds form one term of the power of their sum.
Consequently that term is at most the entire power. Pair-swapping the
two physical species squares (5), which gives exactly (6) after the
external entropy normalization.

For an exact test choose integers `a,b>0` with `a^2<2b^2`. At binary
precision `p`, set

```math
N=\sum_i\left\lceil 2^p\sqrt{r_i c_i}\right\rceil,
\qquad Q=2^p.                                        \tag{7}
```

Integer square roots compute every summand. Since `S<=N/Q`, the strict
integer comparison

```math
N^{2b}<Q^{2b}d^a                                     \tag{8}
```

implies `rho<a/(2b)<1/sqrt(2)`. There are no logarithmic rounding errors
in this gate. The default rational exponent is `a/b=140/99`;
`140^2=19600<19602=2*99^2`.

For comparison the tool also records `M=max_i r_i c_i` and the number
`L` of retained pieces. The simpler test

```math
(M L^2)^b<d^a                                        \tag{9}
```

uses `S<=L sqrt(M)`. The sharper test can pass when (9) fails. Both tests
are only sufficient: failure is not an impossibility statement.

## 4. Implementation and exact controls

The implementation is
[`finite_typed_adaptive_certificate.py`](../tools/finite_typed_adaptive_certificate.py).
Its `finite_certificate` function accepts the already certified
`GaugeProfile` catalogues. `gauge_from_heights` constructs their mean
shifts directly from ramified letter heights and integer type counts.

The independent checker recomputes complete profiles by block-Toeplitz
ranks, constructs full typed powers over `F_3`, reconstructs the actual
opposite-flag pieces, and verifies their endpoint sum and both ranks.
It covers two equal-alphabet types and both orientations of unequal
two-by-three physical alphabets. Its exact totals are:

| Control | Count |
|---|---:|
| Complete local Smith profiles | 24 |
| Local Smith channels | 180 |
| Actual typed Hasse-jet ranks | 78 |
| Actual reconstructed-piece cut ranks | 52 |
| Strict endpoint entries | 1,960 |
| Agreements with the older equal-alphabet finite routine | 2 |
| Integer square-root upper-bound checks | 4,004 |
| Invalid-input rejections | 4 |

The artificial table has positive rank products `100,1`. Its exact
square-root sum is `11`; with `d=40`, the stronger integer gate passes
and the coarse `M L^2` test fails. This checks only arithmetic. Every
fully reconstructed tensor block in this checker fails the threshold.

```sh
python3 tools/check_finite_typed_adaptive_certificate.py \
  --output evidence/finite-typed-adaptive-certificate-controls.json
```

The checker and replay record are
[`check_finite_typed_adaptive_certificate.py`](../tools/check_finite_typed_adaptive_certificate.py)
and
[`finite-typed-adaptive-certificate-controls.json`](../evidence/finite-typed-adaptive-certificate-controls.json).
