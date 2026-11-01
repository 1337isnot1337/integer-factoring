# Mismatched mode weights: a noncommutative Weyl family and a spectral obstruction

This note tests characteristic two with mode weights `1,3,9,27`, rather than
powers of the characteristic.  The family is an exact crossing construction.
The matched-weight additive law does not apply.  The new finite result below
excludes **every polynomial recombination of its constant matrix at energy 63**,
including all spectral projectors.  It does not exclude the algebra generated
by all coefficients, their opposite orientations, or their transposes.
No improved factoring bound follows from this experiment.

## 1. The curve does not require characteristic-power weights

For arbitrary positive integer weights `w_j`, let

\[
\mathcal L_h=\{\alpha\in\mathbb Z_{\ge0}^{r}:\sum_jw_j\alpha_j=h\},
\qquad d=|\mathcal L_h|.
\]

Over a field of characteristic two, define

\[
P_h(z)(\alpha,y,u,\beta)=
 [\alpha+u=\beta+y]\,[t^\alpha]
 \prod_j(1+t_j)^{\beta_j}
       (1+(1+z)^{w_j}t_j)^{y_j}.                 \tag{1}
\]

At `z=1`, only `\alpha\le\beta` can contribute.  Equal weighted energies and
positive weights force `\alpha=\beta`, then `y=u`.  Thus `P_h(1)=E` exactly.
Its degree is at most `h`.

Row multiplication means the matrix with rows `(alpha,y)` and columns
`(beta,u)`.  The row blocks have fixed integer difference `alpha-y`.
The other cut has rows `(alpha,u)` and columns `(y,beta)`, with fixed sum
`alpha+u`.  These block decompositions use only the conservation condition
in (1); they remain valid for arbitrary weights.

For `w_j=3^j`, the first useful small shells are:

| Energy | Labels | Independent coefficients | Raw coefficient rank sum | Target `d^sqrt(2)` |
|---:|---:|---:|---:|---:|
| 57 | 105 | 58 | 118078.626 | 721.761 |
| 63 | 132 | 64 | 215867.266 | 997.572 |

These are exact binary ranks with decimal square roots only in the final sum.
They are failures of the raw coefficient decomposition, not universal bounds.

The normalized curve `B(z)=P(1+z)` has coefficients `B_i=C_iD_i`, but its
global Gauss coefficients no longer collapse to the characteristic-power
formula.  Noncommutativity already occurs at energy 3, with labels `(3,0)` and
`(0,1)`.  On the difference-zero block,

\[
B_1=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
B_3=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\]

so `B_1B_3\ne B_3B_1`.  This distinguishes the family from the matched binary
additive algebra.

## 2. A general nilpotent-kernel detector

Let `A` be a row matrix supported on the same integer-difference blocks as
the crossing, with the crossing matrix ordered as `I`.  Suppose:

1. `ker(A^q)` has dimension `n` for some positive integer `q`;
2. exactly zero rows of `A` include row pairs `(alpha,y)` in `m` distinct
   integer-sum classes `alpha+y`.

For every polynomial `f` over any extension field with `f(0)\ne0`,

\[
\operatorname{rank}_{\rm row} f(A)\ge n,\qquad
\operatorname{rank}_{\rm col} f(A)\ge m.       \tag{2}
\]

For the first inequality, `A` is nilpotent on `ker(A^q)`.  A polynomial with
nonzero constant coefficient is invertible on that invariant subspace, by
the finite geometric series for `f(0)I` plus a nilpotent operator.

For the second, a zero row of `A` is zero in every positive power.  That row
of `f(A)` is therefore `f(0)` times the corresponding row of `I`.  Its diagonal
entry becomes a nonzero entry of the other flattening, in the sum block
`alpha+y`.  Selecting one such entry per distinct sum block gives at least
one rank in each block.  Row products preserve integer-difference support,
so the other flattening remains a direct sum over those sum blocks.

The functional `f(A)\mapsto f(0)` is intrinsic whenever there is a zero row:
that diagonal entry recovers it.  Consequently every finite decomposition
`I=\sum_i f_i(A)` has a summand with `f_i(0)\ne0`, and that same summand
satisfies both inequalities in (2).  Nothing here assumes a perfect or finite
coefficient field.

## 3. Exact certificate at energy 63

For weights `(1,3,9,27)`, energy 63, and `A=P_63(0)`, exact computation gives

\[
\begin{aligned}
d&=132, & d^2&=17424,\\
\operatorname{rank}A&=11924,
&\operatorname{rank}A^{32}=\operatorname{rank}A^{64}&=7280.
\end{aligned}
\]

There are 2531 zero rows of `A`, occupying 496 distinct sum blocks.  Hence
`n=17424-7280=10144`, and every polynomial decomposition of the crossing
contains a piece with

\[
r\ge10144,\qquad c\ge496,\qquad
rc\ge5031424>132^3>132^{2\sqrt2}.              \tag{3}
\]

Thus no such decomposition can meet the finite-seed criterion
`\sum_i\sqrt{r_i c_i}<132^{\sqrt2}`.  The result includes arbitrary polynomial
degrees and all spectral-projector recombinations, over every field of
characteristic two.

For orientation only, the full minimal polynomial of `A` has degree 215.
Its irreducible primary degrees and multiplicities are
`(1,19),(1,18),(2,4),(3,2),(3,2),(4,2),(4,2),(6,2),(8,2),(57,2)`.
The ten primary projectors sum to `I`, but their actual rank sum is about
23338.842.  The certificate (3) is stronger than checking those projectors
and does not depend on polynomial factorization.

## 4. Replay and limits

Run:

```sh
python tools/check_mismatched_weyl_spectrum.py
```

The standard-library driver compiles `explore_mismatched_weyl.cpp`, regenerates
the energy-63 block matrices, checks all 632732 supported polynomial entries
at the crossing endpoint, computes the three power ranks, and reconstructs
the 496 protected sum blocks.  It checks the integer inequality in (3) and
reports source, dump, and witness hashes.  `--prefix PATH` accepts existing
`PATH.row.bin` dumps for a quick replay.

The exploratory `explore_mismatched_spectrum.py --projectors` routine also
factors the exact minimal polynomial by binary Berlekamp elimination and
evaluates the primary projectors.  That larger exploratory calculation is
not needed for the obstruction proof or its replay gate.

The full coefficient algebra is already noncommutative, and the proof above
does not constrain words using two different coefficients.  Those words,
opposite orientations, and the analogous mismatch of characteristic three
with binary weights remain distinct constructive possibilities.
