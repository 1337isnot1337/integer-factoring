# Fixed-gauge floors for up to eight pure corrections

## Theorem

Let `T_5` be the retained five-letter crossing. Use its tensor square
on four 25-dimensional physical boundaries, with the tensor-product
natural physical heights. Consider a correction

\[
 \Delta(z)=\sum_{j=1}^r \lambda_j(z)\,
     a_j(z)\otimes b_j(z)\otimes c_j(z)\otimes d_j(z),
 \qquad a_j,b_j,c_j,d_j\in\mathbb F_4(z)^{25}.
                                                        \tag{1}
\]

Each physical vector may correlate its two copy indices arbitrarily.
Assume the resulting tensor is a polynomial strict crossing after
the permitted common Laurent monomial translation. There is no
entry-floor restriction on the correction.

The fixed-natural-gauge paired-type rate satisfies

\[
 \boxed{r\le5\ \Longrightarrow\ \rho\ge\frac57,\qquad
        r\le8\ \Longrightarrow\ \rho\ge\frac{17}{24}.}   \tag{2}
\]

Both lower bounds exceed `1/sqrt(2)`. The second comparison is the
integer inequality `2*17^2=578>576=24^2`. Thus adding at most eight
physically pure rational terms cannot cross the factoring threshold
at these fixed physical heights.

The theorem concerns this reference, these heights, and coefficient
field `F_4` or `F_2`. It does not cover larger coefficient fields,
different heights, nine or more pure terms, or arbitrary tensors
in the aligned Gaussian source module. It is a limitation of this
construction and pressure calculation, not a lower bound on factoring.

## 1. Two different matrix ranks

At either local place, the two-copy crossing is a `625 x 625`
matrix. Its regular Smith frame gives image vectors of length 625.
Each such vector can also be reshaped as a `25 x 25` matrix on the
two physical boundaries of that image space.

The rank distances below refer to this **physical reshape of an
image vector**, not the rank of the entire crossing matrix. For a
pure term in (1), its image vector has physical reshape rank one.
Any linear combination of `r` such image vectors has physical
reshape rank at most `r`.

The [single-pure-term theorem](pure-two-copy-flag-floor.md) supplies
the complete one-copy profile and the corresponding constant image
flags. If `S_i` is the five-by-five reshape of the `i`th image vector
and `a_i` is its relative order, the two-copy flag through order `b`
is

\[
 F_b=\operatorname{span}_{\mathbb F_4}
     \{S_i\otimes S_j:a_i+a_j\le b\},                   \tag{3}
\]

where the row-copy coordinates and column-copy coordinates are
grouped together. The two image flags and two coimage flags agree
at every one-copy prefix in the exact retained frames. The checker
verifies all one hundred prefix comparisons.

## 2. Exact finite rank distances

The first rank frontiers in the one-copy image flags are:

| Rank allowed | First flag dimension containing such a nonzero matrix | First relative order | Projective column spaces checked |
|---:|---:|---:|---:|
| 1 | 12 | 10 | 341 |
| 2 | 6 | 6 | 5,797 |
| 3 | 3 | 3 | 5,797 |
| 4 | 2 | 1 | 341 |

For each rank `r`, the checker enumerates the complete Grassmannian
of possible `r`-dimensional column spaces in `F_4^5`, represented
by canonical reduced row-echelon annihilators. A flag matrix has
column space in that subspace exactly when the annihilator kills
every column. This is a linear system in the flag coefficients.
The checker proves every earlier flag intersection is zero and
retains a literal matrix at the first successful flag.

The two-copy flags have substantially larger rank distances:

| Relative order cutoff | Flag dimension | Projective `F_4` vectors exhausted | Minimum physical matrix rank |
|---:|---:|---:|---:|
| 4 | 10 | 349,525 | 12 |
| 5 | 14 | 89,478,485 | 12 |

These are complete finite enumerations. A projective representative
has its first nonzero coefficient equal to one, so a space of
dimension `d` has exactly `(4^d-1)/3` representatives. The native
checker computes the exact rank of every representative, stopping
an elimination early only when its rank already meets the smallest
rank previously witnessed. That shortcut cannot hide a smaller rank.
Literal rank-twelve witnesses are retained for both spaces.

In particular, the relative-order-at-most-five flag contains no
nonzero matrix of rank at most eleven. It is not enough to tensor
the one-copy rank-distance bound to obtain this fact; the displayed
enumeration proves it directly.

At relative order six, there is a rank-nine positive control:
`S_2 tensor S_2`, with zero-based index `2`. Each factor has rank
three and relative order three. This is a low-rank image vector,
not a construction of a successful nine-term four-leg correction.

## 3. The higher-jet-safe rank-r lemma

Work over a discrete valuation field `K((t))`. Normalize a reference
matrix so its minimum order is zero, and let

\[
 LMR=\operatorname{diag}(d_i),\qquad \nu(d_i)=a_i.
\]

Only the finite jet through the orders used below is needed. Both
frames are regular and invertible at zero. Write `F_b` for the span
of the constant image vectors `L(0)^(-1)e_i` with `a_i<=b`.

Suppose a correction has cut rank at most `r`, and the constant
reductions of its saturated image space meet `F_b` only at zero.
Then, for every integer `ell<=b+1`,

\[
 R_{M+\Delta}(\ell)\ge R_M(\ell),\qquad
 R_a(\ell)=\sum_i(\ell-a_i)_+.                         \tag{4}
\]

Here `R` is the rank of the coefficient Toeplitz matrix at that
precision, after both matrices have been made integral.

To prove (4), choose a regular basis for the saturated image space.
Its reductions are independent. In the reference frame, projection
of those reductions onto the coordinates with `a_i>b` is injective.
Some square minor on high-order rows is therefore a unit. Regular
row operations clear the entire correction from every other row,
leaving it supported only on those high-order rows. The same
operations change each original diagonal row only by combinations
of high-order diagonal entries, which vanish modulo `t^ell`.
All original nonzero truncated rows survive. The correction can
only add rows at this precision, so rank cannot decrease.

### Why pure terms satisfy the saturated-space hypothesis

At each cut, the image of a sum of `r` pure terms lies in the
rational span of `r` physically rank-one matrices. Every vector
in that rational span has physical matrix rank at most `r`.
Normalize a nonzero vector to have a nonzero constant coefficient.
All its `(r+1)`-minors vanish before reduction and hence after
reduction. Its constant coefficient still has rank at most `r`.

The saturated image lattice is free over the discrete valuation
ring and admits a basis with independent constant reductions.
Every vector in its reduction has the same rank-at-most-`r`
property. If the reference flag has minimum rank greater than `r`,
the two reduced spaces have zero intersection. This proves the
hypothesis used in (4), including combinations whose leading jets
cancel. Keeping only the individual leading coefficients of the
original pure factors would not justify this step.

For arbitrary negative Laurent orders, multiply both matrices by
`t^L` until integral. The flag cutoff becomes `L+b` without changing
the constant flag. Translating the Toeplitz inequality back gives

\[
 R_{\beta}(x)\ge R_{\alpha}(x)\qquad(x\le b+1)           \tag{5}
\]

for the original relative order lists. Integer orders allow linear
interpolation, and below the common lower endpoint both functions
are zero. No bound on the negative orders is assumed.

## 4. From rank distance to a complete partition bound

Rank-`r` valuation interlacing gives

\[
 N_\beta(x)-N_\alpha(x)\ge-r,
\]

where `N` counts orders at most `x`. It follows by intersecting
the valuation min-max subspaces with the kernel of the correction,
whose codimension is at most `r`.

Combine this with (5), setting `c=b+1`. The difference
`D=R_beta-R_alpha` satisfies `D(x)>=0` for `x<=c` and
`D(x)>=-r(x-c)` for `x>c`. Therefore, for every `tau>0`,

\[
 \begin{aligned}
 \sum_i e^{-\tau\beta_i}-\sum_i e^{-\tau\alpha_i}
 &=\tau^2\int_{-\infty}^{\infty}e^{-\tau x}D(x)\,dx\\
 &\ge-r e^{-c\tau}.                                    \tag{6}
 \end{aligned}
\]

Let `P(s)` be the complete relative one-copy partition polynomial
from the single-term note, and let

\[
 Q(s)=1+s^2+s^6+s^8+s^{10}.
\]

Each two-copy reference local partition is `P(s)^2`. The opposite
reference minima cancel in the product of the two places. The
rank-twelve flag through order five consequently gives

\[
 B_{\rm new}(s)\ge\bigl(P(s)^2-rs^6\bigr)^2
 \quad\text{for }r\le11,
                                                        \tag{7}
\]

whenever the displayed lower partition is positive. The cases
`r=5` and `r=8` used below are positive throughout `[0,1]`.

## 5. Exact uniform pressure certificates

For a proposed floor `rho=m/n`, put

\[
 A(s)=P(s)^2-rs^6.
\]

The required lower-pressure inequality reduces to

\[
 A(t^{2m})^n\ge Q(t^n)^{4m},\qquad 0\le t\le1.          \tag{8}
\]

The two principal certificates are:

| Pure terms `r` | Rate floor `m/n` | Exact dyadic intervals on `[1/2,1]` | Largest denominator |
|---:|---:|---:|---:|
| 5 | 5/7 | 149 | `2^13` |
| 8 | 17/24 | 260 | `2^16` |

For `r=5`, `A` has nonnegative coefficients. For `r=8`, its sole
negative coefficient is `-3s^6`; the preceding coefficient is
`4s^5`. Their sum and derivative are nonnegative on `[0,1]`,
because `4-3s>=1` and `20-18s>=2`. In both cases `A` is increasing
and `A(s)>=1+2s`.

On `[0,1/2]`, write `u=t^(2n)`. The exact integer comparison

\[
 16m(2^{2n}+4)^{4m-1}
  <(16m+1)(2^{2n})^{4m-1}
\]

bounds the derivative of `(1+4u)^(4m)` by `16m+1`.
Consequently

\[
 A(t^{2m})^n\ge1+2n t^{2m},\qquad
 Q(t^n)^{4m}\le1+(16m+1)t^{2n}.
\]

The checker verifies `2n*2^(2(n-m))>16m+1`, proving (8)
on this interval. On each retained dyadic interval `[a,b]`, it
checks `A(a^(2m))^n>Q(b^n)^(4m)` with integer numerators and
denominators. Monotonicity then proves (8) on the entire interval.
The exact interval cover has no gaps or overlaps.

The same record also retains the earlier two-term cutoff-six and
three-term cutoff-five `5/7` certificates. No floating-point
optimizer is used in any of these uniform inequalities.

## 6. Reproduction and controls

```sh
python3 tools/check_low_rank_two_copy_flags.py
```

The [Python checker](../tools/check_low_rank_two_copy_flags.py)
reconstructs the reference frames, checks all constant prefixes,
enumerates the one-copy Grassmannians, independently verifies
small complete rank-space controls, and expands every retained
rank witness as a literal matrix. A binary two-by-two pencil
whose rank drops only over `F_4` tests extension-field arithmetic;
a dependent-basis example tests detection of rank zero.

The [native enumerator](../tools/low_rank_copy_flag_enumeration.cpp)
performs the two complete projective enumerations. The
[certificate](../evidence/low-rank-two-copy-flag-floors.json)
contains its source hash, exact input hashes, all basis pairs,
minimum-rank witnesses, the one-copy rank frontiers, and every
dyadic pressure interval. The 89-million-point enumeration took
about three minutes in the recorded run.

The factoring exponent claim is unchanged. The theorem directs
future correction searches toward different heights, larger fields,
or more than eight pure terms.

This research and its implementation used substantial OpenAI assistant support.
