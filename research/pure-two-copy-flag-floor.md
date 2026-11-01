# A 5/7 floor for pure corrections to the two-copy reference

## Theorem and scope

Let `T_5(z)` be the retained degree-18 five-letter strict crossing.
Take its tensor square, with each physical label now an ordered pair
of five-letter labels. Use the tensor-product natural physical heights.
Let

\[
 \Delta(z)=\lambda(z)\,a_X(z)\otimes b_Y(z)
                    \otimes c_U(z)\otimes d_V(z),       \tag{1}
\]

where all factors have coefficients in `F_4(z)`. Each physical vector
has dimension 25 and may correlate its two copy indices arbitrarily.
Assume the resulting crossing is polynomial after the permitted common
Laurent monomial translation and has the strict endpoint.

**The fixed-natural-gauge paired-type pressure of
`T_5 tensor T_5 + Delta` has rate parameter**

\[
 \boxed{\rho\ge\frac57>\frac1{\sqrt2}.}                \tag{2}
\]

Consequently this family cannot cross the threshold needed for an
improved factoring exponent using these physical heights. The theorem
covers arbitrary Laurent orders of the correction; preserving the
reference entry floors is not required. The correction is physically
pure across the four 25-dimensional boundaries, even though each vector
may be entangled between its two copy indices.

The statement is specific to one pure correction, this reference,
these physical heights, and coefficient field `F_4` or its subfield
`F_2`. It does not address several correlated pure corrections,
different physical heights, larger coefficient fields, or arbitrary
perturbations in the aligned source module. It is a pressure-method
obstruction, not a lower bound on integer factoring.

## 1. Complete one-copy data

At the two local places, the exact one-copy profiles are

\[
 \begin{aligned}
 \alpha_\infty={}&(-8,-7,-5,-4,-3,-2,-1,0,0,1,1,2,2,3,3,\
                  &4,4,4,5,5,6,7,7,8,12),\\
 \alpha_0={}&\alpha_\infty+16.
 \end{aligned}                                         \tag{3}
\]

After subtracting the minimum, both lists have partition polynomial

\[
 \begin{aligned}
 P(s)={}&1+s+s^3+s^4+s^5+s^6+s^7
       +2s^8+2s^9+2s^{10}+2s^{11}\\
       &+3s^{12}+2s^{13}+s^{14}+2s^{15}+s^{16}+s^{20}.
 \end{aligned}                                         \tag{4}
\]

The five physical energies are `(0,10,6,2,8)`, so their partition
polynomial is

\[
 Q(s)=1+s^2+s^6+s^8+s^{10}.                             \tag{5}
\]

For two copies, the local profiles consist of all pairwise sums.
Their minima are `-16` and `16`; both normalized partitions are
`P(s)^2`. The physical alphabet partitions are `Q(s)^2` on both pairs.
This includes all physical types on the 25-letter alphabets.

## 2. The constant image flag contains no pure vector

Use a regular frame at either place,

\[
 L(t)M(t)R(t)=\operatorname{diag}(d_i(t))
               \pmod {t^{64}},
 \qquad \nu(d_i)=\alpha_i+\text{clearing shift}.
\]

Both `L(0)` and `R(0)` are invertible. The physical image vectors
are the **columns of `L(0)^(-1)`**. They are not the rows of `L(0)`.
The corresponding physical coimage vectors are the rows of
`R(0)^(-1)`.

At relative orders at most nine, there are eleven channels. The two
image flags and the two coimage flags have the same eleven-dimensional
constant span `U` after identifying the two physical coordinates of
each cut. Reshape its vectors as five-by-five matrices.

**Exact finite fact.** `U` contains no nonzero rank-one matrix over
`F_4`.

The checker proves this exhaustively without enumerating all matrices
in `U`. It constructs fourteen independent binary annihilating forms
for `U`. For each of the 341 projective nonzero vectors `a` in `F_4^5`,
the condition `a b^T in U` is a homogeneous system of fourteen equations
in the five coordinates of `b`. Every system has rank five. The retained
record specifies five independent constraint rows for every `a`.
Every nonzero left factor has a unique projective representative, so
these checks exclude every nonzero rank-one matrix.

This cutoff is real: the twelve-dimensional flag at relative order
ten does contain a rank-one matrix. A literal pair of factors is kept
as a positive control. The minimum-order coefficient has an even
simpler form. As a cut matrix it is `S S^T`, where

\[
 S=17246593
\]

encodes a 25-entry binary vector. Its five-by-five reshape has rank
five, with support

\[
 (0,0),(1,2),(1,3),(2,1),(2,3),(3,1),(3,2),(3,3),(4,4).
\]

Thus a physically pure correction cannot be proportional to the
complete minimum-order coefficient.

### Tensoring does not create a pure vector

If a matrix space `U` contains no nonzero rank-one matrix over a field,
then `U tensor U`, viewed with the two row spaces grouped together and
the two column spaces grouped together, also contains none.

To see this, suppose a nonzero element were `a b^T` in the grouped
spaces. Contract the first row-copy coordinate and first column-copy
coordinate independently. Some such choice gives nonzero slices of
both `a` and `b`; the contraction is their nonzero outer product. It
also belongs to `U`, a contradiction.

The relative-order-at-most-nine two-copy image flag lies inside
`U tensor U`: a pair of nonnegative one-copy orders whose sum is at
most nine has both orders at most nine. Therefore the two-copy low
image flags contain no physically rank-one vector, including vectors
with arbitrary correlations between the copies inside a boundary.

## 3. A necessary condition that retains every higher jet

Let `M` and `M'` be integral matrices over a discrete valuation field,
with a regular diagonalizing frame for `M` through order `b`. Put

\[
 F_b(M)=\operatorname{span}
  \{L(0)^{-1}e_i:\nu(d_i)\le b\}.
\]

Suppose `M'=M+lambda u v^T`, and normalize `u`
so that it is regular with nonzero leading coefficient `u(0)`.

**Lemma.** If `u(0)` is outside `F_b(M)`, then at every positive
integer precision `ell<=b+1`, the truncated coefficient Toeplitz
matrix of `M'` has rank at least that of `M`.

**Proof.** In the regular frame, put `p=L u`. Some coordinate `p_j`
is a unit with `nu(d_j)>b`. Use this coordinate to subtract
`p_i/p_j` times row `j` from each other row. This moves the entire
rank-one correction into row `j`. The same operation changes each
other diagonal row only by a multiple of `d_j`, invisible modulo
`t^ell`. The original row `j` is itself zero at this precision.
All original nonzero truncated rows therefore remain, and the update
only supplies an additional row. Rank cannot decrease. These are
regular invertible operations, so the conclusion transfers back
through the frame. No higher coefficient of `u` or `v` was discarded.

For an order list `alpha`, the Toeplitz rank at precision `ell` is

\[
 R_\alpha(\ell)=\sum_i(\ell-\alpha_i)_+.
\]

At either two-copy place, the leading coefficient of the local left
factor in (1) is physically rank one. The flag exclusion and the lemma
therefore imply, after normalizing the reference minimum to zero,

\[
 R_\beta(x)\ge R_\alpha(x)\qquad(x\le10).               \tag{6}
\]

Here `beta` is the updated profile. For initially integral matrices,
the integer inequalities imply the real ones by linear interpolation
between integer orders.

### Negative Laurent orders are allowed

If the update has negative orders, multiply both matrices by `t^L`
until both are integral. The reference flag cutoff becomes `L+9`,
and its constant image space stays unchanged. The preceding row
operation gives the rank inequality through precision `L+10`.
Translate by `x=ell-L`. Below `-L` both rank functions vanish;
interpolation proves (6) for every real `x<=10`, including negative
orders. This is why an entry-floor assumption is unnecessary.

Inequality (6) excludes a lexicographic improvement beginning within
these low orders. It permits a later improvement accompanied by an
earlier deterioration. The next step accounts for those trades.

## 4. Rank-one interlacing controls the remaining tail

Rank-one updates obey the usual valuation interlacing inequalities

\[
 \beta_i\le\alpha_{i+1},\qquad
 \alpha_i\le\beta_{i+1}.
\]

A short proof uses the valuation-norm min-max description of the
Smith orders. Intersect a subspace attaining one order with the
kernel of the update; its codimension increases by at most one, and
the two matrices agree on that intersection. Applying this argument
in both directions gives the displayed inequalities. The min-max
description follows directly for a diagonal matrix, and regular
invertible frames preserve the valuation norm.

Consequently the counting functions satisfy
`N_beta(x)-N_alpha(x)>=-1`. Combine this with (6), and define
`D(x)=R_beta(x)-R_alpha(x)`. Then

\[
 D(x)\ge0\quad(x\le10),\qquad
 D(x)\ge-(x-10)\quad(x>10).
\]

For every `tau>0`, integration over the entire real line gives

\[
 \begin{aligned}
 \sum_i e^{-\tau\beta_i}-\sum_i e^{-\tau\alpha_i}
 &=\tau^2\int_{-\infty}^{\infty}e^{-\tau x}D(x)\,dx\\
 &\ge-e^{-10\tau}.                                    \tag{7}
 \end{aligned}
\]

The integral is finite because the order lists are finite and have
finite lower endpoints. Thus the complete two-copy local partition
can save at most one channel of relative order ten:

\[
 Z_{\rm new}(s)\ge P(s)^2-s^{10},\qquad 0<s<1.          \tag{8}
\]

The scalar shifts from the two reference minima `-16` and `16`
cancel in the product of the two local partitions. Applying (8) at
both places gives the exact pressure lower bound

\[
 B_{\rm new}(s)\ge\bigl(P(s)^2-s^{10}\bigr)^2.          \tag{9}
\]

## 5. A uniform integer certificate for rho>=5/7

At rate `rho=5/7`, the paired-type pressure parameter is `r=2rho=10/7`.
Using the physical partition `Q(s)^2` on each pair, the required
lower pressure inequality is

\[
 P(s)^2-s^{10}\ge Q(s^{7/10})^{20/7}.                  \tag{10}
\]

Put `s=t^10` and define

\[
 A(t)=P(t^{10})^2-t^{100}.
\]

The coefficients of `P(s)^2-s^10` are nonnegative. Therefore `A`
is positive and increasing on `[0,1]`. Raising (10) to the seventh
power reduces the claim to a single integer-polynomial inequality,

\[
 \boxed{A(t)^7\ge Q(t^7)^{20}\qquad(0\le t\le1).}       \tag{11}
\]

On `0<=t<=1/2`, use `A(t)>=1+2t^10` and
`Q(t^7)<=1+4t^14`. The integer inequality

\[
 80\cdot4097^{19}<81\cdot4096^{19}
\]

bounds the derivative of `(1+4u)^20` by 81 for
`0<=u<=2^-14`. Hence

\[
 A(t)^7\ge1+14t^{10},\qquad
 Q(t^7)^{20}\le1+81t^{14}.
\]

Since `81t^4<=81/16<14`, this proves (11), strictly when `t>0`.

The remaining interval `[1/2,1]` has an exact cover by sixty dyadic
intervals, with maximum denominator `2^10`. For every leaf `[a,b]`,
the checker verifies the integer comparison

\[
 A(a)^7>Q(b^7)^{20}.
\]

Monotonicity proves (11) throughout that leaf. The checker verifies
the cover has no gaps or overlaps and performs all endpoint
comparisons with integer numerators and denominators. There is no
floating-point minimization in this certificate.

Inequality (11), together with (9), proves (2). The threshold comparison
is also exact: `2*5^2=50>49=7^2`.

## 6. Reproduction

```sh
python3 tools/check_pure_two_copy_flag_floor.py
```

The [checker](../tools/check_pure_two_copy_flag_floor.py) and
[complete finite certificate](../evidence/pure-two-copy-flag-floor.json)
retain the exact reference digest, complete one-copy profiles, four
constant flag comparisons, all 341 projective rank tests, the
fourteen annihilating forms, a rank-one positive control at the next
flag, the tensor-square inclusion, and every dyadic pressure interval.
The regular-frame identities are checked on the literal polynomial
matrices. A proposed stronger `3/4` rate floor is rejected at an exact
rational evaluation point.

This theorem and certificate used substantial OpenAI assistant support.
