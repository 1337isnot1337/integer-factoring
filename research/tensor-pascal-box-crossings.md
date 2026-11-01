# Tensor-product Pascal boxes with summed mode energies

## Construction and search result

There is a constant-flag crossing on a Cartesian box of modes whose
energy is the sum of its coordinates. It has explicit complete Smith
profiles at both crossing cuts. This changes the physical energy
grading from the single interval in
[the one-species construction](schur-one-species-interval.md).

The new grading does not give a free improvement. For a box with two
coordinates, the displayed centered gauge has pressure at least that
of the one-coordinate interval of the same side length, at a rescaled
evaluation point. The interval's all-scale obstruction therefore
excludes a factoring score above `1/2` for this box and gauge.

The smallest two-coordinate box also admits nontrivial coherent
physical restrictions inside its repeated energy levels. An exact
finite screen of 150 such restrictions has pressure strictly above
`3/4` at every positive scale, at the inherited gauges. These cases
include all 30 ways to retract one of its five two-dimensional
occupation-and-energy blocks to a line. The remaining 120 cases are
a declared deterministic sample with every block retracted; they do
not exhaust that larger family.

The constructor and both spectrum formulas are ordinary proofs for
every permitted box size. The stochastic comparison concerns the
centered two-coordinate gauge. The coherent calculation concerns
the stated finite families and heights. No improved factoring
exponent is claimed.

## 1. Polynomial kernel

Work in characteristic two. Let `k` and `M` be powers of two, put

\[
 L=kM,\qquad n=M^k,\qquad
 \mathcal A=\{0,\ldots,M-1\}^k,
\]

and index the `n` modes by `a in A`. Write

\[
 w(a)=\sum_{j=1}^k a_j,\qquad
 E(a)=2w(a)+k,\qquad T=M-1.
\]

The energies satisfy `k<=E(a)<=2L-k`. Define

\[
 \begin{aligned}
 c&=1+u,\quad z=u^2,\quad q=1+z,\\
 P_{ab}&=\prod_{j=1}^k[b_j\subseteq a_j],\\
 Q_{aa}&=q^{w(a)},\\
 H_{ab}&=z^{kT-\sum_j(a_j\mathbin{|}b_j)}.
 \end{aligned}                                                     \tag{1}
\]

Here inclusion refers to binary digits. The Pascal matrix is the
ordinary tensor product `P_M` to the `k`th power and satisfies `P^2=I`.
Use the polynomial matrix

\[
 \mathcal D(u)=
 \begin{pmatrix}
 c^k QP^{\mathsf T}&u^{2L}P\\
 I&c^k H
 \end{pmatrix}.                                                    \tag{2}
\]

Rows have boundary labels `(V,Y)` and columns have labels `(B,X)`.
Its minors form the one-species Gaussian physical tensor.

The one-coordinate Pascal identity, multiplied over all coordinates,
gives

\[
 H=P^{\mathsf T}\mathrm{diag}(q^{kT-w(a)})P.
\]

Consequently the product of the upper-left and lower-right blocks
is `c^(2L)P`. Since `L` is a power of two,

\[
 \det\mathcal D
 =\det\bigl((u^{2L}+c^{2L})P\bigr)=1.                \tag{3}
\]

At `u=1`, both diagonal blocks vanish and

\[
 \mathcal D(1)=\begin{pmatrix}0&P\\I&0\end{pmatrix}.  \tag{4}
\]

This is an instance of the
[constant-flag constructor](constant-flag-crossings.md), with
exponents `E(a)` and outer scale `L`. The strict inequalities
`0<E(a)<2L` hold even when `M=1`.

A nonzero term of a Pascal minor sends each coordinate to a binary
superset. Every strict inclusion raises `w(a)`. On tensor powers,
prescribing a common physical word type on the four boundaries fixes
the same total mode energy. Every surviving endpoint term must
therefore use equality in every inclusion. The diagonal entries are
one, so (4) gives the exact crossing on each such common type.

## 2. Both complete spectra

In boundary order `(B,X,V,Y)`, assign the zero heights

\[
 h_0(a)=(0,E(a)-2L,0,E(a)).                         \tag{5}
\]

The upper-left internal block is regular at zero. A lower-right
entry has weighted valuation

\[
 E(a)+E(b)-2L+2\left(kT-\sum_j(a_j\mathbin{|}b_j)\right)
 =2\sum_j(a_j\mathbin{\&}b_j)\ge0.                 \tag{6}
\]

The two cross blocks are a constant unit matrix times diagonal
monomials with exponents `E(a)`. Thus the zero row-cut profile has
two copies of every positive energy `E(a)`.

At infinity, choose an integer center `C` and set

\[
 h_\infty(a)=(E(a)-C,2L-C,E(a)-C,-C).               \tag{7}
\]

Rational centers can instead be implemented after a common
ramification. The infinity cut separates all plus modes from all
minus modes, so its complete Fock profile is the exterior compound
of the following one-particle profile:

\[
 \boxed{\quad E(a)-2C,\qquad E(a)+2L-2C
                \quad(a\in\mathcal A).\quad}       \tag{8}
\]

To prove (8), write `xi=1/u`, `t=xi^2`, and now put

\[
 D_{aa}=t^{w(a)},\quad Q_{aa}=(1+t)^{w(a)},\quad
 O_{ab}=t^{\sum_j(a_j\mathbin{|}b_j)}.
\]

Remove the common factor `xi^(k-2C)`. Multiplication of the first
block of rows by `(1+xi)^(-k)` and of the second block of columns by
`(1+xi)^k` gives

\[
 F(t)=\begin{pmatrix}QP^{\mathsf T}D&DP\\
                      D&(1+t)^kO\end{pmatrix}.       \tag{9}
\]

The tensor-product Pascal identity is

\[
 QP^{\mathsf T}O=(1+t)^{k(M-1)}DP.                  \tag{10}
\]

Every row `a` of `O` is divisible by `t^(w(a))`. The regular row
and column operations used for the one-coordinate Schur reduction
therefore turn (9) into

\[
 \mathrm{diag}(D,S),\qquad
 S=\bigl(1+(1+t)^{kM}\bigr)DP=t^LDP.                \tag{11}
\]

Its orders in `t` are `w(a)` and `L+w(a)`. Returning to `xi` and
restoring the common factor proves (8).

The total restored physical height is

\[
 G(a)=4(E(a)-C).                                    \tag{12}
\]

At the centered choice `C=L`, the complement
`a_j -> M-1-a_j` sends `E(a)` to `2L-E(a)`. Thus (8) is precisely
the signed energy list `-E(a),E(a)`. Since `sum_a E(a)=nL`, the
product of the two exterior profile sums is

\[
 \mathcal B_{k,M}(s)
 =s^{nL}\prod_{a\in\mathcal A}(1+s^{-E(a)})^4.       \tag{13}
\]

For `R=4rho`, the exact common-type pressure criterion uses

\[
 \mathcal Z_{k,M,R}(s)
 =\prod_a(1+s^{-4(E(a)-L)/R}),\qquad
 \mathcal B_{k,M}(s)<\mathcal Z_{k,M,R}(s)^R.         \tag{14}
\]

The local tensor and all profiles here belong to one common
polynomial kernel. No specialization in other parameter variables
is used.

## 3. A two-coordinate stochastic comparison

Let `M>=2`, put `x=L log(s)>0`, and let `A` be uniformly distributed
on the mode box. Define

\[
 Z=E(A)/L-1,\qquad f(y)=\log(1+e^{-y}).
\]

The distribution of `Z` is symmetric. The logarithmic difference
in (14), divided by `n`, equals the expectation of the even function

\[
 g_R(z)=x-2x|z|+2f((1+z)x)+2f((1-z)x)
                         -R f(4x|z|/R).             \tag{15}
\]

For `0<R<=4`, this function decreases as `|z|` increases on `[0,1]`.
Indeed, for `z>=0`, put `A=x/2`, `B=zx/2`. Differentiating gives

\[
 \frac{g'_R(z)}x
 =\tanh(A+B)-\tanh(A-B)-2\tanh(4B/R).               \tag{16}
\]

Since `A>=B>=0`, the first difference is at most `tanh(2B)`.
Also

\[
 \tanh(2B)\le2\tanh B\le2\tanh(4B/R).
\]

This proves the monotonicity, with strict inequality when `z>0`.

For one coordinate, `Z` is uniform on

\[
 \left\{\frac{2j+1-M}{M}:0\le j<M\right\}.
\]

For two coordinates it is `(Z_1+Z_2)/2`, where `Z_1,Z_2` are
independent copies of that variable. The averaged variable is more
concentrated around zero in the exact sense

\[
 \Pr\left(\left|\frac{Z_1+Z_2}{2}\right|\le a\right)
       \ge \Pr(|Z_1|\le a)\qquad(a\ge0).             \tag{17}
\]

For a direct discrete proof, it suffices to check the positive
grid thresholds `a=(2r-1)/M`, `1<=r<=M/2`. At such a threshold,
the two probabilities are

\[
 1-\frac{(M-2r)(M-2r+1)}{M^2},\qquad \frac{2r}{M}.
\]

Their difference is `(M-2r)(2r-1)/M^2>=0`. Between consecutive
grid thresholds the second probability is constant, whereas the
first is nondecreasing. This proves (17) everywhere.

Combining (15)--(17), at the same value of `x`, gives

\[
 \frac{\log\mathcal B_{2,M}(e^{x/(2M)})
       -R\log\mathcal Z_{2,M,R}(e^{x/(2M)})}{M^2}
 \ge
 \frac{\log\mathcal B_{1,M}(e^{x/M})
       -R\log\mathcal Z_{1,M,R}(e^{x/M})}{M}.        \tag{18}
\]

The right side is strictly positive at `R=2sqrt(2)` by the
one-species interval theorem. Thus the centered two-coordinate
box cannot cross the half-score boundary at any positive scale,
for any permitted side length. Formula (18) does not assert an
obstruction for other centers, other physical gauges, or coherent
restrictions of the physical alphabet.

## 4. Small coherent restrictions and exact finite certificates

At `k=M=2`, the four mode energies are `(2,4,4,6)` and `L=4`.
A physical symbol is a subset of these modes, encoded by its
four-bit mask. Five occupation-and-energy spaces have dimension
two:

\[
 (2,4),\quad(3,5),\quad(6,9),\quad(10,12),\quad(11,13).
                                                               \tag{19}
\]

Each pair in (19) lists two physical masks. All remaining
occupation-and-energy spaces have dimension one.

For one such pair choose nonzero binary row vectors `ell,r` with
`ell*r^T=1`. There are six choices. Apply `ell` on boundaries
`B,X` and `r` on `V,Y`, retaining all other states. The two maps
are a split projection. They preserve every height in (5) and (7)
because they mix states of the same occupation and energy.

At the endpoint the nonidentity Pascal terms strictly increase
energy. The split projection preserves this filtration and its
identity diagonal. The same common-type argument as in Section 1
therefore proves the exact typed crossing after each restriction.

The screen includes all `5*6=30` single-pair restrictions and
120 deterministic samples from the `6^5` choices which restrict
every pair. It constructs the actual four-leg tensor by contraction
and computes both complete physical Smith profiles.

Let `nu,mu` be these physical profiles and `G` the total physical
height list. At `R=3`, substitute `s=u^3` in the pressure difference
and clear its smallest Laurent power. The resulting integer
polynomial is

\[
 P(u)=u^b\left[
   \left(\sum_i u^{-3\nu_i}\right)
   \left(\sum_i u^{-3\mu_i}\right)
   -\left(\sum_j u^{-G_j}\right)^3\right].          \tag{20}
\]

The integer `b` makes its smallest exponent zero. For every retained
case, every coefficient of `P(1+v)` is nonnegative and at least one
is positive. This exact integer certificate proves `P(u)>0` for
every `u>1`. Hence each case has pressure strictly above `3/4`
at every positive scale, at the inherited heights.

This finite result is a stronger rejection test than a sampled
scale optimization. It does not cover the other `6^5-120` choices,
other height assignments, or extension-field coherent maps.

## 5. Reproduction and limits

The implementation and checker are
[tensor_pascal_box_crossings.py](../tools/tensor_pascal_box_crossings.py)
and
[check_tensor_pascal_box_crossings.py](../tools/check_tensor_pascal_box_crossings.py).
The optional finite screen is
[search_tensor_pascal_box_codes.py](../tools/search_tensor_pascal_box_codes.py).

The constructor uses only polynomial matrices and finite exterior
compounds. The theorem permits arbitrary fixed `k,M` satisfying the
power-of-two conditions. The displayed pressure calculation fixes
the centered gauge, and the all-size comparison fixes `k=2`.

This work supplies an explicit additional crossing family and
locates two unsuccessful routes within it. It supplies no local
seed below the relevant transfer threshold and leaves the complete
factoring score at `1/2`. The derivation and controls used
substantial OpenAI assistant support; no Lean formalization or
external mathematical review is claimed.
