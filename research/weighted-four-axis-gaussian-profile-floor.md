# The four-axis floor survives additive mode degrees and common physical types

## Result and scope

The [four-axis profile floor](four-axis-gaussian-profile-floor.md) extends
to any number of coupled Gaussian modes per physical vertex, any common
physical subalphabet and type, and additive degree bounds that can depend
on both the vertex and the mode. For the common-scale criterion below,

```math
\rho_{\mathrm{profile}}(\pi,t)>\frac12.
```

Individual degree heights may be negative. The antipodal endpoint forces
their averages to be nonnegative, which is the condition used in the proof.
The entries need not be homogeneous: upper degree bounds suffice.

This is a restriction on the stated profile estimate. It leaves open
nonadditive physical degree bounds, coefficient-specific estimates,
non-Gaussian tensors, and other methods of controlling actual ranks. It
does not prove a factoring lower bound or improve the achieved exponent.

The profiles in this estimate are those of the full ambient Gaussian
cuts. A common subalphabet changes the entropy and degree charge in the
formula. Computing smaller actual profiles after restricting the physical
tensor is a different criterion and is not excluded by this theorem.

## 1. Weighted Gaussian setup

Let each vertex `v in {0,1}^4` have `r>=1` modes. Work over a field of
characteristic two. Let `K` be an alternating polynomial matrix indexed
by `(v,j)`, `1<=j<=r`, with endpoint

```math
K_{(v,i),(w,j)}(1,1,1,1)=[w=\bar v][i=j].        \tag{1}
```

Choose real heights `h_(v,j)` such that every nonzero entry satisfies

```math
\deg K_{(v,i),(w,j)}\le h_{v,i}+h_{w,j}.         \tag{2}
```

Total polynomial degree is meant here. Rational heights are sufficient
for an effective construction; the inequality below does not require
their integrality. Define the averaged heights

```math
\bar h_j=\frac1{16}\sum_v h_{v,j}.              \tag{3}
```

The endpoint entry `K_((v,j),(bar v,j))` is a nonzero polynomial, so
`h_(v,j)+h_(bar v,j)>=0`. Summing over the eight antipodal pairs gives
`bar h_j>=0`, even when some individual heights are negative.

Let `C` be any subset of the physical alphabet `{0,1}^r`, and let `pi`
be a probability distribution on `C` with positive entropy

```math
H(\pi)=-\sum_{b\in C}\pi_b\log\pi_b>0.
```

Every physical leg uses this same subalphabet and type. Write

```math
E(\pi)=\sum_{b\in C}\pi_b\sum_{j=1}^r\bar h_j b_j.
\tag{4}
```

The physical tensor is `G_b=Pf K[S(b),S(b)]`. Every matching term in this
Pfaffian has degree at most the sum of heights at its occupied modes.
Consequently a word tensor of length `L` on the common type has total
degree at most `16L E(pi)`. This remains true when individual heights are
negative: a nonzero polynomial cannot have a negative degree bound.

At each axis `a`, the one-particle crossblock has size `8r`. Its endpoint
is a permutation matrix. Let its generic axis valuations be

```math
0\le s_{a,1}\le\cdots\le s_{a,8r},\qquad
\Delta_j(a)=\sum_{i=1}^j s_{a,i},\qquad
Z_a(t)=\prod_{i=1}^{8r}(1+t^{s_{a,i}}).
```

The other variables are treated as rational-function coefficients, without
specializing them to one. The Gaussian exterior-power profile and the
additive degree cap give the numerical common-scale expression

```math
\rho_{\mathrm{profile}}(\pi,t)=
\frac{\sum_{a=0}^3\log Z_a(t)-16E(\pi)\log t}
     {32H(\pi)},\qquad0<t<1.                    \tag{5}
```

If an entry has degree smaller than its cap, (5) still uses the cap. A
coefficient-dependent argument that exploits the smaller degree can be a
different criterion.

## 2. Two complementary minors bound every mode subset

Fix axes `a!=b` and a mode subset `J` of size `k`, `1<=k<=r`. For
`sigma in {0,1}`, use the `4k` rows

```math
R_\sigma=\{(v,j):v_a=0,\ v_b=\sigma,\ j\in J\}
```

and their antipodal columns `(bar v,j)`, in the same order. Denote the
determinant by `F_sigma`. Its endpoint is one. It is a minor of the
`a`-cut. At the `b`-cut it is a minor when `sigma=0`, and its transpose
is a minor when `sigma=1`; alternating symmetry preserves its determinant
valuation.

By (2), the degree of `F_sigma` is at most

```math
D_\sigma=\sum_{(v,j)\in R_\sigma}
                 (h_{v,j}+h_{\bar v,j}).
```

Every determinant term has this same upper bound, regardless of its
permutation. Since `F_sigma` is nonzero,

```math
\Delta_{4k}(a)+\Delta_{4k}(b)
\le v_a(F_\sigma)+v_b(F_\sigma)
\le D_\sigma.
```

The two row/antipodal-column sets partition all vertices for the chosen
modes, so `D_0+D_1=16 sum_(j in J) bar h_j`. Taking the smaller bound gives

```math
\Delta_{4k}(a)+\Delta_{4k}(b)
\le8\sum_{j\in J}\bar h_j.
```

Sum over the six unordered axis pairs. Each axis is counted three times:

```math
\sum_{a=0}^3\Delta_{4k}(a)
\le16\sum_{j\in J}\bar h_j.                    \tag{6}
```

No assumption on the couplings between different modes was used.

## 3. Ordered partial sums retain the individual mode heights

Relabel the modes so that `bar h_1<=...<=bar h_r`. In (6), choose the
first `k` modes, for every `k`. Put

```math
u_i=\frac14\sum_a s_{a,i}\quad(1\le i\le4r),
\qquad v_i=\bar h_{\lceil i/4\rceil}.
```

Both sequences are nondecreasing, and (6) says

```math
\sum_{i=1}^{4k}u_i\le4\sum_{j=1}^k\bar h_j
=\sum_{i=1}^{4k}v_i.
```

The partial sums of a nondecreasing sequence form a convex polygonal
function of the index. Between consecutive multiples of four, its values
lie below the line between its endpoint values. The `v` partial sums are
exactly linear on each such interval. Hence

```math
\sum_{i=1}^j u_i\le\sum_{i=1}^j v_i
\qquad(1\le j\le4r).                            \tag{7}
```

For completeness, (7) implies the required decreasing-convex inequality
directly. Let `phi` be differentiable, convex, and decreasing. Its tangent
bound gives
`sum(phi(u_i)-phi(v_i))>=sum phi'(v_i)(u_i-v_i)`. Write the partial-sum
differences as `E_j<=0` and put `a_i=phi'(v_i)`. Then `a_i` is
nondecreasing and nonpositive, and summation by parts gives

```math
\sum_i a_i(u_i-v_i)
=a_{4r}E_{4r}+\sum_{j<4r}(a_j-a_{j+1})E_j\ge0.
```

Apply this to `phi(s)=log(1+t^s)`, and apply Jensen across the four cuts
at each fixed index `i`. The first half of the profiles satisfies

```math
\sum_a\sum_{i=1}^{4r}\phi(s_{a,i})
\ge4\sum_i\phi(u_i)
\ge4\sum_i\phi(v_i)
=16\sum_{j=1}^r\log(1+t^{\bar h_j}).             \tag{8}
```

All remaining profile terms are strictly positive. Thus the left side
of (8), with all `8r` orders retained at each cut, is strictly larger
than its final right side.

## 4. Joint types and subalphabets cannot reverse the inequality

The entropy variational inequality for the distribution `pi` gives

```math
H(\pi)+E(\pi)\log t
\le\log\sum_{b\in C}t^{\sum_j\bar h_j b_j}
\le\log\prod_{j=1}^r(1+t^{\bar h_j}).            \tag{9}
```

The final product is the partition function for the full binary alphabet,
used only as an upper bound. The actual alphabet may be a proper subset,
and `pi` may have arbitrary correlations between its modes. Combining (8)–(9) with the
strictly positive second half of the profiles proves

```math
\rho_{\mathrm{profile}}(\pi,t)
>\frac{16H(\pi)}{32H(\pi)}=\frac12.              \tag{10}
```

Thus additive mode degree caps and arbitrary common types, including
multiple-mode subalphabets, cannot provide a common-scale four-axis
profile value below one half. A different criterion must use information
not accounted for by (2) and (5).

## Exact controls

Run:

```sh
python3 tools/check_weighted_four_axis_gaussian_profile_floor.py \
  --output evidence/weighted-four-axis-gaussian-profile-floor.json
```

The checker constructs small polynomial kernels with several mode heights,
including negative individual heights and nonhomogeneous entries. It checks
both orientations of each shared minor, their endpoint, degree caps, and
coordinate valuations exactly over `F2`. It separately checks the
block-partial-sum implication, its partition-product comparison, and the
cleared-denominator entropy inequality for correlated types and restricted
physical alphabets.

The finite controls do not replace the general proof. This note extends the
shared-minor argument of the companion four-axis note. It used substantial
OpenAI assistant support. No Lean verification or new achieved factoring
exponent is claimed.
