# Shared endpoint minors through eight Gaussian axes

## Results and scope

The shared-minor argument has a uniform form in every dimension. For a
polynomial Gaussian kernel on the vertices of a `d`-dimensional cube,
every set of `q` axes gives a common endpoint minor. Its degree bounds the
first `2^(d-q)` Smith orders at all those axes simultaneously. The complete
family of these constraints excludes the common-scale, unweighted
profile criterion through **eight axes**.

The following floors hold pointwise for the pressure defined below.

| Axes `d` | Proved floor `rho_d` | Exact comparison with `1/sqrt(d)` |
| --- | --- | --- |
| 2 | `3/4` | `2*3^2 > 4^2` |
| 3 | `3/5` | `3*3^2 > 5^2` |
| 4 | `21/40` | `4*21^2 > 40^2` |
| 5 | `9/20` | `5*9^2 > 20^2` |
| 6 | `5/12` | `6*5^2 > 12^2` |
| 7 | `19/50` | `7*19^2 > 50^2` |
| 8 | `177/500` | `8*177^2 > 500^2` |

These inequalities are ordinary proofs with exact integer certificates.
They require neither Cayley symmetry nor equal profiles at the different
cuts. The three-axis case extends the
[endpoint-minor budget calculation](cube-gaussian-profile-budget.md).
The [four-axis proof](four-axis-gaussian-profile-floor.md) gives a shorter
argument for the weaker floor `1/2` using only pairwise minors.

At **nine axes**, these same degree constraints permit an explicit integer
profile with pressure strictly below `1/3`, at rational physical frequency
`13/100` and rational scale `447/448`. It is a relaxed construction target.
No kernel realizing that profile is supplied. Thus nine is the first
dimension admitted by this particular relaxation, among `d>=2`; it is
not a new factoring exponent.

The theorem concerns one binary mode per physical leg, a common kernel
degree cap, a common physical frequency, and one common profile-count
scale. Homogeneous kernels are included, but homogeneity is unnecessary
when the same cap is charged. It is not a lower bound on actual coefficient
ranks. Weighted physical modes, non-Gaussian tensors, different rank
estimates, and information about actual coefficient support beyond the
common degree cap are outside its scope. The achieved factoring score
remains `1/2`.

## 1. The common-scale criterion

Let `K(x_1,...,x_d)` be an alternating `2^d` by `2^d` polynomial matrix
over a field of characteristic two. Index its rows by `v in {0,1}^d`.
Every entry has total degree at most `D>0`, and

```math
K_{u,v}(1,...,1)=[v=\bar u],
```

where `bar u` is the antipodal vertex. Put `m=2^(d-1)`. At axis `a`,
the one-particle crossblock `B_a` has size `m`. Its endpoint is a
permutation matrix. Its ordered Smith valuations in `x_a`, over the
rational function field in the other variables, are therefore finite:

```math
0\le s_{a,1}\le\cdots\le s_{a,m}.
```

Write `Delta_j(a)=sum_(i=1)^j s_(a,i)` and

```math
Z_a(t)=\prod_{i=1}^m(1+t^{s_{a,i}}),\qquad 0<t<1.
```

The associated Gaussian physical tensor consists of the principal
Pfaffians of `K`. On words with common occupied fraction `p`, its total
coefficient degree at word length `L` is at most `mDpL`. Its unweighted,
common-scale profile expression, charging this cap, is

```math
R_d(p,t)=
\frac{\sum_{a=1}^d\log Z_a(t)-mDp\log t}
     {dmH(p)},\qquad
H(p)=-p\log p-(1-p)\log(1-p).                         \tag{1}
```

All logarithms in this note are natural. Formula (1) is the direct
`d`-axis analogue of the expression used by the
[three-axis construction](cube-homogeneous-gaussian.md). The natural
half-score target for repeated `d`-axis regrouping is `R_d<1/sqrt(d)`.
This note bounds the numerical profile criterion; it does not assert an
effective general-dimensional factoring transfer from arbitrary
multivariate Smith factorizations.

## 2. Every subset of axes supplies a shared minor

Fix a set `Q` of `q` coordinate axes. Select rows whose coordinates in
`Q` are all zero, and select their antipodal columns in the same order.
There are `j=2^(d-q)` rows and columns. At the endpoint the resulting
minor is the identity. Its determinant `F_Q` is a nonzero polynomial
of total degree at most `jD`.

For every `a in Q`, this is a minor of `B_a`. The valuation of any
`j` by `j` minor bounds the first `j` Smith orders. Consequently

```math
\sum_{a\in Q}\Delta_j(a)
\le\sum_{a\in Q}v_{x_a}(F_Q)
\le jD.                                                   \tag{2}
```

The final inequality follows because every monomial of `F_Q` has total
degree at most `jD`. It remains valid when determinant terms cancel:
the surviving polynomial is still nonzero and respects the degree cap.

Summing (2) over all `q`-element axis sets counts each axis
`binom(d-1,q-1)` times. Thus, for every `q=1,...,d`,

```math
\boxed{\quad
\sum_{a=1}^d\Delta_{2^{d-q}}(a)
\le\frac{d}{q}\,2^{d-q}D.\quad}                         \tag{3}
```

This includes the shared-entry constraint at `q=d` and the complete
determinant budget at `q=1`. No structure of the matrix beyond the degree
cap and the endpoint was used.

## 3. One explicit relaxed profile dominates the entire family

Average the normalized ordered valuations across the axes:

```math
z_i=\frac1{dD}\sum_a s_{a,i}.
```

The `z_i` are nonnegative and nondecreasing. Equation (3) gives

```math
\sum_{i=1}^{2^{d-q}}z_i\le\frac{2^{d-q}}q.              \tag{4}
```

Define a nondecreasing reference list `u_1,...,u_m` by the following
values and multiplicities:

```math
\begin{array}{c|c}
\text{value}&\text{multiplicity}\\ \hline
1/d&1\\
2/(q-1)-1/q&2^{d-q}\quad(q=d,d-1,...,2).
\end{array}                                             \tag{5}
```

Its cumulative sums are exactly `2^(d-q)/q` at every index occurring
in (4), and its total is `m`. Between those indices the cumulative
sum is linear. Since cumulative sums of an increasing list are convex,
(4) implies every prefix inequality

```math
\sum_{i=1}^jz_i\le\sum_{i=1}^ju_i\qquad(1\le j\le m).
                                                               \tag{6}
```

For completeness, (6) implies
`sum_i phi(z_i)>=sum_i phi(u_i)` for every twice differentiable decreasing
convex function `phi`. Indeed,

```math
\sum_i(c-z_i)_+=\max_{0\le j\le m}
\left(jc-\sum_{i=1}^jz_i\right)
\ge\sum_i(c-u_i)_+.
```

Choose `M` at least as large as every list entry and use

```math
\phi(s)=\phi(M)+\phi'(M)(s-M)
        +\int_0^M(c-s)_+\phi''(c)\,dc.
```

The prefix inequality proves the integral comparison. The total-sum
inequality and `phi'(M)<=0` prove the linear comparison. This also
handles a strict deficit in the total Smith budget.

Take `phi(s)=log(1+t^(Ds))`. Jensen's inequality across the `d` axes,
followed by (6), now gives

```math
\sum_a\log Z_a(t)
\ge d\sum_{i=1}^m\log(1+t^{Du_i}).                       \tag{7}
```

Let `x=-D log t` and define

```math
A_d(x)=\frac1m\sum_{i=1}^m\log(1+e^{-u_i x}).
```

Equations (1) and (7) imply the uniform relaxation

```math
R_d(p,t)\ge\frac{A_d(x)+(x/d)p}{H(p)}.                  \tag{8}
```

The list (5) satisfies every constraint (4), so (7) is the exact
convex minimum allowed by these aggregated budgets. Realizing the list
by one common polynomial kernel is a separate algebraic problem.

## 4. Exact all-scale certificates through dimension eight

For every floor `rho_d` in the opening table, the checker proves

```math
A_d(x)>\rho_d\log(1+e^{-x/(d\rho_d)})
\qquad(x\ge0).                                         \tag{9}
```

The entropy inequality

```math
\rho H(p)-(x/d)p
\le\rho\log(1+e^{-x/(d\rho)})
```

then combines with (8) to give `R_d(p,t)>rho_d` for every
`0<p<1` and `0<t<1`. At `t=1`, the pressure is at least one directly.
The strict inequalities are pointwise; their infima are at least the
stated rational floors, already strictly above `1/sqrt(d)`.

### The tail is analytic

The smallest value in (5) is `1/d`, with multiplicity one. The elementary
bounds `y/2<=log(1+y)<=y`, valid for `0<=y<=1`, give

```math
A_d(x)\ge\frac1{2m}e^{-x/d},\qquad
\rho\log(1+e^{-x/(d\rho)})\le\rho e^{-x/(d\rho)}.
```

Thus (9) follows whenever

```math
e^{x(1-\rho)/(d\rho)}>2m\rho.                           \tag{10}
```

For each retained `d,rho`, the first 33 positive Taylor terms of the
left side at `x=32` already exceed the right side. This is an exact
rational comparison. Monotonicity proves (10) for every `x>=32`.

### The compact interval uses integer enclosures

The [checker](../tools/check_multiaxis_gaussian_profile_floor.py) uses
interval endpoints at binary scale `S=2^128`. A pair of integers
`[L,U]` means the true number lies in `[L/S,U/S]`.

For `exp(-r)`, it first divides the nonnegative rational `r` by a power
of two until the argument is at most one half. Directed integer
rounding bounds each Taylor term. The even 64th partial sum is an
upper bound, and subtraction of the next term gives a lower bound.
Repeated interval squaring restores the original argument.

For `log(1+y)`, `0<=y<=1`, put `z=y/(2+y)`, so `0<=z<=1/3`, and use

```math
\log(1+y)=2\sum_{j=0}^{63}\frac{z^{2j+1}}{2j+1}
          +\text{a nonnegative remainder smaller than }3^{-128}.
```

Every multiplication and division rounds outwards. The final one-unit
tail allowance is also rounded upwards. These operations establish
the enclosures without a floating-point transcendental routine.

On a candidate interval `[a,b]`, the left side of (9) is bounded below
using its value at `b`; the right side is bounded above using its value
at `a`. A strict integer gap certifies the whole interval. Recursive
bisection produces the following complete covers of `[0,32]`:

| `d` | Intervals | Maximum dyadic depth |
| --- | --- | --- |
| 2 | 65 | 10 |
| 3 | 31 | 8 |
| 4 | 47 | 9 |
| 5 | 25 | 7 |
| 6 | 51 | 9 |
| 7 | 59 | 9 |
| 8 | 195 | 12 |

The checker independently verifies that every consecutive interval
touches the next and that their union is exactly `[0,32]`. The
[retained evidence](../evidence/multiaxis-gaussian-profile-floor.json)
records every interval, the minimum strict integer gap, and a hash of
each cover. Direct rational Taylor sums at higher precision separately
check ten exponential enclosures and six logarithm enclosures.

## 5. An exact nine-axis construction target

For `d=9`, set `D=2520`. The integer list `D*u_i` from (5) is

| Order | Multiplicity |
| --- | --- |
| 280 | 1 |
| 350 | 1 |
| 405 | 2 |
| 480 | 4 |
| 588 | 8 |
| 756 | 16 |
| 1050 | 32 |
| 1680 | 64 |
| 3780 | 128 |

Use this same hypothetical list at all nine cuts. It has 256 orders,
total `256D`, and for every `q=1,...,9` its first `2^(9-q)` orders
sum to `2^(9-q)D/q`. Hence it satisfies even the individual subset
constraints (2), with equality.

Let `Z(t)` be its counting product. At the entirely rational choices

```math
t=447/448,\qquad p=13/100,
```

the same integer logarithm enclosures verify

```math
\log Z(t)-\frac{256D}{9}p\log t
<\frac{256}{3}H(p).                                     \tag{11}
```

The saved strict gap is

```math
\frac{11169004197368338418068673439992410253768}
     {75\cdot2^{128}}>0.
```

Consequently the relaxed value in (8) is strictly below `1/3`.
This does not supply any matrix entries, compatible Smith frames,
physical coefficient decomposition, or effective factoring evaluator.
It says exactly where these shared-minor degree constraints cease to
exclude a favorable profile. Additional polynomial identities or minor
constraints may still rule out the proposed list.

## 6. The relaxation eventually becomes much weaker

The failure at nine axes is consistent with an asymptotic limitation of
these degree budgets. For `d>=16`, put

```math
L=\lceil3\log_2d\rceil\le d,\qquad x=3L\log d,
\qquad p=1/d.
```

For the groups in (5) with `q<=L`, every order is at least `1/L`,
so their total contribution to `A_d(x)` is at most `d^-3`. All
remaining groups, including the first singleton when necessary, have
total normalized multiplicity at most `2^(1-L)<=2d^-3`. Their
logarithms are at most `log 2`. Therefore

```math
A_d(x)<3d^{-3},\qquad H(1/d)\ge(\log d)/d,
```

and the pressure of this relaxed list is at most

```math
\frac{3L}{d}+\frac{3}{d^2\log d}
=O\!\left(\frac{\log d}{d}\right).                      \tag{12}
```

This is eventually smaller than `1/sqrt(d)`. Thus a uniform exclusion
of all dimensions requires more information than the shared coordinate
minors used here. Conversely, a favorable relaxed profile still needs
an actual polynomial construction before it has any algorithmic force.

## Replay and provenance

Run

```sh
python3 tools/check_multiaxis_gaussian_profile_floor.py \
  --output evidence/multiaxis-gaussian-profile-floor.json
```

The checker uses only Python's standard library. Its certificate is
deterministic. It also verifies all 77 exact reference-prefix identities
for dimensions two through twelve. No finite kernel search is being
represented as a completeness proof: the arbitrary-kernel result is
the shared-minor argument in Section 2.

This ordinary proof builds on the repository's three-axis endpoint-minor
argument and the concurrent four-axis proof, with substantial OpenAI
assistant support. No new Lean formalization or external mathematical
review is claimed. The complete factoring exponent remains `1/2`.
