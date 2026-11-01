# Physical subalphabets strengthen graded adjugate compatibility

Applying graded adjugate inequalities to matched physical subalphabets
gives necessary constraints missed by whole-alphabet checks. These exclude
the five-letter integer target retained in the
[graded adjugate note](graded-adjugate-caps.md), including its entire
proposed row-prefix sector. Rectangular subalphabets give an additional
cofactor gap. The statements hold for scalar strict crossings over every
field, with arbitrary degree and independent physical heights.

Three old five- and six-letter parameter points now have exact positive
critical-pressure lower certificates. These are fixed-parameter results,
not a floor over all scales or all energies. The complete factoring score
remains **1/2**. No realizing tensor or stronger factoring algorithm is
claimed.

## 1. Smith interlacing permits restriction

Choose `p` physical letters for the matched `x,v` legs and `b` for the
matched `y,u` legs. Restricting the tensor preserves its strict crossing
endpoint. Each resulting cut is a `pb`-square submatrix of the original
cut. Its individual increasing Smith orders are at least the corresponding
orders of the full cut, at the same place and with the same physical weights.

For completeness, this half of Smith interlacing has a short DVR proof.
Write a full matrix as `U diag(t^alpha_i) V`. Split off its first `i-1`
Smith terms. The remaining matrix is divisible by `t^alpha_i`. After any
row/column selection, the selected matrix is therefore `L+t^alpha_i C`
with `rank L <= i-1`. Diagonalize `L` over the DVR: modulo `t^alpha_i`
the selected matrix has at most `i-1` nonzero columns. Its image over
`O/(t^alpha_i)` needs at most `i-1` generators. A matrix with Smith orders
`beta` needs exactly `#{j: beta_j < alpha_i}` such generators, by its
Smith decomposition and Nakayama's lemma. Thus `beta_i >= alpha_i`.
Scalar shifts handle negative Laurent orders. This is the needed direction
of the classical [invariant-factor interlacing theorem of
Thompson](https://doi.org/10.1016/0024-3795(79)90144-7).

Consequently any necessary subalphabet inequality whose coefficients on
individual **absolute** Smith orders are nonnegative lifts to the full
profiles. Mere prefix dominance would not suffice for every rectangular
constraint below.

Use full-alphabet normalized profiles `r,s`, with `r_1=s_1=0`, normalized
deficit `d`, and sorted paired energies `X,Y` of minimum zero, as in the
[cell-set note](cell-minor-caps.md). The subalphabet constraints use the
first `pb` full orders. Its own profile minima need not equal the full
minima: the absolute-order conversion below handles that distinction.

## 2. Square subalphabets

For every `2 <= p <= q`, take the `p` smallest energies in each physical
list, set `n=p^2` and

```math
C=(n-1)d+p(\sum_{i\le p}X_i+\sum_{j\le p}Y_j).
```

Every inequality (1)--(2) of the graded adjugate note holds with this `C`,
these energy lists, and the first `n` full orders. Include both cut
orientations. For example,

```math
(k-1)s_2+k\sum_{i<n}r_i
 \le kC-TopCell_k(X_{\le p},Y_{\le p}),\quad 1\le k\le n. \tag{1}
```

To justify lifting, use absolute profiles `mu,nu` with paired minima zero,
so `r=mu-mu_1`, `s=nu-nu_1`, `d=-mu_1-nu_1`. Moving the deficit to the
left in (1) yields

```math
k\sum_{i<n}\mu_i+(k-1)\nu_2+[k(n-2)+1]\nu_1. \tag{2}
```

Every coefficient is nonnegative. The shared-line version has
`l sum(mu)-Top_l(mu)` in place of `k sum_{i<n}(mu)` and the analogous
nonnegative coefficients on `nu_1,nu_2`. Interlacing therefore proves
the asserted full-profile conditions even when restricted minima change.

The smallest energy subsets give the strongest lifted conditions. For
example, `k*p*sum(X,Y)-TopCell_k` is nondecreasing in each selected energy:
one energy can occur at most `p` times in the top-cell sum. The shared-line
right sides are minima of expressions with nonnegative coefficients and
are also coordinatewise nondecreasing. Ordered selected energy lists
dominate the corresponding smallest subsets. No enumeration of all
physical subsets is needed for these profile-only bounds.

## 3. Rectangular cofactor gaps

Take alphabet sizes `p,b`, write `n=pb`, `h=min(p,b)^2` and
`g=n-1-h`. Assume `g>0`; the square case uses the stronger rank-one
reshuffle lemma instead. Set

```math
W=b\sum_{i\le p}X_i+p\sum_{j\le b}Y_j,
\qquad C=(n-1)d+W.                                    \tag{3}
```

Then the following necessary constraint holds in both orientations:

```math
g s_2+\sum_{i<n}r_i \le C-X_p-Y_b.                    \tag{4}
```

More generally, for `1 <= l <= max(p,b)`, put

```math
L_l=\max\{TopX_l+lY_b\ (l\le p),\quad
             TopY_l+lX_p\ (l\le b)\}.
```

Include only the terms whose indicated condition holds. Then

```math
l g s_2+l\sum_{i\le n}r_i-Top_l(r)\le lC-L_l.          \tag{5}
```

Here and in (4), energies are the smallest chosen prefixes and profile
indices stop at `n`, rather than the full `q^2`. The `l=1` case is (4).

To prove the gap, normalize the restricted zero and infinity tensors as
in the graded adjugate proof, temporarily using their **own** minima.
If the crossed cut of the infinity tensor `T` has second order `s_2=k`,
a unit pivot gives a rank-one crossed lift modulo `t^k`. As an ordinary
operator on `p`-by-`b` matrices, the lift is

```math
M\longmapsto A M^T B^T,
\qquad A\in O^{p\times b},\quad B\in O^{b\times p}.
```

It has rank at most `h`. Write `T=R+t^k U`. Every `(n-1)`-minor uses at
least `n-1-h=g` factors from the perturbation: determinant multilinearity
and `rank R <= h` kill the other terms over the fraction field and hence
over the DVR. Thus **every entry** of `adj T` is divisible by `t^(gk)`.
Every order of its crossed cut is at least `gk`, including its first.
For `k=0` the conclusion is immediate.

The graded identity for the rectangular row operator is

```math
adj T=t^C D_{-a}\,(adj R_{\rm zero})(1/t)\,D_{-p}.
```

The weight sum in its cofactors is precisely (3). Its opposite complete
profile is `theta_i=sum r-r_(n+1-i)-C`, and its paired energies are
`-X,-Y`. The endpoint shared-line minor using `l` varying physical
letters gives `S_l(eta)+S_l(theta) <= -L_l`. Substituting
`eta_i >= g s_2` proves (5). This minor proof works for unequal physical
sizes: only the length of the chosen line matters.

The absolute-order left side of (4) is

```math
\sum_{i<n}\mu_i+g\nu_2+h\nu_1.                        \tag{6}
```

For (5) it is
`l sum(mu)-Top_l(mu)+l*g*nu_2+l*h*nu_1`.
All individual coefficients are nonnegative, even when rewriting (6)
in prefix sums would produce a negative coefficient. Individual Smith
interlacing therefore lifts all aspect ratios. The right sides are
coordinatewise nondecreasing in the selected energies, so again the
smallest physical subsets suffice. Averaging the two orientations makes
all the same linear conditions valid for `gamma=(r+s)/2`.

## 4. The old integer target is impossible

For `d=6`, energies `(0,0,12,16,18)` in both physical lists, the target

```text
0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,
36,36,37,38,40,41
```

passes all whole-alphabet graded caps. The zero-energy `2`-by-`2` core
has `C=18`, however. Equation (1) with `k=4` requires

```text
3*6 + 4*(0+6+9) = 78 <= 72,
```

which fails. The `2`-by-`5` rectangle has `g=5`, `C=146`. Equation (4)
requires `5*6+118 <= 146-18`, failing by 20. Nonnegative absolute
coefficients also exclude its entire proposed row-prefix sector,
including the degree-42 affine target with absolute bounds
`18+gamma` and `-24+gamma`. Changing frames, coefficient fields or degree
does not repair those fixed paired energies and lower profiles.

This is a stronger later compatibility update. The earlier negative
pressure computation and statement that the target passes the **older**
caps remain correct; neither supplied a realizing tensor.

## 5. Exact positive pressure at three fixed parameter points

For an averaged normalized profile `gamma`, the balanced critical
pressure is bounded below by

```math
d+2\log\sum_i e^{-\gamma_i}
 -\sqrt2(\log\sum_x e^{-X_x/\sqrt2}
                +\log\sum_y e^{-Y_y/\sqrt2}).           \tag{7}
```

For any rational probability vector `pi`, the entropy variational bound
is `2 log sum exp(-gamma) >= 2 H(pi)-2 pi dot gamma`.
The [dual packet](../evidence/subalphabet-adjugate-duals.json) supplies
nonnegative rational multipliers of cell/shared-prefix caps, ordering
constraints and the new subalphabet caps. Their combined coefficients
dominate `2*pi_i` for every `i>=2`; `gamma_1=0`. If their combined right
side is `B`, then `2*pi dot gamma <= B` exactly.

The checker verifies each rational inequality and evaluates
`d+2H(pi)-B-physical` by 90-bit outward dyadic arithmetic:

| Alphabet | Deficit | Both energy lists | Certified pressure lower bound |
| --- | --- | --- | --- |
| 5 | `1` | `(0,0,2,8/3,3)` | `>0.11494194` |
| 5 | `7/6` | `(0,0,19/9,23/9,37/9)` | `>0.12334661` |
| 6 | `8/7` | `(0,0,15/8,27/10,4,4)` | `>0.08859371` |

No numerical optimizer is trusted by the replay. These certificates concern
the listed deficit, energies and critical scale; they do not exclude other
scales along the same ray. They exclude arbitrary feasible complete profiles
at those points, rather than only the old specific optimizers.

## 6. Evidence and limits

Run `make verify-subalphabet-adjugate-caps`. The
[checker](../tools/check_subalphabet_adjugate_caps.py) verifies 4,016 literal
subalphabet inequalities on prime-field tensors and six five-letter
characteristic-two gauges, 96 individual Smith interlacing comparisons,
and 167 independently interpolated rectangular cofactors. Rectangular
controls include `2`-by-`5`, whose gap exceeds the channel rank bound,
and nontrivial rank-one jets of orders one and two. The
[packet](../evidence/subalphabet-adjugate-controls.json.xz) retains every
literal input, profile, exclusion, dual replay and source hash.

A separate registered numerical screen completed 600 parameter cases on
alphabets 5--8. Two optimizer calls failed; among the successful solves
with acceptable constraint residuals none had negative pressure. Its
[raw packet](../evidence/subalphabet-relaxation-random.json.xz) is exploratory
evidence only. It proves no global floor. Unequal energy lists, larger
alphabets, asymmetric supporting criteria and other constructions remain
research directions. Substantial OpenAI assistant support was used; no
new Lean formalization was performed.
