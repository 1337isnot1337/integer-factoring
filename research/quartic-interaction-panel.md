# A fixed four-fermion interaction panel

## Registered mechanism and result

The [registered panel](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338563)
tests an actual change of tensor class: one four-fermion interaction
couples two binary modes inside each of four physical labels. This
breaks the Gaussian Pfaffian relation while preserving the crossing
endpoint and occupied-bit homogeneity.

All six specified perturbations have complete two-cut profiles. Their
paired channel partition functions are **no smaller** than the Gaussian
baseline at every scale 0<t<1. Consequently their tested common-type
profile envelope has the exact all-scale floor

~~~
rho_profile >= 3/4 > 1/sqrt(2).
~~~

No numeric angle or physical-frequency optimization is needed to settle
this panel. This is a negative result for these six tensors under the
stated estimate, not a theorem about all non-Gaussian tensors, all
four-fermion interactions, other physical types or actual coefficient
rank compression. The complete factoring score remains **1/2**.

## 1. Acquire the tensors, not just relaxed profiles

Index the four physical labels by X,Y,U,V, each with copy modes 0,1.
The eight mode indices are `2q+i`. Start with two copies of a homogeneous
four-mode degree-two Gaussian. In one copy, the alternating Cayley
kernel has nonzero group steps

~~~
k_01=y^2+xy,       k_10=x^2+xy,       k_11=xy.
~~~

The eight-mode kernel is the direct sum of these copies. At (x,y)=(1,1)
it is the matching X_i--V_i, Y_i--U_i. Let G_mask be its principal
Pfaffians. The empty coefficient is one and every nonzero coefficient
is homogeneous of degree |mask|.

Use the two copy-mixing supports

~~~
S153={X0,Y1,U0,V1},       S149={X0,Y0,U0,V1}.
~~~

For each support S and one of the three degree-four polynomials

~~~
h=x^4+y^4,       x^4+x^2y^2,       y^4+x^2y^2,
~~~

acquire the complete perturbed state

~~~
T_mask=G_mask+h*G_(mask minus S)       if S subset mask,
       G_mask                         otherwise.       (1)
~~~

This is multiplication of the Grassmann state by `1+h*xi_S`. Its
endpoint is unchanged because h(1,1)=0, and (1) is still homogeneous
of degree |mask|. Vacuum and all two-mode coefficients are unchanged,
but T_S differs from G_S by the nonzero polynomial h. A normalized
Gaussian is determined by those two-mode coefficients through its
Pfaffians; therefore each perturbed state is genuinely non-Gaussian.

The interaction couples the two mode copies at a common word position.
The independent-position-class binary word barrier therefore cannot
simply be invoked without checking its independent-permutation premise.
The exact profile calculation below settles this panel directly.

## 2. Complete profiles with independent hinge checks

Flatten across (X,Y)|(U,V) for the row cut and (X,U)|(Y,V) for the
column cut. Each is a 16-by-16 physical matrix. For row orders set y=1
and use x-adic Smith elimination; for column orders set x=1 and use
y-adic elimination.

This dehomogenization retains the generic valuations. Each entry has
degree equal to row occupation plus column occupation, so the generic
matrix is obtained from the dehomogenized matrix at x/y, or y/x, by
monomial row and column scalings in the other parameter. Those scalings
are units for the chosen valuation.

All 16 pivots resolve at both cuts for all seven tensors, including the
baseline. No infinite-order or unresolved tail is inferred. A separate
dense binary Toeplitz calculation verifies

~~~
rank modulo z^ell = sum_j max(ell-s_j,0)
~~~

at every cutoff through one above the largest reported order. These
checks independently reconstruct the complete hinge profile.

The [receipt](../evidence/quartic-interaction-panel.json) retains all
256 polynomial entries per tensor, both complete profiles and every
Toeplitz rank. The maximum order in the panel is eleven. In total it
contains 1,792 endpoint/homogeneity entries, 14 complete profiles and
131 independent hinge checks.

## 3. An exact certificate for every scale

For a candidate's row and column orders r_i,c_j, form the 256 paired
orders a_(i,j)=r_i+c_j. Their counting polynomial is exactly
`Z_row(t)*Z_column(t)`. Let b_k be the corresponding sorted baseline
paired orders.

For each of the six candidates, the checker verifies all 256 inequalities

~~~
sum_(j<=k) a_j <= sum_(j<=k) b_j.                       (2)
~~~

Every prefix gap is retained as an integer. For a sorted order list s,
its hinge sum is `R_s(ell)=max_k(k*ell-sum_(j<=k)s_j)`. Thus (2)
implies `R_a(ell)>=R_b(ell)` for every ell. The exact generating identity

~~~
sum_j t^(s_j)=(1-t)^2 sum_(ell>=1) R_s(ell)t^(ell-1)
~~~

then proves candidate partition-product domination at every 0<t<1.
This is a complete all-scale certificate, rather than a grid of angles.
The comparisons survive any extension of the coefficient field because
the source tensors and all their minors are defined over F_2.

## 4. Normalize the tested physical type

At word length m, give each of the eight component modes occupied
fraction p. Each physical parent letter contains two such words and
has entropy at most 2mH(p). Homogeneity gives total coefficient degree
8pm. The common-scale profile envelope is consequently

~~~
R(p,t)=(log Z_row(t)+log Z_column(t)-8p log t)/(8H(p)).
~~~

For the baseline, each physical profile has counting polynomial

~~~
Z_base(t)=(1+t)^2(1+t^3)^2.
~~~

The [exact partition inequality](fixed-cube-projection-profile-floor.md#4-an-exact-all-scale-partition-function-certificate)
already proves
`(log Z_base(t)-4p log t)/(4H(p))>=3/4` for every p and t.
The candidate partition product is at least Z_base(t)^2 by (2), while
its degree cap and physical entropy are the same. Its envelope therefore
has the same floor. Nothing here asserts that the profile estimate equals
the actual factor ranks or that other physical restrictions obey this
same expression.

## 5. Scope, reproduction and stopping condition

Run make verify-quartic-interactions. The
[checker](../tools/check_quartic_interaction_panel.py) acquires (1) from
the original Pfaffian implementation and uses two existing independent
Smith/Toeplitz implementations. All 1,536 paired-prefix comparisons pass,
and it replays the exact Bernstein certificate for the baseline floor.
The scalar-field arithmetic is exact throughout.

This meets the registered stopping condition: the specified panel cannot
cross the target under this common-type envelope. No larger field, degree
or support search is justified by this checkpoint alone. A new physical
restriction, interaction mechanism, rank estimate or evaluator would
require its own question and full cost argument.

The Gaussian, polynomial, Smith and Toeplitz machinery retain their
existing contributors' credit. The interaction diagnostic and all-scale
comparison were developed with OpenAI assistant support. There is no
new Lean, stronger complete factoring bound or independent mathematical
acceptance.
