# A two-axis Gaussian profile floor with arbitrarily many modes per label

## Result and scope

The [registered extension question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338160)
asks whether putting multiple Gaussian modes inside each parent label
escapes the earlier one-mode two-axis profile floor. With a common
unweighted physical type, it does not.

Let each of the four labels X,Y,U,V contain L>=1 modes. Let K(x,y) be
an alternating polynomial matrix of order 4L in characteristic two, with
entry degree at most D>0. Its endpoint pairs the X modes with V and the
Y modes with U, by fixed bijections. At each of the cuts

~~~
row: (X,Y) | (U,V),       column: (X,U) | (Y,V),
~~~

take the 2L one-particle Smith orders at its corresponding variable, over
the rational function field in the other variable. Write Z_row(t) and
Z_column(t) for their subset-product generating functions. Then the
common-scale, common-unweighted-type profile envelope satisfies

~~~
R_L(p,t)
 = (log Z_row(t)+log Z_column(t)-2LDp log t)/(4L H(p))
 >= 3/4,              0<p<1, 0<t<1.                    (1)
~~~

Here H(p)=-p log p-(1-p)log(1-p). The floor exceeds the growth-base-two
target 1/sqrt(2). It is uniform in L, D, field, kernel and scale, without
Cayley symmetry or homogeneity.

This extends the [one-mode shared-minor argument](multiaxis-gaussian-profile-floor.md)
to equal block multiplicities. It includes two ordinary Gaussian modes
per superlabel, as in the cube regrouping proposals. It does **not**
cover unequal or weighted physical types, negative global polynomial
powers, arbitrary rational poles, actual coefficient-rank improvements,
joint support beyond the degree cap, non-Gaussian tensors or another
global evaluator. In particular, a large profile-envelope bound is not
a lower bound on an implementation's actual running time. The complete
factoring score remains **1/2**.

## 1. A common L-minor replaces the one-mode scalar minor

Let B_r and B_c be the two one-particle crossblocks, each of order 2L.
Their endpoint matrices are permutations, so their determinants are
nonzero polynomials. The total Smith orders at each cut are therefore
at most 2LD, the degree cap on that determinant:

~~~
sum_j s_(r,j) <= 2LD,       sum_j s_(c,j) <= 2LD.        (2)
~~~

Now take the same submatrix F_0=K[X,V] in both cuts. It has order L,
and its endpoint is a permutation matrix. Thus

~~~
F=det F_0 != 0,       degree(F)<=LD.
~~~

Let Delta_L(r) and Delta_L(c) be the sums of the first L Smith orders.
Every L-minor has valuation at least the corresponding determinantal
prefix. The common minor consequently gives

~~~
Delta_L(r)+Delta_L(c)
 <= v_x(F)+v_y(F)
 <= degree(F)
 <= LD.                                                 (3)
~~~

The middle inequality uses nonnegative polynomial exponents: every
monomial of F has x exponent at least v_x(F) and y exponent at least
v_y(F). Cancellation cannot break the degree bound. No field-size or
generic-coefficient assumption is involved.

## 2. Jensen turns the prefixes into a partition-function lower bound

There are 2L orders in the two first-L prefixes combined, and another
2L orders in the two tails combined. Call their averages A and B.
Equations (2)--(3) say

~~~
0<=A<=D/2,          A+B<=2D.
~~~

For a fixed 0<t<1, f(s)=log(1+t^s) is decreasing and convex. Jensen
on the two pools gives

~~~
log Z_row(t)+log Z_column(t) >= 2L*(f(A)+f(B))
                              >= 2L*(f(A)+f(2D-A)).
~~~

The second inequality uses B<=2D-A and decreasing f. The function
g(A)=f(A)+f(2D-A) decreases for A<=D, since
`g'(A)=f'(A)-f'(2D-A)<=0`. With A<=D/2 this proves

~~~
log Z_row(t)+log Z_column(t)
 >= 2L*(log(1+t^(D/2))+log(1+t^(3D/2))).                (4)
~~~

This is the all-multiplicity step. The number of modes cancels from the
normalized bound; putting more equal Gaussian modes into a physical label
does not weaken this pair of prefix constraints.

## 3. The same exact partition inequality proves the floor

Set r=t^(D/2). Substitution of (4) into (1)'s expression gives

~~~
R_L(p,t)
 >= (log(1+r)+log(1+r^3)-2p log r)/(2H(p)).              (5)
~~~

The [fixed-projection certificate](fixed-cube-projection-profile-floor.md)
proves, for every r>0,

~~~
(1+r)^2(1+r^3)^2 >= (1+r^(4/3))^3.
~~~

Its proof substitutes r=u^3 on (0,1], writes the difference as u^3 Q(u),
and gives all 65 positive degree-64 Bernstein coefficients of Q. The
minimum is `1589593605121/28736878651740`; monomial reconstruction is
checked exactly. Reciprocity handles r>1. This is an interval-wide
algebraic certificate, not a numerical scan.

Taking logarithms and using entropy duality at v=(4/3)log r yields

~~~
2log(1+r)+2log(1+r^3)-4p log r
 >= 3*(log(1+exp(v))-p*v)
 >= 3H(p).
~~~

Equation (5) proves R_L>=3/4. The target comparison is exact:
9/16>1/2. No optimization of p, t, block size or degree can lower this
particular envelope below the target.

## 4. Why the normalization is 4L H(p)

At word length m, put the same occupied fraction p on every one of the
4L physical modes. A parent superlabel consists of L such words and has
entropy at most LmH(p). A Gaussian entry has degree at most D per paired
occupied mode, so the coefficient-degree cap is 2LDpm.

The two-cut width budget divides the sum of logarithmic channel counts
by four times the parent-letter entropy. After cancelling m, these are
precisely the numerator term 2LDp and denominator 4LH(p) in (1).
No factor of L is gained by calling a collection of modes one label.

The expression is a profile estimate from that degree cap. Actual
coefficient support or smaller factor ranks can give extra information;
this theorem does not assert they saturate the expression. A different
physical restriction changes the calculation and needs its own proof.

## 5. Exact controls and next step

The [checker](../tools/check_block_gaussian_profile_floor.py) and
[receipt](../evidence/block-gaussian-two-axis-floor.json) reuse the exact
partition certificate and test block multiplicities L=1,2,3,4. The
polynomial kernels are congruence-mixed copies of a four-mode degree-two
kernel. Separate constant changes of basis within the four groups preserve
the endpoint and cut Smith data while mixing mode coordinates.

All eight finite Laurent Smith transformations reconstruct exactly and
have profiles `(1 repeated L,3 repeated L)`. Their common minor is
`(xy)^L`, saturating (3), and their total orders saturate (2). All 480
one-particle endpoint entries are checked. These controls test the
block indices, normalization and tight budget conventions; the ordinary
minor/Jensen proof supplies every L and every qualifying kernel.
Run make verify-block-gaussian-floor.

The shared-minor technique, Gaussian exterior-power profile construction
and finite Laurent implementation retain the existing contributors' credit.
This extension and its checks used OpenAI assistant support. There is no
new Lean or independent mathematical acceptance.

The common-type block-multiplicity escape is now closed under the stated
profile method. Further work must supply information outside that envelope
or a different representation/evaluator, rather than increase the mode
count or repeat local profile optimization.
