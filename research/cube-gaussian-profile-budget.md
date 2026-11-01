# Endpoint-minor budgets beyond the Cayley kernel class

## Scope and result

The homogeneous [common eight-leg construction](cube-homogeneous-gaussian.md)
has a profile bound near `rho=0.6248`. Replacing its Cayley kernel by an
arbitrary homogeneous eight-by-eight alternating kernel cannot reach the
half-score target through the same common-scale, unweighted Smith-profile
method. In fact, that method always has `rho>3/5`.

This statement concerns Gaussian seeds with one binary mode on each of
the eight physical legs, a common homogeneous kernel degree, a common
constant-weight physical type, and the profile-count bound at the three
coordinate variables. It is not a rank lower bound for the coefficient
pieces themselves. It does not cover coupled multiple modes on each
physical leg, coherent physical restrictions, different decompositions,
or stronger coefficient-rank estimates.

The argument uses endpoint minors shared by two or three cube cuts. It
does not use the Cayley form, commutation, or the reciprocal Smith-order
identity of the preceding construction.

## 1. Three coupled determinantal budgets

Let `K(x,y,z)` be an alternating eight-by-eight matrix over a field of
characteristic two. Assume every nonzero entry is homogeneous of degree
`D>0`, and `K(1,1,1)` is the antipodal matching on cube vertices.
At axis `a`, let the four Smith orders of its one-particle crossblock be

\[
0\le s_{a,1}\le s_{a,2}\le s_{a,3}\le s_{a,4}.
\]

The orders are taken at that axis variable, over the rational function
field in the other two variables. The endpoint makes every crossblock
generically invertible.

**Proposition.** These twelve orders satisfy

\[
\sum_a s_{a,1}\le D,\qquad
\sum_a(s_{a,1}+s_{a,2})\le3D,\qquad
\sum_a\sum_{j=1}^4s_{a,j}\le12D.
\tag{1}
\]

For the first inequality, take the single entry `K_(000,111)`. It is
nonzero and has degree `D`, because its endpoint is one. Its valuations
at `x,y,z` sum to at most `D`. The smallest Smith order at each cut is
at most that entry's valuation.

For the second inequality, fix two axes `a,b`, and let `c` be the third
axis. Take the two antipodal pairs with representatives `000,e_c`.
At both cuts `a,b`, the same two-by-two minor occurs, with rows
`000,e_c` and columns `111,111 xor e_c`. This minor has endpoint one,
so it is a nonzero homogeneous polynomial of degree `2D`. Therefore

\[
(s_{a,1}+s_{a,2})+(s_{b,1}+s_{b,2})\le2D.
\tag{2}
\]

Sum (2) over the three pairs of axes and divide by two. This proves the
middle inequality in (1).

Finally, each complete crossblock determinant is a nonzero homogeneous
polynomial of degree `4D`. Its valuation at its own variable is at most
`4D`. Summing gives the last inequality in (1).

## 2. The most optimistic profile product allowed by the budgets

For `0<t<1`, put `f(s)=log(1+t^s)`. This function is decreasing and
convex. Split the twelve orders into the three smallest orders, the
three second orders, and the six remaining orders. Applying Jensen's
inequality inside these three groups, followed by (1), gives

\[
\sum_{a,j}\log(1+t^{s_{a,j}})
\ge3\log(1+t^{D/3})
  +3\log(1+t^{2D/3})
  +6\log(1+t^{3D/2}).
\tag{3}
\]

One direct way to finish the balancing step is to use the three group
means `u<=v<=w`. Their constraints are `u<=D/3`, `u+v<=D`, and
`u+v+2w<=4D`. Monotonicity permits increasing the last total to equality.
Convexity then permits transferring mass from a larger mean to a
smaller mean until `u=D/3` and `v=2D/3`; the remaining mean is `3D/2`.
Every such transfer decreases the displayed convex sum. This gives (3).

Equivalently, after putting `r=t^(D/6)`, the optimistic profile for one
cut has counting function

\[
Z_*(r)=(1+r^2)(1+r^4)(1+r^9)^2.
\tag{4}
\]

This is a relaxation: there need not be a kernel attaining all these
orders. Allowing it can only make a proposed common-profile bound more
optimistic.

## 3. An all-scale rational separation from the half target

For the common physical frequency `p`, homogeneity fixes the total
coefficient degree to `4Dpm` over `m` sites. The same coefficient-count
calculation as in the constructive note shows that the common-profile
method's proposed contraction is bounded below by

\[
R(p,r)=
\frac{\log_2 Z_*(r)-8p\log_2r}{4H_2(p)}.
\tag{5}
\]

For every `0<r<=1`,

\[
Z_*(r)^5>(1+r^{10/3})^{12}.
\tag{6}
\]

The following exact positive-coefficient certificate proves (6). Put
`r=s^3` and let

\[
P(s)=[(1+s^6)(1+s^{12})(1+s^{27})^2]^5-(1+s^{10})^{12}.
\]

The polynomial `P` vanishes to order six at zero. The quotient has
degree 354, and every coefficient of

\[
(1+u)^{354}\,(P(s)/s^6)\big|_{s=u/(1+u)}
\tag{7}
\]

is a positive integer. There are 355 coefficients, and the smallest is
five. Positivity on `0<s<1` follows; at `s=1`, it follows by direct
evaluation. The checker retains every coefficient, so this step is
independent of floating-point optimization or root finding.

For `rho=3/5`, the binary entropy identity gives

\[
4\rho H_2(p)+8p\log_2r
 \le4\rho\log_2(1+r^{2/\rho})
 =\frac{12}{5}\log_2(1+r^{10/3}).
\]

Combining this with (6) proves

\[
R(p,r)>\frac35>\frac1{\sqrt3}
\quad(0<p<1,\ 0<r\le1).
\tag{8}
\]

Thus the scoped method cannot pass the cube half-score threshold, even
after allowing arbitrary homogeneous non-Cayley kernels.

## 4. Why multiple physical modes remain a different construction problem

With several binary modes on each physical leg, a complete crossblock
has more than four one-particle orders. Physical words can have unequal
mode heights, so their entropy is governed by a larger graded alphabet.
The particular three-group relaxation in section 2 and the binary
entropy denominator in (5) no longer describe that problem.

The endpoint-minor idea still gives useful constraints on such kernels,
but it does not establish the floor (8) for them. This is the reason to
search coupled weighted physical modes instead of continuing an
unrestricted eight-by-eight homogeneous kernel search.

## 5. Checks

The [checker](../tools/check_cube_gaussian_profile_budget.py) and
[retained certificate](../evidence/cube-gaussian-profile-budget.json)
evaluate all crossblock minors of 32 seeded, generally non-Cayley
homogeneous kernels of degrees one through four. They check the shared
endpoint minor for every pair of axes and all three budgets. They also
retain and verify all 355 positive integer coefficients in (7).

Run `python3 tools/check_cube_gaussian_profile_budget.py` to reproduce
the exact controls. The arbitrary-kernel statement follows from the
endpoint-minor argument, not from completeness of the finite sample.
