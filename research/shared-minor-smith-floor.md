# Shared minors strengthen the all-degree floor through six letters

## Result and scope

For the same bare, ungauged polynomial crossing interface as the
[degree-only theorem](all-degree-small-alphabet-smith-floor.md), every
field, arbitrary polynomial degree, and arbitrary fixed auxiliary
dimension satisfy

```math
\sqrt{B(s)}\ge
\min_{x\ge1}\left[x^{q-1}+\frac{q-1}{x}
                         +\frac{q(q-1)}{x^q}\right].       \tag{1}
```

For each physical alphabet size 2<=q<=6, exact rational certificates prove

```math
B(s)>q^{20/7},\qquad
\rho_{\rm bare}>5/7>1/\sqrt2.
```

Thus merely increasing the degree, changing the field, or adding fixed
auxiliary states cannot yield a half-rate improvement on these alphabets
through this bare-profile criterion. This strengthens the earlier
ungauged exclusion q<=4 to q<=6. It is not a factoring lower bound;
independent gauges, coefficient recombinations, cancellation refinements,
and larger local alphabets remain outside the conclusion.

The common-minor principle is used for scalar binary crossings with
independent heights in [PR #39's universal binary theorem](https://github.com/1337isnot1337/integer-factoring/blob/0ac85be63d6eff434eddd92a32eb5c7e4ee75cb2/research/universal-binary-smith-floor.md).
Here it becomes a q-by-q minor, with a rectangular-minor argument for
fixed auxiliary states and a simple one-variable ungauged optimization.
This ordinary proof used substantial OpenAI Codex assistance. No new
Lean proof or stronger factoring exponent is claimed.

## 1. One minor belongs to both cuts

Remove a common monomial as before, and let M be the resulting actual
degree. Let mu and nu be the increasing infinity and zero Smith orders
of the two rectangular cut matrices, regularizing infinity by t^M.
Their minimum orders are zero and their sums are at most nM, where n=q^2.

Contract the auxiliary indices with the fixed boundary vectors that
give the crossing at parameter one. Fix y=u=a at any physical letter a,
and allow x,v to vary. The resulting q-by-q scalar polynomial matrix
is a submatrix of **both** scalar cut flattenings. Its determinant Delta
has Delta(1)=1 and degree at most qM.

For the rectangular auxiliary cuts, the scalar contraction is a constant
linear map on the auxiliary side. Cauchy--Binet expresses every minor
after contraction as a constant linear combination of original minors.
The minimum q-minor valuation therefore cannot exceed the valuation of
this scalar shared minor. The same argument applies at infinity. Hence

```math
\sum_{i=1}^q\nu_i\le\operatorname{ord}_0\Delta,\qquad
\sum_{i=1}^q\mu_i\le qM-\deg\Delta.
```

Because the nonzero polynomial Delta has order at most its degree,

```math
\boxed{\sum_{i=1}^q(\mu_i+\nu_i)\le qM.}                  \tag{2}
```

The two auxiliary indices may be placed on different physical halves in
the two cuts, as long as they remain together within each cut. The fixed
boundary contraction is still a constant map on that cut's auxiliary side.

## 2. Collapse the profile optimization

Assume M>0, put t=M log s, and set a_i=(mu_i+nu_i)/(2M). Then

```math
a_1=0,\qquad \sum_i a_i\le n,\qquad \sum_{i=1}^q a_i\le q/2.
```

Cauchy--Schwarz gives

```math
\sqrt{B(s)}\ge e^{t/2}\sum_{i=1}^n e^{-ta_i}.             \tag{3}
```

Write b=sum_(i=2)^q a_i. Jensen's inequality on the two remaining groups
of sizes q-1 and n-q gives

```math
\sum_i e^{-ta_i}\ge
1+(q-1)e^{-tb/(q-1)}+(n-q)e^{-t(n-b)/(n-q)}.               \tag{4}
```

The right side decreases in b throughout 0<=b<=q/2: its derivative has
the sign of -exp(-tb/(q-1))+exp(-t(n-b)/(n-q)), which is nonpositive
because b<=q/2<n(q-1)/(n-1). Thus substitute b=q/2. Using n=q^2 and
x=exp(t/(2(q-1))) turns (3)--(4) into

```math
\sqrt{B(s)}\ge x^{q-1}+(q-1)x^{-1}+q(q-1)x^{-q}.
```

This proves (1). If M=0, the two cut matrices are constant full-rank
matrices, B=n^2, and the same bound at x=1 is equality.

## 3. Exact global threshold certificate

Define U_q(x) to be the expression in (1). The checker proves U_q(x)>=T_q
for all x>=1 with these rational numbers:

| q | T_q |
| --- | --- |
| 2 | 297/100 |
| 3 | 526/100 |
| 4 | 773/100 |
| 5 | 1029/100 |
| 6 | 12933/1000 |

To check the compact interval [1,3], multiply U_q(x)-T_q by x^q and
clear denominators. This is a polynomial of degree 2q-1. Substitute
x=1+2z and certify its nonnegativity on 0<=z<=1 with exact rational
Bernstein subdivision. For x>=3 the first term x^(q-1) alone exceeds
T_q. There is no uncovered interval or numerical root search.

For T_q=p/r, the exact integer comparison p^14>q^20 r^14 proves
T_q^2>q^(20/7). Finally 2*5^2>7^2 proves 5/7>1/sqrt(2).
The [checker](../tools/check_shared_minor_smith_floor.py) and
[saved record](../evidence/shared-minor-smith-floor.json) reconstruct all
five compact certificates. They need at most eight leaves and depth seven.

The approximate optimized rho floors for q=5 and q=6 are respectively
0.724456775 and 0.714328746. The rational floor 5/7 is slightly weaker
and is used because its strict comparison is entirely exact.

## 4. What the abstraction leaves open

For q=7 the relaxed bound reaches approximately rho=0.706230448, below
1/sqrt(2). This is not a construction: the degree/valuation inequalities
are necessary conditions, and an algebraic tensor need not realize their
optimizers. It identifies the first local size left open by this specific
ungauged argument.

The formal equality pattern at q=7 has degree M=12 and both normalized
Smith profiles (0,7 repeated six times,13 repeated 42 times). Each profile
has 49 entries and total 588=49M, and their seven-entry prefixes saturate
(2). The [inverse-degree compatibility theorem](inverse-degree-crossing-compatibility.md)
subsequently rules out this exact pattern at every integer ramification:
its second infinity order is too large relative to the polynomial inverse
degree. Thus this relaxation does not certify that an algebraic witness
exists at seven letters. Other seven-letter profiles remain outside the
exclusion. A correct column Jordan spectrum alone is insufficient.

Run `python3 tools/check_shared_minor_smith_floor.py` to reproduce the
exact all-scale exclusion through q=6. The next constructive work must
cross an actual structural boundary, rather than extend a polynomial
degree catalogue within this closed certificate class.
