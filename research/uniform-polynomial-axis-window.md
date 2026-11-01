# A universal axis window for strict polynomial crossings

## Result and constructive consequence

The [fixed-gauge theorem](polynomial-crossing-axis-obstruction.md) admits
a stronger uniform statement. At fixed physical alphabet sizes, *no*
strict Laurent-polynomial crossing improves the critical support by
approaching either axis, even if its coefficients, degree, physical
heights and thermal scale all change with the angle.

Put `n=A*B>=4` and

```
X=1024*n^2,        epsilon_n=exp(-2*X).
```

For every coefficient field, every finite strict Laurent-polynomial
crossing on alphabets X,Y,Y,X of sizes A,B,B,A, all real monomial physical
heights, and every positive scale, the critical margin M is **strictly
positive** whenever either positive norm-one supporting weight is at most
`epsilon_n`. The constant is deliberately coarse. Its role is to prove
uniformity, not to set a useful numerical search cutoff.

This closes the previously stated moving-height/tensor escape *at the
axes* for this strict polynomial pressure framework. Interior weights,
other coefficient-rank estimates, and rational paths with additional
finite poles remain outside the conclusion. General type-only endpoints
still require their own matched-minor hypotheses. The complete factoring
score stays `1/2`; no universal factoring lower bound is asserted.

The research question was whether moving heights could defeat a
fixed-height exclusion. The new mechanism is analytic: full endpoint
majorization turns the shared second-prefix deficit into a quantified
partition gap, while the cold-tail bound controls the remaining physical
scale. Thus an axis search cannot be rescued merely by larger height
bounds or by changing the polynomial tensor as its angle tends to zero.

## 1. An abstract vector statement

It suffices to prove the following statement for real increasing lists
r,s,e of length n>=2. Suppose, writing `delta=s_2-s_1`,

```
sum_(i<=k)(r_i+s_1) <= sum_(i<=k)e_i,  1<=k<=n,
sum_(i<=2)(r_i+s_1) <= e_1+e_2-delta.                   (1)
```

For positive norm-one a,b and t>0 define

```
M=a*log(sum exp(-t*r_i/a))+b*log(sum exp(-t*s_i/b))
  -log(sum exp(-t*e_i)).
```

Then M>0 for `b<=exp(-2048*n^2)`. Strict polynomial crossings satisfy
(1), with e the sorted additive physical pair energies, by the credited
endpoint and shared-minor proof. It applies in both orientations, proving
the claim for either small weight. Real heights follow by continuity of
the finite weighted-minor minima.

Absorb t into the lists and subtract e_1 from e and from
`lambda_i=r_i+s_1`. This leaves
`M=a log sum exp(-lambda_i/a)+b log sum exp(-(s_i-s_1)/b)-Phi(1)`
unchanged. Now `e_1=0`; put `x=e_2>=0` and continue to write
delta for the scaled weak gap. Equation (1) is unchanged.
Let `Phi(v)=log sum exp(-v*e_i)` and let H be the entropy of its
unit-scale Gibbs distribution. Then

```
H<=log n,        H<=(n-1)*(1+x)*exp(-x).                 (2)
```

For the second inequality, the ground probability is at least
`1/(1+(n-1)*exp(-x))`. At fixed ground probability, symmetrizing all
other probabilities increases entropy. The resulting entropy decreases
as the ground probability increases above 1/n. It is therefore at most
the entropy of energies `(0,x,...,x)`, namely
`log(1+(n-1)e^-x)+x*(n-1)e^-x/(1+(n-1)e^-x)`. The two elementary
bounds in (2) follow. This also covers x=0 and degenerate grounds.

## 2. Prefix deficit gives a partition deficit

Put `d_i=lambda_i-e_i`. All its prefix sums are nonpositive, and
`d_1+d_2<=-delta`. Hence one of its first two coordinates is at most
`-delta/2`. For decreasing positive weights `w_i=exp(-e_i/a)`, Abel
summation gives `sum w_i*d_i<=0`. The exact identity

```
sum exp(-lambda_i/a)-sum exp(-e_i/a)
 = -(1/a)*sum w_i*d_i
   +sum w_i*(exp(-d_i/a)-1+d_i/a)
```

has nonnegative terms in its second sum. The identified first-two term
has weight at least `exp(-x/a)`. Since `exp(z)-1-z>=z^2/2` for z>=0,

```
log(sum exp(-lambda_i/a))-Phi(1/a)
 >= log(1+exp(-x/a)*delta^2/(8*n*a^2)).                  (3)
```

Here `sum exp(-e_i/a)<=n` was used. No equality of total prefix sums,
simple cut minimum, bounded gauge or bounded tensor degree is needed.

The perspective derivative of `a*Phi(1/a)` is the corresponding Gibbs
entropy. Entropy decreases under cooling, so

```
a*Phi(1/a)-Phi(1) >= -(1-a)*H.
```

Keeping two weak terms and using (3) proves the key lower bound

```
M >= a*log(1+z)+b*log(1+exp(-delta/b))-(1-a)*H,
z=exp(-x/a)*delta^2/(8*n*a^2).                          (4)
```

## 3. One cold reduction

Set `c=a+b` and `D=-lambda_1>=0`. The second prefix in (1) gives
`(r_2-r_1)+(s_2-s_1)<=x+2D`. The two-gap convexity argument in
[the cold-tail theorem](uniform-cold-physical-gap.md) therefore applies
also to the abstract physical list e. With the coarse count n-1, it gives

```
M>0 if x >= c/(c-1)*log(2*(n-1)/c).                    (5)
```

If x=0 or (5) does not apply, then for b<=1/2,

```
x < (3/b)*log(2*n),   unless x=0.
```

Indeed `c<=3/2` and `c-1>=b/2`. With `b<=exp(-2X)`, this implies
`x*b^2<=1`. Also `a>=3/4`, `1/a-1<=b^2`, and `1-a<=b^2`.
These deliberately weak constants suffice below.

## 4. The remaining scales and weak gaps

If z>=1 in (4), its first term exceeds 1/4, while
`(1-a)H<=b^2 log n<1/4`. Thus M>0.
If z<1, `log(1+z)>=z/2`, and `x*b^2<=1` gives

```
a*log(1+z) >= delta^2*exp(-x/a)/(16*n*a)
             >= delta^2*exp(-x)/(48*n).                 (6)
```

Use (2) for the negative term in (4).

- If `x<=X` and `delta<=bX`, the weak term exceeds
  `(b/2)*exp(-X)`, which exceeds `b^2 log n` because
  `b<=exp(-2X)` and `exp(X)>2 log n`.
- If `x<=X` and `delta>bX`, (6) exceeds
  `b^2*X^2*exp(-x)/(48n)`, larger than
  `b^2*(n-1)*(1+x)*exp(-x)`. The exact comparison is
  `X^2>48n(n-1)(1+X)`.
- If `x>X` and `delta<=bx/2`, the weak term is at least
  `(b/2)*exp(-x/2) >= (b/4)*(1+x)*exp(-x)`, larger than
  `b^2*(n-1)*(1+x)*exp(-x)` since `b(n-1)<1/4`.
- If `x>X` and `delta>bx/2`, (6) exceeds
  `b^2*x^2*exp(-x)/(192n)`, larger than that same negative
  bound because `x^2/(1+x)>X/2>192n(n-1)`.

All comparisons hold for `X=1024n^2`, n>=2. Each case makes the
right side of (4) strictly positive. Together with (5) this proves the
abstract vector statement and the uniform crossing theorem.


## Controls and credit

The [checker](../tools/check_uniform_axis_window.py) tests the partition
deficit and perspective formulas on the earlier 408 independent
untruncated crossing orders, replayed from their pinned receipt without
recomputing them. It makes 576 source comparisons. Separately, 512
explicitly unrealized abstract vector cases exercise full majorization,
the shared deficit, 6,144 partition inequalities and 512 exact rational
Abel sums. The [receipt](../evidence/uniform-axis-window-controls.json)
records source hashes and the coarse constants. These are ingredient
controls. The extremely small window is proved by the analytic case
split, not sampled in floating point; abstract vectors are not tensor
constructions. Independent mathematical acceptance remains open.

```
make verify-uniform-axis-window
```

This strengthens the repository's credited shared-minor and fixed-gauge
axis results using the new cold-tail lemma and a quantified Bregman
partition deficit. Directed by Cinna Davis / 1337isnot1337 with
substantial OpenAI assistant support. No new Lean, external novelty
claim, or complete factoring improvement is asserted.

The subsequent [approximate rational strictification](approximate-rational-strictification.md)
also excludes a negative axis margin for scalar rational kernels over finite
fields when every paired-type word is globally Laurent-regular. That
corollary charges and removes their extra poles. Arbitrary rational poles,
selected-type-only constructions and matrix auxiliaries are not assumed
to satisfy those hypotheses.
