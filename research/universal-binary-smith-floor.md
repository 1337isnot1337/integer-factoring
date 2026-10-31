# Every binary polynomial seed has positive two-place Smith pressure

## Result and scope

Let `K` be any field. Let `T(w;x,y,u,v)` be a Laurent-polynomial tensor,
with all four physical letters binary, and with the strict endpoint

\[
T(1;x,y,u,v)=[x=v][y=u].
\tag{1}
\]

Use the cost cuts `(x,y)|(u,v)` at zero and `(x,u)|(y,v)` at infinity.
Allow every independent physical-letter height at both places, every
polynomial degree, and every regular local change of basis. The complete
two-place Smith pressure is strictly positive at the critical parameter
`rho=1/sqrt(2)`, even after optimizing the two physical types independently.

Thus **no binary Laurent-polynomial seed can improve the factoring
exponent through this two-place Smith-profile criterion**. This includes
non-Gaussian seeds, arbitrary coefficient fields, and unbounded regular
row frames. It replaces the degree and flag restrictions of the earlier
binary searches with an all-seed theorem.

The theorem is about this sufficient crossing certificate. It is not a
factoring lower bound or an obstruction to every tensor algorithm.
Larger physical alphabets, auxiliary word modules with stronger rate
bounds, multipoint reconstruction, and rational kernels whose polynomial
endpoint emerges only after type projection are outside its hypotheses.
The established factoring score remains `1/2`.

## 1. Pressure and harmless normalizations

Subtract the height at letter zero on each leg and at each place, so
`h_(e,l)(0)=0`, and write `h_(e,l)=h_(e,l)(1)`. Put

\[
H_l=h_{0,l}+h_{\infty,l},\qquad
E_X=H_x+H_v,\qquad E_Y=H_y+H_u.
\tag{2}
\]

Let the sorted Smith orders be `sigma_0` and `sigma_infinity`, each with
four entries. For `t>0`, define

\[
B(t)=\left(\sum_i e^{-t\sigma_{0,i}}\right)
     \left(\sum_i e^{-t\sigma_{\infty,i}}\right).
\tag{3}
\]

The [independent paired-type criterion](paired-type-crossing-pressure.md)
has critical pressure

\[
G(t)=\log B(t)-\sqrt2\left[
 \log(1+e^{-tE_X/\sqrt2})+
 \log(1+e^{-tE_Y/\sqrt2})\right].
\tag{4}
\]

A post-half certificate would require `G(t)<0`. Using the complete Smith
profiles gives at least as strong a test as replacing them by sufficient
row-lattice bounds.

Adding a constant to a physical height shifts every Smith order at that
place by the same constant. Its contribution to the physical term in
(4) changes by exactly the same amount, so `G` is unchanged. Simultaneously
flipping the letters of a matched pair preserves (1), permutes the
matrices' rows and columns, and changes its energy difference's sign.
Consequently we may assume

\[
E_X=-2a,\qquad E_Y=-2b,\qquad a,b\ge0.
\tag{5}
\]

The individual `H_l` can still have either sign. No ordering or
regularity restriction on the physical heights has been introduced.

We first prove the algebraic inequalities for integer heights. Rational
heights follow by ramification of `w`; real heights follow from continuity
of the determinantal-prefix functions, as explained in Section 6.

## 2. Normalize the entire tensor using its two minimum orders

Let `mu_0` and `mu_infinity` be the actual minimum orders at the two
places. They are the minimum entry valuations, so every nonzero Laurent
entry is supported between its zero and infinity entry bounds. Set

\[
d=\sum_l H_l-\mu_0-\mu_\infty
  =-2(a+b)-\mu_0-\mu_\infty.
\tag{6}
\]

For a binary four-vector `lambda`, complement every bit and define

\[
R_\lambda(w)=
w^{-\mu_0+\sum_l h_{0,l}(1-\lambda_l)}
T_{\mathbf1-\lambda}(w).
\tag{7}
\]

The entry bounds give the exact support statement

\[
R_\lambda\in K[w],\qquad
\deg R_\lambda\le d-\sum_lH_l\lambda_l.
\tag{8}
\]

A negative upper bound in (8) forces that entry to vanish. Since
`R_0000(1)=1`, its degree bound proves `d>=0`. The endpoint of `R` is
again the literal crossing.

At zero use the unweighted matrix `R(w)`. At infinity use the matrix
with entries

\[
z^{d-\sum_lH_l\lambda_l}R_\lambda(z^{-1}).
\tag{9}
\]

Both matrices are polynomial. Their sorted Smith profiles are

\[
r_i=\sigma_{0,i}-\mu_0,\qquad
s_i=\sigma_{\infty,i}-\mu_\infty,
\qquad r_1=s_1=0.
\tag{10}
\]

Every determinant term in either four-by-four matrix uses each physical
letter one exactly twice. Its degree is therefore at most

\[
C=4d-2\sum_lH_l=4d+4(a+b).
\tag{11}
\]

The determinants are nonzero because their values at one are permutation
determinants. A nonzero polynomial's zero order cannot exceed its degree,
so

\[
\sum_i r_i\le C,\qquad \sum_i s_i\le C.
\tag{12}
\]

## 3. A minor shared by both cuts supplies the missing restriction

Fix `y=u=0`, and let `x,v` vary. The resulting two-by-two determinant
`Delta_X(w)` is a minor of both cost flattenings. Its value at one is
one, and (8) gives

\[
\deg\Delta_X\le2d-H_x-H_v=2d+2a.
\tag{13}
\]

At infinity the same minor is
`z^(2d+2a) Delta_X(z^(-1))`. The sum of the orders at the two places
equals

\[
2d+2a+\operatorname{ord}_0\Delta_X-\deg\Delta_X\le2d+2a.
\tag{14}
\]

The sum of the first two Smith orders is the minimum valuation among
all two-by-two minors. Since `r_1=s_1=0`, (14) implies
`r_2+s_2<=2d+2a`. Fixing `x=v=0` gives the analogous bound with `b`.
Together,

\[
\boxed{r_2+s_2\le2d+2\min(a,b).}
\tag{15}
\]

This restriction couples the two cuts. Separate determinant bounds
alone do not supply it. It holds over every field, irrespective of
local flags or coefficient degrees.

For the formerly searched height ray
`h_0=(-1,-1,0,0)`, `h_infinity=(0,0,-1,-1)`, the same argument before
normalization says

\[
\sigma_{0,1}+\sigma_{0,2}
+\sigma_{\infty,1}+\sigma_{\infty,2}\le-6.
\tag{16}
\]

One can see (16) directly by fixing `y=u=1`: the shared determinant
receives height minus three at each place, and its Laurent support has
nonnegative width. The cleared target convention used by the binary
degree-twelve searches turns (16) into a prefix bound of 26.

## 4. Reduce all profile choices to three nonnegative parameters

The log-sum-exponential is convex. Averaging the two normalized profiles
therefore decreases their sum of log partition functions. Let
`c=(r_2+s_2)/2` and `m=min(a,b)`. Equations (12)--(15), followed by
convexity for the last two coordinates, give

\[
\log B(t)\ge [2(a+b)+d]t+
2\log\left(1+e^{-ct}+2e^{-(C-c)t/2}\right).
\tag{17}
\]

Here `0<=c<=d+m<=C/3`. The expression inside the logarithm is
nonincreasing in `c` on this interval, because its derivative is
`t(e^(-(C-c)t/2)-e^(-ct))`. Substituting `c=d+m` is consequently a
valid lower bound.

Put `f(z)=log(1+e^(-z))` and

\[
A=dt,\qquad X=at,\qquad Y=bt.
\]

The leading physical term `2(a+b)t` cancels in (4). Thus

\[
\begin{aligned}
G(t)&\ge\Phi(A,X,Y),\\
\Phi(A,X,Y)
&=A+2\log\left(1+e^{-A-\min(X,Y)}
 +2e^{-[3A+4(X+Y)-\min(X,Y)]/2}\right)\\
&\quad-\sqrt2\,[f(\sqrt2X)+f(\sqrt2Y)].
\end{aligned}
\tag{18}
\]

All three parameters in (18) are nonnegative. No finite target catalogue
or numerical choice of Smith orders remains in this bound.

## 5. Unequal physical energies reduce to the balanced case

Assume `X>=Y`, and write `q=e^(-X)`, `u=e^(-Y)`, so `0<q<=u<=1`.
For `x=e^(A/2)>=1`, equation (18) becomes

\[
\Phi=2\log\left(x+\frac{u}{x}
        +\frac{2u^{3/2}q^2}{x^2}\right)
-\sqrt2\log(1+q^{\sqrt2})-\sqrt2\log(1+u^{\sqrt2}).
\tag{19}
\]

Set `y=u/x^2`, with `0<y<=u`. The quantity inside the logarithm is
`sqrt(u) H_q(y)`, where

\[
H_q(y)=y^{-1/2}+y^{1/2}+2q^2y.
\tag{20}
\]

Its derivative changes sign exactly once. The unique unrestricted
minimizer `y_*(q)` in `(0,1)` solves

\[
y_*+4q^2y_*^{3/2}=1.
\tag{21}
\]

Let `M_q(u)=inf_(A>=0) Phi(A,-log(q),-log(u))`. If `u>=y_*`,
the minimizer is interior. The only `u`-dependent part of `M_q` is
`log u-sqrt(2) log(1+u^sqrt(2))`, whose derivative is

\[
\frac{1-u^{\sqrt2}}{u(1+u^{\sqrt2})}\ge0.
\tag{22}
\]

If `u<=y_*`, the minimizer is `A=0`. Its derivative has the sign of

\[
\begin{aligned}
&(1+3\sqrt u\,q^2)(1+u^{\sqrt2})
-u^{\sqrt2-1}(1+u+2u^{3/2}q^2)\\
&\qquad=1-u^{\sqrt2-1}+3\sqrt u\,q^2
                  +u^{\sqrt2+1/2}q^2\ge0.
\end{aligned}
\tag{23}
\]

The two formulas agree at `u=y_*`. Hence `M_q` is nondecreasing in
`u`, and its minimum on `q<=u<=1` occurs at `u=q`. We have proved

\[
\Phi(A,X,Y)\ge\inf_{A'\ge0}F(A',\max(X,Y)),
\tag{24}
\]

where the balanced function is

\[
F(A,t)=A+2\log\left(1+e^{-A-t}+2e^{-(3A+7t)/2}\right)
             -2\sqrt2\log(1+e^{-\sqrt2t}).
\tag{25}
\]

## 6. Exact positivity of the remaining two-variable function

The following explicit bound holds for every `A,t>=0`:

\[
\boxed{F(A,t)>e^{-t}/28.}
\tag{26}
\]

For the tail, put `u=e^(-t)<=1/4` and `x=e^(A/2)>=1`. Then

\[
F(A,t)=2\log\left(x+u/x+2u^{7/2}/x^2\right)
                  -2\sqrt2\log(1+u^{\sqrt2}).
\tag{27}
\]

The first logarithm's argument is at least `1+u`. Also,
`log(1+u)>=7u/8`, while

\[
\sqrt2\log(1+u^{\sqrt2})<\frac{10}{7}u\,u^{2/5}
<\frac67u.
\]

The last inequality follows from `(5/3)^5<16`. Subtracting gives
`F>2(7/8-6/7)u=u/28`.

For `u>=1/4`, differentiate the argument in (27) with respect to `x`.
When `x>=2`, its derivative is at least
`1-1/4-1/2=1/4`. Its minimum therefore occurs at `1<=x<=2`.
Both `A` and `t` can be restricted to `[0,log 4]`, contained in
`[0,3/2]`. An exact dyadic interval certificate proves

\[
F(A,t)>1/16\qquad\text{throughout }[0,3/2]^2.
\tag{28}
\]

The [certificate](../evidence/binary-fixed-ray-floor.json) contains 421
positive rectangular leaves. Its independent cover audit checks every
split and every lower endpoint using integer arithmetic with 80 binary
fractional digits. The logarithm and exponential bounds use explicit
series remainders. Removing a leaf or modifying a saved lower bound is
rejected. The compact margin `1/16` exceeds `u/28`, proving (26).

For each fixed `t`, the function `F(A,t)` is continuous and tends to
infinity as `A` tends to infinity, so its infimum over `A>=0` is
attained. Combining (18), (24), and (26) yields the quantitative theorem

\[
\boxed{G(t)>\frac{e^{-\max(a,b)t}}{28}>0.}
\tag{29}
\]

For rational physical heights, replace `w` by a fixed integer power and
multiply every height by the same integer. Smith orders and `a,b,d`
scale by that integer, while `t` scales inversely. The proof applies
without change. For real heights, each determinantal prefix is the
minimum of finitely many affine height functions: row and column
weights factor from each fixed minor. It is therefore continuous.
The closed inequalities (10)--(18) extend by continuity, after which
the strict analytic argument (19)--(29) applies directly. Strictness
is not inferred merely by taking a limit of strict inequalities.

## 7. Literal controls and research consequences

The [fixed-ray checker](../tools/check_binary_fixed_ray_floor.py) constructs
the exact compact cover and compares the normalization with full
determinantal profiles of 72 colored polynomial tensors through degree
eight, including Laurent clearing. It checks 2,304 literal endpoint
entries, 288 complete profiles, and 216 outward pressure comparisons.
A negative control shows that deleting the last two rank channels would
make the purported analytic floor false.

The [all-gauge checker](../tools/check_binary_all_gauge_floor.py) independently
expands every minor for 84 further tensors with arbitrary integer leg
heights of both signs. It checks pair flips, full support intervals,
both shared minors, determinant totals, and the general pressure bound.
Its [record](../evidence/binary-all-gauge-floor.json) contains 2,688 endpoint
entries, 504 complete profiles, and 252 outward pressure comparisons.
These finite controls check the transcription; the general proof is
the argument above.

```sh
python3 tools/check_binary_fixed_ray_floor.py --output /tmp/binary-fixed-ray-replay.json
python3 tools/check_binary_all_gauge_floor.py --output /tmp/binary-all-gauge-replay.json
```

This result supersedes the remaining binary higher-jet and formal-frame
searches as routes to a better bare two-place Smith certificate. Those
searches were stopped after the shared-minor restriction was proved;
their unfinished records are not complete catalogues. The constructive
interfaces remain applicable when their hypotheses are changed, while
the active score search should use larger alphabets or a stronger
coefficient-bound mechanism.
