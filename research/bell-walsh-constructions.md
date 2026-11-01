# Quadratic-parity Walsh paths and the full radial Bell obstruction

## Results and scope

There is an explicit Bell crossing for every positive integer energy table
whose parity is a nondegenerate binary quadratic form. The energies need
not come from a code weight. In particular, a degree-two family exists at
every physical dimension `d=2^k`. For `k>=3` this lies outside the
positive-code-weight quotients, whose reduced degree is at least `k`.
The degree-two family itself has poor rank pressure; it is not a factoring
improvement.

There is also a strict constructive improvement of a fixed folded seed.
Moving the five coefficients of the physical-dimension-16 seed to degrees
`(0,22,40,58,80)` preserves its crossing endpoint and reciprocal identity.
Its complete all-length pressure satisfies
`rho_retimed < 64/75 < rho_original`. Section 7 gives exact rational
certificates and the general optimal-retiming formula. This does not lower
the asymptotic radial optimum or improve the established factoring score.
In fact every fixed folded seed with `k>=3` admits some strict pressure
improvement by reciprocal retiming; the dimension-16 schedule is a small
explicit rational certificate of this general statement.

The radial part of this enlarged family has a complete obstruction. Put
`m=2k` and use the common Hamming-coordinate identification specified below.
In every decomposition of a nonzero multiple of the crossing into tensors
from the **full radial Bell space**, some one tensor has both cut ranks at
least

\[
\binom{m}{\lfloor m/4\rfloor}.
\tag{1}
\]

This holds over every field of odd characteristic and in characteristic
zero, with any number of summands and arbitrary coefficient recombinations.
The optimal leading parameter within this whole space is therefore

\[
\rho_{\rm radial}=H_2(1/4)
=2-\tfrac34\log_2 3=0.811278124459\ldots.
\tag{2}
\]

The earlier [folded-cube result](folded-bell-crossings.md) proved this limit
only within its particular `k+1`-coefficient span. The new lower bound permits
the entire `2k+1`-dimensional radial space, mixtures of radial paths,
nonmonotone radial energies, and modular cancellations. The existing folded
construction attains (2), so the enlarged radial space has the same optimum.

Nonradial energy tables, different common physical identifications, and
general Bell paths remain outside this obstruction. Finite searches below
have not improved the established factoring score `1/2`. The results are
ordinary proofs and exact computational controls, with substantial OpenAI
assistant support. There is no new Lean formalization or independent review.

## 1. Arbitrary spectral energies preserve a strict crossing

Let `V=F_2^(2k)`, `d=2^k`, and write an index as `u=(a,b)` with `a,b` in
`F_2^k`. Fix a binary vector `g`, and define

\[
q_g(a,b)=a\cdot b+g\cdot(a,b)\pmod2.
\]

Choose any integer energy function satisfying

\[
e(0)=0,\qquad e(u)>0\ (u\ne0),\qquad e(u)\equiv q_g(u)\pmod2.
\tag{3}
\]

Set `h=max e`, and let `H` be the unnormalized Walsh matrix on `V`, with
entries `H_(v,u)=(-1)^(v dot u)`. Define the polynomial vector

\[
f_v(z)=2^{-h}\sum_{u\in V}(-1)^{v\cdot u}
 (1+z)^{h-e(u)}(1-z)^{e(u)}.
\tag{4}
\]

This uses only powers of two as denominators. In any odd characteristic,
these are units. At `z=1` all terms except `u=0` vanish, so `f_v(1)=1`
for every `v`.

For `g=0`, the regular Bell transform from
[the Bell construction](bell-diagonal-crossings.md) is

\[
T=d^{-2}H\,\mathrm{diag}((-1)^{a\cdot b})\,H.
\tag{5}
\]

Indeed, the Walsh transform of the binary hyperbolic quadratic phase is
`d` times the same phase. An extra linear phase in the middle of (5) is
translation by `g` after the transform. Thus, writing `P_g` for translation,

\[
S=P_gT=TP_g
=d^{-2}H\,\mathrm{diag}((-1)^{q_g(u)})\,H.
\tag{6}
\]

Each scalar polynomial in (4) has reciprocal identity

\[
z^h(1+z^{-1})^{h-e}(1-z^{-1})^e
=(-1)^e(1+z)^{h-e}(1-z)^e.
\]

Equations (3)--(6), together with `H^2=d^2 I`, prove

\[
\boxed{\quad Sf(z)=z^h f(1/z),\qquad f(1)=\boldsymbol1.\quad}
\tag{7}
\]

Use `P_g f(z)` as the actual row-eigenvalue vector of the physical Bell
tensor. Its column-eigenvalue vector is `T P_g f=S f`. Translation changes
neither the row orders nor their multiplicities. Consequently (7) gives a
literal strict crossing, with the same reciprocal row/column profiles.

This construction contains the folded code weights, but it does not require
an embedding into a positive-weight binary code. Any energy table with a
nondegenerate quadratic parity can first be put in the form (3) by binary
symplectic elimination. The implementation verifies the resulting quadratic
parity on the whole supplied table; it does not infer quadraticity from a
few pair tests.

## 2. Exact normalization and cost

The common gcd of the spectral polynomials in (4) is one: the zero-energy
term is `(1+z)^h`, while a maximum-energy term is `(1-z)^h`. These are
coprime in odd characteristic. The invertible Walsh transform preserves
the common component gcd. Also `f(0)` is supported only at zero and is
nonzero there. Hence no common scalar polynomial or initial monomial is
being silently charged in (4). Equation (7) proves that its actual degree
is exactly `h`.

Let `n_j` count the components whose order at zero is `j`. For every
rational `x>=1`, the complete two-cut pressure is

\[
Q(x)=x^h\left(\sum_{j=0}^h n_jx^{-j}\right)^2.
\tag{8}
\]

At tensor length `t`, every coefficient piece has rank product at most
`Q(x)^t`, by the order/degree argument in
[the polynomial synthesis note](bell-polynomial-synthesis.md). Thus a
fixed seed meeting the required strict pressure inequality would give an
effective family. No such seed is supplied here.

For fixed `k,h`, all coefficients and physical factors are finite explicit
objects, so the usual tensor-power generation has `2^O(lambda)` bit cost.
For a varying family with `h=poly(k)` and an energy table generated in
`2^O(k)` bit operations, a Walsh transform at each degree constructs (4) in
`2^O(k)` bit operations as well. Integer coefficient bit lengths are
`O(h+k)`. Listing all physical entries and Bell factors retains the same
exponential scale. The rank-width inequality must still be proved uniformly;
the explicit construction does not supply it automatically.

## 3. A strict family outside positive code weights

Take `g=0` and

\[
e(u)=\begin{cases}
0,&u=0,\\
1,&a\cdot b=1,\\
2,&u\ne0\text{ and }a\cdot b=0.
\end{cases}
\tag{9}
\]

The degree is two for every `k`. In contrast, the positive-code-weight
construction in [code quotients](code-quotient-crossings.md) has degree
at least half its binary coordinate dimension, namely `k`. Thus (9) is
outside that class for `k>=3`.

The exact profile also shows its limitation. Put `q=d^2` and work over
`F_3`. Away from the zero Walsh coordinate, direct expansion of (4) gives

\[
f_v(z)=z(1-d\eta_v/4)+z^2d\eta_v/4,
\qquad \eta_v=(-1)^{a_v\cdot b_v}.
\tag{10}
\]

The zero coordinate has nonzero constant term. Therefore the valuation
histogram `(n_0,n_1,n_2)` is

\[
\begin{cases}
(1,(q-d)/2,(q+d)/2-1),&k\text{ even},\\
(1,(q+d)/2-1,(q-d)/2),&k\text{ odd}.
\end{cases}
\tag{11}
\]

Its pressure square root is `x+n_1+n_2/x`. Since `n_1` is asymptotic
to `d^2/2`, the associated pressure parameter tends to one. Lower degree
alone therefore does not imply smaller rank cost.

## 4. The entire radial space has a common invariant covector

For the radial theorem, use `m=2k` information coordinates for the even
subspace of `F_2^(m+1)`. A vector of information weight `w` has total even
weight `w+(w mod 2)`, so its quadratic phase is

\[
\epsilon_w=(-1)^{\lceil w/2\rceil}.
\tag{12}
\]

The polar form is `u dot v+(sum u)(sum v)`. It is nondegenerate when `m`
is even, and a fixed symplectic coordinate change identifies its transform
with a translated physical Bell transform. This is the same type of
identification used in the folded-cube construction.

A radial row vector is specified by `a_0,...,a_m`; level `w` occurs on
exactly `binom(m,w)` Bell eigenvalues. Its row rank is consequently

\[
r(a)=\sum_{w:a_w\ne0}\binom mw.
\tag{13}
\]

Let `K` be the radial Walsh matrix,

\[
K_{r,w}=\sum_j(-1)^j\binom rj\binom{m-r}{w-j}.
\tag{14}
\]

It satisfies `K^2=2^m I`. The other-cut vector is

\[
Ta=2^{-m}K\,\mathrm{diag}(\epsilon_w)\,Ka.
\tag{15}
\]

This preserves the whole radial space, not just the coefficients of one
path. Its column rank is (13) applied to `Ta`.

Represent a linear functional on this space by a homogeneous binary form

\[
P(A,B)=\sum_{w=0}^m p_w A^{m-w}B^w,
\qquad L_P(a)=\sum_w p_w a_w.
\tag{16}
\]

The dual action of `K` is substitution `P(A+B,A-B)`, by the generating
identity for (14). In an algebraic extension with `i^2=-1`,

\[
\epsilon_w=\frac{1+i}{2}i^w+\frac{1-i}{2}(-i)^w.
\]

Combining these substitutions in (15), its dual action is

\[
T^*P=\frac{1+i}{2}P(R_+(A,B))
     +\frac{1-i}{2}P(R_-(A,B)),
\tag{17}
\]

where `R_+` fixes `X=A+B` and multiplies `Y=A-B` by `i`, while `R_-`
fixes `X` and multiplies `Y` by `-i`. All factors of `2^m` in (15)
are accounted for by the two half-normalized variables in `R_+` and `R_-`.

Now set `a=floor(m/4)` and

\[
P_m(A,B)=\bigl(AB(A^2+B^2)\bigr)^a(A+B)^{m-4a}.
\tag{18}
\]

The identity

\[
AB(A^2+B^2)=\frac{X^4-Y^4}{8}
\]

shows that both `R_+` and `R_-` fix `P_m`. Consequently `T^*P_m=P_m`.
Although (17) was derived in an extension, this is an identity over the
original odd-characteristic field: (18) has integer coefficients and all
division used powers of two. The argument also applies in characteristic
zero.

Every monomial of `P_m` has its `B` exponent between `a` and `m-a`.
At the crossing endpoint,

\[
L_{P_m}(\boldsymbol1)=P_m(1,1)=2^{m-3a}\ne0.
\tag{19}
\]

No binomial multiplicity has been divided by, and none is reduced modulo
the coefficient characteristic when counting a matrix rank.

## 5. Both rank bounds detect the same summand

Suppose radial tensors with row-eigenvalue vectors `a^(j)` sum to a
nonzero scalar `c` times the crossing. Their vectors sum to `c*1`.
By (19), at least one summand satisfies `L_(P_m)(a^(j))!=0`.

For this summand some nonzero row level lies between `a` and `m-a`.
Equation (13) gives `r_j>=binom(m,a)`. The same functional is invariant,
so `L_(P_m)(T a^(j))!=0` as well; the identical argument gives
`c_j>=binom(m,a)` for this **same** summand. This proves (1).

The physical binary label width is `lambda=k=m/2`. Stirling's formula
therefore gives

\[
\log_2(r_jc_j)\ge
2mH_2(1/4)-O(\log(m+2))
=4\lambda H_2(1/4)-O(\log(\lambda+2)).
\]

Any uniform maximum-product parameter is at least `H_2(1/4)`. The folded
construction supplies the matching upper bound within this space, proving
(2). In particular `3^5<2^8` proves `H_2(1/4)>4/5>1/sqrt(2)` exactly.

The hypothesis of one common radial identification matters. Arbitrary
nonradial coefficients, arbitrary code projections, or mixing tensors that
are radial in incompatible identifications are not covered by this proof.

## 6. Nonradial constructions and completed finite checks

The implementation also constructs a genuinely coupled family from a
quadratic norm over `F_(2^k)`. Choose `c` with absolute trace one and put

\[
N(x,y)=x^2+xy+cy^2,\qquad e(x,y)=\mathrm{wt}_{\rm binary}N(x,y).
\]

The norm vanishes only at `(0,0)`. Every nonzero scalar trace of `N` has
nondegenerate polar form, so symplectic relabeling puts its energy parity
in (3). This gives an explicit degree-`k` path, including arbitrary
invertible binary changes of output coordinates. The retained controls at
`k=1,...,6` and characteristics `3,5,7` do not improve the pressure. At the
tested default bases, the optimum of (8) is already one for `k>=3`.
This is a finite observation about the stated bases, not an exclusion of
all norm bases or all nonradial energies.

Four complete exploratory walks used 5,000 proposed moves each at physical
dimensions eight and sixteen, starting from the folded energies and their
threefold dilations. They evaluated respectively 264, 452, 61, and 148
valid proposed tables, with 159, 323, 42, and 129 distinct profile lists.
They retained the starting rates, approximately `0.8623168` and `0.8534553`.
The walks are not exhaustive, and no global optimum is inferred from them.
Every saved table is rebuilt and its rational pressure checked separately.

Run the exact controls with

```sh
python3 tools/check_bell_walsh_constructions.py \
  evidence/bell-walsh-energy-search-d8x1.json \
  evidence/bell-walsh-energy-search-d8x3.json \
  evidence/bell-walsh-energy-search-d16x1.json \
  evidence/bell-walsh-energy-search-d16x3.json \
  --output evidence/bell-walsh-constructions.json
```

The controls compare literal binomial products and Walsh sums, actual
physical matrix ranks, independently compressed Krawtchouk histograms,
every invariant-covector column through `m=32`, and literal two-copy
coefficient ranks. They verify every finite-field modulus used by the norm
examples. The all-length radial claim is the proof in Sections 4--5,
not an extrapolation from those checks.

## 7. Reciprocal coefficient retiming strictly improves a fixed seed

Let a strict Bell path have coefficients `a_0,...,a_D`, with
`S a_i=a_(D-i)` and `sum_i a_i=1`. Choose strictly increasing integers

\[
0=d_0<d_1<\cdots<d_D=H,
\qquad d_i+d_{D-i}=H.
\tag{20}
\]

The new path `sum_i a_i z^(d_i)` has the same endpoint and satisfies the
same reciprocal identity with degree `H`. If an eigenvalue first appears
in coefficient `i`, its new order is exactly `d_i`. In particular, if its
old valuation histogram is `n_0,...,n_D`, the new pressure square root is

\[
F=\sum_{i=0}^D n_i\exp(\xi_i),
\qquad \xi_i=(H/2-d_i)\log x,
\qquad \xi_{D-i}=-\xi_i.
\tag{21}
\]

This transformation preserves every coefficient vector. It changes the
relative order spacing, which a single monomial substitution cannot do.

Assume the positive profile entries satisfy that
`n_(D-i)/n_i` decreases strictly with `i`. Pairwise arithmetic-geometric
mean gives the infimum over real reciprocal schedules:

\[
\boxed{\quad
F_* = \sum_{i=0}^D\sqrt{n_i n_{D-i}},
\qquad
\xi_i^*=\tfrac12\log\frac{n_{D-i}}{n_i}.
\quad}
\tag{22}
\]

The ratio hypothesis makes these optimal exponents strictly decreasing,
so they respect (20). Rational scales approaching one and integer
reciprocal schedules approximate them arbitrarily closely. Formula (22)
therefore gives the infimum for these actual integer retimings as well.
It is not an assertion about arbitrary changes to the coefficient vectors.

For the folded seed `D=k` and `n_i=binom(2k+1,i)`. The adjacent ratios
have an especially useful gap. Write `K=k+1`, `t=i+1`, and
`R_i=n_(k-i)/n_i`. Direct cancellation of binomial coefficients gives

\[
\frac{R_i}{R_{i+1}}
=1+\frac{2K^2}{Kt-t^2}\ge9.
\tag{23}
\]

Consequently all independently optimal exponents in (22) are compatible.
They also give a strictly smaller infimum than the equally spaced folded
profile whenever `k>=3`. Indeed, equality would require
`R_i=x^(k-2i)` for one `x`, so every adjacent ratio in (23) would be
the same. Its values at `t=1` and `t=2` differ for `K>=4`. The paired
objective is strictly convex, so the equally spaced minimum is strictly
larger. Rational integer schedules approach (22), and therefore some
finite schedule attains a strict improvement for each such fixed `k`.
This argument does not claim that every rational approximation improves
the original seed.

There is an explicit rational schedule with polynomial degree and no
real-number oracle. Set `x=1+1/k^2`. For `i<k/2`, take the nonnegative
integer `b_i` nearest to `log(R_i)/(2 log x)`, equivalently the least integer
satisfying `x^(2b_i+1)>=R_i`. Put `b_(k-i)=-b_i` and put `b_(k/2)=0`
when `k` is even. Finally set `d_i=b_0-b_i`.

Equation (23) makes this schedule strictly increasing; the case `k=1`
is also direct. Each exponent has error at most `log(x)/2`, giving

\[
F\le\sqrt{x}\,F_*.
\tag{24}
\]

Since `log R_0=O(k)` and `log x` has order `k^(-2)`, the degree
`H=2b_0` is `O(k^3)`. Every choice is obtained by comparisons of rational
integer powers. Retiming preserves the coefficient bit lengths. Together
with the construction-cost bounds in Section 2, this is a uniform explicit
family at the same `2^O(k)` generation scale. The radial lower bound still
applies, so its limiting parameter remains `H_2(1/4)`.

For a concrete strict improvement take `k=4`. The old profile is
`(1,9,36,84,126)`. The preceding algorithm gives

\[
x=17/16,\qquad
(b_0,b_1,b_2,b_3,b_4)=(40,18,0,-18,-40),
\]

and hence the new coefficient degrees `(0,22,40,58,80)`. Its pressure
square root is the exact rational

\[
F_{\rm new}
=(17/16)^{40}+9(17/16)^{18}+36
 +84(17/16)^{-18}+126(17/16)^{-40}.
\tag{25}
\]

Direct integer arithmetic verifies `F_new^150<2^1024`. Thus
`rho_new=log_2(F_new)/8<64/75`, approximately `0.8532539341`.

For comparison, the original seed's pressure square root at any `y>=1` is

\[
F_{\rm old}(y)=y^2+9y+36+84/y+126/y^2.
\]

Its second derivative is at least two for all positive `y`. The quadratic
Taylor lower bound at `y_0=13/4` therefore proves, for all positive `y`,

\[
F_{\rm old}(y)\ge
F_{\rm old}(y_0)-\tfrac14F_{\rm old}'(y_0)^2
=\frac{2192857863}{19307236}>\frac{2271}{20}.
\tag{26}
\]

The exact integer inequality `(2271/20)^75>2^512` then gives
`rho_old>64/75`. Equations (25)--(26) establish a strict improvement
without trusting a floating-point minimizer. The old numerical optimum
is approximately `0.8534553058`. Both parameters are above the required
factoring threshold `1/sqrt(2)`; this finite-seed improvement earns no
new competition score.

The replay constructs all coefficients of these retimings through `k=4`,
checks both physical endpoint and reciprocal identities, and recomputes
their pressure from the actual new valuation histograms. It also checks
both rational inequalities separating the dimension-16 rates.
