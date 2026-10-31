# Exact global fixed-type minima for unrestricted Gaussian gauges

## Result and scope

The [Clifford profile reduction](clifford-gaussian-smith-profiles.md) makes
it possible to certify global gauge minima with a small rational dual.
The certificates here fix explicit physical bit marginals and cover
**every real independent mode height at both places and every positive
profile scale**. Internal Gaussian factors may have arbitrary poles.

For the listed marginals, the natural reflection heights attain the
global minimum of the complete paired-type Smith pressure. The result is
not a statement about every physical type. The joint optimization over
all types remains open, and the achieved factoring exponent is unchanged.

The first certificates use reflection interval lengths `M=4,8,16,32`.
Their exact rational marginals are stored in the
[certificate file](../evidence/reflection-gaussian-fixed-type-duals.json).
The two physical species have equal marginal vectors in each certificate.
For orientation, the certified optima are approximately:

| Interval `M` | `exp(-scale)` | Number of selected dual minors | Global fixed-type optimum `rho` |
|---:|---:|---:|---:|
| 4 | `18/25` | 23 | 0.7310811530 |
| 8 | `29/36` | 48 | 0.7166538701 |
| 16 | `69/79` | 95 | 0.7106044209 |
| 32 | `12/13` | 192 | 0.7083210521 |

The decimal values are only displays. The checker proves strictly
positive gaps above `1/sqrt(2)` with outward integer logarithms. It also
recomputes all selected determinant valuations over the rational
function field, without using a truncated power-series window.

## 1. The unrestricted fixed-type objective

Use the notation and exact `2M`-dimensional Clifford matrices from the
preceding note. Let `h` be the `4M` physical mode heights at both places,
and let `ell>0` be the profile scale. The two physical crossing-pair
species have `M/2` modes each. Their total mode energies are

\[
w_{X,i}=h_{0,x,i}+h_{0,v,i}+h_{\infty,x,i}+h_{\infty,v,i},
\]

and the same expression with `(x,v)` replaced by `(y,u)` for `w_Y`.
Fix probabilities `p_{X,i},p_{Y,i}` strictly between zero and one.

Write the positive Clifford Smith orders at place `s` in decreasing
order as `a_{s,1}>=...>=a_{s,M}>=0`. The full logarithmic profile product
is

\[
\log B=-M^2\ell-\frac12\sum_i h_i
       +\sum_{s=0,\infty}\sum_{i=1}^M
          \phi(a_{s,i}),
\qquad \phi(a)=\log(e^{a/2}+e^{-a/2}).
\tag{1}
\]

The fixed-type numerator is

\[
\mathcal J_p(h,\ell)=\log B+\sum_{X,Y}\sum_i p_{X,Y,i}w_{X,Y,i}.
\tag{2}
\]

If `H(p)` is the sum of the `M` binary marginal entropies, the associated
rate is `J_p/(2H(p))`. This normalization agrees with the positive root
`R=2rho` in the existing paired-type theorem.

The following finite dual certificate proves a global minimum of (2).
It is a specialized exact form of the earlier
[fixed-type convex formulation](rectangular-convex-pressure.md), with
only linearly many Clifford spectral variables.

## 2. A finite dual theorem

At either place, a `k` by `k` nonzero Clifford minor has weighted
valuation

\[
V_m(h,\ell)=d_m\ell+c_m\cdot h,
\tag{3}
\]

where `d_m` is an integer and `c_m` is an integer vector. Creation
coordinates contribute their physical heights, and annihilation
coordinates contribute their negatives. Since the `k` smallest signed
Clifford orders are `-a_{s,1},...,-a_{s,k}`, every such minor satisfies

\[
\sum_{i=1}^k a_{s,i}+V_m(h,\ell)\ge0.
\tag{4}
\]

Fix a proposed optimum `(h*,ell*)`. Assume its positive orders are
`a*_{s,i}` and define

\[
c_{s,i}=\phi'(a^*_{s,i})
       =\tfrac12\tanh(a^*_{s,i}/2),\qquad c_{s,M+1}=0.
\tag{5}
\]

Choose finitely many minors at each size `1<=k<=M`, with nonnegative
multipliers `lambda_m`. The certificate consists of these conditions:

1. Each selected minor is active at the proposed optimum: equality
   holds in (4).
2. For every place and size, the total multiplier is
   `sum_{m at (s,k)}lambda_m=c_{s,k}-c_{s,k+1}`.
3. For each physical height coordinate `h_i`,
   `sum_m lambda_m c_{m,i}=p_{type(i)}-1/2`.
4. The scale equation is `sum_m lambda_m d_m=-M^2`.

**Dual theorem.** These conditions imply

\[
\boxed{\mathcal J_p(h,\ell)\ge\mathcal J_p(h^*,\ell^*)}
\tag{6}
\]

for all real independent heights and every `ell>0`.

**Proof.** Convexity of `phi` gives

\[
\phi(a)\ge\phi(a^*)+\phi'(a^*)(a-a^*).
\]

For either place, the spectral mass equations telescope:

\[
\sum_i c_{s,i}a_{s,i}
=\sum_k(c_{s,k}-c_{s,k+1})\sum_{i\le k}a_{s,i}
\ge-\sum_{m\text{ at }s}\lambda_m V_m(h,\ell),
\tag{7}
\]

where the last inequality is (4) and nonnegativity of every multiplier.
Insert the tangent bounds and (7) in (1)--(2). The height and scale
equations cancel every term depending on `(h,ell)`. All selected
inequalities are equalities at `(h*,ell*)`, so the remaining constant
is exactly `J_p(h*,ell*)`. This proves (6). No optimization stopping
criterion or finite gauge box enters the proof.

## 3. Rational certificates at natural heights

Let `v` be a rational number in `(0,1)` and put `ell*=-log v`. Take the
natural integer reflection heights multiplied by `ell*`. The positive
Clifford orders are

\[
a^*_{s,i}=(4M-2-4(i-1))\ell^*,\qquad 1\le i\le M,
\tag{8}
\]

at both places. Therefore every quantity in the dual equations is
rational:

\[
c_{s,i}=\frac{1-v^{4M-2-4(i-1)}}{2(1+v^{4M-2-4(i-1)})}.
\tag{9}
\]

Activity also becomes an integer equality after canceling the positive
common factor `ell*`. The fixed-type objective at this point is

\[
\mathcal J_p(h^*,\ell^*)
=2\sum_{i=0}^{M-1}\log(1+v^{4i+2})
 -\log v\sum_{X,Y}\sum_{j=0}^{M/2-1}4(M-1-2j)p_{X,Y,j}.
\tag{10}
\]

The dual equations additionally imply the exact useful check

\[
\sum_{X,Y,j}4(M-1-2j)p_{X,Y,j}
=2\sum_{i=0}^{M-1}\frac{(4i+2)v^{4i+2}}{1+v^{4i+2}}.
\tag{11}
\]

Every logarithm in (10) and in `H(p)` is the logarithm of a positive
rational. The checker reduces its argument by an integer power of two
to `[1,2)`, then uses the established outward integer `log(1+x)` kernel.
It proves

\[
\mathcal J_p(h^*,\ell^*)-\sqrt2 H(p)>0
\]

and encloses the actual global minimum ratio. A half-precision or
floating logarithm is never used to prove this sign.

The theorem also covers correlated physical alphabet distributions with
these same bit marginals. Their expected height costs agree with (2),
while entropy subadditivity bounds their total entropy by `H(p)`.
Product Bernoulli distributions attain that entropy, so they give the
optimal rate among distributions with the displayed marginals.

## 4. Discovery versus verification

The optional [finder](../tools/find_reflection_gaussian_fixed_type_dual.py)
samples active prefix minors near the natural
height vector and solves a finite linear feasibility problem. An
infeasible sampled system does not prove that no certificate exists.
It uses an elastic linear-program dual to choose a new gauge direction,
asks the exact field oracle for additional active minors, and repeats.
This targeted step found certificates that random samples alone missed.
After feasibility, it continues with the target-distance dual to obtain
physical marginals closer to the natural Gibbs distribution. Closeness
is a search preference; it is not required by the proof.

After a floating support is found, rational row reduction reconstructs
the multipliers and marginals. Every dual equation and nonnegativity
condition is then checked exactly. A failure stops generation.

The independent checker recomputes each selected determinant valuation
with rational DVR elimination, checks the complete natural Clifford
profiles, and verifies every rational dual equation. It rejects a
changed multiplier, a removed minor, and an altered physical type.
The [replay record](../evidence/reflection-gaussian-fixed-type-dual-check.json)
contains the exact counts and outward intervals.

```sh
python3 tools/check_reflection_gaussian_fixed_type_duals.py \
  --output /tmp/reflection-gaussian-fixed-type-dual-check.json
python3 tools/find_reflection_gaussian_fixed_type_dual.py \
  --interval 8 --exp-negative-scale 29/36 --trials 100 \
  --output /tmp/reflection-gaussian-fixed-type-M8.json
```

The proof and exact certificates used substantial OpenAI assistant
support. No Lean source was changed. These results strengthen the
search's global verification boundary for explicit types; they do not
close the remaining optimization over all physical marginals.
