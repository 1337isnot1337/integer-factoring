# A sharper rational type for the existing cube kernel

## Result and scope

The degree-six homogeneous Gaussian cube construction admits the exact
strict interface bound

\[
\rho < \frac{312359221}{500000000}=0.624718442.
\]

This strengthens the previously certified `278/445` bound in
[cube-gaussian-tight-pressure.md](cube-gaussian-tight-pressure.md).
It changes only the common physical frequency and the pressure scale.
It does not change the kernel or improve the complete factoring exponent,
which remains `1/2`. The cube target is `rho < 1/sqrt(3)`.

The construction and its finite Laurent factorizations are those of
[cube-homogeneous-gaussian.md](cube-homogeneous-gaussian.md). All three
one-particle cut profiles remain `(2,4,8,10)`.

## Exact certificate

Write

\[
Z(r)=(1+r)(1+r^2)(1+r^4)(1+r^5),\qquad
H(p)=-p\log p-(1-p)\log(1-p).
\]

For every fixed rational frequency `0<p<1`, the original coefficient
argument gives

\[
\rho_{p,r}=\frac{\log Z(r)-4p\log r}{4H(p)}.
\]

Choose

\[
p=\frac{208173}{1000000},\qquad
r=\frac{434044514}{1000000000}.
\]

The [checker](../tools/certify_cube_optimized_type.py) proves

\[
4\left(\frac{312359221}{500000000}\right)H(p)
-\log Z(r)+4p\log r>0.
\]

It uses outward dyadic intervals for logarithms of positive rationals.
Range reduction and the positive atanh series include an explicit tail
bound. Every sign decision uses integers and fractions, not floating-point
logarithms. The certificate passes at both 192 bits with 96 series terms
and 224 bits with 112 terms. A smaller bound `0.624718441` has a strictly
negative gap at these same parameters. That negative control does not
prove a lower bound over other parameters or constructions.

The [receipt](../evidence/cube-optimized-type.json) contains the exact gap
intervals, both precision settings, source hashes, and the three replayed
one-particle profiles. The checker does not repeat the complete Gaussian
tensor factorization; that identity is inherited from the cited construction.

## Uniform word construction

Take word lengths divisible by one million, with occupied fraction `p`.
Choose the smallest such length `m` for which
`binom(m,pm) >= 2^lambda`, and retain any `2^lambda` words of that weight.
Standard type counting and minimality give

\[
mH_2(p)=\lambda+O(\log(\lambda+2)),\qquad m=O(\lambda+1).
\]

The large fixed spacing affects constants, not the asymptotic rate.
Every eight-leg word entry has total degree `24pm`, and the same
coefficient tensors are used at all three cuts. The original finite
Laurent coefficient extraction therefore gives total width at most
`12 rho lambda + O(log(lambda+2))`, with each individual width at most
`4 lambda`. The uniform effectiveness argument is unchanged.

This is a small numerical refinement, not a new tensor family or a
certificate of the global optimum of the pressure problem.
