# A sharper rational certificate for the common cube family

The [homogeneous three-cut Gaussian construction](cube-homogeneous-gaussian.md)
admits the uniform contraction

\[
\rho=\frac{278}{445}<\frac58.
\tag{1}
\]

This improves the certified local constant, using the same fixed kernel
and the same three profiles `(2,4,8,10)`. The kernel, physical endpoint,
coefficient pieces, and finite-alphabet evaluation argument are unchanged.
The resulting cube-transfer exponent is `log_3(445/278)`, still below the
established factoring exponent one half. No new Big O Speedrun score is
claimed.

Use the physical frequency and deflated pressure scale

\[
p=\frac5{24},\qquad r=\frac{23}{53},\qquad
Z_0(r)=(1+r)(1+r^2)(1+r^4)(1+r^5).
\tag{2}
\]

For arbitrary wire width `lambda`, take `m` divisible by 24, with
`binom(m,5m/24)>=2^lambda`, and use a common set of `2^lambda` words
of weight `5m/24` on all eight legs. The proof in the preceding note
applies with `H_2(5/24)` in place of `H_2(1/5)`. In particular, it gives

\[
\sum_{a=x,y,z} r_k^a
\le12\lambda\,
\frac{\log_2 Z_0(r)-4p\log_2 r}{4H_2(p)}
 +O(\log(\lambda+2)).
\tag{3}
\]

The pressure fraction in (3) is approximately `0.624718611117`, whereas
`278/445` is approximately `0.624719101124`. The strict comparison has
the following exact integer certificate. Put

\[
A=(53+23)(53^2+23^2)(53^4+23^4)(53^5+23^5).
\]

Clearing the logarithm and entropy denominators gives the equivalent
inequality

\[
A^{10680}\,5^{5560}\,19^{21128}
 <53^{119260}\,23^{8900}\,24^{26688}.
\tag{4}
\]

The [checker](../tools/check_cube_gaussian_tight_pressure.py) evaluates
both integers in (4) exactly. It also replays all three finite Smith
factorizations of the fixed kernel. The
[receipt](../evidence/cube-gaussian-tight-pressure.json) retains the
bases and exponents, bit lengths, hashes of both integers, and the source
hashes. Floating-point calculations are display values only.

The number of pieces remains polynomial in `lambda`, each separate cut
has width at most `4 lambda`, and all factor entries remain constructible
in `2^O(lambda)` bit operations. Thus (1) is an effective uniform
certificate, rather than a finite-word numerical observation.

Run `python3 tools/check_cube_gaussian_tight_pressure.py` to reproduce
the exact comparison.
