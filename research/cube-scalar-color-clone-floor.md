# A uniform floor for scalar Gaussian color clones

For every integer \(M\ge33\), replicating any homogeneous scalar Gaussian cube kernel over \(M\) independent colors, then retaining the vacuum and single-particle letters, has common-type, common-scale, unweighted profile pressure

\[
\boxed{\rho_{\mathrm{prof}}>3/4>1/\sqrt3.}
\]

This includes non-Cayley scalar kernels. It does not cover interactions between colors or other physical code spaces. It explains why scalar replication cannot realize the large-alphabet target from the [graded minor budget](cube-graded-minor-budget.md).

## 1. Repeated exact-color blocks

Work in characteristic two. Let \(K(x,y,z)\) be an alternating eight-by-eight scalar matrix with homogeneous entries of degree \(D>0\) and antipodal-matching endpoint. Replicate it as \(K\otimes I_M\). The physical vacuum-plus-one-particle space has \(q=M+1\) letters and is closed under every particle deletion. Thus the within-cut Gaussian triangular factors descend to invertible polynomial factors. The complete cut orders equal those of the restricted cross compound matrix.

That matrix splits by the exact multiplicity of every color. A block with \(r\le4\) distinct colors, each occurring once, has dimension \(4!/(4-r)!\), corresponding to their injective placements among four vertices. There are \(\binom Mr\) identical blocks. Other blocks contribute nonnegative terms and can be discarded below.

Let \(\mu_{a,r}\) be the smallest order of this block at axis \(a\), with \(\mu_{a,0}=0\). Choose any fixed \(r\) antipodal pairs and give them distinct colors. One entry of the block at every cut is the same polynomial

\[
P_r=\prod_{i=1}^r K_{v_i,\bar v_i}.
\]

It has degree \(rD\) and endpoint one. The first Smith order equals the minimum entry valuation. Hence

\[
\sum_a\mu_{a,r}\le\sum_a v_a(P_r)\le rD.
\]

For the complete cut partition functions, convexity of the logarithmic exponential sum consequently gives, for \(0<t<1\),

\[
\frac13\sum_a\log Z_a(t)\ge\log Z_M(u),\qquad
u=t^{D/3},\qquad Z_M(u)=\sum_{r=0}^4\binom Mr u^r.
\]

The three cut minima come from one common scalar kernel.

## 2. An exact all-size polynomial certificate

For every integer \(M\ge33\),

\[
Z_M(u)>(1+M u^{4/3})^3\qquad(0<u<1).
\]

Set \(P_M(s)=Z_M(s^3)-(1+M s^4)^3\) and \(Q_M(s)=P_M(s)/s^3\). The quotient has degree nine. The coefficient \(C_j(M)\) of \(v^j\) in \((1+v)^9Q_M(v/(1+v))\) is a polynomial in \(M\) of degree at most four. Its exact Newton expansion is
\(C_j(M)=\sum_{k=0}^4 a_{j,k}\binom{M-33}{k}\), with:

| \(j\) | \(a_{j,0}\) | \(a_{j,1}\) | \(a_{j,2}\) | \(a_{j,3}\) | \(a_{j,4}\) |
|---:|---:|---:|---:|---:|---:|
| 0 | 33 | 1 | 0 | 0 | 0 |
| 1 | 198 | 6 | 0 | 0 | 0 |
| 2 | 396 | 12 | 0 | 0 | 0 |
| 3 | 528 | 33 | 1 | 0 | 0 |
| 4 | 1782 | 156 | 6 | 0 | 0 |
| 5 | 1881 | 210 | 9 | 0 | 0 |
| 6 | 176 | 300 | 29 | 1 | 0 |
| 7 | 3102 | 825 | 78 | 3 | 0 |
| 8 | 5973 | 963 | 81 | 3 | 0 |
| 9 | 7634 | 2447 | 352 | 28 | 1 |

Every constant is positive and every other coefficient is nonnegative. Thus every \(C_j(M)>0\) for integer \(M\ge33\). Substituting \(s=v/(1+v)\) proves the asserted inequality on the whole open unit interval. This finite coefficient identity proves the all-size statement directly.

## 3. Entropy conversion

Let a common physical type have total nonzero probability \(p\) and alphabet entropy \(H>0\). Its color distribution can be arbitrary, and
\(H\le H_2(p)+p\log_2 M\). The Gaussian tensor has degree \(D/2\) per occupied leg, so the common coefficient total in a word tensor power is \(4Dmp\). Therefore

\[
\rho_{\mathrm{prof}}\ge
\frac{\log_2 Z_M(u)-4p\log_2u}{4H}.
\]

The entropy identity gives

\[
H_2(p)+p\log_2M+\frac43p\log_2u
\le\log_2(1+M u^{4/3}).
\]

The polynomial certificate makes the numerator strictly larger than \(3[H_2(p)+p\log_2M]\ge3H\). This proves the claimed floor. A constructive escape must alter the repeated distinct-color structure, for example through genuine interactions between colors.

## 4. Exact replay

```sh
python3 tools/check_cube_scalar_color_clone_floor.py \
  --output /tmp/cube-scalar-color-clone-floor.json
```

The checker reconstructs all fifty Newton coefficients, verifies additional exact evaluations, and checks the complete color-block multiplicities. The [receipt](../evidence/cube-scalar-color-clone-floor.json) retains the certificate. The factoring exponent remains \(1/2\).
