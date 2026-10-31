# Four boundary parameters and a scalar Gaussian surface

## Result

Allow all four extreme-mode boundary vectors in the
[even-core construction](two-parameter-even-core.md) to vary. Put
parameters \((p,L,r,M)\) on the legs \((x,y,u,v)\), respectively:
each input is \(|\varnothing\rangle+a|\{0,3\}\rangle\), and each
output is \(\langle\varnothing|+a\langle\{0,3\}|\), with its stated
parameter \(a\). Keep the middle two modes as the physical alphabet.

This gives a tensor \(K_{p,L,r,M}\) that is polynomial and of degree
at most one in each of the four parameters. Its endpoint is

\[
K_{p,L,r,M}(1)=
\bigl(1+L(1)r(1)+p(1)M(1)\bigr)E
\tag{1}
\]

when all parameter functions are regular at one. In particular, every
choice of four unit-coefficient Laurent monomials \(w^{a_x},w^{a_y},
w^{a_u},w^{a_v}\) is a strict crossing after clearing a common power
of \(w\).

Over any characteristic-two rational function field, the complete
Gaussian locus consists of the four coordinate hyperplanes
\(pLrM=0\), together with

\[
Lr=w^8,\qquad Mr=1.
\tag{2}
\]

There is no new Gaussian family to optimize on the nonzero part of
this surface. Every such tensor is exactly a scalar multiple of the
exceptional tensor already identified in the two-output family:

\[
K_{p,w^8/r,r,1/r}
=\frac p r\,K_{1,w^8,1,1}.
\tag{3}
\]

Here Gaussian contractions include degenerate specializations with a
zero vacuum coefficient. Equation (3) has nonzero vacuum when \(p,r\)
are nonzero. The note gives an algebraic classification and a family
for further searches; it asserts no improved pressure bound.

## Polynomial construction and endpoint

Before contraction, the tensor is \(G(1+Q/c^2)\), with
\(Q=x_0y_3u_0v_3\) and \(c=1+w\). Recording the extreme occupation
on each leg separately gives sixteen coefficient tensors. Their
coefficientwise polynomiality follows from the same index-energy
bound as the two-output family. The finite constructor verifies all
4096 divisions by \(c^2\) exactly.

At the crossing endpoint, precisely three of these coefficient tensors
survive: the constant coefficient, the coefficient of \(Lr\), and the
coefficient of \(pM\). Each equals \(E\). This proves (1) as a
polynomial identity in the four parameters. Setting \(p=r=1\)
recovers the previously checked two-output construction.

General rational parameter substitutions remain subject to denominator
control if used for the polynomial crossing method. In contrast, a
common monomial clearing factor for Laurent-polynomial parameters is
one at \(w=1\), so it preserves the endpoint in (1).

## Necessity and sufficiency of the Gaussian conditions

Use the occupation order
\((x_1,x_2,y_1,y_2,u_1,u_2,v_1,v_2)\) and the Wick residual convention
of the two-output note. Exact expansion gives

\[
\begin{aligned}
\mathcal W_{\{x_1,x_2,y_1,u_1\}}(K)
 &=pLrM^2w^2c^4(w^8+Lr),\\
\mathcal W_{\{x_1,x_2,u_1,v_1\}}(K)
 &=pL^2rMw^2c^7(1+rM).
\end{aligned}
\tag{4}
\]

If all four parameters are nonzero, these two necessary Gaussian
identities force (2). Conversely, substituting (2) and multiplying by
\(r^2\) gives the polynomial identity

\[
r^2K_{p,w^8/r,r,1/r}=prK_{1,w^8,1,1}.
\tag{5}
\]

The checker compares all 2048 coefficients of \(p\), \(r\), and the
physical tensor labels in (5). The known eight-edge exterior
factorization of \(K_{1,w^8,1,1}\) therefore proves sufficiency.

If any one parameter is zero, its boundary vector selects vacuum on
both extreme modes of that leg. The quartic \(Q\) requires an occupied
extreme mode on every leg, so its contribution vanishes. What remains
is a contraction of the original Gaussian with Gaussian two-mode
boundary states. This proves sufficiency on all four coordinate
hyperplanes and completes the classification.

## Reproduction

```sh
python tools/check_four_parameter_even_core.py
```

The [checker](../tools/check_four_parameter_even_core.py) uses exact
binary polynomial arithmetic. Its
[evidence](../evidence/four-parameter-even-core.json) records the full
coefficient construction, all endpoint coefficients, recombination
with the two-output family, the two four-variable Wick identities,
and the scalar identity on the Gaussian surface. Optional
`--seed-directory` exports the complete sixteen-tensor basis.
