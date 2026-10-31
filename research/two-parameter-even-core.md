# Two output parameters and the exact Gaussian locus

## Result

The four-letter [even-core postselection](even-core-postselection.md)
admits independent output parameters. Over any characteristic-two field,
its polynomial tensor has the form

\[
K_{L,M}=A+L B_y+M B_v+LM C.
\tag{1}
\]

The four coefficient tensors have degrees \(17,24,32,14\), respectively.
Their endpoint values are

\[
A(1)=B_y(1)=B_v(1)=E,\qquad C(1)=0.
\tag{2}
\]

Consequently polynomial parameters satisfying \(L(1)=M(1)\) give a
strict local crossing. Independent Laurent monomials \(L=w^a\) and
\(M=w^b\) are also valid after clearing a common power of \(w\).
All these monomials equal one at the endpoint.

There is an exact classification of the Gaussian members over the
rational function field: the two coordinate axes \(LM=0\), together
with the single point

\[
(L,M)=(w^8,1).
\tag{3}
\]

On the axes the tensor is a Gaussian contraction; this includes
degenerate specializations whose vacuum coefficient vanishes. At (3)
the vacuum coefficient is \(w^8\), and an explicit ordinary Gaussian
factorization is given below. Every other member with both parameters
nonzero violates a quartic Wick identity. No pressure or factoring
exponent improvement is claimed by this classification.

## 1. Construction and polynomiality

Use the full four-mode seed, with \(c=1+w\), and the rational tensor
\(G(1+x_0y_3u_0v_3/c^2)\) from the one-parameter construction.
On both input legs take the extreme-mode state
\(|\varnothing\rangle+|\{0,3\}\rangle\). On output \(y\) take
\(\langle\varnothing|+L\langle\{0,3\}|\), and on output \(v\)
take \(\langle\varnothing|+M\langle\{0,3\}|\).
The two middle modes remain free on each leg.

Keeping track of the two output occupations separately gives (1).
The charge identity used in the one-parameter proof is unchanged:
the physical index gauge has exponent
\((\Delta n_2-\Delta n_1)/2\), independent of the selected extreme
occupations. Thus the same ramification argument applies coefficient
by coefficient in the two formal parameters. Each of the 1024
coefficientwise divisions by \(c^2\) is exact.

The endpoint calculation gives (2), hence

\[
K_{L,M}(1)=(1+L(1)+M(1))E
\tag{4}
\]

whenever the parameter functions are regular there. The special case
\(L=M=\lambda\) recovers the earlier family: its linear coefficient
is \(B_y+B_v\), and its quadratic coefficient is \(C\).

For arbitrary rational parameter functions, the algebraic classification
below still applies over the function field. Their use as polynomial
crossing paths may require additional denominator control. Equation (4)
alone does not remove those denominators.

## 2. Two Wick identities exclude every other point

Order the occupation variables as

\[
(x_1,x_2,y_1,y_2,u_1,u_2,v_1,v_2).
\]

For a four-element set \(\{a,b,c,d\}\), define

\[
\mathcal W_{abcd}(K)=K_\varnothing K_{abcd}
+K_{ab}K_{cd}+K_{ac}K_{bd}+K_{ad}K_{bc}.
\tag{5}
\]

Exact polynomial expansion of (1) gives

\[
\begin{aligned}
\mathcal W_{\{x_1,x_2,y_1,u_1\}}(K_{L,M})
 &=LM^2w^2c^4(w^8+L),\\
\mathcal W_{\{x_1,x_2,u_1,v_1\}}(K_{L,M})
 &=L^2Mw^2c^7(1+M).
\end{aligned}
\tag{6}
\]

Every Gaussian tensor satisfies both identities. In a rational function
field, \(w\) and \(c\) are nonzero. If \(L,M\ne0\), the first
identity forces \(L=w^8\), and the second forces \(M=1\). This proves
that (3) is the only possible Gaussian point off the axes, over every
extension of the coefficient field as well.

If either output parameter is zero, that output projects the extreme
modes to vacuum. The quartic term requires the top mode on both outputs,
so it vanishes. The remaining tensor is a contraction of the Gaussian
\(G\) with Gaussian two-mode boundary states. It is therefore Gaussian.
The checker additionally verifies all 99 even Pfaffian identities
symbolically along each entire axis.

## 3. The exceptional point has an explicit Gaussian factorization

Let \(\theta_0,\ldots,\theta_7\) denote the variables in the order
above. At (3), the full exterior generating polynomial is

\[
K_{w^8,1}
=w^8\prod_{(i,j,p)\in\mathcal E}(1+p\theta_i\theta_j),
\tag{7}
\]

where the eight edges are

| \(i\) | \(j\) | \(p\) |
|---:|---:|---|
| 0 | 3 | \(c^5\) |
| 0 | 6 | \(w^8\) |
| 1 | 2 | \(c^3\) |
| 1 | 7 | \(w^8\) |
| 2 | 4 | \(1\) |
| 3 | 5 | \(1\) |
| 4 | 7 | \(c^5\) |
| 5 | 6 | \(c^3\) |

Literal exterior multiplication reproduces all 256 coefficients of
(1). In particular its vacuum is the nonzero monomial \(w^8\).
This proves sufficiency at the exceptional point and completes the
classification.

The independent-output family therefore has a Gaussian member that is
absent from the equal-parameter slice. This distinction matters when
screening nearby non-Gaussian paths: an improvement over a restricted
parameter slice need not itself use a non-Gaussian tensor.

## 4. Reproduction

Run from the repository root:

```sh
python tools/check_two_parameter_even_core.py
```

The [checker](../tools/check_two_parameter_even_core.py) uses only exact
binary polynomial arithmetic. The
[evidence](../evidence/two-parameter-even-core.json) records 1024 exact
divisions, 1024 endpoint coefficient checks, 768 comparisons with the
equal-parameter family, both bivariate Wick identities, 198 formal
axis Pfaffian checks, and the complete exceptional-point factorization.
The `--seed-directory` option exports the four coefficient tensors and
the exceptional Gaussian seed for subsequent searches.
