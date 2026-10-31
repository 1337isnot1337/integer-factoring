# An exact entropy floor for the sparse Gaussian bulk model

## Scope

This note isolates a limit model from the
[Boolean down-set crossing](schur-downset-crossing.md). It permits
independent constant heights for every boundary and species at one common
species density, removes only regular internal Gaussian factors, and allows Laurent cross
blocks. Its best two-cut pressure still satisfies
\[
\rho\ge\frac34>\frac1{\sqrt2}.
\]

The proof is a short entropy inequality with an exact rational dual
certificate. It is a statement about the displayed model. It does not
exclude index-dependent gauges with unbounded slopes, arbitrary
four-letter physical gauges, or gains from Schur energy projection.
In particular, no universal obstruction for sparse Schur crossings is
claimed.

The calculations were developed with substantial OpenAI assistant
support. The exact source kernels and complete Schur references are in
the linked down-set note.

## 1. The normalized model

Order the four boundaries as \(B,X,V,Y\). Use the row cut
\((X,Y)\mid(V,B)\) at zero and the column cut
\((X,B)\mid(V,Y)\) at infinity. Divide valuations by the ambient
scale \(M\). In the sparse bulk model, the two-by-two cross blocks
have the following entry valuations and determinant valuations:

~~~math
\begin{array}{c|c|c}
\text{block}&\text{entry valuations}&\text{determinant valuation}\\\hline
R^+_0&\begin{pmatrix}1&1\\0&0\end{pmatrix}&2\\
R^-_0&\begin{pmatrix}1&0\\2&0\end{pmatrix}&1\\
C^+_\infty&\begin{pmatrix}-1&0\\-1&0\end{pmatrix}&0\\
C^-_\infty&\begin{pmatrix}-1&0\\-1&0\end{pmatrix}&-1 .
\end{array}                                               \tag{1}
~~~

The determinant entries in (1) include cancellations that entrywise
minimum matching would miss. They arise from the exact holonomies
\(I+qGH\) and \(z^MI+qGH\), not from a generic matrix with
those entry valuations.

Let \(\ell=\log t\ge0\), where \(t\) is the pressure parameter
on this normalized scale. Write zero-place row and column heights as
\(\alpha_i,\beta_i\), and infinity-place heights as
\(\gamma_i,\delta_i\). These heights already include the factor
\(\ell\). Thus all variables below have the same logarithmic units.

For each of the four blocks let \(k_j\) be its lower Smith order in
these units. Its upper order is \(D_j-k_j\), where

~~~math
\begin{aligned}
D_1&=2\ell+\alpha_X+\alpha_Y+\beta_V+\beta_B,\\
D_2&=\ell+\alpha_V+\alpha_B+\beta_X+\beta_Y,\\
D_3&=\gamma_X+\gamma_B+\delta_V+\delta_Y,\\
D_4&=-\ell+\gamma_V+\gamma_Y+\delta_X+\delta_B.
\end{aligned}                                             \tag{2}
~~~

Set

~~~math
\boldsymbol\kappa=(k_1,D_1-k_1,k_2,D_2-k_2,
                    k_3,D_3-k_3,k_4,D_4-k_4),\qquad
G=\sum_i(\alpha_i+\beta_i+\gamma_i+\delta_i).               \tag{3}
~~~

At common species density \(0<p<1\), the Gaussian pressure numerator
per mode is

~~~math
\mathcal P(p)=pG+\sum_{j=1}^{8}\log(1+e^{-\kappa_j}).       \tag{4}
~~~

The physical label entropy per mode is \(2H(p)\), with
\(H(p)=-p\log p-(1-p)\log(1-p)\). Accordingly the rank-width
exponent certified by this model is

~~~math
\rho=\frac{\mathcal P(p)}{8H(p)}.                         \tag{5}
~~~

## 2. Only six inequalities are needed

The lower Smith order of \(C^+_\infty\) cannot exceed either of
its diagonal entry valuations. Hence

~~~math
\gamma_X+\delta_V-\ell-k_3\ge0,\qquad
\gamma_B+\delta_Y-k_3\ge0.                               \tag{6}
~~~

Regularity of the four internal infinity edges requires

~~~math
\begin{aligned}
\gamma_B+\delta_X-2\ell&\ge0,\\
\gamma_X+\delta_B-\ell&\ge0,\\
\gamma_V+\delta_Y&\ge0,\\
\gamma_Y+\delta_V-\ell&\ge0.
\end{aligned}                                             \tag{7}
~~~

No assumption that the cross blocks themselves are regular is used.

Take one half of each left side in (6)--(7), and call their sum
\(\mathcal L\). Direct substitution from (2)--(3) gives

~~~math
\mathcal L
 =G-\left(D_1+D_2+k_3+\frac{D_4}{2}\right)\ge0.            \tag{8}
~~~

All zero-place heights cancel in this identity. This is the complete
linear dual certificate; the remaining regularity and Smith
inequalities can only further restrict the model.

## 3. Entropy proves the floor

For every real \(x\) and \(0\le\theta\le1\),

~~~math
\log(1+e^{-x})\ge H(\theta)-\theta x.                     \tag{9}
~~~

Use the eight occupancies

~~~math
\boldsymbol\theta=(p,p,p,p,p,0,p/2,p/2).                  \tag{10}
~~~

By (3), their weighted Smith sum is
\(p(D_1+D_2+k_3+D_4/2)\). Equations (4), (8), and (9)
therefore imply

~~~math
\mathcal P(p)\ge5H(p)+2H(p/2)+p\mathcal L
             \ge5H(p)+2H(p/2).                          \tag{11}
~~~

Entropy is concave and \(H(0)=0\), so
\(H(p/2)\ge H(p)/2\). Thus

~~~math
\boxed{\qquad
\mathcal P(p)\ge6H(p),\qquad \rho\ge\frac34 .
\qquad}                                                   \tag{12}
~~~

This is an exact analytic inequality. Numerical optimization is not
part of its proof.

## 4. What a finite-family comparison would require

A useful stability fact explains both the scope and the remaining
gap. Let \(A\) be an invertible square matrix over a Laurent series
field, with sorted Smith orders \(\boldsymbol\nu\). Multiplying
its rows and columns by monomials of integer orders
\(r_i,c_j\) changes those orders to \(\boldsymbol\nu'\), with

~~~math
\|\boldsymbol\nu'-\boldsymbol\nu\|_1
 \le\sum_i|r_i|+\sum_j|c_j|.                              \tag{13}
~~~

For nonnegative shifts, the sorted Smith orders increase coordinatewise.
Their total increase is exactly the determinant-order change. Splitting
signed shifts into positive and negative parts and using the triangle
inequality proves (13). The derivative of
\(\log(1+t^{-x/M})\) has absolute value at most
\((\log t)/M\), so (13) also gives

~~~math
\left|\sum_j\log(1+t^{-\nu'_j/M})
      -\sum_j\log(1+t^{-\nu_j/M})\right|
\le\frac{\log t}{M}
       \left(\sum_i|r_i|+\sum_j|c_j|\right).              \tag{14}
~~~

For index-linear gauges with fixed bounded slopes, the right side is
controlled by the average normalized index. On the down-set of
\(b\)-bit masks of weight at most \(k\), this average is at most
\(k/b\). For fixed pressure parameter and density away from zero and
one, the corresponding normalized error is therefore \(O(k/b)\).

This estimate is not uniform over slopes growing with the ambient
parameters, or over pressure and density parameters approaching
singular limits. Nor does diagonal scaling alone prove that every
matrix orientation converges to (1). A finite-family exclusion requires
those additional controls. The exact offset interval is separately
handled by its [explicit kernel](schur-offset-interval-reduction.md).

## 5. Exact reproduction

~~~bash
python tools/check_sparse_bulk_species_floor.py \
  --output evidence/sparse-bulk-species-gauge-floor.json
~~~

The checker verifies the dual identity (8) over rational coefficient
vectors, all eight occupancy coefficients in (10), and the exact
comparison \(3/4>1/\sqrt2\). No floating optimization result is used
as a certificate.
