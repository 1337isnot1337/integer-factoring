# Pascal crossings in truncated symmetric algebras

## Result and scope

In characteristic \(p\), replace each binary exterior mode by a
commuting variable with \(x^p=0\). This gives an ordinary tensor-power
version of the Pascal construction, with \(p\) states per mode and no
graded signs. Its one-particle Smith order \(\kappa\) contributes the
factor \(1+s^{-\kappa}+\cdots+s^{-(p-1)\kappa}\).

The endpoint, both rank-profile interfaces, and the finite construction
are proved below. This opens \(p\)-level non-Gaussian perturbations.
It does not supply a pressure witness below \(1/\sqrt2\), and the
interval bulk asymptotic still has that same critical rate.

## 1. A functor in characteristic \(p\)

Let \(K\) have characteristic \(p\), and define

~~~math
\mathcal A_r=K[x_1,\ldots,x_r]/(x_1^p,\ldots,x_r^p).
                                                               \tag{1}
~~~

Its basis consists of \(x^\alpha\), with
\(\alpha\in\{0,\ldots,p-1\}^r\). Every linear substitution
preserves the ideal in (1), since

~~~math
\left(\sum_i c_i x_i\right)^p=\sum_i c_i^p x_i^p=0.
~~~

For a matrix \(D\), let \(S_p(D)\) be its induced matrix on this
monomial basis: column \(\beta\) is the coefficient vector of
\(\prod_j(\sum_iD_{ij}x_i)^{\beta_j}\). Direct substitution gives

~~~math
S_p(DE)=S_p(D)S_p(E),\qquad S_p(I)=I.                 \tag{2}
~~~

Thus invertible linear substitutions induce invertible basis changes.
This remains true over a power-series ring and its field of fractions.

There is also a useful kernel form. Write

~~~math
\mathrm{Exp}_{<p}(t)=\sum_{d=0}^{p-1}\frac{t^d}{d!},
\qquad
G_D(x,y)=\prod_{i,j}\mathrm{Exp}_{<p}(D_{ij}x_i y_j).
                                                               \tag{3}
~~~

Every displayed factorial is a unit in \(K\). For factors sharing a
row variable, truncated exponential addition is valid: all missing
terms have a factor \(x_i^p\). Regrouping by column instead gives

~~~math
[x^\alpha y^\beta]G_D
=\frac{S_p(D)_{\alpha,\beta}}{\prod_j\beta_j!}.        \tag{4}
~~~

The factorial normalization is a constant diagonal change of basis.
It can be made separately on the physical input legs.

### Smith spectrum

If the nonzero one-particle Smith orders are
\(\kappa_1,\ldots,\kappa_r\), functoriality transports the regular
Smith basis changes through (2). The resulting diagonal entries are

~~~math
z^{\sum_i d_i\kappa_i},\qquad 0\le d_i<p.             \tag{5}
~~~

For a singular one-particle matrix, each zero direction permits only
digit zero in a nonzero induced entry. In particular, one-particle
rank \(r\) induces rank \(p^r\).

## 2. The global Pascal endpoint

Use the polynomial one-particle matrix \(\mathcal D_S(w)\) from
[the all-field Pascal construction](all-field-pascal-crossings.md).
Its row boundaries are \((V,Y)\), and its column boundaries are
\((U,X)\). For \(m\) copies, replace every entry by that entry times
\(I_m\), then apply \(S_p\).

A physical label is a word of \(m\) vectors in
\(\{0,\ldots,p-1\}^{S}\). Fix the total exponent of each mode on
\((X,V)\), and independently on \((Y,U)\). The resulting ordinary
four-boundary tensor has endpoint

~~~math
\Gamma_m(-1;x,y,u,v)=[x=v][y=u].                      \tag{6}
~~~

**Proof.** At \(w=-1\), the upper-left block of the one-particle
matrix vanishes. Every selected \(V\) exponent must come from an
\(X\) variable. The total \(V\) and \(X\) exponents agree, so no
\(X\) exponent can flow to a \(Y\) variable. The identity block
therefore forces \(y=u\), coordinate by coordinate.

The remaining triangular Pascal block can only increase a mode index.
The matching mode totals equalize the sum of the indices, counted with
their exponents. Every increase is zero, forcing \(x=v\). The
diagonal linear coefficients are one, so the surviving monomial
coefficient in \(S_p\) is one. This proves (6). \(\square\)

There are no determinant or Pfaffian signs. Direct sums of one-particle
matrices give ordinary tensor products under (2). Hence the global
tensor is an ordinary power of the local tensor before the physical
type restriction. Full word types are also allowed, since they fix
all mode totals required by the endpoint argument.

## 3. Both cuts and their gauges

At infinity, put \((X,U)\) on one side and \((V,Y)\) on the other.
Equation (4) identifies the cut matrix with \(S_p(\mathcal D_S)\),
up to constant diagonal factorial changes. The all-field one-particle
Smith theorem therefore applies through (5).

At zero, use \((X,Y)\mid(V,U)\). In the kernel (3), separate factors
whose two variables lie on one side from factors crossing the cut.
Each internal factor is a polynomial with constant term one in a
nilpotent ideal. If its gauged coefficients are regular, multiplication
by it is an invertible regular map on that side. Removing all internal
factors leaves the truncated symmetric kernel of the two cross blocks.
Equation (5) again gives the complete orders. All factorial changes
are constant units and preserve these orders.

The same ramification \(w=-z^2\) and natural physical gauges as in
the all-field note are valid. If \(S=\{a_1,\ldots,a_n\}\), the
zero one-particle orders are two copies of \(2a+1\). The infinity
orders are \(2\nu+1-2n\), where \(\nu\) ranges through \(S\)
and its greedy Pascal complement over \(K\). Replace every binary
profile factor by

~~~math
L_p(t)=\sum_{d=0}^{p-1}t^d.                           \tag{7}
~~~

At independent physical matching types, use the same paired pressure
criterion with \(R=2\rho\). For a mode of total physical energy
\(g\), its physical partition factor is \(L_p(s^{-g/R})\).
For a common type the corresponding parameter is \(R=4\rho\).
The coefficient-rank proof and type optimization are the ordinary
tensor-power special case of the
[general transfer argument](graded-even-power-transfer.md), with all
physical parities set to zero.

## 4. Non-Gaussian perturbations and effectivity

In the finite algebra on all four physical boundaries, let \(H_I\)
be any nonconstant monomial. No even-degree condition is required.
For polynomials \(f_I(w)\) satisfying \(f_I(-1)=0\), set

~~~math
\widetilde G(w)=
\left(1+\sum_I f_I(w)H_I\right)G_{\mathcal D_S(w)}.     \tag{8}
~~~

After the same input-factorial normalization, (8) has the same typed
crossing endpoint. It has ordinary tensor powers in every prime
characteristic. Such tensors need not be induced by any one-particle
matrix, so their two cuts require direct finite profile calculations.

For fixed \(p,S\), and fixed perturbation polynomials, the local
alphabet and polynomial degree are fixed. The \(m\)-copy degree is
\(O(m)\), entries are products of explicitly computable local entries,
and all full tables and ordinary rank factorizations take
\(2^{O(m)}\) operations. Rational type probabilities and rational
gauges are realizable after fixed divisibility and ramification steps.
The [fixed-prime transfer audit](fixed-prime-transfer-source-audit.md)
provides the separate factoring interface. A strict pressure improvement
is still required to improve the score.

## 5. The interval bulk still has the same critical rate

For an initial interval, put
\(f_p(x)=\log L_p(e^{-x})\). Its integral is

~~~math
I_p=\int_0^\infty f_p(x)\,dx
=\frac{\pi^2}{6}\left(1-\frac1p\right).               \tag{9}
~~~

Indeed, \(L_p(e^{-x})=(1-e^{-px})/(1-e^{-x})\); expanding the two
logarithms into absolutely integrable series gives (9).
In the common-type normalization, the leading interval pressure is

~~~math
\frac{(8-R^2)I_p}{4\ell},\qquad \ell=\log s.          \tag{10}
~~~

Thus its bulk cancellation still occurs at \(R=2\sqrt2\), or
\(\rho=1/\sqrt2\). Replacing binary digits by \(p\)-ary digits
alone does not change this critical value. The new search opportunity
is the larger local algebra and its non-Gaussian corrections.

The accompanying exact checker verifies functoriality, the factorial
kernel identity, typed endpoints, and the digit-sum Smith spectrum
in small prime fields. These finite controls support the implementation;
the all-size construction is the proof above.

Run [the checker](../tools/check_truncated_symmetric_crossings.py):

```sh
python tools/check_truncated_symmetric_crossings.py \
  --output evidence/truncated-symmetric-crossings.json
```

The [recorded evidence](../evidence/truncated-symmetric-crossings.json)
contains 8,706 functorial coefficient equalities, 288,148 independently
typed endpoint entries, and ten complete actual cut profiles. The
controls include \(p=2,3,5\), three copies, and a noninterval mode
set in characteristic three.
