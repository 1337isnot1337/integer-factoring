# Fixed field extensions preserve the crossing exponent

## Result

A coefficient-rank certificate over a fixed finite extension field can be
converted to the base field without enlarging its physical alphabet. The
coefficient-cost sum increases by at most the extension degree, once for
the whole tensor power. This constant does not change the width exponent
or the resulting factoring exponent.

This permits searches over a fixed field such as \(\mathbb F_4\). It does
not itself improve any pressure inequality. The explicit local descent in
[the conjugate-code construction](conjugate-residual-phase-codes.md) remains
useful because it supplies a binary local seed and additional choices of
physical weights.

## 1. A rank bound for a linear coefficient projection

Let \(L/F\) be a finite field extension of degree \(d\). Choose an
\(F\)-basis \(e_0=1,e_1,\ldots,e_{d-1}\), and let
\(\lambda:L\to F\) be the projection onto the \(e_0\) coefficient.
In particular, \(\lambda(1)=1\), and \(\lambda\) fixes every element of
\(F\).

For every matrix \(A\) over \(L\), entrywise projection satisfies

\[
 \operatorname{rank}_F\lambda(A)
 \le d\operatorname{rank}_L A. \tag{1}
\]

To prove this, write \(A=UV\), with inner dimension
\(r=\operatorname{rank}_L A\), and expand every entry of \(U,V\) in
the chosen basis. Define

\[
 M_{ab}=\lambda(e_ae_b).
\]

The expanded factors give the exact factorization

\[
 \lambda(UV)=U_{\rm coord}(I_r\otimes M)V_{\rm coord}, \tag{2}
\]

whose inner dimension is \(dr\). This proves (1) and also supplies an
effective base-field rank factorization. The multiplicative factor can
be sharp: the rank-one matrix \((e_ae_b)_{a,b}\) projects to the
nondegenerate matrix \(M\), which has rank \(d\). Nondegeneracy follows
because, for nonzero \(a\in L\), choosing \(b=a^{-1}\) gives
\(\lambda(ab)=1\).

## 2. Apply the projection after taking the tensor power

Suppose a crossing construction over \(L\), including any common-type
restriction, has polynomial coefficient tensors

\[
 K_m(w)=\sum_k A_{m,k}w^k,\qquad
 \sum_k A_{m,k}=E_{D_m}, \tag{3}
\]

where \(E_D(x,y,u,v)=[x=v][y=u]\). Define the base-field coefficients
\(\widehat A_{m,k}=\lambda(A_{m,k})\). Their sum is still the crossing,
because \(\lambda\) is linear and fixes zero and one.

Write \(r_L,c_L\) for the two coefficient flattening ranks over \(L\),
and \(r_F,c_F\) for the ranks after projection. Applying (1) separately
to the two flattenings gives

\[
 \begin{aligned}
 r_F(\widehat A_{m,k})&\le d\,r_L(A_{m,k}),\\
 c_F(\widehat A_{m,k})&\le d\,c_L(A_{m,k}),\\
 \sum_k\sqrt{r_F(\widehat A_{m,k})c_F(\widehat A_{m,k})}
 &\le d\sum_k\sqrt{r_L(A_{m,k})c_L(A_{m,k})}. \tag{4}
 \end{aligned}
\]

The same bounds apply after a physical type restriction because entrywise
projection commutes with selecting entries. They also apply to a rational
local seed once the common-type product has been proved polynomial.

The order of operations matters. First form the coefficient tensor of the
whole typed tensor power over \(L\); then apply \(\lambda\). Applying
\(\lambda\) independently to the local factors produces a different
construction and can incur a factor \(d^m\). In general, \(\lambda\)
is not multiplicative.

## 3. Cost and exponent

For fixed \(L/F\), arithmetic, coefficient projection, and the expanded
factorizations (2) have constant overhead per field operation or matrix
entry. Thus the existing effective \(2^{O(\ell)}\) construction cost
remains \(2^{O(\ell)}\). The extra factor \(d\) in (4) is independent
of the tensor-power length and of the external alphabet size \(D_m\).
It changes only a constant in the coefficient-cost bound.

Consequently, a valid fixed-extension pressure certificate with width
parameter \(\rho\) transfers to the base field with the same \(\rho\).
The [effective crossing transfer](pieces-to-factoring.md) then gives the
same exponent \(\alpha=\log_2(1/\rho)\). The extension must be fixed
for this statement; a growing extension degree needs a separate cost
analysis.

One must choose a linear projection with \(\lambda(1)=1\). The usual
field trace is not always suitable: the trace of one from
\(\mathbb F_4\) to \(\mathbb F_2\) is zero. A basis containing one
always supplies the required projection.

## Exact controls

```sh
python tools/check_fixed_field_descent.py --output evidence/fixed-field-descent.json
```

The checker verifies (2) for deterministic factored matrices over
\(\mathbb F_4\) and \(\mathbb F_8\), including sharp rank-one examples.
It independently builds polynomial crossing tensor powers through length
three, projects their full coefficient tensors, and checks both flattening
ranks and every endpoint entry. It also checks the nonmultiplicativity and
trace counterexamples. These controls accompany the proof for arbitrary
fixed finite extensions.
