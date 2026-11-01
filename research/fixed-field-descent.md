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

Section 4 extends the cost analysis to uniformly constructed extensions
of polynomial degree. A fixed coefficient field is not necessary under
those additional hypotheses.

## 1. A rank bound for a linear coefficient projection

Let \(L/F\) be a finite field extension of degree \(d\). Choose an
\(F\)-basis \(e_0=1,e_1,\ldots,e_{d-1}\), and let
\(\lambda:L\to F\) be the projection onto the \(e_0\) coefficient.
In particular, \(\lambda(1)=1\), and \(\lambda\) fixes every element of
\(F\).

For every matrix \(A\) over \(L\), entrywise projection satisfies

\[
 \mathrm{rank}_F\lambda(A)
 \le d\mathrm{rank}_L A. \tag{1}
\]

To prove this, write \(A=UV\), with inner dimension
\(r=\mathrm{rank}_L A\), and expand every entry of \(U,V\) in
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
For a prime base field, the
[scale-one modular transfer](scale-one-modular-transfer.md) then gives the
same exponent \(\alpha=\log_2(1/\rho)\). The extension must be fixed
for this statement; a growing extension degree needs a separate cost
analysis.

One must choose a linear projection with \(\lambda(1)=1\). The usual
field trace is not always suitable: the trace of one from
\(\mathbb F_4\) to \(\mathbb F_2\) is zero. A basis containing one
always supplies the required projection.

## 4. Polynomial-degree extensions with uniform presentations

Fix a prime `p`. At binary label width `lambda`, allow a different field

```math
L_\lambda=\mathbb F_p[X]/(f_\lambda),\qquad
e_\lambda=\deg f_\lambda\le(2\lambda)^A,
```

where `A` is fixed. Require a deterministic `2^O(lambda)` algorithm for
the monic irreducible modulus and for the proposed crossing pieces and
their two factorizations. Field elements use the polynomial basis starting
with one. Assume polynomially many pieces, scale one, and extension-field
channel widths with sum at most

```math
4\rho\lambda+C\log_2(2\lambda).
```

Then there is an effective family **over the one fixed field `Fp`** with
the same leading `rho`. Indeed, let `ell` select the constant coordinate.
For every supplied factorization `A=UV` with `r` inner indices, the explicit
descent is

```math
\ell(A)_{ij}=
\sum_{k=1}^{r}\sum_{a=0}^{e_\lambda-1}
 U_{ik,a}\,\ell(X^a V_{kj}).                            \tag{5}
```

This uses `e_lambda*r` channels, once for the whole piece. Consequently
each width increases by at most `ceil(log2(e_lambda))`, and the sum
increases by at most `2A log2(2lambda)+2`. If a width exceeds `2lambda`,
use the direct factorization through its two physical input labels as
in the modular interface. That cap only decreases the width.

Arithmetic in the supplied polynomial presentation takes a polynomial
number of bit operations in `e_lambda` for fixed `p`, using ordinary
polynomial multiplication, reduction and Euclidean inversion. The
coordinate factors in (5) have the same polynomial overhead. This
preserves the full `2^O(lambda)` construction bound, the polynomial piece
count, and the exact endpoint because `ell(1)=1`.

Thus the hypotheses of the existing fixed-prime modular transfer hold
after descent; the field used by the factoring machine does not change
with its input. A claimed variable-field family still has to provide its
uniform modulus constructor and all rank factors. Polynomial extension
degree alone does not establish those algorithms. A varying characteristic
is outside this argument.

The [growing clock construction](ordered-charge-area-crossings.md) gives
a concrete uniformly generated example with `e_lambda=O(log lambda)`.
Its loss is only `O(log log lambda)` channel bits. It supplies no improved
uniform rank bound. The accompanying
[checker](../tools/check_ordered_charge_area.py) also replays (5), including
sharp rank increases, in extension degrees two through eight.

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
