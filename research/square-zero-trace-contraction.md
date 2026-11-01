# Aggregate trace constraints in a square-zero coefficient algebra

## Result

The [row-packet construction](trace-constraint-contraction.md) handles
r independent binary equations by at most t^r synchronized products.
There is a stronger way to share their work: keep the subset of equations
used in a square-zero coefficient algebra and perform **one** joint carry
contraction. Its deterministic cost is

```math
2^{O(r)}(S+rt+2)^{O(vw/t)}\operatorname{poly}(L,d,w,v,m).
```

This removes log(t) from the leading cost of the constraint rank. The
base field, synchronized rational product, domains, degree S, and affine
trace equations are exactly those in the preceding note. Field construction
and the construction of a complete arithmetic representation are still
separate. The complete factoring score stays **1/2**.

The later [divided-power lift](divided-power-trace-lifting.md) retains the
moments needed for exact 0/1 trace predicates modulo a growing power of two.
Its coefficient dimension is (2H)^r and is explicitly charged; it recovers
integer counts at sufficient precision for that represented input class.
The field-valued construction and its original scope below remain intact.

Square-zero variables and subset coefficient extraction are standard
algebraic tools. This application uses the existing cyclic carry theorem
and trace row packets, with OpenAI Codex assistance. It is an ordinary
argument with exact finite controls, not a new Lean theorem or a claim of
independent mathematical acceptance.

## 1. The coefficient algebra

Row-reduce the equations over F_2 as before; an inconsistent system gives
zero. For r independent equations introduce

```math
A=E[z_1,\ldots,z_r]/(z_1^2,\ldots,z_r^2).
```

It has the explicit E-basis z_J=product_(j in J) z_j indexed by subsets
J of {1,...,r}. Addition is coordinatewise. Multiplication is

```math
z_I z_J=\begin{cases}z_{I\cup J},&I\cap J=\varnothing,\\0,&\text{otherwise}.
\end{cases}
```

Thus represented operations have cost 2^O(r) times a polynomial in d;
even direct multiplication of all pairs of the 2^r coordinates suffices.
The algebra is commutative but is **not a field**. Its nonzero nilpotent
elements do not have inverses and do not obey the field's Frobenius period.

## 2. One product records every row choice

For the kth equation, use the earlier explicit packets Q_(k,rho), so its
Boolean indicator I_k(U) equals sum_(rho=0)^(t-1) Q_(k,rho)(U). Define

```math
H(U)=\prod_{k=1}^r\prod_{\rho=0}^{t-1}
                  (1+z_k Q_{k,\rho}(U)).
```

Because z_k^2=0, all terms choosing two packets for the same equation
vanish. Therefore

```math
H(U)=\prod_{k=1}^r(1+z_k I_k(U)),\qquad
[z_1\cdots z_r]H(U)=\prod_{k=1}^r I_k(U).
```

Each factor 1+z_k Q_(k,rho) is a synchronized row-rho polynomial over A
with formal degree at most one. Append these at most rt factors to the
base product F(U), considered in A by the scalar embedding E -> A.
Coefficient extraction is E-linear, hence the desired constrained sum is

```math
[z_1\cdots z_r]\sum_U F(U)H(U).
```

This identity is pointwise before summation. It applies to coupled
constraints, redundant original systems, punctured variable domains,
and arbitrary weights in the original scalar product. It computes the
specified field-valued sum; cancellation is not an existence test.

## 3. Why the carry engine still applies

The cyclic-carry proof groups explicit polynomial coefficients, adds and
multiplies them, and uses scalar power-sum identities in the subfield.
Those power-sum identities remain true after embedding E into A. They
require no division in the coefficient algebra. The unique-carry argument
concerns nonnegative integer exponents and is unchanged. Therefore the
polynomial-product part of the theorem works over A verbatim.

Denominators require a separate check. All denominator factors here are
the **original E-valued factors** of F. At every allowed assignment, the
identity iota(x)=x^(2^(d+1)-3) is valid for their scalar values in E,
including x=0. Its image in A remains valid. Expanding these denominators
into their d Frobenius factors consequently costs the same S as before.
No square-zero factor is inverted, and the new numerator coefficients
are never assumed to have Frobenius period d.

One must not extend this argument to arbitrary A-valued denominators.
For example, z_1 is nonzero with square zero; the field inverse formula
does not make it invertible. The implementation exposes no inverse on A
and retains the original scalar denominator list unchanged.

After expanding the original inverses, the product has formal degree
at most S+rt. The old carry-state bound is therefore
(S+rt+2)^O(vw/t). Multiplying by the represented cost of A operations
proves the announced bound. This conservative estimate does not require
a finer degree filtration or cancellation-based state pruning.

## 4. Arithmetic consequences and limits

At t=w, the exponential part is

```math
2^{O(r+v\log(S+rw+2))}.
```

In particular, r=O(n/(log n)^beta) independent constraints no longer
incur an additional log w factor. A sufficient scalar-evaluator budget is

```math
r+(vw/t)\log(S+rt+2)=O(n/(\log n)^\beta),
```

together with all input construction, field arithmetic, and reconstruction
costs at that same scale. This is a wider admissible interface, not an
arithmetic reduction satisfying it.

For the [operand-prefix construction](frobenius-prefix-contraction.md),
substituting complete fixed words leaves r<=h(w-1). For a fixed number h
of operands, t=w, v=O(n/w), and S polynomial in n, prefix and interval
queries then cost

```math
2^{O(w+(n/w)\log n)}\operatorname{poly}(L,d,n).
```

Choosing w=Theta(sqrt(n log n)) gives 2^O(sqrt(n log n)) under these
representation hypotheses. This improves the prefix interface bound;
it does not improve factoring. The complete integer phases and frequency
sum still lack the required representation. A naive circuit introducing
Theta(n) independent bit variables still pays 2^O(n) here.

## 5. Reproduction and stopping condition

Run `make verify-square-zero-traces`. The
[implementation](../tools/square_zero_trace_contraction.py) packs coefficients
in the squarefree monomial basis, uses the existing carry engine once,
and extracts the top coefficient. The
[checker](../tools/check_square_zero_trace_contraction.py) compares with
literal constrained field assignments and the previous packet algorithm.
It includes full-rank systems, contradiction, mixed domains, a proper
subfield, scalar zero denominators, and explicit nilpotent controls.
The [receipt](../evidence/square-zero-trace-contraction.json) records the
actual finite checks and source hashes.

The mechanism is now explicit and its cost charged. Further work must
address a concrete arithmetic phase representation with sublinear
constraint rank; larger generic constraint tests do not establish that.
