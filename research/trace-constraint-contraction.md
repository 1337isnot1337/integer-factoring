# Keep a small binary constraint interface outside the Frobenius product

## Research question and result

The [joint carry theorem](joint-frobenius-product-contraction.md) loses its
small relative dimension if a bit trace is retained as one synchronized
factor. Can a small number of bit constraints be handled separately?
Yes: Gaussian elimination followed by a sum of synchronized products gives
an exact parameterized algorithm. It changes the representation assumption
behind that limitation, without supplying a complete factoring reduction.
The complete factoring score remains **1/2**.

The later [square-zero aggregation](square-zero-trace-contraction.md)
strengthens the bound below to `2^O(r) (S+rt+2)^O(vw/t)` by sharing the
row choices in a coefficient algebra. The packet implementation here
remains a separate exact comparison.

This is an ordinary argument and implementation developed with OpenAI Codex
assistance. It uses the preceding cyclic-carry theorem, elementary binary
linear algebra, and finite-field trace duality; no novelty claim is made
for those underlying tools. No Lean formalization is added.

## 1. Explicit input and cost

Use the earlier theorem's explicit field E=F_(2^d), subfield F_(2^w),
row modulus t dividing w, relative dimension D=w/t, and v variables.
Let F be its synchronized rational product, with zero-preserving inverses
and charged numerator/inverse degree S. Each variable ranges over either
the full subfield or its nonzero elements. Additionally impose an explicit
system of affine binary equations

```math
\sum_{i=1}^v \operatorname{Tr}_{F_{2^w}/F_2}(a_{ki}U_i)=b_k,
\qquad a_{ki}\in F_{2^w},\quad b_k\in\{0,1\}.
```

Let r be the binary rank of this system's coefficient functionals. An
inconsistent system gives zero. Otherwise the constrained sum of F is
computable in

```math
t^r (S+r+2)^{O(vw/t)}\operatorname{poly}(L,d,w,v,m)
```

field operations, with polynomial bit overhead. Here m is the number of
written constraints and L includes their coefficients and the original
product. Field construction and production of this representation remain
separate costs. The bound charges the **rank**, not the number of redundant
constraints. It preserves arbitrary correlations already present in F.

## 2. Independent equations suffice, including on punctured domains

The trace pairing (a,U) -> Tr(aU) on a finite field is nondegenerate.
Indeed Tr(X)=sum_(e=0)^(w-1) X^(2^e) is a nonzero polynomial of degree
less than 2^w, so its function is nonzero; scaling by a nonzero a is a
bijection. Consequently the rank of the functionals is exactly the binary
rank of their coefficient vectors (a_k1,...,a_kv).

Expand these vectors in an explicit binary field basis and row-reduce the
augmented system. Extension-field coordinates are also valid: the subfield
is a binary subspace of E, and this inclusion preserves rank. Binary row
operations preserve the solution set pointwise, so restricting variables
to nonzero elements causes no difficulty. Inconsistency is detected as a
zero coefficient row with right side one. Keep the r independent rows.

This preprocessing does not require enumerating field elements or computing
the trace-dual basis. XOR of encoded coefficient vectors implements the row
operations directly.

## 3. Split each indicator by its row

For a retained equation let beta(U)=sum_i Tr(a_i U_i). The indicator of
beta(U)=b, as an element of E, is 1+b+beta(U). Write it as

```math
1+b+\beta(U)=\sum_{\rho=0}^{t-1} Q_\rho(U),
```

where

```math
Q_\rho(U)=[\rho=0](1+b)+
 \sum_{i=1}^v\sum_{j=0}^{D-1}
 a_i^{2^{\rho+tj}} U_i^{2^{\rho+tj}}.
```

Each Q_rho is one legal row-rho factor with formal degree at most one
in the vD Frobenius symbols. Expand the product of the r indicators by
choosing one row per equation. There are at most t^r choices. Each choice
appends r legal factors to F, increasing S by at most r and leaving D
unchanged. Apply the existing joint carry theorem and add the answers.
This proves the bound, including empty systems, zero coefficients,
inconsistent equations, mixed summation domains, and zero denominators.

These are sums of **values in characteristic two**. They compute the
specified constrained weighted sum, not an integer count or an existence
test. Even numbers of equal contributions can cancel. An arithmetic
application must supply the appropriate scalar and reconstruction argument.

## 4. A genuine parameter escape, with a charged boundary

At t=w the relative dimension is one, and the exponential part is

```math
2^{\,r\log_2 w+O(v\log_2(S+r+2))}.
```

Thus a single bit equation need not force the earlier vw-dimensional carry
table. For example, for n=vw free input bits, w=Theta((log n)^2),
r=O(n/(log n)^2), and S polynomial in n, the stated cost is
2^O(n/log n). This is an evaluator parameter example, not a construction
of a factoring query at score one. Arbitrary synchronized products in this
regime can genuinely be evaluated with the stated extra constraints.

More generally a score-beta factoring transfer requires an exact arithmetic
reduction with all representation, field, query, and recovery costs at most
2^O(n/(log n)^beta), and

```math
r\log_2 t+(vw/t)\log_2(S+r+2)=O(n/(\log n)^\beta).
```

The constants in the carry bound are absolute. One can inspect the divisors
of w and choose an allowed modulus by this cost estimate. A base product
synchronized at t0 remains synchronized at each divisor of t0, after
reindexing its formal columns; an arbitrary larger modulus is not allowed.
Binary prefix restrictions may have small rank, but a multiplication or
phase evaluation cannot simply be treated as an affine binary equation.

Stop here unless an actual arithmetic representation meets this budget.
Expanding Theta(n) independent bit equations at t=w gives a charge
Theta(n log w) in this bound. This does not prove an algorithmic lower bound;
it explains why a naive complete bit circuit supplies no improvement here.

## 5. Exact controls and next question

The [implementation](../tools/trace_constraint_contraction.py) calls the
existing cyclic carry engine after row reduction and row splitting.
The [checker](../tools/check_trace_constraint_contraction.py) compares it
with literal field assignments and direct trace predicates, including
proper subfields, mixed domains, redundant and contradictory equations,
coupled equations, and zero denominators. It also compares alternate row
moduli and a literal t=1 encoding of the unsplit indicator.

Run `make verify-trace-constraints`; the
[receipt](../evidence/trace-constraint-contraction.json) records finite
controls and source hashes. Those controls do not prove the asymptotic
transfer or count as independent mathematical review.

The next concrete arithmetic question is whether the bounded operand domains
in [carry Fourier elimination](carry-fourier-elimination.md) can be expressed
with a small-r interface while its bilinear phases remain synchronized.
The unresolved phase representation is still essential. Merely selecting
larger field widths or checking more examples does not resolve it.
