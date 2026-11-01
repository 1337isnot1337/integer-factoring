# From rational-function crossings to one fixed prime field

## Result and hypotheses

A future low-rank construction may use rational functions in a fixed number
of auxiliary parameters over `F_p`. Under explicit degree and generation
bounds, this gives the same factoring exponent as a construction directly
over `F_p`. The auxiliary field need not embed into one common finite field
at every label width: specialize and descend each width before recursion.

Fix a prime `p`, a constant number `d` of parameters, and
`1/2<=rho<1`. Suppose a scale-one crossing family over
`K=F_p(t_1,...,t_d)` has, for every binary width `lambda`:

1. At most `(2lambda)^C` pieces on four `lambda`-bit labels, and two exact
   factorizations, with separate
   binary channel widths at most `2lambda` and total width at most
   `4rho lambda+C log2(2lambda)`.
2. Uniform algorithms producing the index set, widths, and every rational
   piece/factor value in
   `2^(C lambda)` bit operations. Numerators and nonzero denominators
   are explicit coefficient arrays of total degree at most `2^(C lambda)`.
   The parameter count `d` is fixed, so these arrays have `2^O(lambda)` size.
3. Exact crossing and factorization identities over the rational-function
   field. These are hypotheses to prove, not conclusions of specialization.

**Transfer theorem.** These hypotheses give an effective scale-one crossing
family over the single fixed field `F_p`, with the same leading width
parameter `rho`, polynomially many pieces, and `2^O(lambda)` complete
construction cost. The [fixed-characteristic transfer](scale-one-modular-transfer.md)
therefore gives deterministic factoring time

```math
exp(O(n/(log n)^alpha)),\qquad alpha=log2(1/rho).
```

This is a sufficient interface for symbolic rank constructions, not such a
construction. It proves no exponent above one half. In particular, generic
ranks of the local polynomial tensors in the Hamming notes do not by
themselves supply the assumed crossing decomposition. The argument is
ordinary mathematics, with substantial OpenAI assistant support and no
new Lean formalization.

The degree bounds are essential. An exponential-size arithmetic circuit
can have doubly exponential degree. Its size alone does not meet the
explicit hypotheses above. Nor does this theorem cover arbitrary infinite
fields or changes in characteristic.

## 1. One finite extension suffices at each width

Generate the complete source factor tables at a fixed width `lambda`.
Each factor has two input labels and a channel string of at most
`2lambda` bits, so all entries of all factors and pieces number
`2^O(lambda)`. Their total generation cost remains `2^O(lambda)`.

List the nonconstant denominators and let `L` be the sum of their total
degrees. Then `L<=2^O(lambda)`. Choose the least positive integer `k` with
`q=p^k>L`. Thus `k=O(lambda)` and `q=2^O(lambda)`.

A nonzero degree-`e` polynomial in `d` variables over a field with `q`
elements has at most `e q^(d-1)` zeros. The usual induction on the last
variable proves this from the one-variable root bound: separate the leading
coefficient's zero set and count the remaining fibers. Applying this bound
to each denominator, their union occupies at most `L q^(d-1)<q^d` points.
Some point in `F_q^d` therefore avoids every denominator.

Enumerate all `q^d=2^O(lambda)` points and test them on the generated
coefficient arrays. Each test takes `2^O(lambda)` bit operations, so the
complete search does too. Choosing the first successful point is a uniform
algorithm. Evaluation at that point is a homomorphism from the localization
containing all displayed rational values to `F_q`. It preserves the exact
crossing sum and both factorizations.

No input-table value is inverted. All inversions here are of explicitly
checked nonzero constants in a finite field, during family generation.

## 2. Construct the finite field within the budget

An explicit construction is enough; a fast irreducible-polynomial algorithm
is unnecessary. Enumerate monic degree-`k` polynomials over `F_p` and test
irreducibility, choosing the first one. There are `p^k=2^O(lambda)` choices.

The test uses polynomial gcds and repeated Frobenius powers modulo the
candidate. Check `X^(p^k)=X` modulo the candidate, and for each prime divisor
`ell` of `k` check

```math
gcd(f,X^(p^(k/ell))-X)=1.
```

The first condition says that every irreducible factor has degree dividing
`k` and has multiplicity one. If a proper factor has degree `a<k`, then
`a` divides `k/ell` for some prime divisor of `k/a`, and a gcd detects it.
An irreducible degree-`k` polynomial passes every test. Existence of finite
fields supplies such a polynomial. Factoring the integer `k` by trial division
and doing these polynomial operations have polynomial cost in `k` for fixed
`p`.

Represent `F_q` as `F_p[X]/(f)`. Schoolbook polynomial multiplication,
reduction, inversion, and evaluation have polynomial bit cost in `k`.
Thus construction and all specialization tests fit `2^O(lambda)`, with no
unproved fast factoring subroutine hidden in the field setup.

## 3. Descend each factorization, not each factor separately

Use the basis `e_a=X^a`, `0<=a<k`, and let `L:F_q->F_p` select the constant
coordinate. Then `L` is linear and `L(1)=1`. Replace each specialized tensor
`A_i` by `A'_i=L(A_i)`. Their sum is the crossing.

For any one factorization `A_i(x,y,u,v)=sum_j U_j(x,y)V_j(u,v)`, write
`U_j=sum_a U_(j,a)e_a`, with prime-field coordinates `U_(j,a)`. Linearity
gives the explicit identity

```math
A'_i(x,y,u,v)=sum_(j,a) U_(j,a)(x,y) L(e_a V_j(u,v)).
\tag{1}
```

This uses `k` times as many channels. It is not the generally false rule
of multiplying `L(U_j)` by `L(V_j)`. Apply (1) independently to the row and
column arrangements.

Padding the `k` component choices by zeros increases each binary width by
at most `ceil(log2 k)=O(log lambda)`. Their sum remains
`4rho lambda+O(log lambda)`. Replace a separate width above `2lambda` by
the explicit trivial input-label factorization, which only decreases it.

All coordinate extractions and products with the fixed basis elements
have polynomial cost in `k`. The factor tables have only `2^O(lambda)`
entries. Hence complete descended factor generation and every metadata/query
algorithm retain the required cost. The piece count is unchanged.

The final family is over `F_p` at every width. Its later base change to a
commutative `F_p`-algebra uses only those prime-field values and no inversions
in that algebra. The recursion therefore needs no common embedding of the
width-dependent extension fields used during construction.

## 4. A separate varying-extension corollary

More generally, suppose a family is already given over
`F_p[X]/(f_lambda)`, where `f_lambda` is explicitly irreducible of degree
`k_lambda<=lambda^C`, its presentation is uniformly generated in
`2^O(lambda)` steps, and its factor algorithms meet the same bit-cost bounds.
The descent in Section 3 applies directly. Its width penalty is still
`O(log lambda)`, so it preserves the same factoring exponent.

Here generation of the field presentation is an explicit hypothesis.
The brute-force field construction in Section 2 was budgeted for
`k=O(lambda)`; it must not silently be used at `k=lambda^C` for `C>1`.
This corollary extends the fixed-extension observation in
[hamming-sectors.md](hamming-sectors.md#5-fixed-extension-fields-are-compatible-with-the-search-target).

## 5. Finite arithmetic controls

```sh
python3 tools/check_extension_descent.py
python3 tools/check_extension_descent.py --output evidence/extension-descent.json
make verify
```

The standard-library checker constructs binary extensions of degrees one
through eight and verifies the inverse of every nonzero element. It checks
4,800 entries in 240 explicit factorization descents, using the coordinate
and multiplication-pairing formula (1), and exercises specialization outside
unions of denominator zero sets. A negative control shows why ordinary field
trace cannot replace a functional with `L(1)=1` when the characteristic
divides the extension degree.

These computations check finite arithmetic and the explicit descent rule.
The uniform specialization budget and factoring implication follow from
the proof above; no finite test stands in for those all-size arguments.
