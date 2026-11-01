# Integer prefixes and intervals in a synchronized field representation

## Structural result

The [binary-constraint contraction](trace-constraint-contraction.md) pays
for the rank of affine bit restrictions. For a binary operand prefix,
most fixed bits need not be paid for that way: substitute every completely
fixed word before applying the evaluator. Only one partially fixed word
per operand remains. This gives a uniform closure theorem for prefix and
interval queries in an explicitly specified word encoding.

It does **not** represent integer multiplication phases as synchronized
field products. That is still a separate, unresolved hypothesis. The
complete factoring score remains **1/2**. This ordinary argument and its
implementation were developed with OpenAI Codex assistance, using the
existing joint carry evaluator; there is no new Lean formalization.

## 1. Word encoding and the exact query

Fix an explicitly represented field E=F_(2^d), w dividing d, and an
explicit binary basis e_0,...,e_(w-1) of its subfield F_(2^w). Encode an
integer digit x in [0,2^w) by

```math
\operatorname{enc}(x)=\sum_{j=0}^{w-1}\operatorname{bit}_j(x)e_j.
```

This is a **binary linear encoding**, not a ring map from integer digits.
Integer products and powers of characters do not become field products
merely by using this notation.

There are h operands. Operand a has v_a words, so it occupies n_a=v_a*w
bits; its integer value is sum_i x_(a,i) 2^(wi), with low words first.
Let v=sum_a v_a. The input F is an explicit synchronized rational product
in the v encoded word variables, with row modulus t|w and charged degree
S as in the [joint theorem](joint-frobenius-product-contraction.md).
Each original word can have full or punctured subfield domain.

Specify one high binary prefix for each operand, of any length between
zero and n_a. The task is the same field-valued sum of F restricted to
those prefixes. A complete prefix specifies one value, and an empty
prefix imposes no restriction.

## 2. Substitute fully determined words

For a prefix of length k, the top floor(k/w) words are fixed, and the
next word has k mod w high bits fixed. Substitute each full word's
encoded value into every input factor. In a term with row rho, a fixed
variable contributes the explicit scalar

```math
u^{\sum_{j=0}^{w/t-1}\mu_j2^{\rho+tj}}.
```

Compute it by binary powering and absorb it into the term coefficient.
The bit length of its exponent is O(w+log(S+1)); the operation has
polynomial cost in the written parameters. Merge equal remaining terms.
The written term count cannot increase, the formal degree cannot increase,
and the row modulus is unchanged. This applies equally to denominator
factors, including those becoming zero under the zero-preserving convention.

If a fixed value violates a punctured domain, the answer is zero. If no
variables remain, evaluate the resulting constants directly. Otherwise let
v' be the number of remaining variables.

## 3. Charge only partial words

Construct the trace-dual basis e_j^* by inverting the binary Gram matrix
Tr(e_i e_j). Its nonsingularity follows from trace nondegeneracy. Then

```math
\operatorname{bit}_j(x)=\operatorname{Tr}(e_j^*\operatorname{enc}(x)).
```

The remaining prefix bits therefore give

```math
r=\sum_{a=1}^h(k_a\bmod w)\le h(w-1)
```

independent binary equations. Independence holds because they are distinct
coordinates within each partial word, and different operands use different
word variables. The [square-zero aggregation](square-zero-trace-contraction.md)
of these constraints gives total cost

```math
2^{O(r)}(S+rt+2)^{O(v'w/t)}\operatorname{poly}(L,d,w,v,h).
```

At t=w a uniform upper bound for every prefix query is therefore

```math
2^{O(hw+v\log(S+hw^2+2))}\operatorname{poly}(L,d,w,v,h).
```

In particular, a long known prefix does not require a proportional number
of trace factors. Its complete words decrease the number of variables.
This is an exact closure statement for the given field representation,
including mixed domains and coefficients in a larger extension field.

## 4. Arbitrary integer intervals

A half-open integer interval [A,B) inside [0,2^n) is a disjoint union of
at most 2n+1 binary prefix cylinders. One construction takes the largest
dyadic interval starting at A that fits before B, advances A, and repeats.
Equivalently, the interval is covered by maximal nodes of the binary trie;
at each depth only the two boundary paths can produce such nodes.

For fixed h operands, decompose each interval this way and distribute
the product of the covers. There are at most product_a(2n_a+1) prefix
cells. Each integer tuple occurs once, so adding their field-valued sums
gives the interval sum without overlap corrections. The overhead is
polynomial for fixed h (in particular, h=2). For variable h this product
must be charged; it is not hidden in a polynomial with an absolute degree.
Empty and singleton intervals are included.

This covers arbitrary top-word limits and operand prefixes *if* the
arithmetic scalar already has the stipulated synchronized representation.
It does not assume a full-field sum equals a bounded integer sum.

## 5. Score implications and remaining arithmetic gap

For h=2, t=w, v=O(n/w), and S polynomial in n, prefix and interval
queries cost

```math
2^{O(w+(n/w)\log n)}\operatorname{poly}(L,d,n).
```

For example w=Theta(sqrt(n log n)) makes this 2^O(sqrt(n log n)), provided
the base representation genuinely has t=w and polynomial written size
and field degree. Constructing the field by exhaustive irreducible search
at degree d=O(w) also fits 2^O(w) times polynomial overhead. More general
field degrees and field construction must be charged separately.
These conditional parameter values would beat every fixed logarithmic
saving score; no such representation of the complete factoring query is
provided here.

The retained [Fourier-dual formula](carry-fourier-elimination.md#7-the-exact-remaining-fourier-dual-problem)
contains phases zeta^(sxy), where x and y are **integer** digits and s is
an integer frequency. The encoding above does not turn this phase into a
small synchronized product. Neither the frequency sum nor the carry
weight is automatically included. Resolving those representations and
preserving deterministic recovery remain necessary for a factoring bound.

The useful structural change is narrower: once a suitable global arithmetic
representation exists, two operand prefixes and intervals can be retained
at the displayed cost. Do not revisit the generic full-rank trace encoding
for this boundary operation. The next research task concerns the phase
and frequency representation, not a larger prefix sweep.

## 6. Reproduction

Run `make verify-frobenius-prefix`. The
[implementation](../tools/frobenius_prefix_contraction.py) substitutes full
words and calls the prior constrained evaluator. The
[checker](../tools/check_frobenius_prefix_contraction.py) independently
enumerates integer digit assignments, applies their literal prefix tests,
encodes the surviving digits, and evaluates the original factors.
It includes alternate binary bases, a proper subfield, mixed domains,
all-fixed and empty-prefix cases, zero denominators, and multiple operands.
The [receipt](../evidence/frobenius-prefix-contraction.json) records exact
finite counts and source hashes, not an improved factoring theorem.
