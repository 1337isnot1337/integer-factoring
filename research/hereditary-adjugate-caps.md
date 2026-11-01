# Adjugate inequalities on every smaller physical alphabet

The existing graded adjugate inequalities can be applied after **any
square physical restriction**, then transferred to the complete profiles
of the original scalar crossing. This gives additional necessary prefix
inequalities over every coefficient field. The strongest subset choice
always uses the smallest paired energies, so no subset enumeration is
needed when compiling the constraints.

The [zero-energy instance](zero-energy-hereditary-adjugate.md) rules out
the formerly retained five-letter target. A new **integer six-letter
point of the stronger relaxation** survives every compiled inequality and
has rigorously negative pressure. It is not a realizing tensor or a
factoring algorithm. The established score remains `1/2`.

## 1. The inherited constraints

Let `r,s` be normalized complete profiles of a scalar strict crossing on
`q` physical symbols, with deficit `d` and sorted nonnegative paired
energy lists `X,Y`, both starting at zero. Write `S_k` for prefix sums.
For `2<=m<=q`, take the first `m` entries of each energy list and put

```
n=m^2,
C_m=(n-1)d + m(sum X[:m] + sum Y[:m]).
```

Let `TopCell_k` be the sum of the `k` largest energies in the truncated
`m`-by-`m` cell grid. Let `TopX_l,TopY_l` and `Xmax,Ymax` refer to those
same truncated lists. Then

```
k S_(n-1)(r) + (k-1) S2(s)
    <= k C_m - TopCell_k,                         1<=k<=n,

(l-1) S_n(r) + S_(n-l)(r) + (l-1) S2(s)
    <= l C_m - max(TopX_l+l Ymax, TopY_l+l Xmax), 1<=l<=m.
```

Both also hold with `r,s` exchanged. At `m=q` these recover the old
full-alphabet inequalities. Proper subalphabets add new constraints.
For `gamma=(r+s)/2`, averaging the two orientations gives the same
formulas with every profile argument replaced by `gamma`.

These are scalar constraints under the hypotheses of the
[graded adjugate theorem](https://github.com/1337isnot1337/integer-factoring/blob/3670b7df8db2bca84bc1dc1ea492c7044fa93e60/research/graded-adjugate-caps.md).
Its cofactor and rank-one arguments are inputs, with their original
credit retained. The new work is restriction, restoration of minima,
monotone prefix transfer and the subset reduction. Matrix-valued auxiliary
interfaces and endpoints obtained only after type projection are outside
this statement.

## 2. Why core normalization does not lose the inequality

Choose any `m` physical labels on each paired alphabet. Use the original
normalization on the restricted profiles, calling them `u,v`, with minima
`a,b>=0`. Every core minor is an original minor, so its prefixes dominate
the full prefixes. If the subset energy minima are `x0,y0`, its normalized
deficit is `d*=d+x0+y0-a-b` and its paired energies subtract those minima.

In the input theorem its cofactor shift becomes

```
C*=C_m-(x0+y0)-(n-1)(a+b).
```

Its top-cell sum subtracts `k(x0+y0)`, and its top-line expression subtracts
`l(x0+y0)`. Thus the energy minima cancel. Expanding the normalized
prefixes cancels `a` and leaves these **stronger** inequalities on the
core profiles:

```
k S_(n-1)(u)+(k-1)S2(v)+[k(n-3)+2]b <= k C_m-TopCell_k,

(l-1)S_n(u)+S_(n-l)(u)+(l-1)S2(v)+[l(n-3)+2]b
    <= l C_m-max(TopX_l+l Ymax,TopY_l+l Xmax).
```

The extra minimum terms are nonnegative. Every remaining prefix
coefficient is nonnegative too. Dropping the extra term and comparing
core minors with full minors proves the stated inherited constraints.
This applies to integer heights and, by ramification, rational heights.

## 3. Why the smallest energies suffice

For fixed `m,k`, increase one chosen energy by `delta>=0`. The main term
of the cell right side increases by `km delta`. Its top-cell term can
increase by at most `min(k,m) delta`. The right side is nondecreasing.

In the shared-line right side the main term increases by `lm delta`;
the maximum subtracted expression increases by at most `l delta`.
This right side is also nondecreasing in every chosen energy. Replacing
any subset by the `m` smallest entries therefore gives the strongest
bound for those profile indices.

The [compiler](../tools/hereditary_adjugate_caps.py) constructs all these
constraints exactly using rational arithmetic.

## 4. A surviving integer six-letter target

Take deficit six, both paired lists `(0,0,12,16,18,24)`, and both target
profiles equal to

```
0,6,7,11,14,15,21,22,22,22,22,22,
24,24,24,24,30,30,30,30,30,34,34,36,
36,38,40,40,42,42,42,46,46,47,49,49.
```

Every original cell/shared cap and all 220 oriented inherited adjugate
inequalities hold exactly. At scale `1/6`, 100-bit outward dyadic
arithmetic proves critical pressure **less than `-1/250`**. Its numerical
value is approximately `-0.00467945`.

This is a point of necessary inequalities, not a tensor. A balanced
splitting of the paired energies suggests a degree-54 coefficient window
for an attempted realization. Neither feasibility of that realization
nor global optimality of this relaxed point is claimed.

Its saturated four-prefix and third prefix 13 still trigger the
[elementary parity-core obstruction](zero-energy-parity-core.md) if the
two zero-energy labels have opposite parity and their core preserves
parity. Same-parity zero labels or parity mixing in that core remain
available for this six-letter target.

The [discovery driver](../tools/search_integer_hereditary_q6.py) uses a
bounded integer convex search. Its solver status does not establish the
mathematical claim: the [checker](../tools/check_hereditary_adjugate_caps.py)
independently verifies every integer cap and the outward pressure bound.
The [target receipt](../evidence/hereditary-q6-integer-target.json) retains
the exact point and discovery source hashes.

## 5. Exact replay and limits

The [controls](../evidence/hereditary-adjugate-controls.json) include twelve
literal Laurent tensors over F10007 with independently split heights:
48 complete profiles and 992 cap checks. Proper core restrictions verify
96 normalization identities exactly. Another 9,320 rational comparisons
test that the smallest-energy subsets minimize the right sides. These
finite controls support the implementation; the restriction proof gives
the all-field theorem.

Run `make verify-hereditary-adjugate`. To repeat the optional discovery:

```
python tools/search_integer_hereditary_q6.py --output /tmp/hereditary-q6.json
```

The six-letter point requires a new crossing construction. A later
[clean full replay](../evidence/hereditary-full-verification.json) at exact
`d31ffec` passed defaults and all nine named joint/hereditary targets:
350 Python command blocks, 66 unit tests, five crossing tests and
984 Markdown files. It adds no Lean result or independent review claim.
