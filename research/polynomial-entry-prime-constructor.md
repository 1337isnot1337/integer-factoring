# Polynomial-time entries for the adaptive prime constructor

## Result

The [width-sensitive prime schedule](width-sensitive-prime-constructor.md)
has an elementary polynomial-time algorithm for each unweighted local
table entry. Exact closest-vector search is unnecessary for that
schedule. The algorithm uses only bounded digit enumeration, modular
Vandermonde interpolation, and the forced internal arithmetic checks.

For each fixed rational `0<epsilon<=1`, the resulting constructor has

\[
 P=\Theta_\epsilon(\sqrt{n/\log n}),\quad
 \lambda=\Theta_\epsilon(\log n),\quad
 P^2\lambda<(12+\epsilon)n,\quad B_{\rm tab}=1,            \tag{1}
\]

and writes all tables in

\[
 \operatorname{poly}_\epsilon(n)\,
                 2^{(P-1)\lambda}.                       \tag{2}
\]

The polynomial factor includes prime search, root construction over the
integer modulus, and every table-entry computation. With local algebra
weights, each entry requires a polynomial number of coefficient
operations as well as polynomial bit work on the integer labels; the
cost of the chosen coefficient algebra is charged separately.

The same elementary decoder gives `exp(O((log n)^4))` entry time in the
earlier [explicit narrow-prime schedule](prime-transform-boundary-constructor.md),
and a finite written-table bound stated in Section 3. Thus that schedule
also has a complete alternative which does not invoke closest-vector
algorithms.

## 1. Enumerate only the unobserved high digits

Use the common-index grouping and boundary labels of the
[observable constructor](observable-boundary-constructor.md). A group
contains `k` first-transform indices and has `A-k` external recipients.
At a fixed internal residue index, let the two operand heights be
`a_p,a_q`. The operand-length construction gives

\[
 a_p+a_q\le A,\qquad 0\le a_p,a_q\le A-1.                 \tag{3}
\]

An operand polynomial of height `a` has
`u=min(a,A-k)` coefficients determined by any selected `u` external
observations once its highest `a-u` coefficients are fixed. Enumerate
those high digits in their actual bounded domains, subtract their
observed contributions, and solve the square Vandermonde system.
Its determinant is a unit: the evaluation points are distinct powers
of the principal `A`th root modulo the prime `Q`.

Because `k<=A/2`, at most one of `a_p,a_q` can exceed `A-k`. Its excess
is at most `k-1`, by (3). Summing over the `k` internal residue indices
leaves at most

\[
 k(k-1)\le k^2                                           \tag{4}
\]

free radix digits. Each lies in `[0,2^w)`. The entire local search costs

\[
 2^{k(k-1)w}\operatorname{poly}(n)
       \le2^{k^2w}\operatorname{poly}(n).                 \tag{5}
\]

The polynomial factor includes exact modular inversion, interpolation,
domain checks, and downstream evaluation. There is no lattice oracle.
The smaller operand, shortened final digits, leading-bit constraints,
and fixed prefix can only reduce the search.

Every remaining internal variable is forced from these inputs and the
boundary: the middle and inverse transforms are deterministic, and each
internal carry run has its incoming carry supplied by its boundary or
by the initial zero. Reject out-of-range residue codes, nonzero padding,
violated extra observations, invalid digit domains, excessive outgoing
carries, and failed target digits. This is exactly the complete decoder
already proved in the observable constructor.

Uniqueness is separate from the enumeration bound. The cyclotomic norm
argument in the prime constructor gives observation-lattice minimum
at least `sqrt(Q/2)>2^w A/2`. Therefore a contiguous hidden block of size
at most `A/2` admits at most one bounded input polynomial of height
at most `A-1`. All accepted internal settings are consequently unique,
and the unweighted entry is zero or one. The argument uses the norm
bound, but no algorithm for closest vectors.

For weighted tables, return zero or the product of the weights on that
unique setting. The digit enumeration and all acceptance decisions are
independent of the coefficient values, so this is a uniformly generated
arithmetic circuit over any of the allowed commutative algebras.

## 2. Constant block size makes entries polynomial

In the width-sensitive schedule the block length is the fixed power
of two `K>=128/epsilon`, and `w=Theta_epsilon(log n)`. Thus

\[
 k=K,
 \qquad k(k-1)w=O_\epsilon(\log n).                       \tag{6}
\]

Equation (5) is polynomial in `n`, with an exponent depending only on
the chosen fixed `epsilon`. All arithmetic and no-aliasing hypotheses,
the prime search, and the volume bound are unchanged. There are exactly
`P 2^((P-1)lambda)` padded boundary labels to write. Multiplying their
count by the polynomial entry cost proves (2), and hence the full
`exp(O_epsilon(sqrt(n log n)))` bound.

This strengthens the local entry algorithm in the previous theorem.
It does not shrink the number of labels, and does not change the
crossing exponent used after construction.

## 3. A finite bound for the earlier explicit schedule

For the schedule valid for every `n>=2^40`, put

\[
 u=\lceil\log_2(n+2)\rceil,\qquad
 P=\left\lfloor\frac{\sqrt n}{16u^2}\right\rfloor,
 \qquad A=\operatorname{ceilpow2}(Pu).
\]

Let `k=ceil(A/P)`. Since `A<2Pu`,

\[
 k\le2u.                                                \tag{7}
\]

One has `P>=sqrt(n)/(32u^2)`, `A>=Pu`, and `A>=1024`. The rounded
digit width therefore satisfies

\[
 w<\frac{n+1}{A^2-A-1}+2<2050u^2.                        \tag{8}
\]

This proves `k^2w=O(u^4)` and the stated quasipolynomial entry time.
For an explicit uniform estimate, put `F=k^2w`. Since
`n/P>=16u^2 sqrt(n)`,

\[
 \frac{FP}{n}<\frac{1025}{2}\frac{u^2}{\sqrt n}<1.
                                                               \tag{9}
\]

Here is a direct check of the last inequality over the entire range.
If `u=41`, use `sqrt(n)>=2^20`. If `u>=42`, use
`n>=2^(u-2)` and the fact that `u^2/2^((u-2)/2)` decreases for
integer `u>=42`. Its value at 42 gives the rational upper bound

\[
 \frac{FP}{n}<\frac{904050}{1048576}<1.                   \tag{10}
\]

The same monotonic comparisons prove the floor estimate for `P` used
above: `sqrt(n)/(16u^2)>=2`, so its floor is at least half its value.

There is a sharper full-table estimate than (9) alone supplies. Since
`f>2w` and `lambda>=6k^2 f`, one always has `F<lambda/12`. Thus
`(P-1)lambda+F<P lambda`. The existing finite volume bound is
`P^2lambda<13n`, so all tables can be written, after prime preparation, in

\[
 \operatorname{poly}(n)\,2^{(P-1)\lambda+F}
       <\operatorname{poly}(n)\,2^{13n/P}.                \tag{11}
\]

Prime preparation retains its separately charged `2^O(u^2)` cost.
No exact closest-vector routine is used in (11).

## 4. Exact replay

The companion checker
[`check_polynomial_entry_prime_constructor.py`](../tools/check_polynomial_entry_prime_constructor.py)
compares interpolated entries with complete local assignment enumeration
over small prime-modulus networks, including unbalanced operand lengths,
mutated boundary labels, and invalid padding. It separately checks the
free-digit count and the exact finite schedule inequalities. Run:

```sh
python tools/check_polynomial_entry_prime_constructor.py --output evidence/polynomial-entry-prime-constructor.json
```

The committed [evidence](../evidence/polynomial-entry-prime-constructor.json)
records 36 genuine entries and 72 changed-boundary entries checked against
55,902 complete local assignments, including 12 cases with a positive
free-digit count. All valid entries are one. It also records 93,406
paired-height controls, 5,121 explicit-schedule checks, 119 constant-block
checks, and rejection of 36 invalid padding and 36 invalid residue labels.
All comparisons passed.
