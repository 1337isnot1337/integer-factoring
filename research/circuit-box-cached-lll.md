# Exact cached LLL preserves every small-root certificate

## Result and scope

The reference [quarter-bit endpoint](circuit-box-ratio-small-roots.md)
recomputes rational Gram--Schmidt after every lattice row reduction and
adjacent swap. The new [cached engine](../tools/circuit_box_cached_lll.py)
updates those same rational values exactly. It uses the same nearest-integer
rule, row order, Lovasz comparison, and early stopping test.

Consequently, on every valid reference query, both implementations return
the same short polynomial, integer lattice combination, reduction count,
and swap count. The [centered-window wrapper](circuit-box-centered-windows.md)
also returns the same complete result. The change removes repeated full
orthogonalization; it does not change a small-root threshold or the complete
factoring exponent.

These are classical LLL update identities. The underlying source is
A. K. Lenstra, H. W. Lenstra, and L. Lovasz,
[*Factoring polynomials with rational coefficients*](https://ir.cwi.nl/pub/9304),
1982, Section 1, especially the updates in (1.22) and the cost analysis
in (1.26). The proof below derives the exact identities needed for
this implementation and its compatibility with the retained certificates.
No novelty claim is made for maintaining Gram--Schmidt data during LLL.

## 1. The cache invariant

Let the independent integer rows be
`b_0,...,b_(d-1)` in `Z^w`, where `d<=w`. Write their exact orthogonalization
as

```math
b_i=b_i^*+\sum_{j<i}\mu_{i,j}b_j^*,\qquad
B_i=\langle b_i^*,b_i^*\rangle>0.                 \tag{1}
```

The cache holds all the rational `mu_(i,j)` and positive rational `B_i`.
Its invariant is equality with a fresh exact computation of (1), rather
than an approximate or merely bounded value. The integer basis and its
integer transformation matrix are maintained separately.

At initialization, the existing rational Gram--Schmidt routine constructs
the cache. The following two local identities preserve its invariant.

### A size reduction changes one coefficient row

For `j<k`, replace `b_k` by `b_k-q b_j`, with integer q. Every prefix span
is unchanged. Thus all orthogonal vectors and all their squared lengths
remain unchanged. Expanding the replacement in (1) gives

```math
\begin{aligned}
\mu'_{k,\ell}&=\mu_{k,\ell}-q\mu_{j,\ell},&&\ell<j,\\
\mu'_{k,j}&=\mu_{k,j}-q.
\end{aligned}                                             \tag{2}
```

Every other coefficient is unchanged. Later basis rows are unchanged and
their orthogonal coordinate vectors are unchanged, so their coefficients
need no repair.

The implementation uses exactly the reference quotient

```math
q=\left\lfloor\frac{2a+b}{2b}\right\rfloor
\quad\text{when}\quad\mu_{k,j}=a/b,\ b>0.                  \tag{3}
```

This rounds a half-integer toward the greater integer, including a negative
half-integer. Keeping that tie convention is necessary for exact output
equivalence.

### An adjacent swap changes two orthogonal directions

Now interchange `b_(k-1)` and `b_k`. Put

```math
\mu=\mu_{k,k-1},\qquad A=B_{k-1},\qquad C=B_k,
\qquad D=C+\mu^2A>0.                                      \tag{4}
```

After projecting away all earlier directions, the old two rows are
`u=b_(k-1)^*` and `v+mu u`, where `v=b_k^*` and `u` is orthogonal to `v`.
The new first direction is `u'=v+mu u`. Its squared length is D. The new
second direction is

```math
v'=u-\nu u',\qquad \nu=\frac{\mu A}{D}.                    \tag{5}
```

It is orthogonal to `u'`. Direct calculation gives

```math
B'_{k-1}=D,\qquad B'_k=\frac{AC}{D},\qquad
\mu'_{k,k-1}=\nu.                                         \tag{6}
```

The coefficients in columns less than `k-1` are exchanged between the two
rows. For a later row `i>k`, let
`alpha=mu_(i,k-1)` and `beta=mu_(i,k)`. Since

```math
\alpha u+\beta v
 =\bigl(\beta+\nu(\alpha-\mu\beta)\bigr)u'
   +(\alpha-\mu\beta)v',                                  \tag{7}
```

the required later-row updates are

```math
\begin{aligned}
\mu'_{i,k}&=\alpha-\mu\beta,\\
\mu'_{i,k-1}&=\beta+\nu\mu'_{i,k}.
\end{aligned}                                             \tag{8}
```

All other entries stay fixed. In particular, later orthogonal vectors
depend on the same preceding span, which the adjacent swap preserves.
Equations (2) and (6)--(8) prove the cache invariant at every operation.
Their denominators are nonzero because the input rows are independent.

## 2. Identical decisions imply identical certificates

Both engines begin with the same integer rows and identity transformation
matrix. Suppose their states agree before an operation. The cache invariant
shows that the rational coefficient queried by the next reduction equals
the reference recomputation. Equation (3) then gives the same integer
quotient. Both engines perform the same row operation on the basis and
transformation matrix.

After reducing the current row, both compare exactly

```math
B_k\ \ge\ \left(\frac34-\mu_{k,k-1}^2\right)B_{k-1}.       \tag{9}
```

They therefore make the same swap decision and the same index update.
Induction proves equality of the entire integer basis, transformation
matrix, loop index, and both operation counters after every iteration.

For a supplied strict first-row norm threshold T, the cached routine tests
`sum(abs(b_0))<T` at the same position in the loop as the reference routine.
It does not move that test before the first iteration. Both therefore stop
on the same row and return the same integer combination of the original
rows. If the threshold is not attained after full reduction, the routine
rejects the query instead of treating a longer row as a valid polynomial.

The reference small-root builder then performs identical coefficient
unscaling, and the unchanged exact Sturm routine sees the same polynomial
and interval. All divisor tests and all returned counters consequently
agree. This includes the centered coordinate change and its negative root
intervals.

The generic cached routine also accepts `threshold=None`. It then completes
ordinary exact LLL, and returns the full reduced basis and integer
transformation. This mode is used for independent finite controls.

## 3. Charged operations and unchanged bit bounds

A fresh Gram--Schmidt computation on d rows of width w uses
`O(d^2 w)` rational arithmetic operations. The reference engine makes
`1+S+R` such computations, where S and R are its swap and nonzero reduction
counts. With auditing disabled, the cached engine makes one full computation.

Each later size reduction updates at most d cached coefficients. Each swap
updates at most `2d-1` cached entries, including its two squared lengths.
The integer basis-row and transformation-row operations cost `O(w+d)`
integer operations per reduction. These costs remain charged; the cache
does not make arbitrary-precision arithmetic unit cost.

Thus the repeated orthogonalization work is replaced by

```math
O\bigl(d^2w+(S+R)d\bigr)                                  \tag{10}
```

rational operations. Integer row operations take
`O((S+R)(w+d))` integer operations, and the full final transformation check
takes `O(d^2 w)` integer operations. Storage for the basis, transformation, and cache
is `O(dw+d^2)` entries. A final direct verification checks every transformed
row against the original basis, including all w coordinates.

The maintained fractions are exactly the ordinary Gram--Schmidt fractions
of the current LLL basis. Their numerator and denominator lengths therefore
have the same polynomial bounds as exact LLL, obtainable from Gram minors
and the usual basis-length bounds. The new recurrence does not introduce
an additional sequence of approximate values or larger unreduced formal
fractions: Python's `Fraction` normalizes each rational result.

This is a finite-operation improvement within a classical polynomial-time
subroutine. The full centered factor search still has
`2^(n/4) poly(n)` cost. The repository's stronger complete score bound is
unchanged, and collective inverse-window aggregation remains unresolved.

## 4. Isolated adapters and complete replay

The module exports the same small-root and centered-window query entry
points with cached lattice reduction. It clones the existing function code
objects into private namespace dictionaries and changes the lattice solver
dependency there. It does not mutate either reference module, the reference
checkers, or their retained evidence.

The [checker](../tools/check_circuit_box_cached_lll.py) performs three kinds
of comparison:

- On complete small LLL runs, it independently recomputes every changed
  Gram--Schmidt state, checks the final LLL inequalities, and verifies that
  the integer transformation has determinant plus or minus one. Rectangular
  row sets and all positive/negative half-integer tie cases are included.
- It repeats thresholded queries with the original reference solver and
  requires equality of the complete short-vector tuple, including the
  integer combination and both counters.
- It runs both original complete small-root checker code objects with
  cached query dependencies, and requires equality of the entire resulting
  dictionaries with the previously retained JSON receipts. All ordinary
  literal-factor, parameter, integer-root, large local, and complete small
  factor controls in those suites are included.

The [deterministic receipt](../evidence/circuit-box-cached-lll.json) records
the exact source and reference-receipt hashes. Native CPU and wall timings
are printed separately; they are not part of the reproducible result and
are not interpreted as a uniform bit-complexity bound.

In the retained 64-bit centered local query, the basis has dimension 19,
with 819 swaps and 2,550 reductions. Its count of full Gram--Schmidt
computations changes from **3,370 to one**, while the entire result stays
identical. The arithmetic updates in (2) and (6)--(8) replace the other
3,369 computations.

```sh
python3 tools/check_circuit_box_cached_lll.py \
  --output evidence/circuit-box-cached-lll.json
```

This implementation and analysis used substantial OpenAI assistant support.
