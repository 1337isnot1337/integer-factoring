# Matrix rank and affine recovery with polylogarithmic workspace

## Result and scope

An entry oracle for a matrix over a prime field supports deterministic
rank, affine consistency, and individual solution-coordinate queries in
quasipolynomial time and polylogarithmic auxiliary workspace. No array of
pivot indices or pivot-selection mask is needed. For a matrix with at
most `N` rows and columns over `F_q`, put `N>=2`. The bounds are

```math
T=N^{O(\log N)}\operatorname{poly}(\log q)(1+T_{\rm entry}),\qquad
S=S_{\rm entry}+O\bigl(\log N\,(\log N+\log q)\bigr).       (1)
```

The input is read only through its oracle. Output coordinates can be
streamed; an explicitly retained output vector would require its own
storage. For prescribed free coordinates, their oracle cost is included
in the same manner as the matrix and right-hand-side oracles.

This strengthens the `O(N)` pivot-mask workspace term in the earlier
[scalar Berkowitz note](small-space-berkowitz-rank.md). For
`N=O(n/log n)` and `q=n^O(1)`, the local query bound becomes

```math
T=\exp(O(\log^2 n))(1+T_{\rm entry}),\qquad
S=S_{\rm entry}+O(\log^2 n).                              (2)
```

The [digit-box](hankel-digit-box-contraction.md) and
[carry-lift](padic-digit-carry-transfer.md) constructions still have
their outer frequency and digit enumeration costs. Improving this
linear-algebra subroutine does not establish a higher factoring score.

The underlying polynomial-weighting idea belongs to the classical
parallel rank literature, notably Ketan Mulmuley, *A fast parallel
algorithm to compute the rank of a matrix over an arbitrary field*,
Combinatorica **7** (1987), 101--104,
[DOI 10.1007/BF02579205](https://doi.org/10.1007/BF02579205).
The determinant evaluator is the classical Berkowitz algorithm cited
in the earlier note. We give a self-contained version of the weighting
argument and explicit sequential entry-oracle, finite-field, and
coordinate-recovery accounting. This is an adaptation of those ideas,
not a claim to a new general rank algorithm. No new Lean formalization
is used.

## 1. A low-degree polynomial detects the full rank

Let `A` be an `m` by `n` matrix over an arbitrary field `F`. Use row and
column indices beginning at one and define diagonal polynomial matrices

```math
D_m(t)=\operatorname{diag}(t,t^2,\ldots,t^m),\qquad
D_n(t)=\operatorname{diag}(t,t^2,\ldots,t^n).
```

Set

```math
B(t)=D_m(t) A D_n(t) A^{\mathsf T}.                       (3)
```

Write the characteristic polynomial in the form

```math
\det(XI_m-B(t))=\sum_{k=0}^m c_k(t)X^{m-k}.
```

The principal-minor formula and Cauchy--Binet give

```math
c_k(t)=(-1)^k
 \sum_{\substack{I\subseteq[m],\ J\subseteq[n]\\|I|=|J|=k}}
       \det(A[I,J])^2\,t^{\sum_{i\in I}i+\sum_{j\in J}j}. (4)
```

Both identities are polynomial identities over the integers, so (4)
holds in every characteristic. In particular, it does not rely on
positivity of squares. Every coefficient has degree at most

```math
D=\frac{m(m+1)}2+\frac{n(n+1)}2\le N(N+1).               (5)
```

Let `r=rank A`. Then `c_k=0` for `k>r`, because all the corresponding
minors vanish. We next prove that `c_r` is not the zero polynomial.
The case `r=0` uses `c_0=1` and needs no separate detection step.

### 1.1 The greedy basis has a unique minimum weight

Order a spanning list of vectors as `v_1,...,v_s`. Scan it in that
order and retain a vector precisely when it increases the span. Let
the retained indices be `g_1<...<g_r`.

For any other basis with indices `b_1<...<b_r`, its first `j` vectors
are independent and lie in the prefix ending at `b_j`. That prefix
has rank at least `j`, so the greedy procedure has already selected
its `j`th vector by then. Hence

```math
g_j\le b_j\quad(1\le j\le r).
```

If the bases differ, at least one inequality is strict. Thus the
greedy basis is the unique basis minimizing the sum of its indices.
Arbitrary subsets can have equal sums; this argument uses the basis
property and proves uniqueness at the minimum.

Apply this fact separately to the row vectors and column vectors of
`A`. Let their unique minimum-weight bases be `I_*` and `J_*`.

### 1.2 Every full row-basis and column-basis pair gives a nonzero minor

Choose any rank factorization `A=UV`, where `U` is `m` by `r` of full
column rank and `V` is `r` by `n` of full row rank. The row bases of
`A` are exactly the nonsingular row minors `U[I,:]`; its column
bases are exactly the nonsingular column minors `V[:,J]`. Therefore

```math
\det A[I,J]=\det U[I,:]\,\det V[:,J]\ne0                 (6)
```

for every such pair. Conversely a nonzero `r` by `r` minor requires
both a full row basis and a full column basis.

Among the nonzero terms of (4) for `k=r`, the exponent is consequently
minimized uniquely at `(I_*,J_*)`. Its coefficient is
`(-1)^r det(A[I_*,J_*])^2`, which is nonzero over every field. There
can be no cancellation of that unique least-degree term. We have proved

```math
\boxed{\operatorname{rank}A
       =\max\{k:c_k(t)\not\equiv0\}.}                    (7)
```

The rank factorization and greedy bases occur only in the proof of
(7). Computing them is not a prerequisite for the rank algorithm.

## 2. Deterministic evaluation over a sufficiently large finite field

Now assume `F=F_q`, with `q` prime. A nonzero univariate polynomial
of degree at most `D` cannot vanish at `D+1` distinct field elements.
Choose an explicitly represented extension `K/F_q` with `|K|>D`.
Then (7) becomes the finite test

```math
\operatorname{rank}A
 =\max\{k:c_k(a)\ne0\text{ for some }a\in\mathcal A\},    (8)
```

where `mathcal A` is any list of `D+1` distinct elements of `K`.
All these elements can be generated by a counter. There is no
probability of failure and no random polynomial identity test.

If `q>D`, use `K=F_q`. Otherwise let `e` be the smallest integer
with `q^e>D`, and construct `K=F_q[T]/f(T)` for a monic irreducible
polynomial of degree `e`. Here

```math
q^e\le qD\le D^2,\qquad e\log_2q=O(\log N).              (9)
```

The harmless boundary case `D=0` is handled directly; nonempty
matrices have positive `D`. One can find `f` by enumerating all
monic degree-`e` polynomials and rejecting a candidate if any monic
polynomial of degree at most `floor(e/2)` divides it. Every reducible
polynomial has such a factor. Irreducibles of every positive degree
exist over a finite field, so this search terminates. For completeness,
the identity `X^(q^e)-X=prod_(d|e) prod_(f irreducible, deg f=d) f`
and the degree count show there is a degree-`e` irreducible. Applying
the identity at degree `d` bounds `d` times the number of degree-`d`
irreducibles by `q^d`. The proper divisor degrees are at most `e/2`,
and the total degree of all monic
irreducibles of degree at most `e/2` is at most
`sum_(d<=e/2) q^d < q^e` when `e>=2`; degree one is immediate.

There are only polynomially many candidates and trial divisors in
`N`, by (9). Polynomial division needs a constant number of arrays
of `O(e)` base-field coefficients. This costs `O(e log q)` bits,
and the construction time is `N^O(1) poly(log q)`. The promise that
the input coefficient field is the prime field `F_q` is part of the
interface; no slow trial-division primality test is included.

In either case a `K` element occupies

```math
b=O(\log N+\log q)                                      (10)
```

bits. Arithmetic uses a polynomial basis and ordinary modular
arithmetic, in `poly(b)` time and `O(b)` temporary space. A fixed
number of polynomial arrays of length `e` is still `O(b)` bits.
We can use one extension large enough for the maximum dimensions
throughout all the smaller prefix queries below.

### An extension can be necessary

For the all-ones `2` by `2` matrix over `F_2`, the rank is one but
the only potentially nonzero lower characteristic coefficient of
(3) is

```math
c_1(t)=t^2+t^4.
```

It vanishes at both elements of the base field. A test confined to
`F_2` would return the wrong rank. A field with more points than the
degree bound prevents this error. Likewise one must inspect the
characteristic coefficients, rather than assume that the rank of
one evaluated `B(a)` determines the answer; an evaluated matrix can
be nonzero and nilpotent.

## 3. Scalar characteristic coefficients give the space bound

For a fixed evaluation point `a`, a requested matrix entry is

```math
B(a)_{ij}=a^i\sum_{k=1}^n A_{ik}a^k A_{jk}.             (11)
```

The sum and powers are computed with counters, a constant number of
field elements, and sequential calls to the original entry oracle.
An entry costs `N^O(1) poly(b)(1+T_entry)` time and
`S_entry+O(log N+b)` workspace. No copy of `B(a)` is stored.

Use the balanced scalar Berkowitz recursion from
[small-space-berkowitz-rank.md](small-space-berkowitz-rank.md) to
compute one `c_k(a)`. It uses `O(log N)` frames, each holding a
constant number of indices and field elements. Its time is
`N^O(log N) poly(b)(1+T_entry)`, and its workspace is
`S_entry+O(log N(log N+b))`. The exponentiation, matrix-power, and
Toeplitz-product stacks are nested a fixed number of times; their
depths add.

Scan all `D+1` evaluation points and all `k<=min(m,n)`, retaining
only the largest successful `k`. This adds polynomial time and a
constant number of counters. Substituting (10) gives (1). Zero rows,
zero columns, and the zero matrix are handled by the same formula
with the empty characteristic coefficient or by direct base cases.

## 4. Recover basis indices when requested

Define the row-prefix ranks

```math
r_i=\operatorname{rank} A[1:i,1:n],\qquad r_0=0.
```

The greedy row basis consists precisely of the indices where
`r_i>r_(i-1)`. Its `j`th index is

```math
I(j)=\min\{i:r_i\ge j\}.                               (12)
```

Each requested prefix rank is computed by Section 3 with an oracle
that checks the prefix limit and then calls the original `A` oracle.
Scanning for the first successful `i` needs only a counter. It is
also possible to binary-search the monotone ranks, but a linear scan
already fits the time bound. Define the column selector `J(j)` in
the same manner using column-prefix ranks.

These selectors reconstruct the two bases independently. By (6),
their intersection

```math
U_{ij}=A_{I(i),J(j)}\qquad(1\le i,j\le r)                (13)
```

is nonsingular. An entry query for `U` computes its two selected
indices, reads the original matrix entry, and discards the indices.
No `N`-bit mask and no list of `r` indices is retained.

## 5. Affine consistency and requested solution coordinates

Test `Ax=b` by comparing the ranks of `A` and `[A|b]`; the maximum
dimension increases by at most one. Suppose they agree. Set every
nonpivot coordinate to zero, and solve

```math
U x_J=b_I.                                             (14)
```

The pivot rows span the full row space, so these equations together
with consistency imply all the original equations. For a requested
coordinate `l`, compare column-prefix ranks before and after `l`.
If no rank increase occurs, return its prescribed free value. If
the rank increases to `j`, then `l=J(j)`. Cramer's rule gives its
value as the determinant of `U` with column `j` replaced by `b_I`,
divided by `det U`.

Both determinants use scalar Berkowitz over the original prime
field. Their virtual entry oracles call the selectors (12). A selector
calls the rank procedure of Section 3, which calls only the original
entry oracle. There is no recursive pivot search inside itself:

The selected-minor determinant calls a prefix-rank selector. The
selector calls weighted characteristic coefficients, whose entries
come directly from the original matrix. This is a fixed nesting
depth. A product of a fixed number of
`N^O(log N)` time bounds remains `N^O(log N)`; their simultaneously
active scalar recursion stacks add only
`O(log N(log N+log q))` workspace. The nonzero base-field denominator
can be inverted by binary powering to `q-2` with `O(log q)` bits.

If nonpivot values are supplied by a separate oracle, replace `b_I`
by

```math
b_{I(i)}-\sum_{l\notin J}A_{I(i),l}x_l.                 (15)
```

Each membership test uses two prefix ranks. The sum is streamed,
so (15) adds polynomial time, a field accumulator, and the stated
free-coordinate oracle costs. A full solution can be emitted one
coordinate at a time with an additional factor `n` in time.

## Verification and limitations

The companion exact checker
[`check_polylog_space_rank.py`](../tools/check_polylog_space_rank.py)
verifies the weighted-minor identity,
rank detection, the necessity of the extension in small characteristic,
prefix-basis selection, and affine coordinate recovery against ordinary
elimination. The test matrices and reference elimination may use dense
storage. The workspace theorem follows from the explicit recomputation
algorithm above, not a memory measurement of Python reference objects.

Only the matrix-query workspace is reduced. A caller that explicitly
enumerates all digit assignments or frequencies must still charge its
own counters, input access, arithmetic, and total enumeration time.
The unconditional factoring saving exponent remains `1/2`.
