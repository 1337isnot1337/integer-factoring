# Quasipolynomial-time matrix queries with linear auxiliary space

## Result and use

The later [polynomial-weighting construction](polylog-space-rank-recovery.md)
removes the linear pivot-mask term below. It retains quasipolynomial time
and gives `O(log m (log m+log q))` auxiliary bits for rank and requested
affine coordinates. The scalar Berkowitz identities and the simpler
mask-based algorithm in this note remain valid.

An entry oracle for a matrix over a prime field is enough to compute its
rank, retain a nonsingular pivot minor, test affine consistency, and
recover requested solution coordinates in quasipolynomial time and
linear auxiliary space in the matrix order. The characteristic can be
arbitrarily small.

For at most `m` rows and columns over `F_q`, the bounds are

```math
T=m^{O(\log(m+2))}\operatorname{poly}(\log q)(1+T_{\rm entry}),
\qquad
S=S_{\rm entry}+O\bigl(m+\log(m+2)\log(mq+2)\bigr).
\tag{1}
```

Here the input matrix is read-only through its oracle. Its query time
and workspace are charged explicitly. No dense matrix, inverse matrix,
or list of pivot indices is retained. For determinant alone, the `m`
term in the auxiliary space bound is unnecessary.
Take the size parameter `m>=1`; a zero-dimensional minor returns one.

This improves the local `2^O(m)` procedures in the
[digit-box Fourier note](hankel-digit-box-contraction.md) and the
[small-characteristic carry note](padic-digit-carry-transfer.md).
With `m=O(n/log n)` and `q=n^O(1)`, the new query time is
`exp(O(log^2(n+2)))`, with `O(n/log n)` additional working bits.
It does not remove the exponentially many outer Fourier frequencies
proved in those notes, and gives no new factoring score.

The determinant identity is the standard division-free algorithm of
S. J. Berkowitz, *On computing the determinant in small parallel time
using a small number of processors*, Information Processing Letters
18 (1984), 147--150. The same rectangular Toeplitz identity and reference
appear in the inspected [SymPy implementation at commit
52de66e](https://github.com/sympy/sympy/blob/52de66ee854dade90871fb39072c0c2296b2b588/sympy/matrices/determinant.py).
The contribution here is explicit entry-query, bit-space, and greedy
pivot accounting for this factoring research interface. The identities
and costs are proved below; no new determinant algorithm is claimed.

## 1. The rectangular Toeplitz product

For an `s` by `s` matrix, split off its first row and column:

```math
A=\begin{pmatrix}a&R\\C&M\end{pmatrix}.
```

Define `b_0=1`, `b_1=-a`, and

```math
b_j=-R M^{j-2}C\qquad(2\le j\le s).
```

Let `B(A)` be the `(s+1)` by `s` lower-triangular Toeplitz matrix
with entry `b_(i-j)` when `i>=j` and zero otherwise, with indices
starting at zero. If `c(A)` lists the coefficients of
`det(XI-A)` from the leading coefficient to the constant term, then

```math
c(A)=B(A)c(M).                                      \tag{2}
```

Here `c(empty)=(1)`. Thus, for the trailing principal submatrices
`A_i=A[i:m,i:m]`,

```math
c(A)=B(A_0)B(A_1)\cdots B(A_{m-1})(1).              \tag{3}
```

This identity is valid over every commutative ring. One proof first
uses the universal matrix of independent integer indeterminates.
Over its rational function field, the Schur complement identity is

```math
\det(XI-A)=\det(XI-M)
             \bigl(X-a-R(XI-M)^{-1}C\bigr).
```

Expand the inverse formally at infinity as
`sum_(k>=0) M^k X^(-k-1)`. The nonnegative powers of `X` give
exactly (2); only `k<=s-2` can contribute. This is an equality of
integer polynomials in the matrix entries. It therefore specializes
to every commutative ring, including rings with zero divisors.
Neither formal inverses nor division are used by the algorithm.

The determinant is `(-1)^m` times the last coordinate of (3).

## 2. Query a product entry without materializing its factors

A requested entry of `M^k` is computed by balanced recursion:

```math
(M^k)_{uv}=\sum_w
 (M^{\lfloor k/2\rfloor})_{uw}
 (M^{\lceil k/2\rceil})_{wv}.                        \tag{4}
```

The cases `k=0,1` use the identity matrix and the input entry oracle.
For `k<=m`, the depth is `O(log(m+2))`. At a node retain the
indices, the next summation index, one accumulator, and the value
of the first child while recomputing the second. This uses a constant
number of field elements and `O(log m)` index bits per frame.
The entry-query count is at most `(2m)^(O(log(m+2)))`.

An entry of a factor `B(A_i)` is zero, one, `-A_ii`, or a double
sum of input entries times one power entry from (4). It has the
same time scale and recursion-space bound.

Now compute one entry of a consecutive product
`B(A_l)...B(A_(r-1))` by splitting the factor interval in half,
and summing products of the two requested half-product entries.
The intermediate dimension is at most `m+1`. This second recursion
also has depth `O(log(m+2))` and at most
`(2m+2)^(O(log(m+2)))` leaves.

At one leaf of the product recursion, a Toeplitz entry invokes the
power recursion. The two stacks are nested, so their depths add.
There is no retained table of factor entries or power entries.
Each frame stores `O(log m+log q)` bits. Consequently one requested
coefficient of (3), including the determinant, has

```math
T_{\det}=m^{O(\log(m+2))}\operatorname{poly}(\log q)
                    (1+T_{\rm entry}),
\quad
S_{\det}=S_{\rm entry}
             +O(\log(m+2)\log(mq+2)).                \tag{5}
```

The polynomial field-operation overhead can use ordinary binary
addition and multiplication modulo the prime `q`. Each operation
uses `O(log q)` working bits. The determinant computation needs no
inversion, so no condition such as `q>m` appears.

All coefficients can be emitted one at a time with an extra factor
`m+1` in time and no coefficient-vector storage. Formula (5) is an
ordinary deterministic finite-alphabet algorithm, obtained by
recomputation rather than a parallel-machine assumption.

## 3. Greedy pivot minors avoid exhaustive minor enumeration

Maintain two subset masks `I,J`, initially empty, with the same
cardinality `r` and `det A[I,J] != 0`. The empty determinant is one.
Scan every row `i` outside `I` and every column `j` outside `J`.
If

```math
\det A[I\cup\{i\},J\cup\{j\}]\ne0,                \tag{6}
```

add these indices to the masks and restart the scan. If no such pair
exists, return `r` and the masks.

To prove correctness, temporarily order the pivot rows and columns
first and write the matrix in blocks as

```math
\begin{pmatrix}U&V\\W&Z\end{pmatrix},\qquad \det U\ne0.
```

Its rank is `r+rank(Z-WU^(-1)V)`. Every entry of this Schur
complement is, up to a nonzero determinant and a sign, an augmented
minor in (6). If the rank exceeds `r`, some scan entry is nonzero.
If all are zero, the rank is exactly `r`. The inverse occurs only
in this proof; the algorithm computes the determinants directly.

There are at most

```math
\sum_{r=0}^{m-1}(m-r)^2=O(m^3)                       \tag{7}
```

determinant queries, rather than `4^m` minor candidates. Each mask
has `m` bits. To answer an entry query for a selected minor, scan
the masks to find the requested selected row and column. This takes
`O(m)` extra time and `O(log m)` extra workspace, already absorbed
by (1). The sorted mask order is used consistently, so no array of
`m` indices of `log m` bits is necessary.

The same argument works for rectangular matrices, with their respective
row and column masks. For affine consistency, compare the rank of
`C` with that of `[C|b]`; increasing the dimension bound by one is
harmless. Retain a pivot minor for `C` when consistency holds.

## 4. Recover individual coordinates and retain digit constraints

For a fixed choice of nonpivot coordinates in `Cx=b`, a requested
pivot coordinate is a ratio of two pivot-minor determinants by
Cramer's rule. Their order is at most `m`. The denominator is
nonzero by the stored pivot certificate. Inversion in `F_q` can use
binary powering to exponent `q-2`, using `O(log q)` space.

Thus consistency, rank, and coordinate recovery all satisfy (1).
The digit-box contraction still has its two complementary candidate
counts `B^(s-r)` and `q^r`; this note reduces their local linear
algebra overhead. It does not remove the digit-counter or frequency-
counter storage charged in the earlier Fourier analysis.

Likewise, the lifted-digit method uses the rank and pivot system
modulo its base prime, then checks the original higher-precision
equations. This procedure never treats a nonzero nonunit of a
prime-power ring as an invertible field pivot.

## Verification

The standard-library checker
`tools/check_small_space_berkowitz.py` compares every requested
characteristic coefficient with independent principal-minor formulas,
checks determinants and greedy ranks against elimination, and tests
affine consistency and coordinate recovery. It includes small
characteristics and singular matrices with delayed pivots. Its
reference matrices are stored for independent finite verification;
they are not evidence of the stated packed workspace bound. That
bound follows from the explicit recursion and mask accounting above.

No Lean files are changed. This local algorithmic improvement leaves
the established unconditional factoring score at `1/2`.
