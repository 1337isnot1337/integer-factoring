# Uniform construction of the integral Schur crossing

## Result

The integral crossing in [integral-schur-crossing-pr30.md](integral-schur-crossing-pr30.md)
has a deterministic construction using `2^O(lambda)` bit operations for
external label width `lambda`. This includes constructing every piece and
both of its factor tables. Coefficient bit lengths are polynomial in
`lambda`; no oracle for symmetric-group characters or Schur functions is
required.

The construction uses the full integral Schur lattice and has crossing
multiplier one. Reduction to any fixed finite field is therefore effective,
including characteristics dividing the old factorial multiplier. The
declared channel exponent remains `rho=1/sqrt(2)`. This audit supplies an
implementation obligation for the ordinary transfer; it does not improve
the factoring exponent above one half.

The source operator recurrence and crossing identities are the pinned
Appendix B construction cited in the linked note. Its integral normalization
is proved there. The algorithms and cost accounting below are an additional
ordinary proof, developed with substantial OpenAI assistant support.

## 1. A finite ring containing every required Schur conversion

Write `p(h)` for the number of partitions of `h`, with `p(0)=1`, and

```math
N_H=\sum_{j=0}^{H}p(j)=\exp(O(\sqrt H)).                       \tag{1}
```

The estimate is the ordinary partition estimate already used in the
[half-rate construction](half-rate-factoring.md). Enumerate partitions
recursively and store a polynomial in the power sums by its partition
monomials. Work in the finite quotient

```math
R_H=\mathbb Q[p_1,\ldots,p_H]/
       (\text{all monomials of weighted degree greater than }H).
                                                                    \tag{2}
```

There are `N_H` monomials. A direct multiplication, merging two sorted
partition lists and discarding weights above `H`, takes
`O(N_H^2 poly(H))` rational arithmetic operations. Discarding such terms
is valid throughout: every variable has positive weight, so subsequent
addition and multiplication cannot move a discarded term to lower weight.

The complete functions have the explicit expansions

```math
h_j=\sum_{\mu\vdash j}\frac{p_\mu}{z_\mu},\qquad
z_\mu=\prod_k k^{m_k}m_k!,                                  \tag{3}
```

with `h_0=1`, and `h_j=0` for `j<0`. For `j<=H`, every denominator
in (3) divides `H!`. Their tables take `N_H poly(H)` bit operations.

## 2. Schur determinants without enumerating permutations

For a partition `lambda` of weight at most `H`, use its Jacobi--Trudi
matrix `A_(i,j)=h_(lambda_i-i+j)`. Its dimension is at most `H`.
The Leibniz sum would introduce an unnecessary factorial running time.
The following division-free recursion uses polynomially many ring operations.

For a matrix over any commutative ring, partition

```math
A=\begin{pmatrix}a&r\\c&B\end{pmatrix}.
```

Formal inversion of `I-tB` in the power-series ring gives

```math
\det(I-tA)=\det(I-tB)
 \left(1-ta-t^2\sum_{j\ge0}(rB^jc)t^j\right).                \tag{4}
```

This uses only a series with constant term one; it divides by no ring
element. To obtain coefficients through degree `n=dim A`, only the
moments `rB^jc` for `0<=j<=n-2` are needed. If
`det(I-tB)=sum_k b_k t^k`, the coefficient of `t^k` in (4) is

```math
b_k-a b_{k-1}-\sum_{j=0}^{k-2}b_{k-2-j}(rB^jc),              \tag{5}
```

where out-of-range coefficients are zero. Recursively construct the
characteristic coefficients of the trailing submatrix `B`, then obtain
the moments by repeated matrix-vector multiplication. Across all trailing
submatrices this uses `O(n^4)` ring operations. The coefficient of `t^n`,
multiplied by `(-1)^n`, is `det A=s_lambda`.

All weights in the ring are nonnegative, and every term of the desired
Jacobi--Trudi determinant has weight `|lambda|`. Consequently (2) gives
the exact answer. Computing every required Schur function, including all
weights through a fixed multiple of `h`, costs

```math
\operatorname{poly}(h)N_{O(h)}^3=\exp(O(\sqrt h))             \tag{6}
```

arithmetic operations. This deliberately loose bound already suffices.

### Bit lengths

Every intermediate expression in the recursion (4)--(5) is a sum of
products of at most `O(H)` entries of Jacobi--Trudi matrices. A moment
of length `j` has degree `j+2` in matrix entries, and its multiplier in
(5) brings the total degree to at most the current matrix dimension.
The vector-power intermediates satisfy the same `O(H)` bound. A common
denominator `(H!)^O(H)` therefore suffices.

Each input entry has at most `N_H` terms with coefficients of absolute
value at most one. The expanded arithmetic expressions have at most
`exp(poly(H))` terms: matrix paths give at most `H^O(H)` choices per
recursion stage, and there are `O(H)` stages. Both numerator and
denominator hence have polynomially many bits. For example,
`O(H^3 log(H+2))` is a sufficient loose bound. Exact rational arithmetic
therefore adds only a polynomial factor to (6).

## 3. Source modes and both crossing factorizations

The source recursion for the coefficient

```math
[p_\beta]\,p_\mu[|\beta|-|\gamma|]p_\gamma                  \tag{7}
```

removes one part of the parameter `mu` at every step. It either removes
one part of the output `beta`, with an integer binomial multiplier, or
one part of the input `gamma`, with an integer binomial and multiplicity
multiplier. The exact recurrence is implemented independently in
[`check_half_piece_ranks.py`](../tools/check_half_piece_ranks.py).

All crossing factors here need parameter weight at most `h`, external
state weight `h`, and intermediate weight at most `4h`. Memoization has
at most

```math
N_h N_{4h}^2=\exp(O(\sqrt h))                                \tag{8}
```

possible triples. Every update has polynomially many branches and recursion
depth at most `h`. Binomial coefficients and the remaining integer factors
have `poly(h)` bits. Summing over at most `exp(poly(h))` recursion paths
still gives polynomial coefficient bit lengths. Thus (7) is constructible
in `exp(O(sqrt(h)))` bit operations.

Apply the integral transform `T` by its explicit polynomial formulas from
the companion note. Truncated polynomial multiplication over (2) constructs
all transformed parameters in the same time bound. Expand Schur input
states by (3)--(6), apply the memoized modes, and convert output
coefficients back to Schur coordinates using Hall pairing. The latter
conversion is explicit: multiply a power-sum coefficient by
`z_mu [p_mu]s_lambda` and sum. Integrality of these character values
also follows from the usual integral Schur basis; the finite checker
independently verifies them at its tested weights.

This constructs each factor entry exactly over the rationals before
any reduction. The integral normal-symbol proof establishes that the
resulting factor entries are integers. Divisions used in this offline
construction are exact rational operations, not inversions of nonunits
in the target field.

There are `O(h)` pieces, at most `p(h)^2` external index pairs per factor,
and `O(h)N_{4h}` intermediate channels. All tables therefore have
`exp(O(sqrt(h)))` entries. Combining them, or computing smaller actual
ranks by field elimination, preserves this time bound. Reduction to a
fixed finite field adds only a polynomial bit cost.

## 4. Uniform label width and the transfer interface

For a binary label width `lambda>=1`, choose the least `h` such that
`p(h)>=2^lambda`, and take the first `2^lambda` partitions in a fixed
computable order. Partition dynamic programming finds this `h`; the
partition estimate gives `h=Theta(lambda^2)`. Finite exceptions are
absorbed in the fixed construction constant.

The preceding bounds become `exp(O(sqrt(h)))=2^O(lambda)`. The source
partition estimate and both factor lengths give

```math
\log_2 R_s+\log_2 C_s
 \le 2\sqrt2\,\lambda+O(\log(\lambda+1)).                   \tag{9}
```

Thus the construction meets the ordinary effective crossing interface
at the same leading parameter `rho=1/sqrt(2)`. If the interface requests
a rational slack parameter, choose any fixed rational `rho>1/sqrt(2)`
and absorb the logarithmic overhead for sufficiently large widths.
No stronger parameter is proved by changing the integral lattice.

The homogeneous prime-power variant needs a separate width analysis:
restricting weights to powers of the characteristic can overshoot the
least usable weight. This audit makes no assertion that the sparse-angle
variant automatically satisfies (9) at every label width.

## 5. Reproducibility

[`check_schur_construction.py`](../tools/check_schur_construction.py)
implements (4)--(5) in the weighted truncated ring and compares the
result with an independent permutation-based Jacobi--Trudi implementation
at small weights. The recorded run through weight eight checks 67 Schur
polynomials, 919 character entries, and 919 Hall pairings. It also includes
a truncation control and records arithmetic-operation and bit-size counts.

```sh
python3 tools/check_schur_construction.py --max-weight 8 \
  --output evidence/schur-construction.json
```

Those finite comparisons check the implementation. The uniform bounds
come from Sections 1--4 and do not rely on extrapolating the measurements.
