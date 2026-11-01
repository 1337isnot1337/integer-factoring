# Polynomial joint frames of the four-letter reflection tensor

## Exact result

Let `T(z)` be the binary four-letter reflection tensor with physical
heights `(0,3,1,4)`, as in the
[constant joint-conjugation note](joint-height-conjugation-floor.md).
View its row and column physical pairs as sixteen-dimensional spaces.
Let `N` be any constant binary matrix preserving total physical height
and satisfying `N^2=0`. For any polynomial `f` over `F_2` with `f(1)=0`,
put

\[
 U(z)=I+f(z)N,\qquad T'(z)=U(z)T(z)U(z)^{-1}.
 \tag{1}
\]

The inverse is again `I+fN`. The tensor `T'` is polynomial, has the
literal crossing endpoint `T'(1)=swap`, and has the same zero-endpoint
Smith profile with the original height gauge. No commutation with swap
is assumed for `N`.

**Theorem for this seed and this binary family.** If `T'` is regular at
infinity with the original height gauge and scalar shift, then `T'=T`.
Thus this family cannot improve the two-cut pressure while retaining
those endpoint charges. The theorem allows arbitrary polynomial degree
for `f`; its proof reduces the check to degree at most eight.

This is stronger than a failure to find an improving frame on a scale
grid. It is also narrower than a general restriction on polynomial
conjugations. Frames over extension fields, sums of several independent
nilpotent directions, other heights, or frames requiring a new pole
charge are outside the statement.

## 1. The endpoint calculation

Let `h(x)` be the given height. Before the opposite-cut reshuffle, write

\[
 W_{xy,uv}(w)=w^{h(u)+h(v)+4}T_{xy,uv}(1/w).
 \tag{2}
\]

This is a polynomial matrix of degree eight. Height preservation means
that `N` commutes with the diagonal matrices implementing both endpoint
gauges. At zero, `U` and its inverse are regular units, so they leave
the gauged Smith profile unchanged. At infinity,

\[
 W'=W+f(1/w)C(w)+f(1/w)^2P(w),
 \qquad C=NW+WN,\quad P=NWN.
 \tag{3}
\]

Regularity of the physical matrix in (3) is equivalent to regularity
after reshuffling: the reshuffle merely permutes entries. The desired
condition is that every negative coefficient of the correction in (3)
vanish.

If `C=P=0`, every polynomial frame leaves `W`, and therefore `T`,
unchanged. Otherwise let `c=ord(C)` and `p=ord(P)` when the respective
matrix is nonzero. Each finite order is at most eight.

Suppose `f` is nonconstant of degree `L`. The least orders of the two
correction terms are `c-L` and `p-2L`. If they differ, a negative smaller
order cannot cancel. If they agree and are negative, equality forces
`L=p-c`. If both are nonnegative, they separately bound `L` by eight.
If only one correction term is present, its regularity likewise bounds
`L` by eight. Consequently every nontrivial admissible frame has

\[
 1\le L\le8.
 \tag{4}
\]

A constant polynomial with `f(1)=0` is zero. Also, for any positive
degree `L`, the coefficient at `w^{-2L}` in (3) is the leading scalar
coefficient squared times `N W(0) N`. Thus

\[
 N W(0)N=0
 \tag{5}
\]

is necessary before the finite polynomial test.

## 2. A linear check for every polynomial frame

Write `f(z)=sum_{j=0}^8 a_j z^j`. Over `F_2`,

\[
 f(1/w)^2=\sum_{j=0}^8 a_jw^{-2j}.
\]

Therefore all negative coefficients in (3), together with
`sum_j a_j=0`, form a linear system with only nine unknown bits.
The checker constructs its columns from the exact coefficient matrices
`C_l` and `P_l`. It retains the full negative matrix coefficients, not
just their ranks. For every changing `N` that passes (5), this system
has rank nine. Its only solution is the zero polynomial.

The total-height fibers have dimensions

\[
 (1,2,1,2,4,2,1,2,1).
\]

A binary scalar square-zero matrix is zero. There are four square-zero
binary two-by-two matrices and 316 square-zero binary four-by-four
matrices. Taking their independent direct sums enumerates exactly

\[
 316\cdot4^4=80,896
 \tag{6}
\]

height-preserving square-zero matrices. The exact exhaustive result is:

| Class | Count | Outcome |
|---|---:|---|
| Fails the necessary leading-pole condition (5) | 74,496 | No nonzero admissible polynomial |
| Passes (5), with a nonzero correction pair `(C,P)` | 6,368 | Nine-column constraint matrix has full rank |
| Has `C=P=0` | 32 | Every polynomial frame leaves the tensor unchanged |

This proves the theorem for all binary polynomials, using the degree
reduction (4) rather than imposing a discovery cutoff as an assumption.

## 3. Exact replay

The [checker](../tools/check_polynomial_joint_height_rigidity.py) reconstructs
the seed and both endpoint normalizations. It independently matches all
infinity coefficient entries with the restored opposite-cut tensor.
It enumerates every height-preserving square-zero matrix, verifies the
nilpotence and height conditions, and computes all 6,368 linear ranks.
Direct Laurent convolution controls also verify the linear encoding for
explicit polynomials. The 32 unchanged matrices are stored individually.

```sh
python3 tools/check_polynomial_joint_height_rigidity.py \
  --output evidence/polynomial-joint-height-rigidity.json
```

The [evidence](../evidence/polynomial-joint-height-rigidity.json) includes
the counts, exact source hashes, the unchanged matrices, and a digest
of the full enumeration. This packet does not change the scalar
reflection bound or the proved factoring exponent.
