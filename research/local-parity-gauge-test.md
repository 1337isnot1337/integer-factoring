# Exact recognition of constant local even-parity gauges

## Status

This note gives a constructive recognition theorem for polynomial tensors over
fields of characteristic two. It decides whether constant, invertible changes
of basis on individual binary physical legs can make every odd-parity
coefficient vanish. The changes of basis may use a coefficient-field extension;
they do not depend on the polynomial parameters.

The general test is linear in `2n+1` unknowns for `n` binary legs. A product of
pairwise equality endpoints reduces it to `n` unknowns and permits an explicitly
endpoint-preserving output. Thus the eight-leg cube needs 17 unknowns in general,
or eight unknowns with its antipodal equality endpoint.

For an input over a finite field `F_q`, every positive instance has gauges over
`F_(q^2)`. The implementation handles `F_4` inputs and constructs gauges over
`F_16`. It also returns small exact inconsistency certificates for negative
instances. These are recognition and search-filter results. They do not improve
the current factoring exponent.

## 1. General recognition theorem

Let `F` have characteristic two, and let

\[
 0\ne T\in (F^2)^{\otimes n}\otimes_F F[z_1,\ldots,z_d].
\]

Write `t_m` for the polynomial coefficient of physical word
`m in {0,1}^n`, and put `D=diag(0,1)`. The following are equivalent.

1. Over some coefficient-field extension `L/F`, there are constant
   `S_i in GL_2(L)` such that `(tensor_i S_i)T` has only even physical support.
2. There are matrices `A_i in Mat_2(F)`, all of trace one, satisfying

   \[
    \left(\sum_{i=1}^n A_i^{(i)}\right)T=0.                 \tag{1}
   \]

Here `A_i^(i)` acts on physical leg `i` and as the identity on every other leg.
The equality is a polynomial identity, checked on every coefficient.

### Necessity and descent

If the transformed tensor is even, its parity operator
`sum_i D^(i)` annihilates it. Conjugating gives (1) with
`A_i=S_i^(-1) D S_i`, initially over `L`.

The entries of (1), together with `tr(A_i)=1`, form a linear system over `F`.
A linear system over `F` is solvable over an extension if and only if it is
solvable over `F`. Thus the matrices can be chosen over the original field.

### Sufficiency

Set `N=sum_i A_i^(i)`. Operators on distinct legs commute. Cayley–Hamilton
and characteristic two therefore give

\[
 A_i^2+A_i=\det(A_i)I,
 \qquad
 N^2+N=\left(\sum_i\det(A_i)\right)I.
\]

Since `NT=0` and `T` is nonzero,

\[
 \sum_i\det(A_i)=0.                                      \tag{2}
\]

In a finite coefficient-field extension, choose roots

\[
 \lambda_i^2+\lambda_i=\det(A_i).
\]

Their sum satisfies `s^2+s=0`, so it is zero or one. Replacing one root by its
other root, when necessary, makes `sum_i lambda_i=0`. Then

\[
 B_i=A_i+\lambda_i I
\]

is an idempotent of trace one, hence a rank-one idempotent. It has a basis in
which it is `D`. Choose `S_i` with `S_i B_i S_i^(-1)=D`. Because the scalar
shifts sum to zero, `sum_i B_i^(i)` still annihilates `T`. Conjugating proves
that the transformed tensor has only even support.

For a finite field `F_q`, all of the displayed quadratic polynomials split in
`F_(q^2)`. For an arbitrary characteristic-two field, a common finite splitting
extension suffices; a quadratic extension is not asserted in that generality.

## 2. The explicit linear system

Write

\[
 A_i=\begin{pmatrix}a_i&b_i\\c_i&1+a_i\end{pmatrix},
 \qquad A=\sum_i a_i.
\]

Only `A`, rather than the individual `a_i`, enters (1). For every polynomial
monomial `e` and physical word `m`, the equation is

\[
 A\,t_{m,e}
 +\sum_{i:m_i=0}b_i\,t_{m\mathbin\oplus 2^i,e}
 +\sum_{i:m_i=1}c_i\,t_{m\mathbin\oplus 2^i,e}
 = (|m|\bmod 2)\,t_{m,e}.                                \tag{3}
\]

There are `2n+1` unknowns. Any solution is realized by taking `a_1=A` and the
remaining diagonal parameters zero. Gaussian elimination either finds such a
solution or expresses `0=1` as a linear combination of at most `2n+2` original
equations. That certificate remains contradictory over every field extension.

If the tensor has `M` polynomial monomials, directly forming and solving this
system uses `O(M 2^n n^2)` field operations, with only `O(n^2)` elimination
storage apart from the input and the small certificates. No enumeration of
local invertible matrices is needed.

## 3. Equality endpoints reduce the test and preserve normalization

Suppose `n` is even, with a fixed pairing of the physical legs, and

\[
 T(1,\ldots,1)=\bigotimes_{(i,j)}(|00\rangle+|11\rangle).
                                                               \tag{4}
\]

Evaluating (3) at the endpoint forces

\[
 A=0,\qquad b_i=c_j,\quad c_i=b_j
 \quad\text{for every pair }(i,j).                       \tag{5}
\]

For `A=0`, use an aligned physical word. For either off-diagonal equality in
(5), use a word with precisely that pair misaligned. These endpoint constraints
leave two variables per pair, hence `n` variables in total.

For one pair, write these variables as `b,c` and choose

\[
 A_i=\begin{pmatrix}0&b\\c&1\end{pmatrix},
 \qquad A_j=A_i^{\mathsf T}.
\]

This does not change the operator in (1), since its diagonal parameters enter
only through their sum. Both matrices have determinant `bc`. Use the same
quadratic root `lambda` for both partners. Their idempotents are transposes.
If `S_i` diagonalizes the first, choose

\[
 S_j=S_i^{-\mathsf T}.                                  \tag{6}
\]

This diagonalizes the second and preserves their equality tensor exactly.
Repeating for every pair gives an even tensor with the original literal
endpoint (4). Thus, for tensors satisfying (4), arbitrary constant local
parity gauges exist over an extension if and only if paired,
endpoint-preserving ones do.

Equivalently, one could diagonalize a general solution first and repair the
endpoint afterward. The transformed endpoint is a product of invertible
pair matrices. Global even support forces each factor to have purely even or
purely odd support. Diagonal rescaling and, for each odd factor, one physical
bit flip restore equality. There are an even number of odd factors, so the
repair preserves global even parity. Construction (6) performs the repair
directly.

## 4. Consequence for flat cube profile barriers

Constant local changes of physical basis multiply every coordinate flattening
on the left and right by invertible constant matrices. They preserve its
complete generic Smith orders. A coefficient-field extension also preserves
these orders, since the ranks and nonvanishing of minors do not change.
Homogeneous polynomial degree is preserved.

Consequently, the all-field lower bound in
[`cube-flat-parity-floor.md`](cube-flat-parity-floor.md) applies to every
homogeneous degree-six flat cube tensor accepted by the endpoint version of
this test. That family includes tensors whose displayed coefficients have odd
physical support. Such support alone does not establish escape from the parity
barrier. A negative linear certificate does establish escape from this entire
constant-local-gauge family. It does not establish a favorable profile or
exclude parameter-dependent rational gauges.

## 5. A literal reference outside every such gauge

The canonical tensor is retained in
[`local-parity-nongauge-reference.json`](../evidence/local-parity-nongauge-reference.json).
It has degree six over `F_4`, the antipodal equality endpoint, and complete
profiles

\[
\begin{aligned}
 &(0,2,2,3,4,4,5,5,6,6,7,8,9,11,11,13),\\
 &(0,0,2,2,3,3,4,5,5,6,6,8,11,13,13,15),\\
 &(0,0,2,2,3,3,4,5,5,6,6,8,11,13,13,15).
\end{aligned}
\]

Its numerically optimized pressure is approximately `0.77520767586107`.
This is a reference for recentering, not an improvement over the current
`0.67817178087567` seed.

In the endpoint system, order the pairs as `(0,7),(1,6),(2,5),(3,4)` and the
variables as `(b_0,c_0,b_1,c_1,b_2,c_2,b_3,c_3)`. With `alpha^2+alpha+1=0`,
two original coefficient equations already contradict one another:

\[
 \begin{array}{c|c|c}
 \text{monomial}&\text{physical word as an integer}&\text{equation}\\\hline
 z^6&10&\alpha(c_1+c_2)=0\\
 x^2z^4&8&\alpha(c_1+c_2)=1.
 \end{array}
\]

Their sum is `0=1`. The general 17-variable test separately returns a
five-equation inconsistency certificate. Both are replayed against the literal
tensor, so its failure is independent of any sampled field values.

## 6. Implementation and replay

- [`local_parity_gauge.py`](../tools/local_parity_gauge.py) constructs the linear
  systems, returns exact contradiction certificates, and constructs explicit
  constant `F_16` gauges for positive `F_4` inputs.
- [`check_local_parity_gauge.py`](../tools/check_local_parity_gauge.py) applies
  the returned maps directly to all polynomial coefficients. It also checks
  paired endpoint preservation and verifies the negative certificates.
- [`local-parity-gauge.json`](../evidence/local-parity-gauge.json) records the
  deterministic controls and source hashes.

The controls include arbitrary local `F_4` images of an even seed, paired
images preserving equality, an odd-support tensor that becomes even, and an
instance whose unique paired linear solution requires a quadratic extension.
For that last instance the first trace-one factor has determinant `alpha`;
`u^2+u=alpha` has no root in `F_4`, while the constructed `F_16` gauges work.
The retained negative reference also receives three complete generic Smith
certificates, independently verified by polynomial multiplication and unit
constant-term determinants. Deliberately damaged gauges, contradiction
certificates, endpoints, and the excluded zero input are rejected.

```sh
python tools/check_local_parity_gauge.py
```
