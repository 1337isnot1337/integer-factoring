# Three-row Hecke path crossings over F8

## Results and limitations

There is an exact, finite construction of strict five- and nine-letter
crossing tensors from a three-row Hecke path representation over F8.
The representation is not the two-row Temperley–Lieb quotient used by
the existing braid lane.

The local Hecke relations are verified uniformly through a finite alcove
state calculation. The projected block-swap endpoints are verified for
the two retained examples only. No uniform block-swap theorem is assumed.

Both examples fail the bare, unweighted zero/infinity Smith-profile test:
its infimum is exactly `rho=1`. This is an exact consequence of their
complete profiles, not a numerical optimizer result. Other gauges,
coefficient recombinations and other physical subspaces are not excluded.
The complete factoring exponent remains `1/2`.

## 1. A finite-state three-row representation

Work in

\[
K=\mathbb F_8=\mathbb F_2[q]/(q^3+q+1).
\]

The element q has order seven. A path is a word in row labels `0,1,2`.
Its prefix row counts satisfy

\[
\lambda_0\ge\lambda_1\ge\lambda_2,\qquad
\lambda_0-\lambda_2\le4.
\]

There are 15 possible gap states
`(lambda_0-lambda_1, lambda_1-lambda_2)`. Adding a common integer to all
three row counts does not change any coefficient below.

For adjacent steps i and i+1, let d be the difference of their contents,
`content(i+1)-content(i)`, where a box in zero-indexed row r and column c
has content `c-r`. On distinct-row steps, put

\[
a_d=\frac{q+1}{1+q^{-d}},\qquad b_d=q+a_d.
\]

The generator acts by

\[
T_i v_p=a_d v_p+b_d v_{s_i p}
\]

when the swapped path remains admissible. Same-row steps have scalar q.
If a distinct-row swap leaves the alcove, the remaining scalar is one.
Every denominator is nonzero: the relevant content differences are
nonzero modulo seven.

These are seminormal-type Hecke formulas; see the
[standard-tableau representation framework in the Sage documentation](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/hecke_algebras/ariki_koike_specht_modules.html).
The characteristic-two boundary truncation used here is checked directly,
rather than inferred from a characteristic-zero semisimplicity statement.

The [path module](../tools/three_row_hecke_paths.py) enumerates all 147
admissible three-step states over the 15 prefix gap states. It verifies
294 quadratic identities and 147 braid identities:

\[
(T_i+I)(T_i+qI)=0,\qquad
T_iT_{i+1}T_i=T_{i+1}T_iT_{i+1}.
\]

Generators at nonadjacent sites commute because their step swaps preserve
every prefix shape outside their own two-step intervals. Thus the finite
local calculation proves the relations at every path length in this alcove.

For `E_i=T_i+qI`, the Temperley–Lieb relation
`E_i E_(i+1) E_i = q E_i` fails. The retained witness starts from zero
counts with the path `(0,1,2)`; the two scalar outputs are the distinct
F8 elements encoded by 4 and 6.

## 2. Two explicit crossing tensors

For physical path length n, use the block-shuffle word

\[
(i_1,\ldots,i_{n^2})=
\bigl(n+j-1,n+j-2,\ldots,j\bigr)_{j=0}^{n-1}.
\]

Apply the factors `q I + s E_i` in that order on the two-block path
space. Select initial paths `P concatenated Q` and final paths
`R concatenated S`. At `s=1` these factors are the Hecke generators.

The [constructor](../tools/three_row_hecke_crossing.py) checks the entire
selected endpoint matrix before returning a tensor. In the two examples
it is a nonzero scalar times the swap. Dividing by that scalar gives

\[
T_{P,Q,R,S}(1)=[P=S][Q=R].
\]

The retained examples are:

- **Five letters:** all five paths ending at `(2,2,2)`, of length six.
  The two-block space ends at `(4,4,4)` and has 462 paths. The polynomial
  degree is 36 and the unnormalized braid scalar is 3 in the chosen basis.
- **Nine letters:** the union of `(0,1,2) concatenated p` and
  `p concatenated (0,1,2)` for those five paths p. The common path is counted
  once. This is a nine-dimensional subspace of the full 42-path space at
  `(3,3,3)`. The two-block alcove space at `(6,6,6)` has 60,181 paths. The
  polynomial degree is 81 and the braid scalar is 7.

All 625 and 6,561 endpoint entries are checked respectively. Both tensor
families have common scalar polynomial gcd one. The nine-letter example
is not claimed to be the full 42-letter construction.

The implementation stores polynomial coefficients in three binary planes.
It charges every global path and every selected input column. It has an
explicit guard on its three primary buffers; this is not a unit-cost
projection oracle or a practical factoring implementation.

## 3. Independent coefficient verification

The [checker](../tools/check_three_row_hecke_crossing.py) reconstructs the
packed polynomials, then independently propagates ordinary F8 Taylor jets
of each local factor about every point of F8. This second calculation
uses byte-valued field coefficients rather than packed polynomial planes.

For the five-letter tensor, five jets at each of eight distinct points
supply 25,000 scalar identities. For the nine-letter tensor, eleven jets
at each point supply 577,368 scalar identities. These are complete
coefficient checks: the total multiplicities 40 and 88 exceed the respective
degree bounds 36 and 81. A difference polynomial vanishing to those orders
at all eight points must therefore be zero.

The local rank calculations use
[restriction of scalars](binary-extension-smith-restriction.md), with
all extension costs and repeated-order checks explicit. Independent
Toeplitz and determinant controls for that utility are separate from the
path and coefficient checks.

## 4. An exact failure of the bare profile test

Write R for the cut `(P,Q)|(R,S)` and C for `(P,R)|(Q,S)`. Let nu be the
zero orders of one cut and mu the orders of the degree-D reversal of the
other. The bare partition expression is

\[
B(r)=r^{-D}\left(\sum_i r^{\nu_i}\right)
                \left(\sum_j r^{\mu_j}\right),
\qquad 0<r<1.
\]

For zero at C and infinity at R, the complete profiles give:

| Alphabet | D | Maximum zero order | Maximum reversed infinity order |
|---:|---:|---:|---:|
| 5 | 36 | 20 | 12 |
| 9 | 81 | 38 | 13 |

In each case `max(nu)+max(mu)<D`. Every term in B is therefore at least
one. There are `Q^4` terms, where Q is the physical alphabet size, so
`B(r)>=Q^4`. In the other orientation,
both order lists are identically zero and the same conclusion follows.
As r tends to one, B tends to `Q^4`; hence the bare pressure infimum is
exactly one. This argument does not cover independent letter gauges or
new coefficient combinations.

## 5. Reproduction and evidence levels

- `make verify-three-row-hecke` reconstructs and independently checks the
  five-letter case, along with the local algebra and twelve-site relations.
- `make verify-three-row-hecke-full` reconstructs the nine-letter case and
  performs its complete confluent-jet verification. This is a longer optional
  computation and is not part of the quick target.

The [five-letter receipt](../evidence/three-row-hecke-q5-replay.json) and
[nine-letter receipt](../evidence/three-row-hecke-q9-replay.json) retain all
local profiles, endpoint counts, coefficient hashes, jet hashes, relation
controls and source pins. The large coefficient table is regenerated by
the checker rather than stored as an opaque precomputed tensor.

These are ordinary finite construction and verification results. They do
not supply a better asymptotic rank budget, a new Lean theorem, or a new
factoring score.
