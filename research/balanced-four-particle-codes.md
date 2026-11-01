# Five-letter codes from all seventy four-particle states

## Result and scope

The five-state reflection code uses only part of the sixteen
reflection-admissible four-particle states on eight modes. This
checkpoint permits every one of the seventy four-particle states as
a physical letter on all four boundaries.

The completed natural-gauge screen contains 2,089 distinct five-state
alphabets. None records an improvement over the retained code's cost
`0.7263479477044683`. The best tested alphabet containing a nonreflection
state has cost approximately `0.7482523811337818`. These are finite
numerical scale optimizations of exact tensors and complete profiles;
they do not exclude other alphabets or physical gauges.

The construction uses the same strict Gaussian matrix as the
[coefficient-aligned complement construction](aligned-gaussian-complement-perturbations.md).
Here the physical alphabet changes on every boundary, rather than
adding complementary columns on one boundary.

## 1. Arbitrary balanced physical alphabets

Let \(u=1+w\) and define the sixteen-dimensional strict kernel

\[
D(w)=\begin{pmatrix}uG(w)&w^{16}I_8\\I_8&uH(w)\end{pmatrix},
\qquad
G_{ij}(w)=w^{2(i\mathbin{\&}j)},\quad
H_{ij}(w)=w^{2(7-(i\mathbin{|}j))}.
\tag{1}
\]

Choose any five distinct four-element subsets of \(\{0,\ldots,7\}\).
Their masks form a code \(C=(C_0,\ldots,C_4)\). Define

\[
F_{xyuv}(w)=
\det D(w)_{\,C_v\sqcup(8+C_y),\,C_u\sqcup(8+C_x)}.
\tag{2}
\]

At \(w=1\), (1) is the block-swap permutation matrix. Distinctness
of the chosen masks therefore gives

\[
F_{xyuv}(1)=[x=v][y=u].
\tag{3}
\]

No reflection-pair condition is used in this endpoint argument.

Every determinant term in (2) uses equally many upper-left and
lower-right entries, because each selected block row and block column
contains four indices. If that number is \(k\), their factor in
\(u\) is \((1+w)^{2k}=(1+w^2)^k\). All remaining monomial
factors in (1) have even exponents. Thus every coefficient of (2)
is a polynomial in \(z=w^2\).

Write the compressed tensor as \(\widetilde F(z)\), and let
\(g(z)\) be the gcd of all 625 entries. Equation (3) implies
\(g(1)=1\), since some endpoint entries equal one. Consequently

\[
T_C(z)=\widetilde F(z)/g(z)
\tag{4}
\]

is a polynomial strict crossing. The implementation divides every
entry exactly and checks the endpoint again. It removes the complete
common polynomial factor, including any factors beyond a monomial.

The source determinant routine proves its degree bound before using
a finite series representation. In these cases the degree bound in
\(w\) is 128 and the exact calculation retains 256 coefficients.
No unresolved series tail is interpreted as a zero determinant.

## 2. Fixed grade energies

For a physical mask \(a\), define

\[
g(a)=\sum_{i\in a}(2i+1),\qquad
E_j=g(C_j)-\min_k g(C_k).
\tag{5}
\]

These even, nonnegative energies specify the tested gauge. At infinity,
the active heights are \(E/2\) on \(U,V\); at zero they are
\(E/2\) on \(X,Y\). All other heights are zero. Both complete
25-by-25 weighted profiles are recomputed from (4).

For the original code, in increasing-grade order,

\[
C=(15,23,43,51,77),\qquad E=(0,2,6,8,10),
\]

formula (4) agrees coefficient by coefficient with the retained
degree-18 source after the same permutation on all four labels.
The natural energy choice is one explicit gauge, not an assumption
that it minimizes the cost of an arbitrary new alphabet.

## 3. Completed finite schedule

The [search record](../evidence/balanced-four-particle-code-search.json)
contains every requested alphabet, its complete profiles, natural
energies, scalar normalization, tensor digest, and numerical cost.

| Stage | New distinct alphabets |
| --- | ---: |
| Retained reference | 1 |
| All one-symbol replacements by an omitted four-particle state | 325 |
| Every five-subset of the twelve lowest-grade states, after overlaps | 784 |
| Seeded random alphabets and short grade windows, after overlaps | 979 |
| Total | 2,089 |

The complete twelve-state shell has \(\binom{12}{5}=792\) alphabets;
eight already occur in earlier stages. The random stage requests one
thousand choices with seed `2026103208`, then removes duplicates.
Half of those choices use a short window in grade order.

The screen found 1,543 distinct energy/profile combinations and checked
1,305,625 endpoint entries. Its best nonreflection code is
\((15,23,43,51,228)\), with energies \((0,2,6,8,28)\), degree 44,
and cost approximately `0.7482523811337818`. Replacing mask 228 by
216 gives the same normalized literal tensor in this case.

## 4. Independent verification

The [control record](../evidence/balanced-four-particle-code-controls.json)
checks four alphabets, including repeated minimum grades and
nonreflection states. It expands all 2,500 polynomial determinants
directly in Python from (1), without the native selected-minor routine.
It certifies eight complete profiles with exact regular polynomial
frames and verifies all 625 permuted baseline coefficients.

The [saved-search replay](../evidence/balanced-four-particle-code-search-replay.json)
independently checks all 24 retained tensors through 15,000 complete
polynomial determinants and certifies 48 complete Smith profiles.
These profile certificates do not use the native Smith oracle from
the search. Numerical cost minimization remains separate from the
exact polynomial and profile checks.

## 5. Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_balanced_four_particle_codes.py \
  --output /tmp/balanced-four-particle-code-controls.json
OPENBLAS_NUM_THREADS=1 python3 tools/search_balanced_four_particle_codes.py \
  --low-states 12 --random-codes 1000 --seconds 1800 \
  --output /tmp/balanced-four-particle-code-search.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_balanced_four_particle_codes.py \
  --search evidence/balanced-four-particle-code-search.json \
  --output /tmp/balanced-four-particle-code-search-replay.json
```

The [constructor](../tools/balanced_four_particle_codes.py) is also
usable independently of the finite schedule. Its strict endpoint
argument applies over every field of characteristic two.
