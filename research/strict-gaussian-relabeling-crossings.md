# Strict Gaussian endpoints, charge projections, and entrywise monomial changes

## Result and scope

A constant Pascal inverse turns the one-species Gaussian seed into a
polynomial tensor with an exact crossing endpoint on its full alphabet.
This repair preserves both inherited weighted Smith profiles. The strict
endpoint permits three further constructions:

- affine sums of common boundary relabelings;
- projections onto arbitrary physical-label charge conservation;
- independent Laurent monomial factors on individual tensor entries,
  followed by one common monomial that makes the tensor polynomial.

The charge projection modulo three is already binary. It splits both
profile matrices into explicit square charge blocks, with no degree
increase. Entrywise monomial changes can couple different boundary labels;
they need not be sums of the usual four physical height functions.

These are ordinary construction theorems. The finite searches recorded
below have not improved the factoring score `alpha=1/2`. A new bound still
requires a strict pressure certificate with `rho<1/sqrt(2)`.

## 1. A constant repair gives the exact endpoint

Use the one-species polynomial Gaussian from
[the interval construction](schur-one-species-interval.md). Let the
selected indices be `S={s_0<...<s_(n-1)}` and let `M` be an admissible
ambient power of two. A physical letter is an occupation subset of `S`,
so the alphabet has size `q=2^n`.

Put

\[
 P_{ab}=[s_b\subseteq s_a],
\]

where inclusion refers to the binary one-bits of the indices. This is
unit lower triangular in the displayed index order, hence invertible over
the binary field. The one-particle matrix, with rows `(v,y)` and columns
`(u,x)`, has endpoint

\[
 A(1)=\begin{pmatrix}0&P\\I&0\end{pmatrix}.                 \tag{1}
\]

The exterior representation `Gamma(A)` is the scalar Gaussian tensor.
Multiply its `v` boundary by `Gamma(P^(-1))`. Equivalently, first replace
the one-particle matrix by

\[
 A_{\rm strict}(w)=
   \begin{pmatrix}P^{-1}&0\\0&I\end{pmatrix}A(w),
 \qquad C(w)=\Gamma(A_{\rm strict}(w)).                     \tag{2}
\]

Every coefficient is still a binary polynomial. The repair is constant,
so it introduces no denominator and does not increase the maximum entry degree.
Exterior functoriality and (1) give the exact local identity

\[
 \boxed{C_{xyuv}(1)=[x=v][y=u].}                            \tag{3}
\]

Thus the endpoint does not require a type restriction. Either flattening
used in the paired-type theorem is a permutation matrix at `w=1`. Both
flattenings consequently have full generic rank `q^2`; no numerical rank
test is needed to supply that rank to a finite-series profile engine.

### The inherited profiles are unchanged

The natural physical heights on the `v` boundary are

\[
 h_{0,v}(V)=(1-2M)|V|,
 \qquad
 h_{\infty,v}(V)=2\sum_{a\in V}a+(2-2M)|V|.                \tag{4}
\]

A nonzero entry of `Gamma(P^(-1))` from `A` to `V` has `|A|=|V|`
and `sum V>=sum A`. The inverse Pascal matrix has the same inclusion
support property, and a nonzero exterior minor admits a matching of
those inclusions. Equality of the two index sums forces the diagonal
entry.

At zero, the first height in (4) therefore commutes with the repair.
At infinity, its weighted entries acquire powers
`z^(h_(infinity,v)(V)-h_(infinity,v)(A))` with nonnegative exponents.
The weighted matrix is triangular with diagonal one, as is its inverse.
It is a regular unimodular change at that place. Whether `v` is written
as a row or column boundary only transposes the corresponding matrix
product. The two inherited weighted Smith profiles are unchanged.

This statement concerns the heights in (4). An arbitrary later physical
gauge must be evaluated on the repaired seed itself.

## 2. Affine sums and common boundary permutations

For any permutation `pi` of the physical alphabet, define

\[
 C^{\pi}_{xyuv}(w)=C_{\pi(x),\pi(y),\pi(u),\pi(v)}(w).
\]

The two equality tests in (3) are preserved by every permutation. Hence
each relabeling has the same exact endpoint, even if it changes particle
number or the index height. If polynomial coefficients `a_j(w)` satisfy
`sum_j a_j(1)=1`, then

\[
 T(w)=\sum_j a_j(w) C^{\pi_j}(w^d),\qquad d\ge1,           \tag{5}
\]

is another strict polynomial crossing. Its degree is at most
`d*max_degree(C)+max_j deg(a_j)`. The same statement applies after
clearing a common Laurent monomial in the `a_j`.

Two useful special cases are an odd sum over the binary field and

\[
 C(w^d)+\gamma w^k(1+w)
       \bigl(C^{\pi}(w^d)+C(w^d)\bigr).                    \tag{6}
\]

The second summand in (6) vanishes at the endpoint. A fixed extension
coefficient `gamma` is permitted by the
[field-descent theorem](fixed-field-descent.md). The finite engine retains
the full substitution degree and the additional factor in (6).

These sums can be non-Gaussian. For four modes, exchange physical labels
`6={1,2}` and `9={0,3}`, leaving all other labels fixed. Vacuum and
one-particle labels are unchanged. The mixture
`(1+gamma)C+gamma*C^pi` has the same vacuum and quadratic coefficients
as `C`, but a different quartic coefficient at
`x={1,2}, y={2}, u=empty, v={1}`. For nonzero `gamma`, its quartic Wick
residual is `gamma*(w^10+w^11)`. Thus this is an actual change of the
Gaussian tensor, rather than a physical diagonal gauge.

## 3. Arbitrary label-charge projections

Choose any function `f:A->Z/3`. It may depend nonlinearly on the full
physical letter. Let `omega^2+omega+1=0` in the four-element field and
multiply the `x,u` legs by `omega^(j*f)` and the `y,v` legs by their
inverses. The resulting diagonal conjugate has entry

\[
 C^{(j)}_{xyuv}(w)
   =\omega^{j(f(x)+f(u)-f(y)-f(v))}C_{xyuv}(w).
\]

Sum these tensors for `j=0,1,2`. The geometric sum is one when its
exponent is zero modulo three and is zero otherwise. In characteristic
two, three copies of one sum to one. Therefore

\[
 \boxed{C_f(w)_{xyuv}
   =[f(x)+f(u)=f(y)+f(v)\pmod3]\cdot C(w)_{xyuv}.}         \tag{7}
\]

Every nonzero endpoint entry in (3) conserves this charge. Equation (7)
is consequently a strict crossing. Its coefficients are already binary:
the intermediate four-element-field calculation is an exact character
identity, with no averaging denominator and no final extension-field
representation. Its degree is at most the original degree.

More generally, imposing any equality of physical-label charges by
deleting the other entries preserves (3). The modulo-three version has
the convenient three-term affine realization above.

### Exact charge blocks at the two places

At zero, flatten the row as `(x,y)` and the column as `(u,v)`. Define

\[
 r_0(x,y)=f(x)-f(y),\qquad c_0(u,v)=f(v)-f(u).
\]

At infinity, flatten as `(x,u)` versus `(y,v)` and define

\[
 r_\infty(x,u)=f(x)+f(u),\qquad
 c_\infty(y,v)=f(y)+f(v).
\]

In each case (7) retains exactly the entries with equal row and column
charge. Each block is square: its row and column counts are the same
ordered-pair count on the alphabet. The exact endpoint restricts to a
permutation matrix in every nonempty block, so its generic rank is its
full block size. Any diagonal physical height gauge respects these
blocks. The complete weighted profile is the multiset union of the
block profiles. This accounting preserves every physical state; no
entropy or type factor is removed from the pressure denominator.

## 4. Independent entrywise monomial exponents

Let `E(x,y,u,v)` be any integer function on the support of `C`. Form the
Laurent tensor

\[
 L_{xyuv}(w)=w^{E(x,y,u,v)} C_{xyuv}(w).                    \tag{8}
\]

Every monomial equals one at `w=1`, so (8) still satisfies (3). Choose

\[
 h=\max\left(0,-\min_{C_{xyuv}\ne0}
       \{E(x,y,u,v)+v_0(C_{xyuv})\}\right).
\]

Then `T=w^h L` is polynomial and still has endpoint (3). The resulting
degree, including this clearing factor, is retained by the implementation.

The common clearing factor has no effect on the product of profile
partition functions. It adds `h` to every zero-place Smith order and
subtracts `h` from every infinity order. Thus

\[
 W_0(T;s)=s^{-h}W_0(L;s),\qquad
 W_\infty(T;s)=s^hW_\infty(L;s).                           \tag{9}
\]

The two factors cancel exactly. This does not remove the individual
entry shifts, which remain part of both profile calculations.

The familiar four physical gauges correspond to separable choices of
`E`. Equation (8) also permits interactions between different boundary
letters. Examples include a charge-mismatch indicator and

\[
 E(x,y,u,v)=k\,(f(x)-f(v))(g(u)-g(y)).                     \tag{10}
\]

Arbitrary exponents can destroy cancellations that make the Gaussian
profile useful. The search tool
[tie_preserving_crossing_exponents.py](../tools/tie_preserving_crossing_exponents.py)
therefore also constructs linear directions that preserve every tied
leading two-by-two minor at both places. This is a search constraint,
not a theorem about higher minors or a profile bound.

For the two-mode strict seed, the exact constraint matrix has 62 supported
entries, 458 displayed equations, and rank 49. Its nullspace has dimension
13, while physical-leg gauges have rank 12. There is one additional
integer direction modulo those gauges. Its first tested nonzero shifts
worsen the frozen profile; retaining two-by-two ties alone does not
preserve all the useful higher cancellations.

## 5. Finite search status

The [saved screens](../evidence/strict-gaussian-search-screens.json) retain
each tested charge, relabeling, and monomial interaction. The inherited
four-mode Gaussian control is approximately
`rho=0.73108093597359`. The first screens used the exact polynomial
tensors and exact finite-field arithmetic, followed by floating pressure
evaluation. They are not global optimization or exclusion theorems.

The hard charge screen covered 64 functions: 40 linear charges modulo
nonzero scalar equivalence, six quadratic bit charges, parity, the
binomial coefficient of particle count, and 16 single-label indicators.
The total-particle charge leaves the seed unchanged and reproduces the
control. None of the displayed nontrivial projections improves it at the
inherited physical gauges. The best displayed nontrivial result is
approximately `0.8088277481`, for the indicator of label 12.

A separate 25-candidate screen covered equal-height transpositions,
constant and vanishing affine mixtures, degree-four substitutions, single
occupation-bit flips, and odd sums of three bit relabelings. It gave no
improvement at the inherited gauges. Ten transpositions were also tested
with the midpoint of their two relabeled gauges, again without a gain.

The two-mode monomial screen tested 240 displayed quadratic interactions
at `log(s)=1`. It found no improvement over its Gaussian control at that
point. These frozen tests motivate retaining tied cancellations and
optimizing the complete physical-letter gauges. They do not exclude the
construction families in Sections 2--4.

## 6. Implementation and numerical obligations

[relabelled_gaussian_sums.py](../tools/relabelled_gaussian_sums.py) supplies
the exact constructors, both source flattenings, complete physical-letter
gauges, and a profile oracle with a supplied endpoint rank. If a series
window does not resolve all promised pivots, it is enlarged; failure at
the maximum supported precision remains an error.

[search_strict_gaussian_gauges.py](../tools/search_strict_gaussian_gauges.py)
optimizes all nonconstant physical-label weights at both places. Its
gradient was checked against independent finite differences away from
valuation walls. Numerical improvements require a saved exact seed,
rationalized physical weights, complete profile replay, and a strict
rational pressure comparison before changing any factoring claim.

The constructor uses binary or four-element-field polynomial arithmetic.
The optional profile and optimization tools require a C++17 compiler,
NumPy, and SciPy. No new Lean formalization is used.

```sh
python3 tools/check_strict_gaussian_crossings.py \
  --output evidence/strict-gaussian-relabeling-crossings.json
```

The [exact control record](../evidence/strict-gaussian-relabeling-crossings.json)
contains 69,904 base endpoint entries, six complete inherited profiles,
four literal phase-average and charge-block controls, three Laurent
clearing checks, and the exact two-mode tie-constraint calculation.
