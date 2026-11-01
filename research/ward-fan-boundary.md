# The complete lowest Ward boundary vanishes for the 84-label fan

The entire coefficient boundary \([y^3]G(x,y)\) is zero for the
explicit \(P=L=10,k=3\), 84-label fan tensor from
[the sparse fan construction](ward-fan-polynomial-queries.md).
This extends the earlier
[single-block calculation](ward-fan-leading-block.md) to every
logical tensor entry.

Consequently this particular homogeneous tensor is divisible by
\((q+r)^4\), one power more than the three-form Ward theorem
alone supplies. The complete integer rank-sum upper bound drops
from 192,340 to **157,060**, with computed upper-bound exponent
\[
 \frac{\log_2(157060)}{2\log_2(84)}
 =1.3501329122942103\ldots.
\]
The cost remains far above the factoring seed criterion. The
result is an exact finite construction and rank-profile
improvement, not a new factoring exponent or a general theorem
for other Ward codes.

## 1. Exact scope of the calculation

The input vectors, homogeneous preimages, and coordinate outputs
are exactly those recorded in
[ward-fan-leading-block.json](../evidence/ward-fan-leading-block.json).
The new evidence pins their deterministic basis digest. The
specific choice of output left inverse is part of the finite
tensor being verified.

Use centered variables
\[
 p=1+x,\qquad q=1+y,\qquad r=1.
\]
The three-form Ward identity gives
\[
 G(x,y)=y^3\Gamma(x,y).
 \tag{1}
\]
After setting \(y=0\), each physical query uses \(q=r=1\),
\(\eta=0\), and the two-state fan recurrence over
\(\mathbf F_2[x]\). Sparse preimage and input contractions give
the exact entries of \(\Gamma(x,0)\). This polynomial has
degree at most seven.

The 84 labels have thirteen type syndromes. Pair-syndrome
conservation splits the row flattening into 31 nonempty blocks.
Their dimensions sum to \(84^2=7056\), their largest dimension
is 592, and the sum of their squared dimensions is
\[
 1,917,184.
\]
The calculation evaluates every entry in all these blocks. It
obtains
\[
 \boxed{\qquad \Gamma(x,0)=0
       \quad\text{as an entire logical tensor.}\qquad}
 \tag{2}
\]
Every entry outside the blocks is zero by the proved syndrome
selection rule. Thus (2) is exhaustive, rather than a
statement about selected entries or ranks at numeric values
of \(x\).

No separate full column contraction is needed. The row
flattening contains every possibly nonzero tensor entry, and
a zero tensor has zero rank in every flattening. The checker
also maintains the corresponding syndrome-separated column
coefficient matrices while traversing these entries; all their
computed ranks are zero.

## 2. The additional homogeneous factor

Equation (2), together with (1), gives divisibility by \(y^4\)
in the centered affine polynomial. The original tensor is
homogeneous of degree ten in \(p,q,r\). Its three-form Ward
quotient is homogeneous of degree seven.

Indeed, write that homogeneous quotient as \(H(p,q,r)\).
Equation (2) says
\[
 H(1+x,1,1)=0
\]
as a polynomial in \(x\), so \(H(p,1,1)=0\) as a polynomial
in \(p\). Homogeneity now gives \(H(p,r,r)=0\) as a
polynomial identity in \(p,r\). Polynomial division by the
monic linear form \(q+r\) therefore shows
\[
 (q+r)\mid H(p,q,r).
\]
The verified tensor is thus divisible by \((q+r)^4\).

Every coefficient identity is binary and hence remains valid
after extending to any characteristic-two field. The
calculation does not assert that the fourth power is the exact
vanishing order, or that the same extra factor holds for every
output map or every member of the fan family.

## 3. The improved complete upper profile

The previously certified full profile allowed rank 7,056 in
each cut at the five coefficients
\[
 (a,b)=(3,3),(4,3),(5,3),(6,3),(7,3).
\]
Their actual coefficient tensors are now known to be zero.
Deleting their five cost contributions decreases the integer
upper bound by
\[
 5\cdot7056=35280,
 \qquad 192340-35280=157060.
\]
Every other coefficient bound is retained unchanged. The
complete updated profile appears in the evidence, so the
cost comparison accounts for all coefficients.

The single-block bound 191,700 remains a valid intermediate
checkpoint. The exhaustive boundary calculation is strictly
stronger.

## 4. Verification and bounded work

The checker visits the 31 syndrome blocks one at a time and
clears its logical-entry cache after each row. No unblocked
\(7056\times7056\) matrix is assembled. The stored coefficient
matrices retain their syndrome blocks, each of dimension at
most 592.

The computation makes 30,604,288 physical polynomial query
requests. In 30,590,480 of them a collapsed local pair is zero.
Exact caching leaves only 55 distinct nontrivial physical
polynomials to evaluate by the two-state recurrence. The
reported run took approximately 72 seconds.

There are 62 additional direct contractions of the original
Pfaffian kernel, two in every nonempty syndrome block. These
use both full input vectors in the truncated ring
\(\mathbf F_2[x,y]/(x^8,y^4)\), independently of the Ward
factorization and the two-state recurrence. They agree with
the full computed boundary polynomials.

The underlying evaluator and sparse code were already
independently checked in the preceding frozen block package,
including 16 nonzero crossing-endpoint values and 16 zero
endpoint values. Those positive controls apply to the same
implementation and the same input/output maps. The new
package pins the source hashes of that implementation.

## 5. Reproduction

Run

~~~sh
python3 tools/check_ward_fan_boundary.py --output /tmp/ward-fan-boundary.json
~~~

The [checker](../tools/check_ward_fan_boundary.py) stores the
five coefficient matrices in both cuts, their block ranks,
the quotient-polynomial block digests, and the complete
updated upper profile. The
[evidence](../evidence/ward-fan-boundary.json) also records all
source hashes, the exact basis digest, and the previous
profile hash.
