# The one-hot color-class Gaussian misses the 256-letter cube profile

The [scalar color-clone floor](cube-scalar-color-clone-floor.md) rules out
replicating one Gaussian kernel over independent colors and keeping the
vacuum plus single-particle letters. It explicitly leaves open genuine
interactions between colors. This note tests the natural interacting
extension: one vacuum plus \(q-1\) colors per cube vertex, and

\[
 T_q(a)=\prod_{c=1}^{q-1}G_{S_c}(x,y,z),
 \qquad S_c=\{v:a_v=c\},
 \tag{1}
\]

where \(G\) is the homogeneous eight-mode Gaussian cube tensor. This is a
homogeneous tensor with the antipodal endpoint. For \(q=2\) it is exactly
the binary Gaussian cube.

The finite profiles below show that the low-order channels grow much faster
than the 256-letter profile budget allows. The number of orders at most six
is \(\binom{q+2}{3}\) at the tested specializations, while the 256-letter
hull has only \(65{,}536\). At \(q=256\) the observed law gives
\(2{,}829{,}056\). The extrapolation is not proved; the note records a
scoped negative control and the structural question it raises.

## 1. Specialized block Toeplitz ranks

Order one cut axis and write the cut matrix as

\[
 M_q(x)=\sum_{k\ge0}M_{q,k}x^k.
\]

For \(L\ge1\), let \(R_q(L)\) be the rank over \(\mathbb F_{256}\) of the
block lower-triangular Toeplitz matrix with blocks \(M_{q,a-b}\) for
\(0\le b\le a<L\). If the specialized Smith orders are \(s_1,\dots,s_n\),
then

\[
 R_q(L)=\sum_{i=1}^{n}\max(0,L-s_i).
\]

Consequently \(R_q(L)-R_q(L-1)=\#\{i:s_i\le L-1\}\), and the exact
multiplicity at order \(k\) is the second difference

\[
 \#\{s_i=k\}
 =
 \bigl(R_q(k+1)-R_q(k)\bigr)
 -
 \bigl(R_q(k)-R_q(k-1)\bigr).
 \tag{2}
\]

The evaluation uses \(y=2,z=3\) in \(\mathbb F_{256}\) with modulus
\(x^8+x^4+x^3+x+1\). A specialization can only increase Smith orders, so
the resulting small-order counts are lower bounds for the generic tensor.
The \(q=2\) check against the known Gaussian profile validates the
extraction:

\[
 0^1,2^1,4^1,6^1,8^1,10^2,12^2,14^2,16^1,18^1,20^1,22^1,24^1.
 \tag{3}
\]

For \(q=3\), the specialized profile is

\[
 0^1,2^2,4^3,6^4,8^5,10^8,12^9,14^{12},16^{10},18^4,
 20^7,22^4,24^6,28^3,32^2,36^1,
 \tag{4}
\]

with \(81\) orders and order sum \(1224\). The same low-order profile was
obtained at \((y,z)=(3,5),(5,7),(7,11)\). These are still specializations,
not a generic computation.

For \(q=4,5,6\) the first four nonzero layers are

\[
 0^1,\quad 2^{q-1},\quad 4^{\binom q2},\quad 6^{\binom{q+1}3}.
 \tag{5}
\]

The tested values are

| \(q\) | order 2 | order 4 | order 6 |
|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 |
| 3 | 2 | 3 | 4 |
| 4 | 3 | 6 | 10 |
| 5 | 4 | 10 | 20 |
| 6 | 5 | 15 | 35 |

The order-six column is \(1,4,10,20,35=\binom{q+1}{3}\). The cumulative
count through order six is even simpler:

\[
 \#\{s_i\le6\}=4,10,20,35,56=\binom{q+2}{3}.
 \tag{6}
\]

This is the quantity used below. Specialization can only increase orders,
so the specialized count \(\#\{s_i\le6\}\) is a lower bound for the
generic tensor. The order-eight column already departs from the
symmetric-power pattern at \(q=5\) (\(70\) rather than \(35\)), so (5) is
an observed low-order law, not a general free-algebra description.

## 2. Comparison with the 256-letter hull

The [graded minor budget](cube-graded-minor-budget.md) gives the necessary
profile

\[
 0^1,2^{255},4^{255},6^{65025},13^{262650},14^{62985},
 18^{66325500},24^{4228250625},
 \tag{7}
\]

with pressure \(72/125<1/\sqrt3\). The order-six layer is the first large
layer, and the prefix through order six has

\[
 1+255+255+65{,}025=65{,}536=256^2
 \tag{8}
\]

orders.

If the law in (6) persists to \(q=256\), the color-class tensor has
\(\binom{258}{3}=2{,}829{,}056\) orders at most six. The hull also gives
the total order budget

\[
 \sum_i\sigma_i=102{,}676{,}561{,}920
 \tag{9}
\]

per cut. With \(2{,}829{,}056\) orders at most six and total order at
most (9), the smallest possible partition value puts those orders at six
and spreads the remaining total order as evenly as possible. The
resulting extremal profile has

\[
 2{,}829{,}056\text{ orders at }6,\quad
 3{,}940{,}408{,}064\text{ at }24,\quad
 351{,}730{,}176\text{ at }23.
 \tag{10}
\]

At \(t=6/13\) its partition is about \(27{,}386.7\), while the hull
partition is about \(804.8\); the ratio is \(34.0\). An exact integer
comparison certifies that the pressure at this scale is at least
\(72/125\). A fine-grid scan of the same lower bound has minimum about
\(0.724\) near \(t=0.583\), above \(1/\sqrt3\). The scan is a
diagnostic, not a certified minimum; the exact certificate is the single
scale \(t=6/13\). The conclusion is conditional on the prefix lower bound
at \(q=256\).

The structural point does not depend on the exact numerical minimum: the
hull's order-six budget is quadratic in \(q\), while the interacting
color-class construction produces a cubic number of low-order channels.
That is the opposite direction from the near-saturation needed for the
[three-axis crossing](cube-crossing-transfer.md).

## 3. Scope and next question

The finite checks do not prove that (6) holds for every \(q\), and they do
not rule out a different one-hot subspace with fewer order-six channels.
They do show that the most direct genuine-interaction extension of the
Gaussian cube is far below the hull at the first layer where the hull
becomes large.

The next question is structural:

> Is there a homogeneous one-hot cube tensor with fewer than
> \(65{,}536\) orders at most six at \(q=256\)?

If the observed \(\binom{q+1}{3}\) lower bound is general for one-hot
tensors, the one-hot route is excluded and the 256-letter hull requires
multi-particle letters or a different grading. If a one-hot construction
can suppress these channels, it must break the symmetric-color mechanism
that produces (5). Capping color-class sizes is not immediate: the
endpoint already uses color classes of size \(8\), so a naive cap destroys
endpoint words.

The factoring exponent remains \(1/2\); this note adds no score claim.

## 4. Exact replay

```sh
python3 tools/check_color_class_cube_profile.py \
  --output /tmp/color-class-cube-profile.json
```

The checker recomputes the \(q=2\) Gaussian regression, the \(q=3\)
specialized profile and three additional specializations, the \(q=4,5,6\)
low-order layers, and the \(q=256\) order-six comparison. The receipt
records the exact partitions and hashes. The source is
[`check_color_class_cube_profile.py`](../tools/check_color_class_cube_profile.py).
