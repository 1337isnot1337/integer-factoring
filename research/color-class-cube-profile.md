# Specialized profiles of the one-hot color-class Gaussian cube

The [scalar color-clone floor](cube-scalar-color-clone-floor.md) rules out
replicating one Gaussian kernel over independent colors and keeping the
vacuum plus single-particle letters. It leaves open genuine interactions
between colors. This note records a computational control for the natural
interacting extension: one vacuum plus \(q-1\) colors per cube vertex, and

\[
 T_q(a)=\prod_{c=1}^{q-1}G_{S_c}(x,y,z),
 \qquad S_c=\{v:a_v=c\},
 \tag{1}
\]

where \(G\) is the homogeneous eight-mode Gaussian cube tensor. This is a
homogeneous tensor with the antipodal endpoint. For \(q=2\) it is exactly
the binary Gaussian cube.

The physical alphabet of \(T_q\) is **graded**, not equal-height. A letter
is a coloring of the four vertices on one side; its energy is the number of
nonvacuum vertices, so the letter energies are \(0,1,2,3,4\) with
multiplicities \(1,4(q-1),6(q-1)^2,4(q-1)^3,(q-1)^4\). The
[equal-height 256-letter profile](cube-graded-minor-budget.md) therefore
does not apply to this tensor. Comparing the specialized low-order counts
below to that profile would be an error. The correct comparison is the
graded envelope with these letter energies; the construction is below the
graded low-order prefix at the tested sizes, so the finite data do not
exclude it.

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

with \(81\) orders and specialized order sum \(1224\). The same low-order
profile was obtained at \((y,z)=(3,5),(5,7),(7,11)\). These are still
specializations, not a generic computation.

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

The low-order layers follow the monomial counts of degree at most three
in \(q-1\) color variables. The order-eight column already departs from
that pattern at \(q=5\) (\(70\) rather than \(35\)), so this is an
observed low-order law, not a general free-algebra description.

## 2. Graded comparison

The one-hot alphabet has energies \(0,1,2,3,4\), but the tensor degree is
three times the number of nonvacuum vertices. The shared-minor envelope is
therefore evaluated with energies

\[
 0^1,\quad 3^{4(q-1)},\quad 6^{6(q-1)^2},\quad
 9^{4(q-1)^3},\quad 12^{(q-1)^4}.
 \tag{6}
\]

For \(q=3\) and \(q=4\) the envelope's prefix counts through order six are
\(65\) and \(175\). The specialized construction has prefix counts \(10\)
and \(20\). These are counts, while the envelope constrains prefix sums;
the count comparison alone is not a viability test.

The specialized \(q=4\) profile does not satisfy the graded prefix-sum
budget at the tested specializations: its first \(175\) orders sum to
\(2150\), while the envelope's first \(175\) orders sum to \(888\). This
is consistent with specialization increasing orders relative to the
generic tensor, and it is one reason the generic profile is the missing
object. The stable specializations are therefore not automatically the
generic profile, despite agreeing at six tested points. No viability or
obstruction conclusion is drawn.

## 3. Scope and next question

The finite profiles are useful as a regression and as a source of the
monomial pattern in (5). They do not prove a uniform law, and they do not
rule out the one-hot construction. The equal-height 256-letter target is
not the right comparison for this graded alphabet.

The next question is:

> What is the generic low-order profile of the one-hot color-class
> Gaussian, and does it satisfy the graded envelope at \(q=256\)?

A second question is whether a different structured subspace can reduce
the observed monomial growth without breaking the antipodal endpoint.
Capping color-class sizes is not immediate because the endpoint already
uses classes of size \(8\).

The factoring exponent remains \(1/2\); this note adds no score claim.

## 4. Exact replay

```sh
python3 tools/check_color_class_cube_profile.py \
  --output /tmp/color-class-cube-profile.json
```

The checker recomputes the \(q=2\) Gaussian regression, the \(q=3\)
specialized profile and three additional specializations, the \(q=4,5,6\)
low-order layers, and the graded envelope comparison for \(q=3,4\). The
receipt records the exact profiles and hashes. The source is
[`check_color_class_cube_profile.py`](../tools/check_color_class_cube_profile.py).
