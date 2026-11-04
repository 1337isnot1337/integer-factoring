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

The physical alphabet of \(T_q\) is **graded**, not equal-height. A
vertex letter is either the vacuum (energy \(0\)) or one of \(q-1\)
colors (energy \(3\) in tensor-degree units). A side word is a coloring
of the four vertices on one side; its energy is three times the number of
nonvacuum vertices, so the side-word multiplicities are
\(1,4(q-1),6(q-1)^2,4(q-1)^3,(q-1)^4\) at orders
\(0,3,6,9,12\). The
[256-letter abstract target](cube-graded-minor-budget.md) uses exactly
this grading when \(q=256\): one vacuum and 255 height-one letters,
with tensor energy three per occupied vertex. The correct comparison
for each \(q\) is the shared-minor envelope with the **vertex-letter** energies \(0\) and \(3\), whose length-four word
multiplicities are exactly the side-word distribution above. Section 2
performs that comparison; the specialized profile satisfies it.

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
the resulting cumulative counts through any order are lower bounds for
the generic tensor. An individual layer multiplicity has no such
monotonicity guarantee.
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
| 7 | 6 | 21 | 56 |
| 8 | 7 | 28 | 84 |

The low-order layers follow the monomial counts of degree at most three
in \(q-1\) color variables; the cumulative count through order six is
\(\binom{q+2}{3}\). The order-eight column already departs from
that pattern at \(q=5\) (\(70\) rather than \(35\)), so this is an
observed low-order law, not a general free-algebra description.

The checker also verifies the endpoint at \((1,1,1)\): every
antipodal-equal coloring has value \(1\), and every other coloring has
value \(0\), for \(q=3,4\). Endpoint-preserving random perturbations of
the non-endpoint entries at \(q=3\) drive the specialized profile to all
orders zero. That is worse for the partition function (it approaches the
maximum \(n\)), not better. The observed low-order profile is therefore a
special favorable feature of the product tensor, not a generic lower
bound; the sector obstruction is independent of this finite control.

## 2. Graded comparison

The shared-minor theorem of
[cube-graded-minor-budget.md](cube-graded-minor-budget.md) is stated for an
alphabet of physical letters. Here the physical letters are the \(q\)
vertex colors, not the \(q^4\) side words: the vacuum has energy \(0\)
and each of the \(q-1\) colors has energy \(3\). The theorem's
\(S_1,S_2,S_4\) are computed from that letter distribution; its
length-four words are precisely the side words, with multiplicities
\(1,4(q-1),6(q-1)^2,4(q-1)^3,(q-1)^4\).

For \(q=3\) the envelope blocks are

\[
 0^1,2^2,4^2,6^4,\left(\tfrac{25}{2}\right)^{24},18^{32},24^{16},
 \tag{6}
\]

with total multiplicity \(81\) and total order \(1296\). The
specialized profile (4) has order sum \(1224\), and every prefix sum is
at most the envelope's; the maximum excess is \(0\), attained at prefix
\(5\).

For \(q=4\) the envelope blocks are

\[
 0^1,2^3,4^3,6^9,\left(\tfrac{216}{17}\right)^{51},18^{108},24^{81},
 \tag{7}
\]

with total multiplicity \(256\) and total order \(4608\). The
specialized profile of section 1 has \(250\) recorded orders, total
order \(4068\), and first-\(175\) order sum \(2150\), against an
envelope first-\(175\) sum of \(2664\). Every prefix sum is at most the
envelope's; the maximum excess is \(0\), attained at prefix \(7\).

An earlier version of this note reported the envelope's first-\(175\)
sum as \(888\). That was wrong: it passed the side-word multiplicities
to the envelope as if they were letter multiplicities, which overcounts
the one-letter and pair bounds. With the correct vertex-letter
distribution the specialized profile satisfies the necessary envelope at
all tested prefixes. Specialization increases orders, so the generic
profile, if it differs, has smaller orders and also satisfies the
envelope. This is a necessary-condition pass, not a viability proof; the
generic \(q=256\) profile has not been computed. The
[distinct-sector proof](color-distinct-sector-obstruction.md) now excludes
its common-scale unweighted three-axis target without that computation.

## 3. Current scope and stopped family

This specialized contribution is retained from worker 4, exact source
`1772970a7dbaab83eeb0607014a1111e46718d86`, including the correction
history above. The \(q=7,8\) rows and endpoint check are from commits
`657a21e` and `d4fd713`, integrated here at `47aac516`. The finite profiles are regression controls and do not
prove a uniform monomial law. The historical `q4_full_profile` receipt
key records 250 of 256 orders through order 39, rather than an exhaustive
profile.

The subsequent [generic direct-sector obstruction](color-distinct-sector-obstruction.md)
proves that this one-hot construction has common-scale, unweighted
three-axis profile pressure at least \(3/4>1/\sqrt3\) for **every
\(q\ge34\)** and every common physical type. Large instances of this
family are stopped under that criterion. No larger specialized matrix
or complete generic profile is needed to reach this conclusion. The same
proof excludes every tensor with the same homogeneity and antipodal
endpoint that preserves separate color parity, over any field; Gaussian
factorization and color-permutation symmetry are unnecessary.

The earlier statement that the 256-letter target had incompatible grading
was incorrect: its 255 nonvacuum letters all have height one, exactly as
here. The side-word versus vertex-letter correction remains necessary.
An earlier revision also recorded a conditional \(q=256\) pressure
diagnostic from the observed low-order law; the sector obstruction
subsumes it and does not require proving that law.
Smaller alphabets, physical types that differ between vertices, weighted gauges, actual
coefficient ranks and different tensors remain outside the new obstruction.

The factoring exponent remains \(1/2\); this note adds no score claim.

## 4. Exact replay

```sh
python3 tools/check_color_class_cube_profile.py \
  --output /tmp/color-class-cube-profile.json
```

The checker recomputes the \(q=2\) Gaussian regression, the \(q=3\)
specialized profile and three additional specializations, the \(q=4,5,6,7,8\)
low-order layers, the \(q=4\) specialized orders through 39, the graded
envelope comparison for \(q=3,4\), and the endpoint check for \(q=3,4\). The
receipt records the exact profiles and hashes. The source is
[`check_color_class_cube_profile.py`](../tools/check_color_class_cube_profile.py).
