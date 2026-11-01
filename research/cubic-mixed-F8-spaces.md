# Full mixed cubic coefficient spaces over F8

## Result and scope

Five structured pairs of output maps have been searched in their **full combined coefficient spaces** over
\(F=\mathbf F_2[z]/(z^3+z+1)\). Every case fails the strict seven-label seed criterion.
The two cyclic pairs generate the same seven-dimensional space, so the five configurations do not represent five distinct spaces.

This search allows linear combinations of coefficient directions from both maps before either flattening rank is evaluated.
It therefore goes beyond mixing the individually affordable tensors of previously searched spaces.
The conclusion here concerns the recorded spaces over F8. It does not by itself cover larger coefficient fields or tensor powers.

For a tensor \(T\), let
\[
 c(T)=\sqrt{\operatorname{rank}R(T)\operatorname{rank}C(T)}.
\]
A strict seed would decompose the crossing \(E\) with total cost below \(7^{\sqrt2}\).
The exact inequalities
\[
 16^{70}>7^{99},\qquad 99^2>2\cdot70^2
\]
give \(7^{\sqrt2}<16\). Thus every nonzero piece of any strict seed has rank product below 256.

## Output maps and spaces

Use the cubic input maps and complementary output notation of
[the general complementary-map identity](cubic-complementary-output-maps.md).
The fixed input words are \(\Omega e_i\), the output map has \(\Pi\Omega e_i=e_i\), and
\(C=\Pi\Omega f\) specifies the complementary images.
For every output map, retain all five raw polynomial coefficients in the order
\[
 (r^3,\;pr^2,\;p^2r,\;qr^2,\;q^3).
\]
Concatenate the two five-direction lists, remove exact F8 dependencies, and keep the first direction \(E\).
No rank-based pruning is applied to the coefficient directions.

Let \(S_{ij}=1\) when \(j=i+1\pmod7\), and let \(J_{ij}=1\) when \(j=i+1\) without wraparound.
Let
\[
 u=(1,z,z^2,z+1,z^2+z,z^2+z+1,z^2+1),
\]
\[
 v=(1,1,z,z,z^2,z^2,z+1),\qquad R=uv^T.
\]
The last configuration uses elementary matrices in two different output rows.

| Output-map pair | Full dimension | Projective vectors | Rank product below 256 | Affordable-span dimension |
|---|---:|---:|---:|---:|
| \(S,zS\) | 7 | 299,593 | 99 | 6 |
| \(zS,z^2S\) | 7 | 299,593 | 99 | 6 |
| \(J,zJ\) | 8 | 2,396,745 | 283 | 7 |
| \(R,zR\) | 8 | 2,396,745 | 8,786 | 7 |
| \(zE_{01},z^2E_{23}\) | 7 | 299,593 | 1,185 | 6 |
| **Total configurations searched** | | **5,692,269** | **10,452** | |

For each case, the affordable span is precisely the hyperplane whose first selected-basis coordinate is zero.
The crossing has first coordinate one. Consequently, no sum of affordable tensors can equal the crossing.
This proves the finite-field exclusion without a shortest-path assumption or a bound on the number of pieces.

## Scalar-family reduction

The independent audit verifies that each of the first four combined spaces also contains every raw coefficient direction for \(C=0\).
For a fixed matrix \(D\), each raw direction for \(C=\lambda D\) is a polynomial of degree at most two in \(\lambda\):
the two output legs each depend linearly on \(\lambda\).

The three distinct parameters \(0,1,z\) give an invertible Vandermonde matrix over F8.
Therefore the span of their values contains all three formal coefficient matrices of this quadratic polynomial.
After extending scalars, the same space contains the directions for **every** \(\lambda\) in every extension of F8.
The audit checks both flattenings with the same coefficient coordinates.

In particular, the first cyclic space contains every scalar-output family \(C=\lambda S\), and the second cyclic basis is another basis of that same space.
The Jordan and rank-one pairs have the analogous scalar-family property.
This is a reduction of the parameter spaces; the finite F8 rank enumeration alone does not prove a rank bound over field extensions.

## Verification and reproduction

The native evaluator enumerates all \((8^D-1)/7\) projective coefficient vectors, with a unique highest nonzero coefficient equal to one.
It computes the row rank exactly and limits the column elimination only after it has already exceeded the permitted product cap.
Every affordable vector receives two complete ranks.
The saved logs include each affordable projective vector, both ranks, the exact product, and the complete enumeration count.

The Python explorer independently recomputes all 10,452 affordable rank pairs, plus twelve controls per configuration, using a different pivot order.
This gives 10,512 independently checked vectors.
The separate audit verifies source hashes, input hashes, log hashes, projective normalization, all five separating functionals, and full raw-direction reconstruction in both flattenings.
It also verifies 24,010 physical occupation-kernel coefficient entries for the two distinct-row maps and repeats the first complete native enumeration.

Files:

- [Mixed-space explorer](../tools/explore_cubic_mixed_gf8.py)
- [Independent audit](../tools/check_cubic_mixed_gf8.py)
- [Complete finite search archive](../evidence/cubic-mixed-F8.json)
- [Audit certificate](../evidence/cubic-mixed-F8-audit.json)
- [Native evaluator](../tools/screen_cubic_gf8.cpp), shared with the earlier [F8 probe](cubic-complement-F8-probe.md)

Replay the saved rank checks and audit without repeating the complete enumeration:

```bash
python tools/explore_cubic_mixed_gf8.py \
  --replay evidence/cubic-mixed-F8.json --output /tmp/cubic-mixed-F8-replay.json
python tools/check_cubic_mixed_gf8.py --output /tmp/cubic-mixed-F8-audit.json
```

To regenerate all five finite searches:

```bash
c++ -std=c++20 -O3 -Wall -Wextra -Werror tools/screen_cubic_gf8.cpp -o /tmp/cubic-gf8-screen
python tools/explore_cubic_mixed_gf8.py \
  --native /tmp/cubic-gf8-screen --output /tmp/cubic-mixed-F8.json
python tools/check_cubic_mixed_gf8.py \
  --native /tmp/cubic-gf8-screen --output /tmp/cubic-mixed-F8-audit.json
```

The exact finite result rules out these full F8 spaces as strict one-copy seeds. It does not improve the verified factoring exponent.
