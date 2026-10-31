# Broader binary polynomial crossing catalogues

The [quadratic identity-start catalogue](binary-quadratic-catalogue.md)
extends in two independent directions. Neither gives a smaller bare Smith
pressure or a stronger factoring score. The strongest complete factoring
exponent remains 1/2. No Lean files change.

## Arbitrary column-rank-one constant tensors

Let A and B be any nonzero binary two-by-two matrices and put
C(x,y,u,v)=A(x,u)B(y,v). Exhaust all tensors M over F2 in

```math
K(z)=C+zM+z^2(E+C+M).
```

Thus K(0)=C has column-cut rank exactly one and K(1)=E. There are
225 possible C, each with 65,536 middle tensors. Every such path has
bare pressure at least 256/27 for the same ungauged row/infinity and
column/zero cuts as the preceding note. The bound is attained when C
is equivalent to I. It does not cover a zero constant tensor or arbitrary
higher-column-rank constant tensors.

The complete reduction uses constant invertible physical changes

```math
(A,B)\longmapsto(SAT^{-1},TBS^{-1}),\qquad S,T\in\mathrm{GL}_2(\mathbb F_2).
```

They preserve E and act invertibly on each side of both cuts, so they
preserve both Smith profiles and biject the entire middle coefficient
space. Represent two-by-two matrices by bits 2i+j. The six invertible
matrices are 6,7,9,11,13,14, with identity 9. The eleven orbit
representatives (A,B), with their sizes, are

| A | B | Size |
| --- | --- | --- |
| 1 | 1 | 36 |
| 1 | 2 | 18 |
| 1 | 4 | 18 |
| 1 | 6 | 18 |
| 1 | 7 | 36 |
| 1 | 8 | 9 |
| 6 | 1 | 18 |
| 6 | 2 | 36 |
| 6 | 6 | 6 |
| 6 | 7 | 18 |
| 6 | 11 | 12 |

Their sizes sum to 225. Over F2, nonzero outer products have unique
nonzero factors, so there is no extra scalar identification. The exact
[checker](../tools/check_extended_binary_crossings.py) constructs all
orbits from scratch. The [eleven catalogues](../evidence/rank-one-constant-quadratic-catalogues.json)
contain the complete counts and first witnesses for 720,896 representative
paths, covering 14,745,600 paths before quotienting by symmetry.

## Complete parity-preserving cubics and quartics

Keep K(0)=I and K(1)=E, and require every entry with x+y+u+v odd to
vanish identically. Each remaining entry can have arbitrary intermediate
coefficients. There are 65,536 cubic paths and 16,777,216 quartic paths.
The complete catalogues have respectively
[671](../evidence/parity-cubic-crossing-catalogue.json) and
[3,999](../evidence/parity-quartic-crossing-catalogue.json) distinct paired
profiles. All have pressure at least 256/27 at every real scale s>=1.

These spaces include non-Gaussian tensors. For example, taking only
M(0,0,0,0)=1 in a quadratic interpolation makes the even Gaussian
relation ah+bg+cf+de nonzero as a polynomial. They remain small local
alphabets, and the result does not constrain general higher-alphabet
families, gauges, blocks, or combined coefficient ranks.

The parity condition makes each cut a direct sum of two two-by-two
polynomial matrices. In physical entry order 0,3,5,6,9,10,12,15, call
the polynomials a,b,c,d,e,f,g,h. The row blocks are (a,b;g,h) and
(c,d;e,f); the column blocks are (a,c;f,h) and (b,d;e,g). For any such
block, its first Smith order is the minimum entry valuation and the
sum of its two orders is the determinant valuation. At t=1 both blocks
are invertible, so their determinants are nonzero polynomials.

This identity allows a fast exhaustive quartic enumeration. It agrees
byte for byte with the full-minor method on all 65,536 cubics. For the
quartics, the independent checker uses full four-by-four block-Toeplitz
ranks on every saved representative plus 4,096 scheduled tensors; it
does not independently replay all 16,777,216 quartics. This limitation
is distinct from the native enumeration's exhaustive coefficient coverage.

## Exact scale certificates and unstructured probes

The checker removes profiles dominated coordinatewise by another profile.
There are ten maximal cubic profiles and twenty maximal quartic profiles.
For each, it proves 27P(x)-256x^d>=0 on [0,1], where d is the
regularization degree. Exact Bernstein subdivision is used as in the
quadratic note. Whenever necessary it first divides out the exact square
(1-3x^r)^2, with 1<=r<=d. No floating-point minimum is used in this
certificate. Monomially reparametrized elementary paths attain the bound.

A separate [million-path cubic sample](../evidence/cubic-crossing-million-sample.json)
uses all sixteen unrestricted entries in both intermediate coefficients.
It visits masks (i*0x9e3779b1+0x7f4a7c15) modulo 2^32, for
0<=i<1,000,000. The odd multiplier makes all sampled masks distinct.
It records 1,061 paired profiles; none beats the same certified floor.
This is not exhaustive in the 2^32-element cubic space. Every saved
representative and the first 256 scheduled masks receive independent
Toeplitz checks; the sample multiplicities are produced by the native
enumerator and are not independently replayed in full.

The [sampled replay](../evidence/extended-binary-crossing-sampled-replay.json)
separates catalogued counts from independently checked counts. The full
structured replay, when run with --full, independently checks all eleven
quadratic catalogues and all parity cubics, retaining the stated sampled
coverage for quartics and unrestricted cubics.
The completed [full smaller-catalogue replay](../evidence/extended-binary-crossing-full-replay.json)
checks all 786,432 structured tensors. It predates the quartic addition;
the separate sampled record above supplies the 8,092 quartic controls.

## Reproduce

```sh
c++ -O3 -std=c++17 -Wall -Wextra -Werror tools/enumerate_quadratic_crossings.cpp -o /tmp/enumerate-polynomial
/tmp/enumerate-polynomial --constant 1 1
/tmp/enumerate-polynomial --parity-cubic
/tmp/enumerate-polynomial --parity-degree 3
/tmp/enumerate-polynomial --parity-degree 4
/tmp/enumerate-polynomial --cubic-sample 1000000
python3 tools/check_extended_binary_crossings.py
python3 tools/check_extended_binary_crossings.py --full
```

For the eleven constant catalogues, repeat --constant for the table's
representatives. The original no-argument enumerator remains unchanged
in output. These calculations build on the existing polynomial/Smith
framework, with substantial OpenAI Codex assistance. The next unresolved
construction problem is a useful higher-degree path without the parity
restriction, or an effective larger local alphabet; a finite pressure
improvement would still need the complete factoring transfer.
