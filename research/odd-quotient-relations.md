# Smooth-pair relations modulo a certified odd kernel

## Question and result

Must a smooth-pair witness give an exact principal relation before the
[deterministic two-part extraction](arithmetic-torsion-extraction.md) can
use it? No. A nonsingular relation matrix may have nonidentity row residues
of **verified odd order**. A charged odd multiplier converts its determinant
to an annihilator without changing its two-adic valuation.

This permits the [positive relation kernel](positive-relation-kernel.md)
to operate on the computable image `H^M` for a known odd integer `M`, and
replaces exact smooth-pair equality by equality after the Mth-power map.
The image preserves every two-part factoring direction. A uniform
multiplier can remove all small odd primary components without knowing
the class number. The pair-extraction problem is still open; this note
does not improve the complete factoring score **1/2**.

Scope was registered in [Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360508),
base `2e2a126`, with the implementation integrated against `517909e`.
The mechanism changes the equality required of a relation witness. It is
not a new height, orbit or prime-parameter sweep of the stopped
[strict-descent family](prime-power-descent-obstruction.md).

## 1. Odd row residues suffice

Let `g_1,...,g_m` generate a finite abelian group H. Let R be a nonsingular
m-by-m integer matrix, with row residues

```math
h_i=\prod_j g_j^{R_{ij}}.
```

Suppose supplied positive **odd** integers `M_i` pass the verifiable tests
`h_i^(M_i)=1`. Set

```math
M=\operatorname{lcm}_i M_i,\qquad E=M|\det R|.       \tag{1}
```

Every row of MR is an exact relation. The adjugate identity for R implies
that `M det(R) e_j` is a linear combination of those rows. Therefore E
annihilates each generator and H. Since M is odd,

```math
v_2(E)=v_2(\det R).                                 \tag{2}
```

The existing determinant-projection and two-group basis algorithm applies
unchanged. In particular it does not require saturation, a class number,
factoring E, or random relation lifts.

The [certificate implementation](../tools/odd_quotient_relations.py)
checks the row powers using only the supplied group operations. It reuses
the exact Bareiss determinant from the positive-kernel checker. If V is the
encoded matrix length and W is the total bit length of the M_i, row
verification uses `O(V+W)` group operations. Determinant and lcm acquisition
cost polynomially many bit operations in their encoded inputs, and

```math
\operatorname{bits}(E)\le
\operatorname{bits}(|\det R|)+W+1.
```

For `log_2|H|<=L` and two-rank r, the imported basis algorithm then uses
`O(m bits(E)+mL(rL+2^ceil(r/2)))` group operations. Canonical group
identifiers and the bit cost of their operations remain charged in an
arithmetic application. An assertion that a residue has odd order, without
a checked odd annihilator, is not this certificate.

## 2. A computable odd quotient preserves density and two-parts

Fix a known positive odd M and define `pi(h)=h^M`. Its kernel K has odd
order: every element of K has order dividing M. Its image `H^M` has

```math
(H^M)_2=H_2.                                       \tag{3}
```

Indeed powering by M is an automorphism on the two-part. Thus two-torsion
extraction and its map to coprime factorizations retain exactly the same
coverage. This does not repair missing generator/genus directions.

For any subset S of H of density delta, let `Sbar=pi(S)`. Every fiber of
pi has size |K|, so

```math
|\overline S|\ge |S|/|K|,\qquad
\frac{|\overline S|}{|H^M|}\ge\frac{|S|}{|H|}=\delta. \tag{4}
```

If elements of S have words of length at most A in the g_j, their images
have the same word vectors and length bound in the generators `g_j^M`.
The same applies to sample words, including their multiplicities. A
small-bias distribution on H pushes forward with no worse bias, because
characters of the image pull back to nonprincipal characters of H.
More generally one may construct the sample directly on `H^M`; it need
not mix the discarded odd directions at all. Its word vectors still give
known lifts in the original generators.

The predecessor's positive autocorrelation argument therefore supplies,
for every target t, some sample u and original representatives a,b in S
with

```math
\pi(tu)=\pi(ab^{-1}),\quad\text{equivalently}\quad
(tu a^{-1}b)^M=1.                                  \tag{5}
```

It suffices that the nonprincipal image-character bias is at most delta/2.
This is an existence statement. An algorithmic extractor must return the
original smooth representatives and their words, not just two image
identifiers. The image of a recognizable smooth set need not itself be
efficiently recognizable.

The predecessor's principal-word gap, safe-radius injection and event
separation do **not** transfer automatically to this image. As supplied in
[the scoped feedback](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360786),
at discriminant -23 the nonidentity prime form `(2,1,3)` has order three.
The M=3 image kills its nonzero word although its norm two is below 23/4.
Image collisions must therefore be retained; the original p-free boundary
bijection cannot be used unchanged. The positive autocorrelation theorem
is applied anew to the image set and requires no inherited injective ball.

To construct rows, take the usual integer `K_0>W+2A`, where W bounds sample
word lengths. Apply (5) to `t=g_i^(K_0)` and set

```math
r_i=K_0 e_i+w_u-v_a+v_b.
```

The row residue is killed by M. Its correction has l1 norm at most W+2A,
so the matrix is strictly diagonally dominant and nonsingular, exactly as
in the predecessor. Equation (1) now gives the required annihilator.
No exact principality test is substituted for the weaker condition (5).

## 3. Uniformly remove small odd primary components

Given `|H|<=2^L` and a numerical smoothness bound B, acquire

```math
M_{L,B}=\prod_{3\le\ell\le B\atop\ell\text{ prime}}
          \ell^{\lfloor L/\log_2\ell\rfloor}.        \tag{6}
```

Its largest prime powers are computed by integer multiplication and
comparison with `2^L`, without real logarithms or a factorization oracle.
Trial primality tests and school arithmetic give cost polynomial in L,B;
the dependence is on the **numerical** B. Also

```math
\operatorname{bits}(M_{L,B})\le L\pi(B)+1\le LB+1.
```

It kills every B-smooth odd order at most `2^L`. More precisely, the image
is H_2 times the odd primary components whose primes exceed B. The order
of that remaining odd group has at most `L/log_2 B` prime factors counted
with multiplicity when B>=2. This does not give its coordinates, orders,
or a fast smooth-pair extractor; even one large cyclic component may be
hard to handle.

There is one directly certifiable fast branch. Power the supplied generators
by M, then square each result at most L times. If every result reaches the
identity, the image is a two-group; apply the existing two-group basis
algorithm directly, without any relation matrix or pair extractor. This
test succeeds exactly when the odd part of the generated H is B-smooth.
It needs no class number. If a generator fails the test, its image retains
a prime-order component greater than B, and this branch makes no further
claim. The fixed class controls below pass this direct test as well as
the relation-certificate route. A product-group control with residual
odd order thirteen at B=11 correctly fails. This is not a smooth-class-number
theorem for every input.

For the existing arithmetic choice `s=ceil(sqrt(b log_2(b+2)))`,
`B=2^O(s)`, and `L=O(b)`, acquiring M and powering the factor-base
generators cost `2^O(s)` bit operations. With the predecessor's generator
and sample counts, word-bit bounds and trial preprocessing, the resulting
matrix and annihilator also have `2^O(s)` encoded size. All subsequent
two-part extraction fits that scale. The value of M can be enormous;
its **bit length**, all exponentiation and all group operations are paid.

Retain the predecessor's good-multiplier, generator-coverage and smooth-set
density assumptions exactly. If a uniform extractor for (5) on its stated
targets costs F(b), then the full conditional factoring invoice is

```math
(F(b)+1)2^{O(\sqrt{b\log b})}.                       \tag{7}
```

There are only `2^O(s)` targets/sample choices and all other work has that
cost. Thus (7) preserves any fixed positive logarithmic-saving exponent
supplied by a proved extractor. None is provided here. In particular the
quotient does not turn the existing single-form norm bound into a fast
search algorithm, nor remove its analytic/generator assumptions.

## 4. Fixed exact controls and attribution

```sh
make verify-odd-quotient-relations
```

The [checker](../tools/check_odd_quotient_relations.py) validates certificates
on mixed product groups, including nonidentity odd row residues and mixed
generators whose separate two-torsion tops miss a direction. It rejects
even witnesses, false odd witnesses and singular matrices. Fifty-four
fixed multiplier controls exhaust the indicated small odd smooth orders.

The positive-kernel fixture is lifted to `(Z/2)^4 x Z/3`. Its five-set has
density 5/48 and image density 5/16. The lifted six-point sample has full
bias one but image bias 1/3. It misses exact pair witnesses for 32 of the
48 targets, while every target has a quotient pair. The nonsingular row
certificate still recovers all four two-torsion directions. This tests a
genuine weakening of the witness and mixing conditions, not exact equality
with different notation. Finite witnesses are obtained by enumeration.

Two static class groups reuse the boundary forms supplied in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359806)
and its [above-height follow-up](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359827):

| Discriminant | Enumerated classes | Image classes | Smooth classes / image | Target | Smaller smooth preimage |
| --- | ---: | ---: | ---: | --- | --- |
| -7423 | 40 | 8 | 21 / 8 | (29,1,64) | (4,1,464) |
| -26623 | 48 | 16 | 13 / 8 | (13,1,512) | (2,-1,3328) |

Here B=11 and L is the discriminant bit length in (6). Each target/preimage
ratio is a **nonidentity** class killed by M. These are smaller norms in
different classes, so they do not refute the original strict-descent
obstruction. Enumeration supplies the finite class-group generators and
relations; it is not a proposed uniform algorithm. The full two-parts have
orders eight and sixteen, and the verified certificates recover factors
13 of 7423 and 79 of 26623. The
[receipt](../evidence/odd-quotient-relations.json) pins sources and results.

The generic two-group proof, positive Fourier kernel, word construction,
class-group arithmetic, ambiguous-form extraction and boundary examples
retain their linked predecessor credit. The algebraic relaxation and
controls were developed by Sol worker 2 with substantial OpenAI Codex
assistance. No novelty claim is made for localization at the prime two.
No new Lean, independent mathematical review, or stronger complete bound
is claimed. The next required progress is an actual efficient extractor
for (5); enlarging these fixed controls is not that progress.
