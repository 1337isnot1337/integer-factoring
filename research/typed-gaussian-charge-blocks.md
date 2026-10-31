# Exact charge-block compression for typed tensor powers

Integer charge conservation makes large typed coefficient matrices
computable without constructing their full flattenings. For the four-letter
Gaussian seed, the two flattenings split into blocks indexed by a charge
vector. Simultaneous site permutations identify every charge vector with
the same five occupation counts.

The resulting rank formula is exact. It applies to every polynomial
coefficient and every linear combination of coefficients. The implementation
below reaches a physical alphabet of dimension \(1260\), whose flattenings
have dimension \(1\,587\,600\), using 38 distinct blocks, none larger than
\(1260\).

The computed rank costs remain above the post-half target. This is a
reusable compression method and a finite search checkpoint.

## 1. Typed words and charge

Let the physical letters be \(0,+,-,\pm\), encoded as \(0,1,2,3\), with

\[
q(0)=q(\pm)=0,\qquad q(+)=1,\qquad q(-)=-1.
\]

Fix a common word type \(c=(c_0,c_1,c_2,c_3)\), of length
\(m=\sum_i c_i\), on all four physical legs. Its dimension is

\[
D=\frac{m!}{c_0!c_1!c_2!c_3!}.
\]

Suppose a local polynomial tensor \(P_{xyuv}(w)\) satisfies

\[
P_{xyuv}(w)\ne0
\quad\Longrightarrow\quad
q(x)+q(u)=q(y)+q(v). \tag{1}
\]

The Gaussian seed in
[the Schur reduction](schur-gaussian-digit-reduction.md) satisfies (1). Its
ordinary tensor power therefore conserves charge separately at every site.
Write \(A_k\) for the coefficient tensors after restriction to the common
word type.

For the row cut, index the columns by \((v,u)\). Equation (1) becomes

\[
q(x)-q(y)=q(v)-q(u).
\]

For the column cut, the same equation reads

\[
q(x)+q(u)=q(y)+q(v).
\]

Both flattenings are block diagonal. The charge at each site belongs to
\(\{-2,-1,0,1,2\}\).

## 2. Orbit blocks and the rank formula

For a charge vector \(\gamma\in\{-2,-1,0,1,2\}^m\), let
\(\eta_j\) count its entries equal to \(j\). Then

\[
\sum_{j=-2}^{2}\eta_j=m.
\]

In the row cut the charge sum must be zero. In the column cut it must be
\(2(c_1-c_2)\). Thus only signatures satisfying

\[
\sum_{j=-2}^{2}j\eta_j=
\begin{cases}
0,&R,\\
2(c_1-c_2),&C
\end{cases}
\]

can occur.

Simultaneously permuting the sites on all four physical legs preserves the
tensor and the common word type. Consequently all blocks with the same
signature \(\eta\) are identical up to row and column permutations. There
are

\[
N(\eta)=\frac{m!}{\prod_{j=-2}^{2}\eta_j!}
\]

charge vectors with this signature.

For either cut \(\tau\), choose the canonical sorted charge vector of
signature \(\eta\), and let \(M_{\tau,\eta}(B)\) denote its block in a
coefficient combination \(B=\sum_k b_kA_k\). Then

\[
\boxed{\operatorname{rank}_{\tau}(B)=
\sum_{\eta}N(\eta)\operatorname{rank}M_{\tau,\eta}(B).} \tag{2}
\]

The blocks are square because each physical leg uses the same type. Formula
(2) holds over any coefficient field for a tensor satisfying (1). The
implementation uses \(\mathbb F_2\).

There are at most \(\binom{m+4}{4}\) signatures before imposing the charge
sum and feasibility constraints. This is polynomial in \(m\), although an
individual block can still be exponentially large.

### Block sizes

The block sizes can be computed without constructing a tensor. Introduce
four formal variables \(X_a\) and four \(Y_b\). For the row cut, the size for
a canonical charge vector is

\[
[X^cY^c]\prod_{i=1}^{m}
\left(\sum_{q(a)-q(b)=\gamma_i}X_aY_b\right).
\]

For the column cut, replace the difference by \(q(a)+q(b)\). The native
probe computes these coefficients by recursively enumerating only the
admissible ordered word pairs. It also verifies

\[
\sum_\eta N(\eta)\dim M_{\tau,\eta}=D^2
\]

for each cut.

## 3. Coefficients and recombinations

In the chosen Gaussian gauge,

\[
\deg P_{xyuv}\le \operatorname{occ}(x)+3\operatorname{occ}(v),
\qquad
\operatorname{occ}=(0,1,1,2).
\]

The common type therefore gives the exact safe degree bound

\[
\deg P^{\otimes m}\le 4(c_1+c_2+2c_3).
\]

The probe checks the local charge condition, this degree bound, and the
crossing endpoint before using them. It builds each block coefficient over
\(\mathbb F_2\), computes binary ranks, and applies (2).

For small coefficient spans it also enumerates all binary recombinations.
Only even degrees occur in the recorded examples. A bit in code position
\(j\) selects coefficient \(A_{2j}\). Gray-code updates change one
coefficient matrix at a time.

The complete rank table gives a rigorous lower bound for every
\(\mathbb F_2\) coefficient-span decomposition of the crossing. Give code
\(b\) the integer weight

\[
\left\lfloor
\sqrt{\operatorname{rank}_R(B_b)\operatorname{rank}_C(B_b)}
\right\rfloor.
\]

An exact shortest path in the XOR Cayley graph then lower-bounds every
actual rank cost. The checker verifies independence of the coefficient
list, so the crossing is represented by the all-ones code. It also records
an explicit decomposition and a rational upper bound for its cost.

## 4. Recorded results

The complete exact data are in
[typed-gaussian-charge-blocks.json](../evidence/typed-gaussian-charge-blocks.json).

| Type \(c\) | \(D\) | Binary coefficient combinations | Certified cost information |
| --- | ---: | ---: | --- |
| \((1,1,1,0)\) | 6 | 32 | Exact optimum \(26\) |
| \((2,1,1,0)\) | 12 | 32 | Exact optimum \(90\) |
| \((1,1,1,1)\) | 24 | 512 | Exact optimum \(324\) |
| \((2,1,1,1)\) | 60 | 512 | Exact optimum \(1890\) |
| \((6,2,2,0)\) | 1260 | 512 | Lower bound \(291966\); witness cost \(<291967.396\) |

These optima concern binary recombinations of this coefficient span.
They do not optimize over extension-field coefficients or other tensors.

For \(D=60\), the optimum uses codes \(1,254,256\). The two outer pieces
each cost \(60\), and the middle piece has both ranks \(1770\).

For \(D=1260\), the recorded witness uses codes
\(1,256,6,192,8,32,16\). Its coefficient ranks are

\[
\begin{aligned}
R_{0,2,\ldots,16}
={}&(1,98,2024,15118,95116,249110,759360,760086,1587600),\\
C_{0,2,\ldots,16}
={}&\operatorname{reverse}(R_{0,2,\ldots,16}).
\end{aligned}
\]

The evidence also records the higher-density type \((4,3,3,0)\), with
\(D=4200\), 80 orbit blocks, and largest block size \(4200\). Its raw
coefficient cost is approximately \(4.97\cdot10^6\), which is not
competitive.

All five finite-span lower bounds exceed \(D^{17/12}\), checked by integer
exponentiation. Since \(17/12>\sqrt2\), they exceed the post-half crossing
target \(D^{\sqrt2}\).

## 5. Reproduction and independent controls

The native command-line tool is
[probe_typed_gaussian_charge_blocks.cpp](../tools/probe_typed_gaussian_charge_blocks.cpp).
The replay script compiles it, reconstructs the local Gaussian seed, and
checks all recorded cases:

~~~sh
python tools/check_typed_gaussian_charge_blocks.py \
  --output evidence/typed-gaussian-charge-blocks.json
~~~

For the first three types, an independent Python implementation constructs
the full typed flattenings and compares every recombination rank in both
cuts: 576 complete code comparisons in total. The native implementation
uses charge blocks; the independent implementation uses simultaneous site
permutations to build full matrices. The checker also supplies an invalid
charge-breaking seed and verifies that the native program rejects it.

To use the native probe separately, provide the 256 local polynomial bit
masks in lexicographic \((x,y,u,v)\) order:

~~~sh
g++ -O3 -std=c++17 tools/probe_typed_gaussian_charge_blocks.cpp \
  -o /tmp/probe_typed_gaussian_charge_blocks
/tmp/probe_typed_gaussian_charge_blocks seed.txt 6,2,2,0 20000 all
~~~

The block-size argument bounds the largest permitted block. The optional
final argument requests all binary recombinations. Explicit size and
degree guards prevent the packed arithmetic from silently overflowing.
