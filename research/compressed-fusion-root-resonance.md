# Exact compressed fusion and root-of-unity normalization

## Result and scope

The literal six-vertex fusion construction can be evaluated and certified
without expanding its full space on `2m` bits. Two recurrences retain only
weights of symmetric blocks. Every step checks the coefficient identities
that justify that compression. The implementation uses polynomially many
field operations in `m`, and it agrees with all complete literal matrix
controls through the eight-letter alphabet.

This permits exact checks at the order-15 parameter `theta=4` in F16.
For `m=1,...,14`, the divided fusion matrix has degree `m` and a nonzero
scalar identity at the crossing endpoint. For the checked values
`m=15,...,18`, additional common factors lower the degree to 14. Their
endpoint matrices remain invertible but are not scalar identities.
Explicit joint normalization supplies strict crossings on 16 through 19
letters, all of degree 14.

These are seed-construction results. They do not establish an improved
crossing exponent or factoring score. In particular, joint normalization
can change the second cut, so its profiles must be computed anew.

## 1. Fix the existing convention

Use the matrix, word order, and field encoding from
[the literal construction](finite-fused-six-vertex-crossings.md).
For a word `w` on `n` bits, its inversion count counts pairs of positions
`i<j` with `w_i=1` and `w_j=0`. Write

\[
 S_{n,k}=\sum_{|w|=k}\theta^{-\operatorname{inv}(w)}|w\rangle.
 \tag{1}
\]

The word with zeros followed by ones has coefficient one, so these vectors
are independent over every field and at every nonzero value of theta.
No inverse quantum factorial or quantum binomial appears.

The two splitting identities are

\[
 S_{n,k}=\theta^{-k}S_{n-1,k}\otimes|0\rangle
          +S_{n-1,k-1}\otimes|1\rangle,
 \tag{2}
\]

\[
 S_{n,k}=|0\rangle\otimes S_{n-1,k}
          +\theta^{-(n-k)}|1\rangle\otimes S_{n-1,k-1},
 \tag{3}
\]

where out-of-range weights contribute zero. These follow directly by
counting inversions involving the removed last or first bit.

The literal fusion word processes `i=m-1,...,0`. In row `i`, its successive
gates use the parameters

\[
 z\theta^{2i},z\theta^{2(i-1)},\ldots,z\theta^{2(i-m+1)}.
 \tag{4}
\]

That convention is preserved exactly below.

## 2. A checked rectangular row

For fixed `i`, construct the map from one bit and `S_{n,k}` to a symmetric
block of length `n` and one output bit. Start at `n=0`, where the map is the
identity on the carrier bit. For the next step, split the input symmetric
block by (2), apply the already constructed row on the first `n-1` bits,
and apply the final elementary gate to the carrier and the last bit.

The result is expressed in the independent basis

\[
 S_{n-1,c}\otimes|e\rangle\otimes|d\rangle.
\]

For each fixed input pair and output weight `c` and carrier `d`, let
`A_0,A_1` be the coefficients with last output bit zero and one. For
`0<c<n`, membership in the desired symmetric block is exactly

\[
 A_0=\theta^{-c}A_1.
 \tag{5}
\]

At `c=0` and `c=n`, only one term is possible. The implementation checks
(5) as an exact polynomial identity for every input and output. It then
retains `A_1` when `c>0`, or `A_0` when `c=0`. Consequently, successful
execution certifies the entire rectangular row on its full input space.

After `n=m`, divide every row coefficient by its exact common monic
polynomial gcd. This is a scalar division shared by the whole row. It
preserves the just-checked invariant subspace. The original scalar factor
and its degree remain recorded.

For each checked parameter in the evidence, every divided rectangular row
has degree one. The correctness argument itself does not assume that
degree: every elementary product has degree at most `m`, which is enough
for polynomial complexity.

## 3. A checked full fusion word

After `p` complete rows, retain coefficients in

\[
 S_{m-p,a}\otimes S_{m,b}\otimes S_{p,c}.
 \tag{6}
\]

To advance, split the remaining left block by (2). Apply the rectangular
row from Section 2 to its last bit and the middle block. Its output bit
lies immediately before the previously processed right block. To combine
that bit with `S_{p,c}`, use (3).

For an interior new right weight `r`, the exact condition is

\[
 B_1=\theta^{-(p+1-r)}B_0,
 \qquad 0<r\le p,
 \tag{7}
\]

where `B_e` is the coefficient with the new first bit equal to `e`.
The implementation checks every such identity, for every original input
pair and current remaining and middle weights. The two boundary weights
have a unique bit. At the last step, (6) is precisely the desired tensor
product of two symmetric blocks.

Induction on the two recurrences proves the following finite certification
statement: if every check succeeds, the retained polynomial matrix is
exactly the restriction of the original full-bit fusion word, divided only
by the recorded common scalar factors. This remains true at roots of
unity, including cases with vanishing quantum binomial coefficients.

There are at most three integer block weights in each state. Every weight
is at most `m`, the number of stages is `m`, and intermediate polynomial
degrees are at most `m^2`. These bounds give polynomial time and storage
with elementary field arithmetic. The procedure raises an error if a
coefficient identity fails; it does not assume an unchecked all-spin
invariance theorem.

## 4. The order-15 endpoint needs a joint normalization

The parameter encoded by 4 has multiplicative order 15 in the specified
F16. The checked table is:

| theta | m | Alphabet | Removed common degree | Final degree | Raw endpoint |
|---:|---:|---:|---:|---:|---|
| 4 | 1 through 14 | 2 through 15 | `m(m-1)` | `m` | Nonzero scalar identity |
| 4 | 15 | 16 | 211 | 14 | Invertible, not a scalar identity |
| 4 | 16 | 17 | 242 | 14 | Invertible, not a scalar identity |
| 4 | 17 | 18 | 275 | 14 | Invertible, not a scalar identity |
| 4 | 18 | 19 | 310 | 14 | Invertible, not a scalar identity |
| 2 | 3 | 4 | 7 | 2 | Invertible, not a scalar identity |
| 2 | 4 | 5 | 14 | 2 | Invertible, not a scalar identity |

The last two are useful small resonance controls. They were correctly
rejected by the old constructor's scalar-identity endpoint test. The new
operation is an additional explicit normalization, not a change to the
outcome of that test.

For a divided matrix `F(z)` with invertible endpoint `E=F(1)`, set

\[
 \widehat F(z)=F(z)E^{-1},\qquad
 T_{xyuv}(z)=\widehat F_{yx,uv}(z).
 \tag{8}
\]

Then `T(1)=[x=v][y=u]`. Equation (8) does not increase degree. A constant
invertible right factor cannot annihilate the complete leading coefficient
matrix, so the actual matrix degree is unchanged.

The multiplication by `E^{-1}` acts jointly on the two input legs. It
preserves every local Smith profile of the `xy|uv` cut. It generally does
**not** preserve the `xu|yv` cut, because those two input legs lie on
opposite sides of that cut. All subsequent search inputs must use the
complete normalized tensor and recompute that second profile.

## 5. Integer charge keeps the new matrices small

Every elementary gate preserves the number of ones. Therefore

\[
 F_{cd,ab}(z)=0\quad\hbox{unless}\quad a+b=c+d.
 \tag{9}
\]

The endpoint and its inverse split into the same integer-charge blocks.
Their largest block has size `m+1`. Thus the joint normalization can be
performed block by block, without a dense inverse of a matrix with
`(m+1)^2` rows.

After the relabeling in (8), both flattenings still split into square
blocks of size at most `m+1`. The first uses equal total charges. The
second uses `u-x=y-v`. This supports exact local Smith computation through
cached fraction-free elimination, with no enumeration of all minors.
For completeness, the elimination cache has the following algebraic
invariant. A state indexed by selected row and column sets `(I,J)` retains
the minor `det M[I,J]` and every minor obtained by adjoining one remaining
row and column. The next state follows from the bordered-minor identity:
multiply its two diagonal borders, add the two off-diagonal borders, and
divide by the previous selected minor. This division is exact. In
characteristic two, the values depend only on the selected sets, so paths
that reach the same state can share the cache without an ordering sign.

At a local place, multiply rows and columns by their requested powers of
the uniformizer. Choose a residual entry of minimum valuation as pivot.
Every quotient used to clear its row and column then has nonnegative
valuation. The clearing operations are invertible over the valuation
ring, and splitting off the pivot leaves the next local Smith problem.
Its valuation is therefore the next complete Smith order. The cached
bordered-minor quotient is exactly that residual entry before the separate
row and column powers are restored. This proves the evaluator for all
physical weights; a common ramification clears rational denominators.
At infinity, the degree correction adds the homogeneous degree once per
pivot, or `k` times that degree for a `k`-minor. Nonlinear closed places use
the same argument in their residue-field extension.

The independent controls also compare 66 full weighted profiles against
exhaustive minors on the smaller seeds, including nonlinear closed places
and negative rational physical heights.

## 6. Reproduction and independent controls

```sh
python3 tools/check_compressed_fusion.py --output evidence/compressed-fusion.json
python3 tools/check_f16_weighted_local_oracle.py \
  --output evidence/f16-weighted-local-oracle.json
```

The first command includes:

- 15,312 independent full-bit entries for rectangular rows;
- 9,733 complete matrix-entry comparisons with the original literal
  constructor or its separately certified eight-letter caches;
- 54,404 exact rectangular symmetry checks and 556,944 growing-output
  symmetry checks across the recorded cases;
- 563,627 strict endpoint entries after the explicitly recorded scalar or
  joint normalization.

Hashes identify the full normalized matrices, the raw endpoint matrices,
their inverses, and every source file. Full matrix caches for the original
seven- and eight-letter cases are retained separately. They are not
silently replaced by the compressed construction.

For context, higher-spin six-vertex matrices and their explicit formulas
are developed in [Mangazeev, *On the Yang-Baxter equation for the six-vertex
model* (2014)](https://arxiv.org/abs/1401.6494). The finite certification here
uses only the elementary repository matrix and the inversion identities
(2) and (3). No generic-characteristic representation formula or
root-of-unity limiting identity from that literature is imported.
