# A closed fusion row and an all-degree invariance proof

## Result

The bit-by-symmetric-power row in the literal six-vertex constructor has
an explicit degree-one matrix after its complete common scalar factor is
removed. Eight elementary Laurent-polynomial identities prove the formula
for every positive symmetric degree. Three more identities show that
successive rows preserve the symmetric relation on adjacent carrier bits.
Together they prove invariance of the full fusion word for every positive
degree in characteristic two.

This removes the missing all-degree invariant-subspace argument for the
specific repository convention. It does not assert a scalar identity or
even an invertible endpoint for every degree. The endpoint and complete
cut profiles remain separate conditions. The finite root-of-unity cases
and their joint endpoint normalization are recorded in
[the compressed fusion note](compressed-fusion-root-resonance.md).

## 1. The exact row formula

Let the characteristic be two, and let `q` be any nonzero field element
other than one. This section uses `q` for the six-vertex parameter, not
for the physical alphabet. Let `S_{m,b}` be the inversion-weighted symmetric
vector from the preceding note. A row moves one carrier bit across `m`
bits, with consecutive gate parameters

\[
 z,zq^{-2},\ldots,zq^{-2(m-1)}.
 \tag{1}
\]

Write `L_m[c,d;a,b]` for its divided matrix from input bit `a` and input
symmetric weight `b` to output symmetric weight `c` and output bit `d`.
Its only potentially nonzero entries are

\[
\begin{aligned}
 L_m[b,0;0,b]&=q^{b-1}+q^{1-b}z,\\
 L_m[b-1,1;0,b]&=(q+q^{2b-2m-1})z,\\
 L_m[b+1,0;1,b]&=q^m+q^{m-2b-2},\\
 L_m[b,1;1,b]&=q^{m-b-1}+q^{1-m+b}z.
\end{aligned}
 \tag{2}
\]

Omit any entry whose symmetric input or output weight lies outside
`0,...,m`. The complete raw row is

\[
 C_m(z)L_m(z),\qquad
 C_m(z)=\prod_{j=1}^{m-1}\bigl(q^{-1}+q^{1-2j}z\bigr).
 \tag{3}
\]

For literal row index `i`, substitute `z q^{2i}` everywhere in (2) and
(3). This is precisely the original gate schedule, without a change of
basis, variable convention, or endpoint.

The coefficient `L_m[1,0;1,0]=q^{m-1}(q+q^{-1})` is a nonzero constant.
The two diagonal families include nonzero linear coefficients. Thus
`L_m` has polynomial gcd one and actual degree one. Consequently (3),
up to its leading scalar, is the complete common gcd of the raw row.
The implementation normalizes `C_m` to monic and multiplies `L_m` by
the removed leading scalar so their product stays exactly unchanged.

## 2. Eight polynomial identities prove the row induction

The proof uses the last-bit identity

\[
 S_{m,b}=q^{-b}S_{m-1,b}\otimes0+S_{m-1,b-1}\otimes1.
 \tag{4}
\]

Assume (3) for `m-1`. Split the input by (4), apply the first `m-1`
gates, and then apply the final elementary gate at `zq^{2-2m}`.
After extracting the common factor `C_{m-1}`, there are only eight
coefficient identities, indexed by input carrier `a`, output carrier `d`,
and last symmetric output bit `e`.

For a compact exact verification, introduce independent Laurent variables
`Q,N,B,Z`, and subsequently specialize

\[
 Q=q,\qquad N=q^m,\qquad B=q^b,\qquad Z=z.
 \tag{5}
\]

The four entries in (2), viewed as functions of the two carrier bits, are

\[
\begin{aligned}
 l_{00}&=Q^{-1}B+QB^{-1}Z,\\
 l_{01}&=(Q+Q^{-1}B^2N^{-2})Z,\\
 l_{10}&=N+Q^{-2}NB^{-2},\\
 l_{11}&=Q^{-1}NB^{-1}+QN^{-1}BZ.
\end{aligned}
 \tag{6}
\]

In the previous row, replace `N` by `N/Q` and `B` by `B/Q^t`, where `t`
is the last input bit. Multiply by `B^{-1}` for `t=0`, or by one for
`t=1`, and then by the elementary gate coefficient from carrier `h` and
bit `t` to bit `e` and carrier `d`. Conservation fixes `h=d+e-t`.
Sum over the allowed `t`.

For each of the eight triples `(a,d,e)`, direct multiplication of (6)
gives

\[
 \sum_t B^{-(1-t)}
 l_{a,h}(Q,N/Q,B/Q^t,Z)
 G_{e,d;h,t}(Q^2N^{-2}Z)
 =
 (Q^{-1}+Q^3N^{-2}Z)
 (Q^{d-a}B^{-1})^{1-e}
 l_{a,d}(Q,N,B,Z).
 \tag{7}
\]

The first factor on the right is `C_m/C_{m-1}`. The second is exactly
the last-bit coefficient of the target `S_{m,a+b-d}`. This proves the
induction. At `m=1`, (2) is the original elementary gate.

These are identities in the Laurent-polynomial ring
`F2[Q^±1,N^±1,B^±1,Z]`, before any finite-field specialization. The checker
represents a polynomial by its set of nonzero exponent tuples. Addition
is symmetric difference, and multiplication adds tuples and cancels equal
pairs. It stores both sides of all eight identities and their empty
residuals.

At physical boundary weights, terms with an invalid intermediate symmetric
weight are absent. The same formula remains valid: when a forbidden
intermediate weight could contribute to an allowed output, its companion
off-diagonal coefficient is zero. For example, the increasing coefficient
at weight `-1` is `q^m+q^m`, and the decreasing coefficient at weight
`m+1` is `(q+q)z`. The unique extreme output weights have only one
valid last bit. Thus no boundary term or quantum-binomial division is
being suppressed.

## 3. Three identities prove preservation of adjacent carrier pairs

A symmetric two-bit vector of weight one is

\[
 01+q^{-1}10,
 \tag{8}
\]

where the first written bit is the lower-position bit. The other two
weights are `00` and `11`.

Move the higher-position carrier across the symmetric block using
`L_m(q^2z)`, and then the lower-position carrier using `L_m(z)`. For each
of the three input pair weights `k=0,1,2`, let `A_01,A_10` be the two
coefficients of an output pair of weight one. Substitution of (6) gives

\[
 A_{10}=q^{-1}A_{01}.
 \tag{9}
\]

The checker verifies these three identities in the same Laurent ring as
(7). If the first row changes the middle weight by
`input_high-output_high`, the second row uses that exact shifted weight.
The input coefficient of `10` is `q^{-1}`. These two adjustments are
included explicitly in the recorded identities.

The boundary argument again uses the vanishing coefficients at intermediate
weights `-1` and `m+1`. Output pair weights zero and two have a unique bit
pattern and require no relation. Equation (9) therefore proves invariance
of the entire three-dimensional pair subspace for every positive `m`.

## 4. Consequence for the complete fusion word

The span of `S_{m,0},...,S_{m,m}` is exactly the simultaneous solution of
the adjacent-pair relations (8). Within each total weight, adjacent
interchanges connect every word, and each interchange fixes the ratio of
its two coefficients. That ratio produces `q^{-inv(w)}`. This argument
uses only invertibility of `q` and remains valid at roots of unity.

The rectangular-row formula already preserves the middle symmetric
block. For any two adjacent carrier positions, earlier rows act only on
the other carriers and the middle block, so their input pair relation is
unchanged. The two corresponding rows are consecutive and have spectral
parameters in ratio `q^2`. Equation (9) preserves their output pair
relation. Later rows do not act on those two output positions. Every
adjacent pair in the final carrier block therefore satisfies (8).

It follows that the full literal fusion word maps

\[
 \operatorname{span}\{S_{m,b}\}^{\otimes2}
 \quad\hbox{into itself}
 \tag{10}
\]

for every positive `m` in characteristic two. Removing the `m` common
row factors leaves a product of `m` degree-one matrices, so its polynomial
degree is at most `m`. Any further common gcd can be removed exactly.

The finite compressed constructor keeps its intermediate checks even
though (7) and (9) now prove their invariance in all degrees. They provide
an independent implementation control and preserve explicit finite
endpoint evidence.

## 5. Root-of-unity structure visible in the formula

If `q` has odd order `ell>1` and `m=ell`, the two coefficients that would
leave the middle weight range vanish:

\[
 L_{\ell}[0,1;0,1]=0,\qquad
 L_{\ell}[\ell,0;1,\ell-1]=0.
 \tag{11}
\]

Thus the rectangular row preserves the middle space spanned by weights
`1,...,ell-1`. On the quotient represented by weights zero and `ell`,
the carrier acts by the scalar `q^{-1}+qz`. This explains the invariant
middle sector observed in the order-three and order-fifteen finite
resonance cases. It does not by itself determine the joint-normalized
second-cut Smith profiles, which still require a complete calculation.

## 6. Reproduction

```sh
python3 tools/check_explicit_fused_row.py \
  --output evidence/explicit-fused-row.json
```

The exact formal identities are supplemented by 105 complete comparisons
with the separately checked recursive row constructor. They use F16
parameters of orders 3, 5, and 15, degrees 1 through 18, and both extreme
row shifts. All 4,302 admissible matrix entries and every complete removed
common factor agree. Source hashes and the full Laurent identities are
retained in the evidence.

The formula and proof use the existing repository's elementary matrix.
For the broader higher-spin six-vertex setting, see
[Mangazeev (2014)](https://arxiv.org/abs/1401.6494). No identity from that
source is needed to justify a field specialization in this note.
