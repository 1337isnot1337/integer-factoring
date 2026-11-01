# Exact spectral transport through the invertible fusion channel

This note makes a fused spectral slide explicit over F16. A three-strand
projector crosses a fundamental strand by a scalar polynomial times
the constant braid in the vacuum quotient. The ordinary ambient
identity retains a nonzero height-four residual, so discarding it
does not assume an invariant physical subspace.

Projected spectral faces and fusion hierarchies are standard tools;
[Morin-Duchesne, Pearce and Rasmussen, sections3 and6](https://arxiv.org/pdf/1401.7750)
provide the design lead. Their transfer relations are not imported as
characteristic-two theorems here. The identities below are separately
verified in the fixed field, with a faithful local TL representation.
This supplies compatibility for a new spectral-kernel direction;
it does not establish a favorable whole-crossing rank sum.

## 1. Baxterized local gates

Use `q=8 in F16`, of order5, `delta=q+q^-1`, and

```math
B_i=qI+e_i,\qquad B_i^{-1}=q^{-1}I+e_i,\qquad
\mathcal R_i(z)=B_i^{-1}+zB_i.                  \tag{1}
```

The TL relations give the Hecke relation and braid relation. The
spectral gates satisfy the polynomial Yang--Baxter identity

```math
\mathcal R_i(x)\mathcal R_{i+1}(xy)\mathcal R_i(y)
=\mathcal R_{i+1}(y)\mathcal R_i(xy)\mathcal R_{i+1}(x),          \tag{2}
```

and the polynomial inversion identity

```math
\mathcal R_i(z)(B_i+zB_i^{-1})
  =(1+\delta^2z+z^2)I.
```

For example `B_i+B_i^-1=delta I`, so the coefficient of z in the
last product is `B_i^2+B_i^-2=delta^2 I`. No analytic spectral
parameter or nonzero-denominator assumption is needed.

## 2. The full ambient three-by-one identity

On four physical strands put

```math
F_L=f_3\otimes I,\quad F_R=I\otimes f_3,\quad P_4=f_4,
\qquad C=B_0B_1B_2F_L,
```

where lower bit positions are the first strands. Let

```math
W(z)=\mathcal R_0(qz)\mathcal R_1(q^3z)\mathcal R_2(z).
```

The three arguments, in application order, are `z,q^-2 z,q^-4 z`.
Then

```math
\boxed{\quad W(z)F_L=p(z)C+r(z)P_4,\quad}        \tag{3}
```

with

```math
p(z)=q^2(1+z)(1+q^3z)(1+q^4z),
\qquad r(z)=q^2(1+q^3)(1+z)(1+q^3z).            \tag{4}
```

Also `F_R W(z)F_L=W(z)F_L` as an ordinary ambient identity.
At q=8 the coefficient lists of p,r, in increasing degree, are
`(12,5,6,15)` and `(13,6,11,0)`.

These are fixed-field TL identities, not just tests of a few values.
The checker constructs fourteen independent images of TL4 diagram
words in the 16-by-16 physical representation. TL4 has its fourteen
noncrossing diagrams as a basis, so this representation is faithful.
An explicit full-rank minor is recorded. Expanding (3) into its four
coefficient matrices verifies every entry. The same checks run at
all four nontrivial fifth roots in F16.

The distinction between ambient and quotient channels is substantive.
The ordinary ranks of F_L,P_4,F_L-P_4 are8,5,3. The partial quantum
trace of P_4 is zero and P_4 kills turnbacks; the
[negligible-wall argument](quantum-braid-path-quotient.md) therefore
kills its ideal under every vacuum boundary pairing. Modulo that
ideal, (3) becomes

```math
W(z)F_L\equiv p(z)C.                            \tag{5}
```

At z=q, p vanishes while r does not: the ambient operator is a
nonzero rank-five multiple of P_4, whereas the vacuum quotient gives
zero. Conversely, using the same unshifted argument z at all three
gates leaks outside F_R. Spectral shifts cannot be omitted.

## 3. Polynomial normalization and arbitrary-length transport

Every ambient entry of W(z)F_L has the common factor
`q^2(1+z)(1+q^3z)`. Cancel it exactly to define an affine polynomial
operator

```math
\overline W(z)=(1+q^4z)C+(1+q^3)P_4.
```

This is a polynomial definition, valid even where the canceled factor
vanishes. In the quotient it is simply `(1+q^4z)C`.

The output-projector identity permits iteration: after one crossing,
the f_3 block is in exactly the location required for the next.
Embedding the faithful local identity in a larger tangle puts each
residual in the same negligible ideal. Thus a block transported
through m strands, with external parameters `q^j z`, is its constant
block braid times

```math
P_m(z)=\prod_{j=0}^{m-1}p(q^jz),
\qquad\overline P_m(z)=\prod_{j=0}^{m-1}(1+q^{j+4}z)            \tag{6}
```

for the raw and normalized versions respectively. This holds in all
vacuum boundary pairings, by local composition, at every length.
No ordinary invariance of the retained physical frame is asserted.

There is an exact sparse simplification. A complete orbit of the
fifth roots gives `prod_j(1+q^j z)=1+z^5`. Hence for `m=5a+r`,

```math
P_m(z)=(1+z^5)^{3a}\prod_{j=0}^{r-1}p(q^jz),\qquad
\overline P_m(z)=(1+z^5)^a\prod_{j=0}^{r-1}(1+q^{j+4}z).
```

In particular when `m=5*2^s`, the normalized scalar is `1+z^m`
and the raw scalar has exactly the four terms at degrees0,m,2m,3m.
This follows from Frobenius in characteristic two. It is a sparse
transport factor, not a decomposition of the entire crossing.

## 4. Constant-braid normalization

The raw constant braid has nontrivial phases. On the retained
three-by-one channel, the double braid is `q^3 I`; the self-braid
of two f_3 objects on their vacuum channel is `q^2 I`.
Normalize **constant** fundamental crossings by
`B_hat=q^2 B=q^(-1/2)B`. Six crossings then multiply the former
phase by q^12 and nine crossings multiply the latter by q^18;
both normalized phases are one.

This is consistent with the
[invertible fusion object](invertible-fusion-graded-seams.md).
The raw q phases must be kept when using the spectral formula (1);
normalizing a constant braid is not an instruction to silently use
the same spectral parameter in a differently normalized Hecke gate.

## 5. Scope and reproduction

The [implementation](../tools/quantum_braid_spectral_fusion.py) includes
the local physical coefficients, exact normalization, scalar strings
and vacuum-path controls for longer transported blocks. The
[checker](../tools/check_quantum_braid_spectral_fusion.py) verifies the
faithful minor, both spectral identities, the full residual in (3),
exceptional zeros and unshifted negative controls. It also compares
nonzero3/5/7-strand vacuum transports coefficientwise and checks the
sparse scalar strings. The
[receipt](../evidence/fused-spectral-transport.json) pins all sources.

```sh
make verify-braid-spectral-fusion
```

The new positive target is a spectrally compatible fused kernel or
grouping identity with polynomially many pieces and a favorable two-cut
budget. Constant prefix projector placements already miss the small
crossings; (3) supplies a local compatibility relation unavailable to
unshifted spectral gates. It does not supply the missing global
spanning identity, strict whole-family endpoint construction or a
stronger complete factoring score. No new Lean proof is claimed.
