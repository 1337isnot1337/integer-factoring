# Regular joint kernels inside inverse residue blocks

## Constructive result

Let `H=2^L M^-1` be the integral normalized inverse supplied by the
[domain precision certificate](pole-cauchy-domain-precision.md). Let
`K(a,b)` be a regular odd-input arithmetic circuit in two variables,
with the gate contract of the
[joint analytic compiler](pole-cauchy-analytic-kernels.md).

For integer residue classes and an external dyadic denominator `2^r`,
the algorithm computes

\[
 2^{-r}\sum_{\substack{a\in A,\ a=\alpha\ (2^h)\\
                       b\in B,\ b=\beta\ (2^g)}}
       H_{b,a}K(a,b)
 \tag{1}
\]

to absolute precision `p`, with polynomial bit cost in the circuit
description, `p+r`, the numerical residue precisions, and the encoded
matrix input. It does not enumerate matrix nodes. Binary-encoded powers
retain the compiler's logarithmic dependence on their exponents.

This combines the regular joint circuit contract with the nonregular
quadratic residue selectors. It does not accept the singular joint
denominator `N-ab`, and it supplies no Smith-rank or factor-count oracle.

## 1. Pointwise compilation and selection

Set `q=p+r`. The regular circuit compiler produces an integer polynomial

\[
 K_q(a,b)=\sum_{d,e}c_{d,e}a^d b^e
 \tag{2}
\]

of total degree less than `q`, agreeing with `K` modulo `2^q` at every
odd pair. The [Hermite lift](pole-cauchy-residue-blocks.md) supplies
integer-valued rational selectors `S_(alpha,h,q)(a)` and
`S_(beta,g,q)(b)` agreeing with the two indicators modulo `2^q`.

All these values and `H` are integral. Replacing the kernel and the
indicators by their approximations changes the sum by a multiple of
`2^q`. Dividing by `2^r` leaves absolute precision `p`. No dimension
factor is needed in this error bound.

## 2. Polynomial times selector is a proper pole polynomial

Use the quadratic pole `c=alpha-2^h w`. The selector has the form

\[
 S_{\alpha,h,q}(a)=
   \sum_{j=q}^{2q-1}\gamma_j
        \frac{2^{3hj}}{(a-c)^{3j}}.
 \tag{3}
\]

For a monomial of degree `d<q`, expand

\[
 a^d=(c+(a-c))^d
     =\sum_{\ell=0}^d{d\choose\ell}c^{d-\ell}(a-c)^\ell.
 \tag{4}
\]

Every resulting denominator order is `3j-ell`, which is positive.
Its largest possible value is `6q-3`. Thus a polynomial of degree below
`q`, multiplied by (3), is an explicitly written proper rational
polynomial at the single pole `c`.

All coefficient pairs are integers in the basis `1,w`. Their bit lengths
are polynomial in `q`, `h`, and the original coefficient lengths.
The [adapter](../tools/pole_cauchy_restricted_kernels.py) constructs these
coefficients by the binomial identity (4).

Group (2) by row degree:

\[
 K_q(a,b)=\sum_d a^d V_d(b).
\]

For each nonzero row, construct the two proper pole polynomials for
`a^d S_(alpha,h,q)(a)` and `V_d(b)S_(beta,g,q)(b)`. Multiply the row
coefficients by `2^L`, then pass the two arrays to the batched
pole-polynomial evaluator. There are at most `q` such contractions.
Each uses one pair of product jets, regardless of its number of
nonzero coefficients.

The sum before the external division is integral. The implementation
extracts it modulo `2^q` using each batched result's explicit scale,
checking exact divisibility whenever that scale is negative. The final
return scale is `-r`, and its residue precision is `q`, so its absolute
precision is exactly `p`.

## 3. Complexity and scope

The compiler has polynomial cost in its stated regular circuit contract.
Its output has polynomially many coefficients and total degree below
`q`. Equations (3)--(4) construct polynomial-size coefficient arrays,
with numerical pole orders at most `6q-3`. At most `q` batched inverse
contractions are needed. Their exact valuations, coefficient scales,
and short-jet precision are all charged by that evaluator.

The normalization `2^L` is placed in the row coefficients before the
batched precision is selected. This preserves the evaluator's tighter
coefficient-valuation budget. It does not force computation at a
precision enlarged by `L` when the normalized answer needs fewer bits.

The method handles a joint regular numerator with separately imposed
row and column residue conditions. It does not replace a genuinely
singular joint denominator by an uncharged polynomial expansion.

## 4. Exact controls

The [checker](../tools/check_pole_cauchy_restricted_kernels.py) verifies
240 exact polynomial-times-selector identities in the quadratic field.
Thirty complete restricted joint-kernel contractions are compared with
dense exact inverse matrices, including external dyadic divisions.
Six additional singleton controls force nonzero answers at five-bit
absolute precision. They use `N=15`, the kernel `a+2b`, and the positions
`(1,1)`, `(3,5)`, and `(5,3)`, both with and without external division
by two. This checks the scaled output on actual factor positions and on
a nonfactor position, without relying on empty or vanishing blocks.

The implicit example has a 402-bit input and dimension `2^200+13`.
Its joint kernel contains a signed 513-bit exponent. The selected
residue classes are `5 mod 128` and `11 mod 512`. Two-bit and four-bit
computations agree without listing any nodes. The lower-precision
record uses one batched contraction and completes in about 0.11 seconds
in the recorded environment. Timing is a finite control, not an
asymptotic factoring claim.

The [evidence](../evidence/pole-cauchy-restricted-kernels.json) records
the compiler output, contraction budgets, and source hashes.

```sh
python3 tools/check_pole_cauchy_restricted_kernels.py \
  --output evidence/pole-cauchy-restricted-kernels.json
```

The established factoring score remains unchanged.
