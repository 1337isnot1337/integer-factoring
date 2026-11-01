# Inverse residue-block sums at arbitrary precision

## Constructive result

The [quadratic selector algorithm](pole-cauchy-quadratic-selectors.md)
extends from parity to any requested precision `p>=1`. It computes

\[
 B_{r,h;s,g}^{(p)}=
 \sum_{\substack{a\in A,\ a=r\ (2^h)\\
                  b\in B,\ b=s\ (2^g)}}
       (2^L M^{-1})_{b,a}\pmod{2^p}
 \tag{1}
\]

in bit time polynomial in `p`, the numerical residue precisions `h,g`,
the displacement `k`, and the encoded domains and input. Here `L` is
the effective certificate for integrality of the normalized inverse.
The matrix nodes are never enumerated.

The algorithm uses `p^2` quadratic-pole contractions, each with numerical
pole orders at most `6p-3`. A polynomial selector lift replaces any need
to expand a power of degree `2^p`. This supplies a full-precision
residue-block interface, including precisions exceeding `L`.

Computing the Smith multiplicities of an exponentially large matrix
from this interface remains unresolved. There is no new factoring
exponent claim.

## 1. A polynomial lift of a binary residue

For `p>=1`, define

\[
 P_p(X)=X^p\sum_{j=0}^{p-1}
                  {p+j-1\choose j}(1-X)^j.
 \tag{2}
\]

This is an integer polynomial of degree `2p-1`. It is divisible by
`X^p`, and

\[
 P_p(X)\equiv1\pmod{(1-X)^p}.
 \tag{3}
\]

For (3), expand `X^(-p)=(1-(1-X))^(-p)` as a formal power series in
`1-X` and retain its first `p` terms. Multiplication by `X^p` then
gives one modulo `(1-X)^p`. This is an identity of integer polynomials;
no input-dependent inversion is required.

Consequently, for an integral quadratic number `x` whose residue
modulo two is either zero or one,

\[
 P_p(x)\equiv
 \begin{cases}
 0\pmod{2^p},&x\equiv0\pmod2,\\
 1\pmod{2^p},&x\equiv1\pmod2.
 \end{cases}
 \tag{4}
\]

The first case follows from divisibility by `x^p`; the second follows
from divisibility of `1-P_p(x)` by `(1-x)^p`. This argument works in
the full unramified quadratic integer ring, not only for rational
integers.

The only nonzero powers of `P_p` have degrees from `p` to `2p-1`.
Its coefficients can be written

\[
 P_p(X)=\sum_{j=0}^{p-1}\gamma_{p+j}X^{p+j},\qquad
 \gamma_{p+j}=(-1)^j
       {2p-1\choose p+j}{p+j-1\choose j}.
 \tag{5}
\]

Expanding (2) and applying the binomial summation identity gives (5).
Each coefficient has `O(p)` bits. The implementation generates these
`p` coefficients directly. Equivalently, this is the majority
polynomial `sum_(j=p)^(2p-1) binom(2p-1,j) X^j(1-X)^(2p-1-j)`.

## 2. Lifting the quadratic residue selector

The earlier selector is

\[
 \chi_{r,h}(a)=
      \left(\frac{2^h}{a-r+2^h w}\right)^3.
 \tag{6}
\]

It is integral on every integer and has residue zero or one. Its residue
is one exactly on `a=r mod 2^h`. Apply (4) to obtain

\[
 P_p(\chi_{r,h}(a))
    \equiv 1_{a=r\ (2^h)}\pmod{2^p}
 \tag{7}
\]

uniformly over the entire matrix domain.

Put `c=r-2^h w`, using the canonical representative `0<=r<2^h`.
Equations (5)--(6) express this lifted weight as

\[
 P_p(\chi_{r,h}(a))
    =\sum_{j=p}^{2p-1}
       \gamma_j\frac{2^{3hj}}{(a-c)^{3j}}.
 \tag{8}
\]

Each summand before its integer coefficient is the integral quantity
`chi_(r,h)(a)^j`. The corresponding column expansion uses
`d=s-2^g w`.

Because `2^L M^-1` and both lifted weights are integral, replacing the
two residue indicators in (1) by the weights in (7) changes the full sum
by a multiple of `2^p`. The number of entries causes no precision loss.

## 3. Evaluating the block sum

Let `C_(u,v)(c,d)` denote the complete repeated-pole contraction from
the quadratic-jet algorithm. Expanding (8) in both variables gives

\[
 B_{r,h;s,g}^{(p)}\equiv
 \sum_{j=p}^{2p-1}\sum_{\ell=p}^{2p-1}
     \gamma_j\gamma_\ell\,
     2^{L+3hj+3g\ell}
     C_{3j,3\ell}(c,d)
     \pmod{2^p}.
 \tag{9}
\]

There are `p^2` terms, and their largest pole order is `6p-3`. The
quadratic poles satisfy all nonvanishing conditions proved in the
earlier note. For each term, the exact jet algorithm already charges
the precision needed to extract its indicated coefficient.

The [implementation](../tools/pole_cauchy_residue_blocks.py) requests
absolute precision `p` for each contraction in (9). Multiplication by
the displayed nonnegative power of two can only improve that absolute
precision. It then extracts an integral residue modulo `2^p` and forms
the integer linear combination.

The integrality of every term follows directly from

\[
 2^{L+3hj+3g\ell}C_{3j,3\ell}(c,d)
     =\sum_{a,b}(2^L M^{-1})_{b,a}
                      \chi_{r,h}(a)^j\chi_{s,g}(b)^\ell.
\]

If a returned scale is still negative after this normalization, the
program verifies exact divisibility of both quadratic residue
coordinates before dividing. The available residue precision exceeds
the desired `p` bits by the corresponding normalization exponent, so
this division introduces no ambiguity.

All coefficient sizes, pole orders, and required numerical precisions
are polynomial in the stated parameters. The `p^2` calls and their
integer linear combination therefore have polynomial bit cost. Their
cost depends on the encoded interval lengths, not on the numerical
matrix cardinality.

## 4. Scope of higher precision

Modulo two, the normalized inverse can lose all original unit Smith
directions whenever `L>0`. That observation in the earlier note was
specifically about parity. At precision `L+1`, contributions from those
directions need not disappear.

For example, take `N=15`, `A=B={1,3,5}`, and `k=5`. The certificate has
`L=4`. The normalized inverse singleton blocks at `(3,5)` and `(5,3)`
vanish modulo two, but both are `16 mod 32` at five-bit precision. The
new interface computes those higher bits directly.

This is access to entries and structured block sums, not an algorithm
that identifies a Smith basis or counts its unit directions. Arbitrary
matrices can require many such queries, and no sufficiently compressed
rank or Smith-count algorithm has been proved for this Cauchy family.
In particular, a high bit of an individual inverse entry is not asserted
to be a factor-pair test.

## 5. Exact controls

The [checker](../tools/check_pole_cauchy_residue_blocks.py) verifies the
polynomial coefficient identities and both high-order congruences
through precision sixteen. It evaluates the lifts at positive and
negative integers and at quadratic integers with residue zero or one.

Complete block sums are compared with dense exact inversion at one,
two, and three bits. The one-bit outputs are compared with the earlier
parity interface. Three five-bit singleton controls verify the example
from Section 4, including both original factor-pair positions. The
[evidence](../evidence/pole-cauchy-residue-blocks.json) records the results
and source hashes.

```sh
python3 tools/check_pole_cauchy_residue_blocks.py \
  --output evidence/pole-cauchy-residue-blocks.json
```

The large implicit-domain validation of the underlying quadratic-jet
evaluator remains in the preceding bundle. This note introduces the
finite polynomial lift and its explicit precision composition.
