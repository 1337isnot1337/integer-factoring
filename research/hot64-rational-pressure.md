# An exact sixty-four-letter natural pressure bound

## Result

A fixed sixty-four-letter stabilized reflection restriction admits the
uniform balanced crossing bound

```math
\rho<\frac{32}{45}=0.711111\ldots.                    \tag{1}
```

The complete profiles and strict pressure comparison are exact. This
improves the earlier finite certificate `rho<5/7`, but remains above
the factoring-transfer threshold `1/sqrt(2)`. It establishes no new
factoring exponent or score improvement.

## 1. The fixed code and its complete profiles

Order hot masks by `(E(mask),mask)`, where
`E(S)=sum_(a in S)(8a+4)`, and retain the first sixty-four. The verifier
enumerates twelve hot mode positions; all retained positions fit in
eleven. The maximum energy is 84. The interval `M=32` fits the code
and satisfies `84<4M`, so the
[stabilization theorem](stable-hot-reflection-codes.md) applies.
The code is a deletion downset containing the vacuum, and is used on
all four physical legs.

The [occupation-block evaluator](hot-code-charge-profile-constructor.md)
computes all 4,096 common natural Smith orders. The largest block has
order 1,238. Every mask, energy, block, and order appears in the
[evidence](../evidence/hot64-charge-rational-pressure.json).

By [natural-cut duality](reflection-natural-cut-duality.md), the two
actual primitive natural profiles are `a+sigma` and `-a+sigma` for
one scalar shift `a`. Thus their channel-sum product is

```math
Z_0(u)Z_\infty(u)=\left(\sum_i e^{-u\sigma_i}\right)^2. \tag{2}
```

This identity already accounts for the original reciprocal tensor
and its degree contribution. It does not remove a degree by changing
a single flattening.

## 2. Strict integer comparison

Set

```math
z=\frac{459}{460},\quad A=64,\quad B=45,\quad R=A/B,
\quad u=A\log(1/z),
\quad P=\sum_i z^{A\sigma_i},\quad Q=\sum_{x\in C}z^{BE(x)}.
```

The verifier proves

```math
P^{2B}<Q^{2A}.                                       \tag{3}
```

Both sums are evaluated by homogeneous integer arithmetic. After
raising their numerators to the indicated powers and clearing powers
of 460, the resulting two positive integers each have 14,266,197 bits.
Their exact comparison is strict. The evidence stores their bit
lengths and SHA-256 digests, and the digest of their positive difference.
No floating logarithm or approximate root decides (3).

The squared comparison also supports unequal paired codes by replacing
`Q^2` with `Q_X Q_Y`; the current certificate uses equal codes.

## 3. Entropy and uniform crossing cost

Choose the strictly positive rational type

```math
p(x)=z^{BE(x)}/Q.
```

Multiplication by `460^(B max E)` gives positive integer counts.
Their two defining exponents and their exact sum are retained in the
evidence, so a block length clearing all type denominators is effective.
Let `H` be this type's entropy and `E_p` its mean energy. Then

```math
H=\log Q+(u/R)E_p.
```

The opposite scalar shifts in (2) cancel, while the total natural
physical mean remains `2uE_p`. The complete logarithmic rank cost is

```math
K=2\log P+2uE_p<2R\log Q+2uE_p=2RH.
```

Consequently `K/(4H)<R/2=32/45`. The denominator is `4H`: it includes
both paired type entropies and the crossing normalization factor two.
The coefficient prefactors and number of interpolation pieces are
polynomial in the tensor power, and the strict margin absorbs them.
Alternatively, a sufficiently large admissible type block passes the
[finite-block amplification](finite-block-crossing-amplification.md)
certificate for any slightly weaker rational bound.

This uses one fixed finite code and one fixed finite field. The
underlying literal reflection restriction remains an effective strict
crossing. The separate factoring threshold is not passed by (1).

## 4. Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python tools/certify_hot_charge_pressure.py \
  --size 64 --modes 12 --base 459 460 --support 64 45 \
  --output evidence/hot64-charge-rational-pressure.json
```

The [verifier](../tools/certify_hot_charge_pressure.py) reconstructs the
code and every complete occupation-block profile before comparing the
integers. The recorded run took approximately 151 seconds. Source
hashes bind the constructor, unchanged determinant and DVR kernels,
and the literal reflection restriction interface.

The profile engine has independent complete-polynomial, quotient-rank,
endpoint, and retained-profile controls. The certificate does not use
the experimentally observed one-particle odd-order pattern as an
unproved replacement for any complete profile.
