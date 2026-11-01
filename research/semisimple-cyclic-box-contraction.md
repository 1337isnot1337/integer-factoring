# Zero-safe box contraction in the cyclic recovery algebra

## Result

The ordinary cyclic-polynomial recovery algebras already support the
zero-safe Frobenius box circuits. No everywhere-unit weighting or
enumeration of points in a large selector field is required for this
local contraction interface.

Let `E=F_(2^d)` be explicitly represented, let `r` be an odd prime, and
put

\[
 A=E\otimes_{F_2}F_2[X]/(X^r-1),\qquad L=d(r-1).          \tag{1}
\]

The algebra has binary dimension `dr`, may have zero divisors, and
satisfies

\[
 z^{2^L}=z\quad\hbox{for every }z\in A.                  \tag{2}
\]

Consequently the positive-power operation

\[
 \iota_L(z)=z^{2^{L+1}-3}                                \tag{3}
\]

is the inverse on each nonzero field component and zero on each zero
component. For every positive integer `w` and every `z in A`,

\[
 \sum_{j=0}^{2^w-1}z^j
       =(1+z^{2^w})\,\iota_L(1+z).                      \tag{4}
\]

This includes zero-divisor denominators and `z=1`. Width-zero boxes
are the constant one and are handled separately.

## 1. Proof without constructing component fields

Frobenius to exponent `L` fixes `E`, since `d|L`. It fixes the cyclic
generator because `2^L=1 mod r`, by Fermat's theorem and `r-1|L`.
It therefore fixes every element of (1), proving (2).

The polynomial `Y^(2^L)-Y` is squarefree, so a finite algebra satisfying
(2) is reduced and is a product of finite fields whose degrees divide
`L`. This gives the component interpretation of (3), but no component
factorization is needed by the algorithm.

There is also a direct proof of (4). Write `x=1+z` and `B=2^w>=2`.
In characteristic two,

\[
 \sum_{j=0}^{B-1}z^j=x^{B-1},\qquad 1+z^B=x^B.
\]

Equation (2) allows a positive exponent of `x` to be reduced by
`2^L-1`. The right side of (4) is
`x^(B+2^(L+1)-3)`, whose exponent exceeds `B-1` by exactly
`2(2^L-1)`. Both reductions are valid because `B-1>=1`. Thus (4)
holds directly in the algebra, including its zero divisors.

## 2. A prescribed scalar subfield and the circuit hypothesis

When a scalar summation field `F_(2^t)` is needed, use

\[
 A'=E\otimes F_{2^t}\otimes F_2[X]/(X^r-1),\qquad
 L'=dt(r-1).                                             \tag{5}
\]

The same proof fixes all three sets of generators. The binary dimension
is `dtr`, and `L'<=dtr`. Constructing `F_(2^t)`, if it is not already
available, can be charged by the coarse deterministic bound `2^O(t)`.
There is no search through a compositum of degree `L'`.

The algebra in (5) need not be an algebra over the single field
`F_(2^(L'))`. The sufficient hypotheses for the
[zero-safe binomial circuit](frobenius-trace-infinity-correction.md),
Section 4, are instead stated directly: a commutative characteristic-two
algebra with an embedded scalar `F_(2^t)`, `t|L'`, and Frobenius to
exponent `L'` equal to the identity. Its proof is unchanged under these
hypotheses. Expand (3) into its `L'` Frobenius factors, reduce only the
scalar summation-variable exponents modulo `t`, and use the fixed
binary-carry circuit. The same reasoning applies to the
[relative-linearized product algorithm](relative-frobenius-trace-contraction.md).

For the [Fermat-torus star](fermat-torus-binomial-contraction.md), the
root field `E` already supplies the scalar roots. No extra scalar
subfield is needed. Positive-power expansion followed by its
division-free coefficient circuit remains valid over (1).

Each represented algebra operation and the positive-power Frobenius
expansion cost a polynomial in the represented dimension. The binomial
and torus coefficient circuits also have polynomial size in that dimension
and the original factor count. The relative-linearized algorithm retains
its separate `(s+2)^{O(J)}` operation count; this compatibility statement
does not make that count polynomial when `J` varies. In particular, `L'`
is at most the binary dimension, rather than an uncharged exponentially
larger field parameter.

## 3. Compatibility with rough-input recovery

The [sparse cyclic recovery](sparse-divisor-nonzeroness.md) deliberately
uses odd primes `r`, so its rings have exactly the required form.
Under the [rough-observable pretest](rough-observable-factor-constructor.md),
the support bound is `exp(O(P))`. The cyclic prime list and its maximum
prime have size `exp(O(P+log n))`. Once root and scalar fields of binary
degrees polynomial in `n` are represented, constructing the cyclic
extension and performing each represented algebra operation have cost
`exp(O(P+log n))`. The number of local contraction operations remains
separately charged, including the relative-linearized `J` dependence.
Constructing
an additional scalar field retains the separate `2^O(t)` charge from
Section 2. In the adaptive prime schedule, `t=O_epsilon(log n)` and the
[prescribed-subfield root construction](prime-root-field-construction.md)
also has polynomial cost, so preparation fits the same bound there.

Thus these local zero-safe box contractions can be used directly in the
ordinary cyclic recovery family. The separate
[Frobenius-compatible field-point construction](frobenius-compatible-sparse-recovery.md)
still supplies its stronger everywhere-unit properties, but those are
not required for (4) or for the stated positive-power circuits.

This result does not close the full coupled multiplication contraction.
Representation growth under a second elimination, and variable-dependent
exponents from multiplication gates, remain separately charged.

## 4. Why the odd-prime hypothesis matters

In `F_2[X]/(X^2-1)`, the nonzero element `epsilon=1+X` is nilpotent.
Taking `z=X` gives `S_2(z)=epsilon`, whereas `1+z^2=0`. Hence a
zero-safe quotient of the form (4) would give zero incorrectly. Such a
nilpotent algebra does not satisfy (2). The theorem does not authorize
that substitution in arbitrary coefficient algebras.

## 5. Exact replay

The companion checker is
[`check_semisimple_cyclic_box_contraction.py`](../tools/check_semisimple_cyclic_box_contraction.py).
It constructs the tensor/cyclic algebras directly, checks the Frobenius
period, zero-divisor identities, complete scalar-subfield sums, and the
nilpotent countercontrol. Run:

```sh
python tools/check_semisimple_cyclic_box_contraction.py --output evidence/semisimple-cyclic-box-contraction.json
```

The committed [evidence](../evidence/semisimple-cyclic-box-contraction.json)
records 54 Frobenius-period and reflexive-inverse checks, 162 all-point
geometric identities, and 18 complete algebra traces checked against
72 scalar-field evaluations. Six models contain explicit nonzero zero
divisors, and the even-cyclic nilpotent countercontrol fails exactly the
excluded identity. All asserted comparisons passed.
