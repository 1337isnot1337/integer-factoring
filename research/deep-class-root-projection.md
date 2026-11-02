# Always-positive deep class roots suffice for factoring

## Question and conditional result

Can the [high-root certificate contributed in Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360821)
remove relation acquisition itself, rather than only replace odd row
annihilators? Yes, provided a root algorithm handles a specific family
of guaranteed-positive inputs. The [scope](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19362022)
was registered at base `7217aab` before this deduction.

Under the existing [prime-form generation import](polylog-prime-form-generators.md),
an engine costing F(b) on the inputs specified below gives a complete
deterministic smallest-prime-factor algorithm costing

```math
\operatorname{poly}(b)(F(b)+1).                         \tag{1}
```

Here b is the input bit length and F is a uniform nondecreasing bound
for the associated queries, whose discriminants have b+O(1) bits. The
engine is not supplied. No stronger complete proven bound follows;
the score remains **1/2**. This conditional reduction needs no smooth-set
density, bounded smooth words or relation
matrix. It retains the existing external analytic premise for generation.

## 1. An exact projection from any sufficiently deep root

Let G be a finite abelian group with `|G|<=2^L`. For any g in G, form

```math
h=g^{2^L}.
```

The two-primary component of h is the identity. Thus h has odd order,
and powering by any power of two is invertible on its odd-order subgroup.
In particular there exists a z in G with

```math
z^{2^{2L}}=h.                                         \tag{2}
```

Accept **any** z satisfying (2), in the same ambient G, and put

```math
o=z^{2^L},\qquad w=g o^{-1}.                           \tag{3}
```

The element o has odd order because the size bound kills its two-primary
component. Also `o^(2^L)=g^(2^L)`, so `w^(2^L)=1` and w has two-power
order. Equation (3) therefore decomposes g into its unique commuting
two-primary and odd-primary components: w is exactly its two-primary
projection. It is independent of which root z the engine returns.
For a generating list `g_1,...,g_m`, the resulting w_i generate all of G_2.

This specializes the contributed adjugate construction one generator at
a time, with one-row matrix `[2^L]`, determinant valuation L and root
precision `B=2L`. It satisfies `B>=L+v_2(det R)`. Using the tautological
matrix `[2^B]` would not satisfy that inequality and is not the argument.
The [general certificate note](high-root-projection.md#guaranteed-positive-deep-roots-remove-row-selection)
now records and credits this specialization. The all-input engine contract
and complete polynomial invoice are given here.

Each power in this proof is acquired by repeated squaring. Constructing h,
verifying (2) and computing (3) costs O(L) group operations per generator.
The parameter L is charged by value, even though its binary encoding is
short. The powers are evaluated without expanding their values.

## 2. The arithmetic engine contract and all-input transfer

For odd input N, let `Delta=-4N` and `L=bits(4N)`. The existing elementary
reduced-form count gives `|Cl(Delta)|<=4N<2^L`. Canonical reduced forms
have O(b)-bit coefficients and deterministic polynomial-bit-cost
composition, inversion and equality.

The [generation theorem](polylog-prime-form-generators.md) constructs a
polynomial number of prime forms g_i generating this full proper class
group, under its pinned finite-order Hecke zero-free-strip input. Its
uniformity for nonmaximal orders and effective constant are retained.
No factorization of N, the conductor or the class number is computed.

The required engine receives Delta, L and an acquired generator g_i,
with `h_i=g_i^(2^L)` computed explicitly as a reduced form. It must return
a reduced form z_i of the same discriminant satisfying (2). Its uniform
bit cost on these inputs, including output acquisition, is at most F(b).
Every query is positive by Section 1. Neither negative instances nor
roots of arbitrary mixed-order classes are requested. The known shallower
preimage g_i is available; it is not itself the required deeper root.

Verify every output by composition, form w_i by (3), then apply the
two-group phase of [early-stop factor extraction](early-stop-class-factoring.md#3-early-stop-insertion-with-one-cyclic-generator).
That phase needs generators of the full two-part, which the w_i supply;
it does not require them to generate the discarded odd part. The explicit
annihilator `2^L` kills all w_i and has L+1 bits. The phase costs
`O(m L^2)` group operations. Its classification of the identity and the
single possible unhelpful nonidentity ambiguous class applies unchanged.
For N with at least two distinct prime factors it finds a proper divisor.
There is no exponential rank search or preliminary roughness requirement.

Even inputs, primes and perfect powers have the existing polynomial-cost
preprocessing. For other odd inputs, including nonsquarefree ones, the
preceding procedure finds a proper divisor. Recurse on both factors to
obtain the smallest prime factor. There are at most O(b) recursive nodes,
all with bit length at most b. Monotonicity of F, polynomial generation,
root verification, reduction, integer division and primality/root tests
give (1). Thus any fixed logarithmic-saving exponent of a proved engine
would be preserved for the complete problem. This is a conditional
algorithm, not an implementation of the engine or a score improvement.

## 3. Why the inspected square-root routines do not fill this contract

Repeated arbitrary square-root calls do not preserve positivity at the
next call. Already in additive `Z/2`, zero is a fourth power, but a valid
first square root is one, which has no square root. The single deep-root
contract cannot be replaced by an ordinary positive-instance square-root
oracle without a further branch-selection argument.

Bosma and Stevenhagen,
[On the computation of quadratic 2-class groups](https://www.numdam.org/item/JTNB_1996__8_2_283_0.pdf),
pp. 306–309, construct class halves using modular square-root data at
the known prime-power factors of the discriminant. Their stated algorithm
takes that factorization as input. This inspection supplies no algorithm
for our unfactored discriminant and does not exclude other root methods.

The repository's [known-root conic parametrization](class-square-conic.md)
already constructs every rational conic point without factoring Delta.
Its [denominator audit](conic-denominator-descent.md) tracks the original
root class and its possible unhelpful twist unless a displayed gcd already
factors N. Those results on unrestricted conic solubility and the known
shallower root do not supply the deeper class root in (2).
No new search or larger collision-frontier experiment is inferred here.

The higher-dimensional action identity in
[PEGASIS](https://r98inver.github.io/papers/2503_pegasis.pdf), Section 3.1,
starts with ideals b,c already in a supplied class a. Its known
composition uses the principal product `c conjugate(b)` to evaluate a's
oriented action. It does not output a class whose large two-power is a.
Thus this identity does not fill (2). This is an input/output distinction,
not an obstruction to other orientation-preserving constructions; no
implementation or norm-solver performance is imported here.

This ordinary deduction specializes the contributed high-root formula
and retains the predecessor credit for generation, group arithmetic and
early-stop extraction. It was developed by Sol worker 2 with substantial
OpenAI Codex assistance. No new Lean, external-proof rebuild or independent
mathematical review is claimed. The next required step is an arithmetic
root-acquisition mechanism with the exact promise and costs above.
