# Binary product powers with explicit dyadic normalization

## Constructive result

The inverse-contraction engine now accepts weights

\[
 f(a)=\prod_{i=1}^u P_i(a)^{e_i},\qquad
 g(b)=\prod_{j=1}^v Q_j(b)^{d_j},                 \tag{1}
\]

where each base polynomial is explicitly written over the integers and
every exponent is a signed binary integer. The bases must be nonzero on
their respective matrix domains. Their values may be even everywhere.
The total expanded numerator or denominator degree can be exponential
in the input length.

The result is a **scaled** answer. Let

\[
 \mu=\min_{a\in A}v_2(f(a)),\qquad
 \nu=\min_{b\in B}v_2(g(b)),\qquad H=2^L M^{-1},   \tag{2}
\]

where the existing domain certificate makes every entry of `H` integral.
The algorithm computes both minima exactly as binary integers and returns

\[
 C=\sum_{a\in A,b\in B}H_{b,a}
       (2^{-\mu}f(a))(2^{-\nu}g(b))\pmod{2^p}.    \tag{3}
\]

It therefore represents the original contraction by

\[
 \sum_{a,b}(M^{-1})_{b,a}f(a)g(b)
   =2^{\sigma}\bigl(C+O(2^p)\bigr),\qquad
 \sigma=\mu+\nu-L.                              \tag{4}
\]

The running time is polynomial in the explicit base degrees and coefficient
lengths, binary exponent lengths, encoded matrix input, and numerical `p`.
The matrix dimensions and expanded powers are not enumerated. The numerical
value of a potentially enormous `sigma` is not a working-precision parameter.

This extends [regular analytic powers](pole-cauchy-analytic-weights.md),
which require unit bases, and [explicit rational weights](pole-cauchy-rational-weights.md),
which charge expanded polynomial degrees. It does not construct the singular
joint factor-count contraction. The complete factoring exponent remains `1/2`.

Implementation: [pole_cauchy_product_power_weights.py](../tools/pole_cauchy_product_power_weights.py).
Checker: [check_pole_cauchy_product_power_weights.py](../tools/check_pole_cauchy_product_power_weights.py).
Evidence: [pole-cauchy-product-power-weights.json](../evidence/pole-cauchy-product-power-weights.json).
The dependency snapshot is commit `c404b544a538454b9cca8b810ebb224b2c74a0ec`.

## 1. A common valuation partition charges bases once

Apply the [univariate valuation partition](pole-cauchy-rational-weights.md)
to every nonconstant or constant nonzero base. For a polynomial `P_i` of
degree `D_i`, put

\[
 H_i=\mathrm{bits}\left(\sum_r |P_{i,r}|R^r\right),
\]

where `R` is the largest node on the domain. Its partition has
`O(max(1,D_i)(H_i+1))` leaves, and it rejects an actual domain zero.
Every leaf has an exact identity

\[
 P_i(\alpha+2^h Z)=2^{s_i}U_i(Z),\qquad
 U_i(0)\text{ odd},\qquad 2^r\mid[Z^r]U_i.       \tag{5}
\]

The exponents `e_i` are not used to construct these partitions.

Dyadic residue classes are nested or disjoint. Their common refinement
therefore has at most the **sum** of the individual leaf counts. One can
construct it by retaining the most refined classes in the union of all
leaf lists. Every node belongs to a most refined containing leaf, and no
two retained classes overlap. Domain counts verify coverage without listing
nodes.

On a refined class, an older local variable changes by
`Z_i=c+2^r Z` with integers `c,r>=0`. Substitution preserves the coefficient
divisibility in (5), and its constant remains odd. Thus every base has (5)
on every common leaf. The product's exact valuation there is the binary
integer

\[
 v_{\alpha,h}=\sum_i e_i s_i.                    \tag{6}
\]

Taking the minimum of this polynomial-size list gives `mu`. Its maximum
is also available. The bit lengths of these integers are polynomial in
the exponent encodings and the base-value bounds, even when their numerical
values are enormous. The empty product has valuation zero.

## 2. Binary powers of local units

Set `delta=v_(alpha,h)-mu>=0`. On that class the normalized weight is

\[
 F(\alpha+2^hZ):=2^{-\mu}f(\alpha+2^hZ)
      =2^\delta\prod_i U_i(Z)^{e_i}.             \tag{7}
\]

If `delta>=p`, omit the leaf before allocating a power of two. Otherwise
put `q=p-delta` and compute the unit product in

\[
 (\mathbb Z/2^q)[Z]/(Z^q).                       \tag{8}
\]

Inversion uses its odd constant coefficient. Positive and negative powers
use binary powering, so the number of multiplications depends on
`log(1+|e_i|)`. Sums, products and unit inverses preserve coefficient
divisibility by `2^r` in degree `r`. Every discarded degree is therefore
zero modulo `2^q` at integral `Z`.

Multiplying the retained jet by `2^delta` gives an integral polynomial
`C_(alpha,h)(Z)` of degree below `p` that agrees with `F` on that leaf
modulo `2^p`. Its degree-`r` coefficient remains divisible by `2^r`.
Only gaps below `p` are used as literal shift sizes.

## 3. Convert each local jet to an integral global weight

Use the existing quadratic selector, with `w^2+w+1=0`:

\[
 c=\alpha-2^hw,\qquad y=\frac{2^h}{a-c}.
\]

The value `y` is integral at every integer node. It is a unit precisely
on the selected class, and `y^3` is one modulo two there and zero elsewhere.
Let `P_p` be the degree-`2p-1` Hermite selector. The leaf contribution is

\[
 \Phi_{\alpha,h}(a)=P_p(y^3)C_{\alpha,h}(y^{-1}-w). \tag{9}
\]

After expansion, every term has the form `integer_pair * y^t`, with
`2p+1<=t<=6p-3`. Thus (9) is a proper rational polynomial at the single
quadratic pole `c`, with integral coefficient pairs and no external dyadic
denominator. On the selected leaf it agrees with `F` modulo `2^p`.
Off the leaf it vanishes modulo `2^p`, because every displayed positive
power of `y` has sufficient valuation.

The sum over the retained leaves therefore agrees with the entire
normalized weight at every domain node. Omitted leaves have normalized
valuation at least `p`. Summing the errors does not lose 2-adic precision.

Compile the column weight in the same way. For each leaf pair, pass the
two explicit pole arrays to the [batched two-jet evaluator](pole-cauchy-pole-polynomials.md),
placing `2^L` in the row coefficients. The matrix and both normalized
weights are integral, so replacing either weight by its pointwise
approximation preserves the complete contraction modulo `2^p`.
This proves (3).

## 4. Precision and representation boundary

There are polynomially many leaf pairs. Each pole array has order at most
`6p-3`, and each coefficient has polynomial bit length in the leaf depth
and `p`. The batched evaluator's conditioning theorem bounds its internal
precision polynomially in those quantities and the encoded domains.
Neither the numerical exponent values nor `|mu|,|nu|` enlarge this precision.
They only occur in binary valuation arithmetic or binary powering.

Equation (4) is not an absolute-precision claim that suppresses a huge
division. To obtain absolute precision `p0` after the external scale, one
needs `p>=p0-sigma`. If `sigma` is a very large negative integer, that is
a correspondingly large numerical precision request. The normalized
contraction (3) is the polynomial-cost object established here. Cancellation
can increase its actual valuation, so `sigma` is a guaranteed lower scale,
not an asserted exact valuation.

The list of base polynomials is explicit. An exponentially long interval
product does not become an admissible short input by calling each of its
distinct linear factors a base. Likewise the joint denominator `N-ab` is
not separated into a row base and a column base by this theorem.

One consequence is inexpensive high-multiplicity valuation filtering.
For a nonvanishing integer polynomial `P`, let `u=max_a v2(P(a))`.
The normalized weight `2^(uE)/P(a)^E`, with `E=2^(p-1)`, is the indicator
of the maximum-valuation leaves modulo `2^p`. Indeed `E>=p` for `p>=1`,
so smaller valuations are suppressed, while every odd unit to the power
`E` is one modulo `2^p`. The exponent is handled in binary; the base
partition is constructed only once. This selects explicit polynomial
valuation strata and does not identify unknown factor pairs.

## 5. Exact controls

Run:

```sh
python3 tools/check_pole_cauchy_product_power_weights.py \
  --output evidence/pole-cauchy-product-power-weights.json
```

The checker verifies 80 common partitions and 858 literal factor-leaf
identities. It compares 447 compiled rational values with independent scalar
modular powering, and 30 full contractions with dense exact inverse matrices.
Twenty-five small signed-power cases additionally check the original
unnormalized rational sum and the returned binary scale. Another 1,056
controls check the maximum-valuation selector. Five invalid contracts,
including zero bases on the domain, are rejected.

The implicit example uses `N=2^401+45`, dimension `2^200+13`, and

\[
 f(a)=((a-1)^2+4)^{-(2^{512}+3)},\qquad
 g(b)=((b-3)^2+4)^{2^{383}+5}.
\]

Both bases are nonunits. Their valuations vary between two and three.
Normalization discards one row leaf and two column leaves; the extra
refinement enforces the regular coefficient condition in (5). Five-bit
and seven-bit computations use two batched contractions each and agree on reduction.
No matrix nodes or expanded powers are listed. The retained evidence
records the exact residues, binary scale, precision budgets, and hashes.

This is a larger exact input contract for the inverse engine. The
factor-count extraction step remains unresolved.
