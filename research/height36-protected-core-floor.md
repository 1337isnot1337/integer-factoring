# A protected low-core barrier through physical height 36

## Theorem

Consider the homogeneous, deletion-closed coherent physical retractions
on twelve hot modes from
[the coherent construction](coherent-hot-f4-retractions.md).
Work over any field of characteristic two. Assume:

- Both maps retain, literally, every coordinate label of height below 24.
  For these labels the embedding and retraction are the coordinate vector
  and its coordinate dual.
- Every physical label is supported on occupation words of height at most 36.
- The paired maps satisfy `LP=I`, preserve height and occupation, and satisfy
  both deletion-closure conditions required by the construction.

Then their **balanced natural pressure satisfies `rho >= 17/24`**.
Since `17/24 > 1/sqrt(2)`, this family cannot cross the threshold for a
stronger factoring exponent through the balanced natural interface.

The field, coefficients, rank and retained subspaces above the protected
core are unrestricted subject to these hypotheses. This is broader than
changing only the height-25 boundary. It includes the registered six-pair
F8 differential-kernel trials with caps 25 and 27.

This does not exclude changing lower labels, using larger physical heights,
other physical gauges, or asymmetric interfaces. It also does not exclude
small improvements among the retained small-alphabet crossing bounds.
The complete factoring exponent remains `1/2`.

## 1. A binary submatrix that survives every allowed change

The height is `h(w)=sum(2i+1)` over the occupied bit positions.
There are **81**, not 86, coordinate words of height below 24 among twelve
modes. Their paired channel matrices are literal submatrices of both cuts
of every crossing in the theorem.

The exact low-order engine computes the first 3,359 Smith orders of each
core matrix, all at most 64. The two lists agree. Each core matrix has
6,561 channels, so **3,202 higher orders per cut remain unknown** in this
certificate. They are not assigned zero or treated as a complete profile.

Let `Z(t)` be the sum of `t^nu` over only the 3,359 known core orders.
The determinantal-prefix inequality for a submatrix, followed by the
hinge-sum argument for decreasing exponentials, shows that each full
channel partition sum is at least `Z(t)` for `0<t<1`.
Discarding the unknown core tail only weakens this lower bound.
The core is binary, so its Smith orders are unchanged by extension to
any characteristic-two field.

There are 294 coordinate words of height at most 36. Define

```math
Q(t)=\sum_{h(w)\le36}t^{h(w)}.
```

At each height, independent physical embeddings cannot outnumber the
ambient coordinate words. Consequently the physical partition sum of
any allowed retraction is at most `Q(t)`. An ambient reflection interval
of size 64 is sufficient for this fixed height range; the stabilized hot
matrices do not depend on choosing a larger interval.

It therefore suffices to prove

```math
24\log Z(z^{17})-34\log Q(z^{24})\ge0
\qquad(0<z<1).                                      \tag{1}
```

Indeed, set the natural scale to `17 log(1/z)` and test `rho=17/24`.
The two channel lower bounds and physical upper bound then give the
claimed pressure obstruction at every positive scale.

## 2. Exact global sign certificate

Write `Z(t)=sum c_i t^i` and `Q(t)=sum d_j t^j`. The derivative of the
left side of (1) has the sign of the integer polynomial

```math
G(z)=\sum_{i,j}(i-2j)c_i d_j z^{17i+24j-17},
```

omitting the zero `(i,j)=(0,0)` contribution. All exponents are nonnegative.
Its degree is 1,935. The coefficient signs of
`(1+x)^1935 G(x/(1+x))` have exactly three variations.
Exact rational endpoint evaluations isolate three distinct roots in:

```text
[9081026427,9081026428] / 10^10,  signs +,-
[9876550127,9876550128] / 10^10,  signs -,+
[9975603842,9975603843] / 10^10,  signs +,-.
```

Descartes' rule proves that these are all roots in `(0,1)`, including
multiplicity. Only the middle root is an interior local minimum.
At that root, outward 192-bit logarithm bounds give a lower gap numerator

```text
600557868083443333992346191126445124265864003164297843266
```

in units of `2^-192`, which is positive. The limiting gap at zero is zero.
At one, positivity follows from the exact integer inequality
`3359^24 > 294^34`. This proves (1) globally, without a sampled scale grid.

## 3. Reproduction and evidence limits

The compact arithmetic check is:

```sh
python3 tools/check_height36_core_floor.py
```

To rebuild both native low-order core lists as well:

```sh
python3 tools/check_height36_core_floor.py --replay-core
```

The [low-core receipt](../evidence/height24-core-low-orders.json) records
the exact words, both retained lists, all explicit tail counts and source
hashes. The [inequality receipt](../evidence/height36-core-floor.json)
records root intervals, integer evaluation hashes and the positive gap.
The default check validates the stored core receipt and recomputes the
integer/logarithmic proof. Only `--replay-core` recomputes the native matrices.
Neither command claims to compute the missing 3,202 orders per cut.

The low-order construction and tail interpretation are those of
[the certified partial-profile engine](coherent-hot-low-order-certificates.md).
The proof uses its binary specialization. This is an ordinary finite
certificate and argument, not a new Lean formalization or a claim of
independent mathematical acceptance.
