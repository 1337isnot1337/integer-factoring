# Full-degree polynomial realization in an eight-wire packet

## Exact finite result

For one Cauchy packet with `L=8`, fix the row tensor to be the indicator
of the zero state. Every polynomial in `Z_(2)[T]` of degree at most 120
is exactly realized by a suitable column tensor over `Z_(2)`. Here
`Z_(2)` is the ring of rational numbers with odd denominators.

The same fixed row tensor realizes the tapered-Hermite polynomials at the
full wire scale `M=8` for each `u=1,2,3`. The degrees are 64, 80, and 68.
The latter two contain powers of two in some global coefficient
denominators; their membership is checked separately.

This is a finite computer-assisted theorem. It extends the degree range of
the [all-size affine realization](cauchy-packet-realization.md) at one fixed
packet size. It supplies neither an all-size full-degree construction nor
canonical small-precision column entries. The separate
[alternating binary family](canonical-alternating-packet.md) settles a
different point: it controls creation precision, at the lower asymptotic
local-content scale `Theta(L log L)` and only at odd centers.

The source deformation and creation-precision convention are from
[OpenAI's quantitative contraction proof](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex),
revision `4eaee87548250b2682affda7d05d947c05554bae`. The tapered family is
defined in the
[PR #31 witness note](https://github.com/1337isnot1337/integer-factoring/blob/b1ba452fa4a3d19fe6fea3822c53613f3b582104/research/tapered-hermite-precision-witness.md).

## 1. A small exact jet matrix

With

```math
C(T)=\prod_{q=-8}^{7}(T-q),\qquad
W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)},
```

the polynomial for a column state `X` and the fixed zero row is

```math
F_X(T)=\prod_{x=0}^{7}\left(1-\sum_{y=0}^{7}W_{xy}(T)X_y\right).
\tag{1}
```

Every `F_X` has degree at most `8*15=120`. The complete contraction with
column tensor `B` is `sum_X B(X) F_X(T)`.

Let `V=Q[T]_(<=120)`. Record only eight local coefficients at each grid
point:

```math
J(f)=\bigl([z^j]f(t+16z)\bigr)_{-8\le t\le7,\ 0\le j<8}.
\tag{2}
```

This map from `V` to `Q^128` is injective. If all these coefficients vanish,
each of the sixteen grid points is a root of multiplicity at least eight.
The polynomial would have at least 128 roots counted with multiplicity,
unless it were zero.

Form the `128` by `256` matrix whose column indexed by `X` is `J(F_X)`.
Every entry lies in `Z_(2)`. The checker builds (1) by literal truncated
polynomial products, using the rational root product for each `W_xy`.

## 2. Why a finite modular calculation proves exact membership

The following finite-precision lemma is useful here.

**Lemma.** Let `A` be a matrix over `Z_(2)` of rational rank `r`, and let
`b` lie in its rational column span. Run minimum-valuation elimination
modulo `2^P`. Suppose it finds `r` pivots of valuations `v_i<P` and, at
each pivot, the transformed right side is divisible by `2^v_i`. Then
`b` belongs to the column module of `A` over `Z_(2)`.

To prove this, lift the same pivot sequence to exact rational operations.
Each pivot is `2^v_i` times an odd unit. Invert that unit and clear the
remaining entries in its row and column. Minimum valuation makes every
clearing coefficient integral. These operations therefore belong to
invertible matrices over `Z_(2)`.

The modular residues agree with the exact operations at every successful
step. Division of a pivot-column entry by `2^v_i` loses `v_i` bits, but
the right-side pivot is itself divisible by `2^v_i`. Thus the ensuing
product is still determined modulo `2^P`. The same observation applies to
unit normalization. Since `v_i<P`, the stated divisibility test is exact.

After all `r` pivots, the remaining exact matrix is zero. The remaining
exact right side is also zero, because `b` lies in the rational column
span. Solve the resulting diagonal system and undo the column operations.
The solution lies in `Z_(2)`. The argument uses finitely many rational
operations; it does not merely produce an infinite two-adic solution.

For the matrix (2), exact arithmetic modulo `2^32` finds 121 nonzero
pivots, with maximum valuation 26. These pivots prove rational rank at
least 121. The degree bound proves rank at most 121, so the rank is exact.
Consequently its rational image is all of `J(V)`.

The checker applies the right-side divisibility test to `J(T^k)` for every
`0<=k<=120`; all 121 pass. The lemma gives a realizing column tensor for
each monomial. Their `Z_(2)`-linear combinations prove the full claim for
`Z_(2)[T]` of degree at most 120.

## 3. Full wire-scale tapered examples

For `u=1,2,3`, set `h=floor(8/u)` and

```math
f_u(T)=\prod_{\substack{-8\le q\le7\\q\ne0}}
       \left(1-\frac Tq\right)^{h\max(u-v_2(q),0)}.
\tag{3}
```

Each has value one at zero and degree at most 120. The checker constructs
the local coefficients from the integer root products and their exact
denominator; it verifies that every recorded coefficient is two-adically
integral. Each of the three additional right sides passes the same module
test. The lemma therefore supplies exact column tensors for (3), with
`B(0)=1` because `F_X(0)=1_(X=0)`.

This closes the actual-packet realization question for these three finite
examples over `Z_(2)`. It leaves the source's canonical creation precision
obligation open. Reducing a suitable rational column tensor modulo `2^p`
can change its more precise local polynomials. The realization theorem
does not bound the height or precision of a favorable canonical lift.

## 4. A nonmember shows the restriction is real

The fixed-row module is smaller than the space of all polynomials that are
integral on the local disks. Consider

```math
g(T)=\frac{\prod_{j=0}^{119}(T-j)}{2^{112}},\qquad
112=60+30+15+7.
\tag{4}
```

On every disk `t+16z`, each linear factor has content
`min(v_2(t-j),4)`. Counting complete residue classes modulo `2,4,8,16`
shows that the numerator has content at least 112. Thus all coefficients
of every local polynomial in (4) are integral.

Nevertheless (4) is outside the fixed-row module. The global leading
coefficient of any column polynomial (1) has denominator valuation at
most

```math
8v_2(C'(0))=8\cdot11=88.
```

Linear combinations with coefficients in `Z_(2)` preserve this bound.
The leading coefficient of (4) has valuation `-112`, so it cannot occur.
The modular test independently rejects (4), at its first failed pivot.
Multiplying (4) by `2^26` passes the module test; this is included as a
positive control for the finite-precision argument.

## 5. Reproduction

```sh
python3 tools/check_cauchy_packet_lattice.py \
  --output /tmp/cauchy-packet-lattice.json
```

The checker uses only the Python standard library. It records the exact
rank, pivot valuation multiplicities, all successful target classes, the
first obstruction for (4), and source hashes in
[cauchy-packet-lattice.json](../evidence/cauchy-packet-lattice.json).
No floating-point rank calculation or assumed numerical optimum enters
this certificate. The factoring score remains `1/2`.
