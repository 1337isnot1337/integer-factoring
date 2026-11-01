# Exact eight-wire visibility for signed affine-code columns

## Finite result and scope

For the eight-wire Boolean Cauchy packet, fix any row tensor to the
indicator of a Boolean state `Y`. Let the column be any signed character
on an affine binary subspace containing `Y`. If its value at `Y` is one,
at least one coefficient of the complete odd-center packet is nonzero
modulo `2^32`.

There are exactly **2,013,762,304** such row and column pairs. The exact
checker covers all of them. The alternating point-row, point-column
example has minimum content 31, so 32 bits is the sharp visibility
threshold for this finite class.

This extends the finite exploration from length four to a structured
length-eight class. It does not cover arbitrary binary or ternary
columns, arbitrary row tensors, or differences of two arbitrary affine
columns. In particular, it does not supply a general decoder or a better
factoring algorithm. The complete factoring score remains `1/2`.

The source polynomial convention is the one used in the
[canonical decoder frontier](canonical-packet-decoder-frontier.md) and
[packet realization note](cauchy-packet-realization.md). The lower example
is the existing [alternating family](canonical-alternating-packet.md).

## 1. The exact input family

For `L=8`, set

```math
 C(T)=\prod_{q=-L}^{L-1}(T-q),\qquad
 W_{xy}(T)=\frac{C(T)}{C'(0)(T-y+x)}.
```

For each Boolean column state `X`, the point-row polynomial is

```math
 F_{Y,X}(T)=
 \prod_{Y_x=1}(W(T)X)_x
 \prod_{Y_x=0}\bigl(1-(W(T)X)_x\bigr).
 \tag{1}
```

It has degree at most `L(2L-1)=120`, and `F_(Y,X)(0)=[X=Y]`.
The local matrix `J_Y` contains every coefficient of
`F_(Y,X)(t+16z)` at the eight odd centers `t=-7,-5,...,7`.
Thus it has 968 rows and 256 columns. These coefficients are dyadically
integral, so reduction modulo `2^q` is defined using odd denominators.

Take any linear subspace `H<=F_2^8`, with its unique reduced-row-echelon
basis `v_1,...,v_k`, and any `alpha in F_2^k`. The column is

```math
 B\left(Y+\sum_i c_iv_i\right)=(-1)^{\sum_i\alpha_i c_i},
 \qquad c\in\mathbf F_2^k,
 \tag{2}
```

and is zero outside `Y+H`. Its entries are in `{-1,0,1}` and
`B(Y)=1`, so its original contraction against the point row is one.
Multiplying the whole column by minus one preserves visibility and
covers the alternative normalization `B(Y)=-1`.

The number of linear subspaces by dimension is

| Dimension | Number of subspaces | Characters per subspace |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 255 | 2 |
| 2 | 10,795 | 4 |
| 3 | 97,155 | 8 |
| 4 | 200,787 | 16 |
| 5 | 97,155 | 32 |
| 6 | 10,795 | 64 |
| 7 | 255 | 128 |
| 8 | 1 | 256 |

Their sum is 417,199. Including the characters gives 7,866,259 columns
for each fixed point row. Multiplication by the 256 row choices gives
the stated total. The checker computes these counts independently from
the Gaussian-binomial product formula.

## 2. Complete enumeration with necessary projection filters

For each `Y`, form three fixed integer linear combinations of all rows
of `J_Y`, modulo `2^32`. Call the resulting column vectors `h_1,h_2,h_3`.
If every original local coefficient of `J_Y B` vanishes, all three
projected sums also vanish. The projections therefore cannot discard a
true witness. Their pseudorandom seed affects runtime only; correctness
does not depend on a probability estimate.

Enumerate every reduced-row-echelon basis exactly once. The pivot
positions determine which entries can be free, and every assignment
to those free entries is visited. For each basis, list its `2^k`
subset sums in binary coefficient order. A Walsh transform computes
all `2^k` signed sums in (2) simultaneously. Addition and subtraction
are ordinary integer operations modulo `2^32`.

Test the first projected sum for every character, then the remaining
projections for its survivors. Every survivor of all three projections
is checked against **all 968 original local coefficients**. At precision
32 there are no projection survivors at any of the 256 row states.
Consequently there are no complete-packet witnesses in this class.

Unsigned 64-bit arithmetic computes the projected matrix products.
Overflow is reduction modulo `2^64`, which preserves the required
residue modulo `2^32`. The proof uses this divisibility relation;
it does not interpret wrapped integers as exact unbounded sums.
The polynomial construction itself also uses only integer residues.

## 3. Positive controls and independent polynomial checks

At `Y=85`, whose occupied positions alternate, the length-eight
precision-31 search returns three witnesses. Each is replayed against
every original coefficient. The point-column witness alone establishes
the necessary lower bound, by the alternating-family theorem.

At length four, the same program covers all 4,912 signed affine-code
pairs at precision 12 and finds none. The precision-11 alternating
control returns three exact witnesses. These agree with the stronger
existing length-four result for all binary rows and all ternary columns.

The wrapper also reconstructs 44 local polynomials independently using
exact rational arithmetic in the original variable `T`, then substitutes
the local center and dilation. It compares every coefficient with the
modular implementation. The controls include lengths two, four, and
eight. A separate small enumeration checks that the reduced-row-echelon
routine gives all 67 distinct subspaces of `F_2^4`.

The saved receipt records every row state, coverage counts, original
matrix hashes, the positive witnesses, and source hashes. Hashes are
provenance checks. Re-running the finite enumeration is the verification
procedure.

## 4. Reproduction

NumPy and a C++17 compiler are required. No floating-point rank or
optimization result enters the certificate.

```sh
python3 tools/check_canonical_affine_codes.py --full \
  --output /tmp/canonical-affine-codes.json
cmp /tmp/canonical-affine-codes.json evidence/canonical-affine-codes.json
```

Without `--full`, the command runs the independent controls, the entire
length-four class, and the length-eight positive control. The full run
additionally covers all 2,013,762,304 length-eight pairs.

Files: [modular jet constructor](../tools/canonical_affine_code_jets.py),
[replay wrapper](../tools/check_canonical_affine_codes.py),
[native enumeration](../tools/check_canonical_signed_affine_codes.cpp), and
[receipt](../evidence/canonical-affine-codes.json).

The finite result narrows one explicit class in the canonical precision
question. An all-size efficient reconstruction method remains a separate
obligation. No new Lean formalization is claimed.
