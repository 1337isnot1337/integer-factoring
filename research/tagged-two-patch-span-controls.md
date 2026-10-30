# Complete small controls for a tagged trace coordinate

The [two-patch span calculation](two-patch-span-recombination.md) suggests
assigning degree two to one trace factor. This note completes the bounded
experiment for every such choice at widths three and four. No tested
binary coefficient span improves the one-piece crossing cost `d^2`.

## 1. Family and exact endpoint

The coefficient field is F2; the extension-field models for the labels
and the trace coordinates are specified below.

For the two-patch variables `B_0=S(S+w)` and `B_1=(S+1)(S+w)`, choose
one patch `ell in {0,1}` and one trace coordinate `j_0`. Use

```math
G_{\ell,j_0}(z;S,w)=
 (1+z^2\tau_{j_0}(B_\ell))
 \prod_{(i,j)\ne(\ell,j_0)}(1+z\tau_j(B_i)).
```

As before, `tau_j(a)=Tr(alpha^j a)` and the four-leg tensor is obtained
by multiplying the kernel by `[x+y=u+v]`, with `S=x+y,w=x+u`.
The endpoint at one is unchanged and equals the crossing.

The degree is exactly `2t`. Its maximum possible degree `2t+1` would
require the impossible all-one trace pair. Missing just one untagged
coordinate realizes degree `2t`. Any distinct pair `(B_0,B_1)` is
realized by `e=B_0+B_1`, `S=B_0/e`, and `w=S+e`, so this last degree
is attained. There are at most `2t+1` coefficient tensors.

## 2. Complete finite result

| `t` | Field polynomial | Tag choices completed | Span dimension for each | Exact minimum cost for each |
| ---: | ---: | ---: | ---: | ---: |
| 3 | `X^3+X+1` | 6 | 7 | 64 |
| 4 | `X^4+X+1` | 8 | 9 | 256 |

Every choice of patch and trace coordinate was tested. The checker reduces
the coefficient tensors to an independent basis before enumerating its
span, records the original coefficients in that basis, and recovers the
target mask by exact binary elimination. It does not assume independence.

For each span, all cut ranks are exact. The integer-floor XOR-potential
certificate from the companion note gives the lower bound `d^2`; the
single crossing tensor attains it. Arbitrarily many pieces and arbitrary
binary recombinations within each displayed span are included.

The completed evidence contains:

* all 14 profiles and 4,864 tensor evaluations across their coefficient spans;
* every row/column rank pair and every lower potential value;
* 2,190,592 checked potential inequalities;
* 1,596 comparisons with independently constructed full four-leg cuts.

For tag `(0,0)`, every tensor at both widths receives a full dense check.
The other tags receive full checks for their original coefficient tensors,
the crossing, zero, and selected recombinations. Every tensor in every
profile is ranked by the exact convolution-block calculation.

The result concerns only these binary spans with these field models and
trace coordinates. It does not exclude larger widths, different weights,
other trace bases, or recombinations with extension-field coefficients.

## 3. A bivariate extension and one finite control

A larger coefficient space retains two variables instead of merging
the tagged coordinate into the scalar degree. Set

```math
H(z,u;S,w)=
 (1+u\tau_{j_0}(B_\ell))
 \prod_{(i,j)\ne(\ell,j_0)}(1+z\tau_j(B_i)).
```

Its four-leg endpoint at `(z,u)=(1,1)` is exactly `E`. Its bidegree is at
most `(2t-1,1)`, and the coefficient of `z^(2t-1)u` vanishes. Therefore
there are at most `4t-1` nonzero coefficient tensors. The
[bivariate dimension proof](bivariate-two-patch-span-control.md) shows
that these are independent, giving dimension exactly `4t-1` at every width.
Its span contains
both previous scalar spans: substitute `u=z` for the untagged curve and
`u=z^2` for the tagged curve. The separate scalar-span certificates do
not exclude decompositions mixing tensors from this larger space.

Write the bivariate tensor as `sum_j z^j A_j + u sum_j z^j B_j`. The
identity and rank budget to seek are

```math
E=\sum_i Q_i,\qquad
Q_i=\sum_j(a_{ij}A_j+b_{ij}B_j),\qquad
\sum_i\sqrt{r(Q_i)c(Q_i)}<d^{\sqrt2}.
```

The displayed rank budget is a target, not a claimed bound. It retains the higher Boolean
degrees needed for the endpoint and uses no ordinary-kernel composition,
so the earlier low-degree and composition-cancellation proofs do not by
themselves exclude it. The coefficient count remains linear in `t`.

The [subsequent bivariate control](bivariate-two-patch-span-control.md)
completes the small diagnostic `t=3`, tag `(0,0)`. Its dimension is exactly
11, and all 2,048 binary span tensors have been checked. The exact optimum
is again `d^2=64`. No other bivariate tag, width, or trace basis is excluded
by that calculation. The stronger displayed rank budget remains the
requirement for a successful larger construction.

## Reproduction and bounded completion

[Exact checker](../tools/check_tagged_two_patch_span.py) and
[complete fourteen-profile certificate](../evidence/tagged-two-patch-span-verification.json).

```sh
python3 tools/check_tagged_two_patch_span.py \
  --output /tmp/tagged-two-patch-span-verification.json
```

The checker uses the standard-library field and span tools. Its default
time limit is 110 seconds; it atomically saves each completed profile and
labels an interrupted record as partial, with exit status two. The committed evidence has status
`passed` and contains all 14 expected profiles. Its finite claims do not
depend on any unfinished computation or timing extrapolation.
