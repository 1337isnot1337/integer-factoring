# Protected lower states exclude every homogeneous split through height 25

## Result

The [height-16 protected-submatrix obstruction](coherent-boundary-protected-floor.md)
extends substantially without calculating another large Smith profile.
The same fixed 28-letter coordinate space already forces

\[
\inf\rho\ \ge\ \frac{17}{24}\ >\ \frac1{\sqrt2}
\tag{1}
\]

for every constant, height-homogeneous physical split supported on occupation
words of height at most **25**, provided it retains all words of height
below 16. The larger map spaces need not be deletion closed. The two
paired physical alphabets may differ, and their types may be optimized
independently. The scale is common to the two cuts, and the heights are
the inherited natural ones.

This excludes the entire 96-letter height-25 boundary family, including
its nonbinary hyperoval construction, from a strict exponent improvement
through this pressure criterion. It also excludes arbitrary homogeneous
changes at several of the heights 16 through 25, rather than only a single
boundary retraction. All fields of characteristic two are covered.

The theorem is an obstruction for the stated method. It does not cover
changes to the protected lower states, mixtures between unequal heights,
nonconstant interpolation-parameter maps, independent cut scales, other
physical gauges, or direct bounds on actual coefficient ranks. No stronger
factoring exponent is claimed.

## 1. Protect the actual polynomial tensor

Fix a sufficiently large reflection interval and use the regular zero
and infinity tensors in the variable `t=w^2`. Their entries are polynomials.
Their inherited physical height is

\[
h(a)=\sum_{i\in a}(2i+1).
\]

Let `U` be the span of all 28 words with `h(a)<16`. On each paired physical
role, choose constant maps with `LP=I` and the block form

\[
P=I_U\oplus P_{\rm high},\qquad
L=I_U\oplus L_{\rm high}.
\tag{2}
\]

Every paired column and row has one common height. All high supports in
(2) have heights between 16 and 25. Separate choices are permitted for
the `X,V` and `Y,U` physical pairs.

Because the maps preserve height, they commute with the natural diagonal
normalizations of both local places. Thus their actual regular cut matrices
are obtained by applying the same constant maps to the ambient regular
polynomial tensors. Restrict every physical index of these compressed
matrices to the coordinate labels in `U`. The maps then act as identities,
so both cuts contain exactly the protected literal tensor restrictions.

Only this fixed space `U` needs the deletion-closure argument. It is a
height downset. Its two literal restricted profiles therefore equal the
complete 784-order compound profile established in the predecessor note.
Write its counting polynomial as

\[
Z(r)=\sum_{i=1}^{784}r^{n_i}.
\]

It has degree 82 and `Z(0)=1`, with coefficient of `r` equal to one.
Both larger regular matrices have this protected polynomial submatrix,
regardless of whether their high map spaces are deletion closed. The
Toeplitz-rank lemma from the predecessor gives

\[
Z_0(r)Z_\infty(r)\ge Z(r)^2\qquad(0<r<1).
\tag{3}
\]

This is why the extension does not assume that the larger profiles equal
their compressed compound profiles. It protects the literal regular
tensor matrices directly. All candidates and their protected restrictions
use the same ambient regularization. In particular, the common scalar
shifts at zero and infinity are the same before restriction and cancel
in the two-cut product. The protected vacuum entry is one in the regular
zero tensor, hence a pure monomial in the raw tensor. Any removable common
scalar polynomial therefore divides that monomial and is itself a
monomial, up to a nonzero constant. Removing it gives opposite valuation
shifts, which again cancel in the product.

The fixed protected matrices have coefficients in `F_2`. Their Smith
orders remain unchanged after every characteristic-two field extension,
so (3) holds over every such field.

## 2. The full height-25 entropy budget

Every word of height at most 25 uses only the first thirteen hot modes.
The largest possible physical height polynomial is the truncation

\[
A_{104}(t)=\left[\prod_{j=0}^{12}(1+t^{2j+1})\right]_{\deg\le25}.
\tag{4}
\]

Its coefficient list, starting at degree zero, is

```text
1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3,
3, 3, 4, 5, 5, 5, 6, 7, 8, 8, 9, 11, 12.
```

The sum is 104. The rank condition `LP=I` bounds the number of physical
labels at each height by the ambient dimension of that height cell.
Consequently either paired alphabet has `A_X(t), A_Y(t) <= A_104(t)`
for `0<t<=1`.

The final height-25 coefficient includes the single mode of height 25,
the five-particle word 31, and ten three-particle words. The 12-hot-mode
96-letter construction omits that single mode. Thus (4) is a larger
entropy budget than that construction requires; omitting the extra
state does not weaken the exclusion.

The exact scalar certificate proves

\[
Z(s^{17})^{12}>A_{104}(s^{24})^{17}
\quad\text{for every }0<s\le1.
\tag{5}
\]

Set `rho_0=17/24` and `r=s^17`. From (3) and (5),

\[
\log Z_0(r)+\log Z_\infty(r)
>4\rho_0\log A_{104}(r^{1/\rho_0}).
\tag{6}
\]

Let the independently chosen paired types have entropies `H_X,H_Y`
and expected inherited heights `E_X,E_Y`. Their Gibbs inequalities give

\[
H_X+H_Y+\frac{\log r}{\rho_0}(E_X+E_Y)
\le\log A_X(r^{1/\rho_0})+\log A_Y(r^{1/\rho_0})
\le2\log A_{104}(r^{1/\rho_0}).
\]

For positive total entropy, (6) therefore proves

\[
\frac{\log Z_0(r)+\log Z_\infty(r)
      -2\log r\,(E_X+E_Y)}{2(H_X+H_Y)}
>\rho_0.
\tag{7}
\]

This is the paired-type natural pressure, with every physical height and
both entropies charged. Its infimum is at least `17/24`. As in the
predecessor, pointwise strictness is not promoted to a strict infimum
bound. The exact comparison `2*17^2-24^2=2` suffices for the strict
separation from the factoring threshold.

## 3. Integer certificate and independent finite input

The [checker](../tools/check_coherent_protected_height25.py) reconstructs
(4) in two independent elementary ways: truncated product multiplication
and direct enumeration of all 8,192 subsets of the thirteen modes. It
checks the same protected profile used by the predecessor when requested.
Those controls compare the original binary and field-valued normalized
compound constructors on both cuts, totaling 3,136 complete orders. Both
constructors use the established native Smith engine; no separate new
Smith implementation is claimed.

For `0<s<=4/5`, set `t=s^24`. The coefficient list in (4) gives

\[
A_{104}(t)\le1+t+102t^3\le1+2t,
\]

because `102*(4/5)^48<1`. The predecessor's unchanged geometric-series
estimate then proves (5):

\[
Z(s^{17})^{12}\ge1+12s^{17}
>1+\frac{34s^{24}}{1-32s^{24}}
\ge A_{104}(s^{24})^{17}.
\]

The strict middle comparison has the same integer margin
`182169346247795596` as before.

On `[4/5,1]`, the checker uses **685 rational intervals**, with maximum
dyadic depth 14. On each interval, the left side of (5) at its lower
endpoint strictly exceeds the right side at its upper endpoint. Both
polynomials are increasing, so these comparisons cover every point.
The checker verifies exact endpoint adjacency and full interval coverage.
The two homogeneous denominator degrees are 16,728 and 10,200; their
difference is 6,528. Every comparison uses integers only.

The retained [evidence](../evidence/coherent-protected-height25.json)
contains every interval, the profile and height coefficients, and the
precise theorem scope. Reproduce it with

```sh
python3 tools/check_coherent_protected_height25.py \
  --recompute-profiles --output evidence/coherent-protected-height25.json
```

Omit `--recompute-profiles` to replay only the standard-library scalar
certificate. This extension builds on the established reflection
normalizations, coherent compound theorem for the fixed protected space,
and the predecessor's Toeplitz-rank monotonicity argument. It had
substantial OpenAI assistant support and adds no Lean formalization.
