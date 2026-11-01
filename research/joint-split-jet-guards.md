# Rank guards and formal rigidity for rectangular joint splits

Two matrix facts narrow the
[joint-compression search](joint-grid-compressions.md).
First, a protected two-letter core and a degree budget force a necessary
rank pattern for any balanced natural improvement: each regular cut must
have exactly one zero order and exactly one order one. Otherwise its
balanced pressure is at least `5/7 > 1/sqrt(2)`.
Second, at the stated coordinate split of the five-letter source, fixing
the first two actual matrix coefficients at both places makes every
formal deformation of the constant joint maps produce the same tensor.

A concrete two-parameter joint fiber retains both leading coefficients
but changes the first opposite residual. Every nonzero parameter pair,
over every characteristic-two field, has balanced `rho >= 3/4`.
These are ordinary proofs about the declared local constructions.
They do not improve the factoring exponent, which remains **1/2**.
Asymmetric supporting weights are outside the pressure bounds here.

## 1. A spectral guard from a small core

Consider a four-by-four crossing with physical heights `(0,1,3,4)`.
Assume its two naturally regular matrices are polynomial, have generic
rank 16 and determinant degree at most 128. Suppose their two-letter
physical subtensor has complete orders `(0,1,3,4)` at both cuts.
The crossing need not come from a joint split for this section.

For either full profile `theta`, sorted increasingly, the subtensor's
minors give the ascending-prefix bounds `0,1,4,8` through prefix four.
Also `sum theta<=128`, since a nonzero polynomial determinant's valuation
is no larger than its degree. For `j>=4`, monotonicity of the remaining
orders gives

\[
\sum_{i=1}^j\theta_i
\le s_4+\frac{j-4}{12}(128-s_4)
\le8+10(j-4),\qquad s_4\le8.
\]

Every ascending prefix is therefore bounded above by that of
`(0,1,3,4,10,...,10)`, with twelve tens. The positive-part integral
argument from [joint conjugations](joint-height-conjugation-floor.md)
implies, for `0<x<1`,

\[
Z(x)=\sum_i x^{\theta_i}\ge Q(x)
 =1+x+x^3+x^4+12x^{10}.                            \tag{1}
\]

If the leading matrix has rank at least two, the first two orders are
zero. Combining this with the same core prefixes and budget gives the
stronger virtual profile `(0,0,4,4,10,...,10)`. Thus

\[
Z(x)\ge P_0(x)=2+2x^4+12x^{10}.                   \tag{2}
\]

If there is exactly one zero order but at least two orders one, use
`(0,1,1,6,10,...,10)` instead:

\[
Z(x)\ge P_1(x)=1+2x+x^6+12x^{10}.                 \tag{3}
\]

Let `F(x)=1+x+x^3+x^4`. For either `i=0,1`, the balanced natural
criterion at `rho=5/7` follows from the integer polynomial inequality

\[
\bigl(Q(x^5)P_i(x^5)\bigr)^7-F(x^7)^{20}>0
\qquad(0<x<1).                                   \tag{4}
\]

The [checker](../tools/check_joint_split_jet_guards.py) proves (4) in
the integer Bernstein basis. After removing the leading monomial,
the two degrees are 700 and 695 and the minimum integer coefficients
are 127 and 21. No numerical optimization enters this argument.
One cut supplies (2) or (3), and the other supplies (1).

The protected core forces at least one zero order. If there is exactly
one, its second prefix forces at least one order one. Consequently a
balanced `rho<5/7` requires exactly one of each at *both* cuts.
This is also necessary for beating `1/sqrt(2)`, since `50>49`.
In matrix terms, the leading matrix and its first DVR residual must
both have rank one. This necessary condition does not fix their row
or column factors.

## 2. The five-letter joint model satisfies the degree budget

Let `A(z)` be the strict 25-by-25 operator of the existing five-letter
source, in row `(X,Y)` and column `(V,U)` order. Its physical heights
are `(0,5,3,1,4)` and its joint heights are their pair sums. Write
`D_H(z)` for the diagonal of joint monomials. The exact source obeys

\[
\operatorname{ord}_z A_{ij}\ge8-H_i,
\qquad \deg A_{ij}\le H_j+8.                       \tag{5}
\]

These 625 entry bounds are checked directly on the polynomial source.
Any constant height-preserving rectangular maps `P,L` to the target
grid retain them. After zero regularization by `z^(h_row-8)`, an entry
has degree at most `h_row+h_col`. The infinity reversal gives the same
bound, with the new opposite row and column heights. The determinant
degree at either cut is at most twice the sum of the sixteen joint
heights, which is 128.

Thus section 1 applies whenever such a split also retains the stated
two-letter physical subtensor. The [F4 interface](../tools/joint_five_letter_splits.py)
checks `LP=I`, homogeneity, every endpoint entry and regularity before
using the full-rank native profile promise. Constant maps in the
spectral variable are assumed; arbitrary moving maps require their
own degree and core checks.

## 3. The first nontrivial affine fiber already loses a residual

Take the coordinate split `P_0,L_0` retaining physical heights
`(0,1,3,4)` on each paired direction. Source joint coordinates use
`X+5Y`, while target coordinates use `4x+y`. Target coordinate ten is
`(2,2)` and embeds as source coordinate twelve. Put

\[
n=e_8+e_{16},\qquad
P_\beta=P_0+\beta n e_{10}^{\mathsf T},\qquad
L_\alpha=L_0+\alpha e_{10}n^{\mathsf T}.            \tag{6}
\]

All three involved source coordinates have joint height six.
Moreover `L_0 n=0`, `n^T P_0=0` and `n^T n=0`, so (6) has `LP=I`
for arbitrary parameters over every characteristic-two field. It fixes
the protected physical core.

Both first zero matrices and the leading infinity matrix are unchanged.
Eliminate the leading infinity pivot at row and column zero. The next
residual consists of the following disjoint blocks after row and column
permutations: an all-ones four-by-four block, two singleton entries
`alpha`, and two singleton entries `beta`. There are no mixed terms.
Its exact rank is

\[
1+2[\alpha\ne0]+2[\beta\ne0].                      \tag{7}
\]

The nonzero entries and original indices are retained symbolically in
the [receipt](../evidence/joint-split-jet-guards.json); the calculation
uses formal powers, without imposing `alpha^2=alpha` or field-size
identities. In particular any nonzero pair introduces at least three
orders one at infinity. With total order at most 128, the remaining
twelve integer orders have exponential sum at least
`7x^10+5x^11`: smoothing integer pairs minimizes the convex sum at the
most balanced distribution of total 125. Therefore

\[
Z_\infty(x)\ge1+3x+7x^{10}+5x^{11}.
\]

Together with (1), an integer Bernstein certificate of degree 71 proves
the balanced floor `rho>=3/4` for every nonzero pair in (6).
All sixteen F4 pairs additionally have complete native profiles,
independently replayed with the Python fraction-free oracle.

## 4. Fixing both first two jets gives formal rigidity

Define the normalized unreshuffled matrix coefficients

\[
Z_j=[z^j]\,z^{-8}D_H(z)A(z),\qquad
W_j=[z^j]\,z^8 A(z^{-1})D_H(z),\qquad j=0,1.
\]

The infinity coefficients are reshuffled only after compression.
Let `P(epsilon),L(epsilon)` be height-preserving matrices over
`K[[epsilon]]`, constant in `z`, with constant terms `P_0,L_0` and
`LP=I`. Suppose each of `L Z_j P` and `L W_j P` is exactly its
baseline matrix, for both `j=0,1`.

**Theorem.** `L(epsilon) A(z) P(epsilon)=L_0 A(z) P_0`.
This holds over every characteristic-two field `K`.

There are 88 allowed first-order entries in variations `(X,Y)` of
`P,L`. Linearizing `LP=I` and the four matrix moments gives a binary
system of rank 81. Its seven-dimensional kernel consists exactly of

\[
X=P_0 C,\qquad Y=C L_0,
\]

where every generator `C` preserves target height and commutes with
the *entire* polynomial `L_0 A(z) P_0`. The checker retains the seven
generators and verifies every polynomial commutator entry. Binary
linear algebra extends this complete kernel description to every `K`.
For comparison, fixing only `Z_0,W_0` gives tangent dimension twelve,
with a two-dimensional projection outside the coordinate image.

For the formal proof, take the first nonzero variation at order
`epsilon^d`. It lies in that kernel. The invertible formal matrix
`U=I+epsilon^d C` preserves height, commutes with the frozen moments
and commutes with the baseline tensor. Replace `P` by `P U^(-1)`
and `L` by `U L`. This cancels that variation while retaining every
constraint. Repeating at increasing orders gives a convergent formal
product of such gauges and reduces the maps to `P_0,L_0`. The product
commutes with the baseline tensor, so the original tensor was already
unchanged. This argument also handles inseparable parameter powers.

The theorem concerns a formal neighborhood of the specified split.
It does not classify disconnected solutions, parameters singular at
the baseline, or maps that preserve Smith orders while rotating the
first matrix jets. Rank-one leading and residual *factors* may still
move; section 1 only constrains their ranks.

## 5. Reproduction and verification scope

```sh
make verify-joint-split-jet-guards
```

The checker reconstructs the rank-81 system and its complete kernel,
the formal residual pattern (7), all three integer pressure certificates,
all 625 source entry bounds, both protected-core profiles, and all
512 orders of the sixteen F4 fiber pairs with independent fraction-free
elimination. The core supplies eight additional order controls.
These are local author checks, with no new Lean work or independent
mathematical review. The full default run at `decce64` predates this
follow-up and is not claimed as a new full run for these files.
