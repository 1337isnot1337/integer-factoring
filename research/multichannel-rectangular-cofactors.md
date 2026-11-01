# Rectangular cofactors must pay for several Smith channels

The existing rectangular cofactor gap uses the second crossed Smith order.
Larger aspect ratios also force a cofactor to use the third, fourth and
later channels. This gives further necessary scalar crossing inequalities
over **every field**, with arbitrary degree and independent physical heights.

An exact eight-letter point passes all previous cell/shared, square
adjugate and rectangular first-jet constraints and has critical pressure
below `-1/500`. The new constraints exclude it. This proves that the new
conditions strengthen that relaxation; it supplies no realizing tensor or
factoring algorithm. The complete factoring score remains **1/2**.

## 1. Integral rank capacity under reshuffling

Work over a DVR `O` with uniformizer `t`. Let `T` be an ordinary operator
on a `p`-by-`b` matrix space, of dimension `n=pb`. Let its crossed matrix
have increasing Smith orders `s_1,...,s_n`, and set `h=min(p,b)^2`.
An integral Smith decomposition writes the crossed matrix as

```
sum_j t^(s_j) a_j b_j^T,        a_j,b_j in O^n.
```

Un-crossing one outer product gives a transpose channel
`M -> A_j M^T B_j^T`. For `p<=b`, the intermediate matrix `A_j M^T`
has size `p`-by-`p`. This is a factorization through `O^h` with integral
factors. For `b<p`, the other association factors through `O^(b^2)`.
Consequently, concatenating all channels gives an integral factorization

```
T = U diag(t^s1 repeated h times, ..., t^sn repeated h times) V.
```

Cauchy--Binet implies that every `k`-minor of `T` has order at least the
sum of the `k` smallest exponents in that repeated multiset. In particular,
every cofactor has order at least

```
L_(p,b)(s) = h*(s1+...+sa) + e*s_(a+1),
a = floor((n-1)/h),        e = (n-1) mod h.             (1)
```

Zero cofactors satisfy the bound with infinite order. Integral Smith
factors, rather than a factorization over the fraction field, are essential
here: no uncharged poles enter the Cauchy--Binet terms. The argument also
works for Laurent orders after a scalar shift.

For `p=2,b=6`, equation (1) is `4s1+4s2+3s3`. With `s1=0`, this refines
the older rank-one entry bound `7s2` when `s3>s2`. The
[first-jet result](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/rectangular-adjugate-first-jet.md)
separately bounds the rank of the first cofactor coefficient. Neither
argument is discarded: both sets of prefix constraints should be imposed.

## 2. Graded caps and physical restrictions

Use the normalized full profiles `r,s`, deficit `d` and paired energy
lists `X,Y` from the [cell-set framework](cell-minor-caps.md). Choose
`p` X labels and `b` Y labels. Put

```
n=pb,    W=b*sum(X_selected)+p*sum(Y_selected),
C=(n-1)d+W.
```

The graded cofactor identity and its negative paired energies are exactly
those of the credited
[graded adjugate proof](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/graded-adjugate-caps.md)
and its [rectangular extension](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/subalphabet-adjugate-caps.md).
The new input is (1): every crossed adjugate entry, hence every one of its
Smith orders, is at least `L_(p,b)(s)`. Substituting the resulting first-`k`
prefix lower bound into the existing cell and shared-line arguments gives

```
k S_(n-1)(r) + k L_(p,b)(s) <= k C - TopCell_k,        1<=k<=n,  (2)
```

```
l S_n(r) - Top_l(r) + l L_(p,b)(s)
    <= l C - TopLine_l,                    1<=l<=max(p,b).     (3)
```

Here `TopLine_l` is the larger of the permitted terms
`TopX_l+l*Ymax` (if `l<=p`) and `TopY_l+l*Xmax` (if `l<=b`). The top
orders in (3) refer to the first `n` orders, before restriction is lifted.
Both cut orientations hold; their average yields the same conditions on
`gamma=(r+s)/2`.

To lift from a restricted tensor, convert to absolute profiles `mu,nu`.
The coefficients in `L` sum to `n-1`. Thus each of the two order lists
has coefficient sum `k(n-1)` in (2), or `l(n-1)` in (3), exactly matching
the coefficient of `d=-mu_1-nu_1`. No residual minimum term remains. All
absolute coefficients are nonnegative. The individual Smith interlacing
proof in the rectangular extension therefore transfers these inequalities
to the full profiles, including restrictions with changed minima.

The right sides increase with every selected physical energy, by the same
top-cell/top-line multiplicity argument as the existing extension. It
suffices to select the smallest `p` X energies and smallest `b` Y energies.
The [compiler](../tools/multichannel_cofactor_caps.py) uses this choice.
Square cases only give the trivial `L=(n-1)s1`, so the implementation
retains rectangular cases and leaves the stronger square adjugate inputs
unchanged. Rational heights follow by ramification.

## 3. Exact separation from the preceding constraints

The [retained rational point](../evidence/multichannel-q8-separation.json)
has eight letters, denominator `100000`, two zero paired energies, and 64
orders. It was discovered numerically, then rounded with smaller orders
and larger deficit/nonzero energies. The rounded point is checked with
exact rational arithmetic, independently of the discovery optimizer.

It satisfies all original cell/shared caps, all 476 oriented square
adjugate constraints and all 1,428 compiled rectangular first-jet caps.
The identical two profiles make the other cut orientation automatic.
100-bit outward arithmetic proves its critical pressure is below `-1/500`.
It violates 236 new multichannel caps; the largest recorded violation is
`525457/6250` on the `2`-by-`8` determinant case of (2).

Thus the point is impossible over every field. It remains evidence that
the earlier necessary inequalities alone missed a further obstruction.
It is not a crossing or a stronger factoring bound.

## 4. Reproduction and limits

Run `make verify-multichannel-cofactors`. The
[checker](../tools/check_multichannel_cofactor_caps.py) independently
interpolates 269 literal cofactors over F10007. The three examples attain
(1) at orders `6,7,10`, against older rank-one bounds `3,5,7`. Another
61 binary cofactors are expanded directly by permutations, with eleven
independently computed source Smith orders. Nine independently graded
strict scalar tensors supply 36 complete profiles, 1,872 lifted cap checks
and their absolute-normalization identities.

The [compressed packet](../evidence/multichannel-cofactor-controls.json.xz)
retains source matrices, cofactors, profiles, exact separation and source
hashes. The first-jet compiler is reused byte-for-byte from source commit
`365f319fbb62c734bb409c93d39323d566b87c49`, preserving attribution and its
earlier proof. The new theorem is the integral capacity argument and its
graded transfer, not a replacement for the earlier cofactor work.

These are necessary scalar constraints. They do not exclude every energy
shape or supply a useful tensor. The earlier clean full-suite receipt at
`d31ffec` predates this package. No new Lean work or independent mathematical
review is claimed. Substantial OpenAI assistant support was used.
