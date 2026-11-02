# Distinct-color sectors obstruct the large one-hot Gaussian cube

## Result and scope

The one-hot color-class Gaussian cube has a generic profile obstruction
at every physical alphabet size `q>=256`. For its **common-scale,
unweighted three-axis Smith-profile criterion with a common physical
type**, the pressure coefficient satisfies

```math
\rho_{\mathrm{prof}}\ge\frac35>\frac1{\sqrt3}.
```

This excludes the proposed `q=256` instance and every larger alphabet
in that method. It does not require the complete generic profile,
a specialization extrapolation or a large matrix. These large-color
instances are stopped. Smaller alphabets, physical types that differ between vertices,
weighted gauges, actual coefficient ranks and other constructions are
outside this theorem. The complete factoring exponent remains `1/2`.

The construction and specialized controls are the original worker-4
contribution, retained from commit
`1772970a7dbaab83eeb0607014a1111e46718d86`, including its
[grading corrections](color-class-cube-profile.md). The new argument
uses exact direct summands and a complementary minor of that construction.

## 1. The tensor and its correct grading

Let `d=q-1`, with colors `1,...,d` and vacuum zero. Use the homogeneous
Gaussian cube kernel over characteristic two from
[cube_homogeneous_gaussian.py](../tools/cube_homogeneous_gaussian.py).
For a labeling of the eight vertices, put

```math
T_q(a)=\prod_{c=1}^d G_{S_c}(x,y,z),
\qquad S_c=\{v:a_v=c\}.
```

Here `G_S` is the principal Pfaffian and `G_empty=1`; odd cardinalities
vanish. Every entry is homogeneous of degree three times its number of
nonvacuum vertices, and the endpoint is antipodal equality. Thus the
physical letter heights are zero and one, with generating function
`W(v)=1+d v`. They are not equal-height side words. The side-word
multiplicities are `1,4d,6d^2,4d^3,d^4`.

For axis `a`, write `Z_a(t)=sum_i t^{s_(a,i)}`, with ascending generic
Smith orders and `0<t<1`. Two exact kernel identities are used:

- Every cross-axis kernel entry has valuation at least two in that axis.
- Every antipodal entry is the monomial `k=x^2 y^2 z^2`.

These follow directly from the displayed augmentation kernel and are
rechecked symbolically in all axes by the replay.

## 2. Four distinct colors give closed 24-dimensional blocks

Fix a set `C` of four distinct nonzero colors. A row assigning those
colors bijectively to the four vertices has a nonzero tensor entry only
if the column uses exactly the same four colors, once each. Indeed,
every color must have even total cardinality. Each row color already
occurs once, and four column slots are needed to complete all four.

Consequently each color set gives a closed `24 by 24` direct summand,
indexed by the permutations of `C`. Its entry is the product of the four
kernel entries matching equal colors. Every entry has axis valuation
at least eight, and the antipodal matching is exactly `k^4=x^8y^8z^8`.
Its smallest Smith order is therefore eight. The block has a permutation
endpoint, so all its orders are finite.

There are `binom(d,4)` such independent summands. Their contribution to
each partition function is at least `binom(d,4)t^8`.

In the complementary summand, fix all but one antipodal pair to vacuum
and allow that pair to carry either vacuum or one color. Its selected
minor is exactly `diag(1,k I_d)`. It lies entirely outside the distinct
four-color summands. Thus the complement has a unit order and the sum
of its first `d+1` orders is at most `2d`. Convexity of `t^s` gives a
complementary contribution at least `1+d t^2`. Adding the disjoint
summands proves the generic bound

```math
\boxed{Z_a(t)\ge L_d(t):=1+d t^2+\binom d4 t^8}
```

for every axis. The block separation is literal in the physical basis;
no coupled rank-preserving change of basis is assumed.

## 3. An exact pressure certificate for all d>=255

For `d=255+r`, `r>=0`, define

```math
A(r)=(r+252)(r+253)(r+254)(r+255)=24\binom d4
```

and the integer polynomial

```math
P(r,s)=(24+24d s^6+A(r)s^24)^5-24^5(1+d s^10)^12.
```

It has degree at most 20 in `r`, degree 120 in `s`, and is divisible by
`s^6`. The checker forms the polynomial

```math
F(r,u)=(1+u)^120 P(r,u/(1+u))/u^6.
```

Its 2415 coefficients, in degrees `0..20` of `r` and `0..114` of `u`,
are all nonnegative integers. Every coefficient of its constant-in-`r`
part is strictly positive; the minimum is `10152345600`. The
[receipt](../evidence/color-distinct-sector.json) records the digest,
and the checker reconstructs every coefficient. Three separate rational
evaluation identities also check both expansions.

Hence `P(r,s)>0` for `r>=0` and `0<s<1`. Setting `t=s^3` gives

```math
L_d(t)^5>(1+d t^{10/3})^{12},
\qquad L_d(t)>W(t^{10/3})^{12/5}.
```

No search over alphabet sizes or floating-point pressure optimization
is used. The nonnegative parameter coefficients prove the entire
range `d>=255` at once.

## 4. Entropy and the three-axis target

Let the common physical type have entropy `H(p)>0` and mean occupancy
`E`. The [graded cube pressure formula](cube-graded-minor-budget.md)
uses

```math
\rho_{\mathrm{prof}}\ge
\frac{\frac13\sum_a\log Z_a(t)-8E\log t}{4H(p)}.
```

All logarithms use the same base. The physical-letter entropy inequality
is `H(p)+E log v<=log W(v)`. Use `v=t^(10/3)` and the lower bound above:

```math
\frac13\sum_a\log Z_a(t)\ge\log L_d(t)
 >\frac{12}{5}\log W(v)\ge\frac{12}{5}H(p)+8E\log t.
```

Thus the pointwise profile expression is greater than `3/5`; taking its
infimum preserves the floor `rho_prof>=3/5`. Finally `3*3^2=27>25=5^2`
proves `3/5>1/sqrt(3)`. The grading is used at the vertex-letter level
throughout, avoiding the earlier equal-height and side-word comparisons.

## Reproduction and review scope

```sh
make verify-color-distinct-sector
```

The [checker](../tools/check_color_distinct_sector.py) verifies 1728
sector endpoint identities, 3816 zero-coupling identities and 108
complementary-minor entries in the three formal variables, then proves
the integer polynomial inequality. These constant-size checks support
the displayed general parity and minor arguments; they do not enumerate
a `256^4` profile.

The inherited source hashes match the original receipt. Its historical
`q4_full_profile` key contains only 250 of 256 orders through order 39;
that key is not a complete-profile certificate. The per-layer counts
need not be generic lower bounds: the specialization guarantee applies
to cumulative counts through each order. Selected independent regression
commands replayed the binary specialized profile, the ternary profile
through order six, and the saved graded comparisons. The larger full
specialized screens were retained with their original scope and were
not rerun. The new generic obstruction does not use their numerical
profile values. No new Lean or independent mathematical review is claimed.
