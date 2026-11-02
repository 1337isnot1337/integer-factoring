# Color parity obstructs large one-hot cube profiles

## Result and scope

Let a polynomial cube tensor over any field have one vacuum and
`d=q-1` nonvacuum colors. Assume:

- Each entry has homogeneous degree three times its number of nonvacuum
  vertex labels.
- At `(1,1,1)` the tensor is antipodal equality.
- An entry vanishes whenever any individual nonzero color occurs an odd
  number of times among its eight labels.

For **every `q>=256`**, its common-scale, unweighted three-axis Smith-profile
criterion with a common physical type has pressure coefficient

```math
\rho_{\mathrm{prof}}\ge\frac35>\frac1{\sqrt3}.
```

This excludes the proposed one-hot color-class Gaussian and every other
tensor satisfying these hypotheses, over every field. It does not need
Gaussian factorization, color-permutation symmetry, the complete generic
profile, or a specialization extrapolation. Large-color constructions
retaining separate color parity are stopped under this criterion.
Smaller alphabets, physical types that differ between vertices, weighted
gauges, actual coefficient ranks and tensors breaking a displayed hypothesis
remain outside this theorem. The complete factoring exponent stays `1/2`.

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

## 2. A shared-sector theorem beyond Gaussian kernels

Fix a set `C` of four distinct nonzero colors. A row assigning those
colors bijectively to the four vertices has a nonzero entry only if the
column uses exactly the same four colors, once each. Every color must
have even total cardinality, and all four column slots are needed to
complete the four row colors. The reverse coupling also vanishes.

Thus each color set gives a closed `24 by 24` direct summand in **each**
axis cut, indexed by the permutations of `C`. Its endpoint is a permutation
matrix, so its generic rank is 24. Let `m_(a,C)` be its smallest Smith
order at axis `a`.

Choose one tensor labeling that assigns the four colors to four antipodal
pairs, and let `f_C(x,y,z)` be that entry. It has degree 24 and endpoint
one, so it is nonzero. The same entry belongs to the `C` sector in every
axis. Since the smallest Smith order is the minimum entry valuation,

```math
\sum_a m_{a,C}\le\sum_a\operatorname{ord}_a f_C\le24.
                                                               \tag{1}
```

The last inequality holds because the product of the three coordinate
powers divides a nonzero polynomial of total degree 24. No equality or
axis symmetry is assumed.

The remaining rows and columns form a common complementary summand.
Fix all but one antipodal pair to vacuum and let that pair be vacuum or
one color. Parity makes the selected `(d+1) by (d+1)` minor diagonal:
`diag(1,f_1,...,f_d)`, where each `f_c` has degree six and endpoint one.
This is the **same physical minor** in all cuts, up to row/column order.
If `r_(a,i)`, `i=0,...`, are ascending complementary Smith orders, its
unit entry gives `r_(a,0)=0`. The minor valuation bound gives

```math
\sum_a\sum_{i=1}^d r_{a,i}
 \le\sum_a\operatorname{ord}_a\prod_{c=1}^d f_c\le6d.
                                                               \tag{2}
```

For a fixed axis, retain the first `d+1` complementary orders and one
smallest order from each of the `binom(d,4)` closed sectors. They are
disjoint channels. Apply convexity of log-sum-exp to the three vectors
of retained orders. Then apply convexity of `t^s` to the `d` averaged
nonunit complementary orders. Equations (1)--(2) yield

```math
\frac13\sum_a\log Z_a(t)
 \ge\log\left(1+\sum_{i=1}^d t^{\bar r_i}
                  +\sum_C t^{\bar m_C}\right)
 \ge\log L_d(t),\qquad
 L_d(t)=1+d t^2+\binom{d}{4}t^8,
                                                               \tag{3}
```

where bars denote the three-axis averages and `0<t<1`. In particular,
`sum_i bar(r_i)<=2d` and each `bar(m_C)<=8`. Direct sums let us retain
these channels without assuming compatible Smith bases for different cuts.

For the **Gaussian specialization**, the bound is stronger axis by axis:
every cross-kernel entry has valuation at least two, while the antipodal
entry is exactly `k=x^2y^2z^2`. A four-color sector entry is a product of
four matching cross entries, with the antipodal one equal to `k^4`.
Hence its smallest order is exactly eight in every axis, and the
complementary minor is `diag(1,k I_d)`. Thus `Z_a(t)>=L_d(t)` separately.
The general theorem needs only the averaged inequality (3).

The asymmetric control `T(a)=x^(3w(a))*delta_antipodal(a)` illustrates
why this distinction matters. It satisfies all general hypotheses, but
at `q=5,t=1/2` its x-cut partition is smaller than `L_d(t)`; the geometric
mean of all three cut partitions still obeys (3). The replay checks this
exactly from the literal four-letter side words.

## 3. An exact pressure certificate for all d>=255

For `d=255+r`, `r>=0`, define

```math
A(r)=(r+252)(r+253)(r+254)(r+255)=24\binom{d}{4}
```

and the integer polynomial

```math
P(r,s)=(24+24d s^6+A(r)s^{24})^5-24^5(1+d s^{10})^{12}.
```

It has degree at most 20 in `r`, degree 120 in `s`, and is divisible by
`s^6`. The checker forms the polynomial

```math
F(r,u)=(1+u)^{120} P(r,u/(1+u))/u^6.
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
the integer polynomial inequality and replays the asymmetric control. These constant-size checks support
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
not rerun. The general parity theorem is the ordinary argument in section 2;
the Gaussian identities are finite controls, not an exhaustive check of
all tensors. The obstruction does not use the specialized profile values. No new Lean or independent mathematical review is claimed.
