# Exact rank envelopes for every nilpotent-seam coefficient

The [four-state seam construction](nilpotent-seam-crossings.md) has an
explicit rank envelope at every bidegree, including its middle coefficients.
An eight-counter integer recurrence accounts for both nilpotent runs and
the exact residual vacuum dimension. It is polynomial time and uses no
floating-point rank or asymptotic extrapolation.

The bounds give, in particular, a uniform rectangle of coefficients
with ranks at most **`(L^(5/4+o(1)),L^(9/8+o(1)))`**. Another family has the bounds
`(L^(1/2+o(1)),L^(4/3+o(1)))`. These are candidate pieces for a graded
decomposition. Their simultaneous upper bounds do not prove that they
are nonzero or supply the rank mass needed to span the crossing.
No stronger complete factoring exponent is established.

## 1. A vacuum-weighted finite-state counter

Write `v_f` for the dimension of the vacuum space on 2f strands.
At the height-three wall,

```math
v_0=v_1=1,\qquad v_{f+2}=3v_{f+1}-v_f.
```

The two-step path matrix on heights zero and two is
`[[1,1],[1,2]]`; this proves the recurrence and
`v_f=O(phi^(2f))`, where `phi=(1+sqrt(5))/2`.

A central degree pattern is a word `d_1,...,d_n` in `{0,1,2}`.
Its degree is `sum d_j`. Let k be its number of degree-two letters.
Every such letter supplies `N tensor e_i`; a degree-zero letter is I.
Four degree-two letters separated only by degree-zero letters vanish,
because `N^4=0`. Treating every degree-one letter as a reset gives a
necessary condition for a nonzero word. Some allowed patterns can
still vanish; no converse is assumed.

Call a word admissible when, after deleting its zeros, every run of
twos between ones has length at most three. Define the integer

```math
W_{n,a}=\sum_{\substack{d\text{ admissible}\\\sum d_j=a}}
             v_{n-\#\{j:d_j=2\}}.              \tag{1}
```

Each degree-two letter forces a central cup. For a fixed pattern,
the mandatory cup projector has exact rank `v_(n-k)`, by the
insertion/removal splitting proved in the parent note. All remaining
central physical factors commute with these disjoint mandatory cups.
The combined physical/auxiliary operator therefore has rank at most
`4v_(n-k)`, even after summing its degree-one choices. Summing patterns
and extracting the outer y-degree b gives

```math
\boxed{\quad r(G_{a,b})\le4h_{n,b}W_{n,a},\qquad
h_{n,b}=\min(b+1,n-b+1).\quad}                  \tag{2}
```

This uses the actual F4 path quotient, so no extension-degree multiplier
is introduced. The trivial bound `L^2` can always be applied as well.

For an exact small counter, keep the mandatory-N run `h in {0,1,2,3}`
and degree a. At each state store the pair of sums
`(v_f,v_(f+1))` over its pattern-dependent free-pair counts f.
The three transitions are:

| Letter | New run | Degree increment | Counter pair |
| --- | --- | ---: | --- |
| 0 | h | 0 | `(x,y) -> (y,3y-x)` |
| 1 | 0 | 1 | `(x,y) -> (y,3y-x)` |
| 2 | h+1, if h<3 | 2 | `(x,y) -> (x,y)` |

Initialize state `(h,a)=(0,0)` with pair `(1,1)` and sum the first
counter over final h. This is (1), by induction on the word length.
It uses eight integer counters per degree, `O(n^2)` arithmetic
operations and `O(n^3)` bit operations: every counter has `O(n)` bits.
The subtraction in the recurrence is exact; all counter values remain
nonnegative sums of vacuum dimensions.

At the maximal allowed degree `a_*=2n-floor(n/4)`, the weights for
n=4,6,8 are respectively4,2,20. The exact counter generally gives a
much smaller finite bound than the earlier stars-and-bars estimate.
For example the eight-site row bound is `80h_(n,b)` before the `L^2`
cap. Each outer half has y-degree at most n/2; counting its coefficient
convolutions gives the sharper h_(n,b), including near b=n.
Actual coefficients can have substantially smaller ranks.

### Fixed-width near-leading bands

The half-power row exponent also holds a fixed number g of degrees
below `a_*`, even if that top allowed coefficient vanishes. Write
`n=4m+s`, `0<=s<4`. A pattern of degree `a_*-g` with z zeros has
`m+g-2z` ones and `n-m-g+z` twos. The bound on the number of twos
implies z<=g. The runs of twos have total slack
`3-s+4g-7z<=4g+3` from their maximum length three.

For each z there are at most `(n+1)^(4g+3)` possible run-length
patterns and at most `(n+1)^g` ways to insert the zeros. There are
at most m+g free pairs. Consequently

```math
W_{n,a_*-g}\le(g+1)(n+1)^{5g+3}v_{m+g}.         \tag{2a}
```

For any fixed g this gives `r<=L^(1/2+o(1))`, uniformly over b.
This is a new bound for the lower degrees; it does not silently
reuse the leading-degree estimate when the leading coefficient is zero.

## 2. Opposite-cut factors at every boundary degree

In degree b, choose the b gates of the boundary seam that supply a
physical e instead of Q. There are `binom(n,b)` choices. The
[vacuum-interface factorization](fusion-projector-seams.md#1-vacuum-interfaces-instead-of-physical-bit-channels)
compresses their joint interface to the 2b-strand vacuum space of
dimension v_b. This replaces the looser physical4^b channel count.

Diagonalize Q over F16. The fixed auxiliary input and output use one
regular two-dimensional copy, giving two incoming and two outgoing
eigenchannels. On the central seam the first-block gates precede the
second-block gates; splitting there leaves four auxiliary states.
Thus there are at most `2*2*4=16` auxiliary channels. After also
fixing one of the v_b vacuum-interface paths, the physical operator
is a product of two block polynomials in x.
The degree-a coefficient is a convolution of at most a+1 products.
Consequently

```math
\boxed{\quad c(G_{a,b})
 \le16(a+1)\binom nbv_b.\quad}                  \tag{3}
```

The factorization can be taken over F16, but the final tensor has F4
entries and the same matrix rank over both fields. Dense F4 factors,
if required, have the allowed `2^O(lambda)` construction cost.

Equations (2),(3), capped by `L^2`, are an explicit joint envelope
for the whole coefficient grid. The
[standard-library counter](../tools/nilpotent_seam_rank_bounds.py)
computes them using integers only.

## 3. An exactly certified middle window

The following window has row rank at most `L^(5/4+o(1))` and column
rank at most `L^(9/8+o(1))`:

```math
\left\lceil\frac{12n}{7}\right\rceil\le a\le a_*,
\qquad 0\le b\le\left\lfloor\frac n8\right\rfloor.             \tag{4}
```

Here is a rational certificate, rather than a numerical saddle point.
The path matrix above gives `v_f<=phi^(2f)<=(21/8)^f`. To bound
the generating sum of (1), put `c=21/8`, `z=256` and use a four-state
nonnegative transfer matrix T:

* a zero letter adds weight c to `h -> h`;
* a one letter adds weight cz to `h ->0`;
* a two letter adds weight `z^2` to `h ->h+1`, when h<3.

Thus `sum_a W_(n,a) z^a <= 1^T T^n e_0`. The positive rational vector

```math
w=(1,11/4,(11/4)^2,(11/4)^3)^T
```

satisfies `Tw<=24000 w`, entry by entry. Since `e_0<=w`,

```math
W_{n,a}\le\frac{2055}{64}\frac{24000^n}{256^a}.                \tag{5}
```

The exact integer inequality

```math
24000^{28}5^{35}<256^{48}8^{35}
```

proves `24000/256^(12/7)<(8/5)^(5/4)<phi^(5/4)`.
Together with (2), this is the asserted row exponent in (4), with
only polynomial prefactors.

For the column bound use `v_b<=(21/8)^b`. The binomial generating function at `t=8/147`
gives

```math
\binom nbv_b\le(8/7)^n(147/8)^b
 \le\bigl((8/7)(147/8)^{1/8}\bigr)^n.
```

The other exact integer inequality

```math
147\,5^9<64\,7^8
```

shows that this base is less than `(8/5)^(9/8)`.
Finally `L=Theta(phi^n)` at the least admissible path length, and
`phi>8/5`. This proves (4) uniformly after any common binary restriction.
The window has a central width growing as n/28 and a boundary width
growing as n/8; it contains middle bidegrees, not just the extreme edge.

For a shoulder example, combine the earlier highest-degree row bound
`L^(1/2+o(1))` with (3) for `b<=n/6`. The choice `t=8/105` gives
the base `(6/5)(105/8)^(1/6)`, and

```math
105\,6^6 5^2<8^9
```

certifies the column exponent at most4/3. Thus this family has the
rank bounds `(L^(1/2+o(1)),L^(4/3+o(1)))`.

## 4. Rank windows are pieces, not a complete crossing

The middle window has polynomially many pieces and rank cost
`sum sqrt(rc)<=L^(19/16+o(1))<L^sqrt(2)`. It must still participate
in an exact decomposition. A collection of pieces with a favorable
upper budget alone supplies no spanning theorem or factoring algorithm.

The [minimum-cut constraint](minimum-cut-rank-budget.md) still requires
minimum-rank mass L. The stronger
[weighted constraints](weighted-cut-rank-budget.md) require intermediate
aspect ratios too. For example, the middle window together with both
original edges and their cut reflections is still insufficient:
at `s=floor(L^(3/8))`, every such piece contributes at most
`L^(5/4+o(1))` to `sum min(r,s^2c)`, whereas that sum must be at
least `sL=L^(11/8+o(1))`. Polynomially many pieces do not close this gap.
This rules out a search restricted to just those few regimes; it does
not exclude a graded selection using the full envelope or new recombinations.

Cut reflection is available without a rank penalty. Interchange Q,R
in a constant tensor; this exchanges the two cuts and fixes SWAP.
Append each reflected candidate beyond the old degree as
`(t^k+t^(k+1))T`. Both strict endpoints are preserved, and only
polynomially many coefficient slots are needed. This permits both
orientations of every certified candidate window.

The open positive task is a compatible graded spanning identity,
using actual coefficients or a new hybrid kernel, with the full
rank-sum target and every necessary cut budget respected. The
[scalar-relaxation analysis](geometric-cut-relaxation-limit.md) already
shows that the necessary budgets alone do not settle the tensor problem.

## 5. Replay and source scope

[The checker](../tools/check_nilpotent_seam_rank_bounds.py) compares
the counter with independent literal ternary-word enumeration, checks
the rational certificates with exact fractions and integer powers,
and verifies every rank pair in the complete4/6/8-site coefficient
grids. It also checks cut reflection and degree boundaries. The
[receipt](../evidence/nilpotent-seam-rank-envelope.json) pins the source
bytes and records the exact finite results.

```sh
make verify-braid-rank-envelope
```

The complete finite matrix check needs NumPy; the rank-bound counter
itself uses only the standard library. This extends the nilpotent seam,
the prior exact cup-image splitting and the existing cut-budget framework.
It does not assert a rank-sum optimum, a complete spanning identity,
independent mathematical acceptance or a score above `1/2`.
