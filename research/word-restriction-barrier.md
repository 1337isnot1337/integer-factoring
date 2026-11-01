# Arbitrary binary-word restrictions of the 23-term profile

## Result

Over any field of characteristic two, no injective restriction to binary
words of subexponential length can improve
the leading contraction of the 23-term Pfaffian coefficient profile. This
holds for arbitrary label subsets, arbitrary collective recombinations of
the complete coefficient symmetry span, and **any number of final pieces**.

Here is an exact finite statement. Let `C` be any set of `D` distinct binary
words of length `m>=2`. Remove its constant words and let `D_*>=D-2` be
the number remaining. Write `H_(m-1)=sum_(j=1)^(m-1) 1/j` and set

```math
L=\left\lceil\frac{D_*}{2(m-1)H_{m-1}}\right\rceil.     (1)
```

Every decomposition `E=sum_i F_i` in the coefficient symmetry span defined
below has a single piece with both cut ranks at least

```math
r_i,c_i\ \ge\ \left\lfloor L^2/4\right\rfloor.          (2)
```

The claim is nontrivial when `L>=2`. The proof selects two disjoint word
families `X,Y` on which the complete polynomial becomes

```math
P(z)|_{X\times Y\times Y\times X}=z^{2m+6}E_{X,Y}.      (3)
```

Every coefficient image and its eight crossing-preserving leg symmetries
therefore restricts to zero or the same full-rank permutation matrix.
Some one piece must have a nonzero restriction in both cuts.

When `D=2^lambda` and `log_2 m=o(lambda)`, (1)--(2) force both widths
to be `2lambda-o(lambda)`. The leading contraction is at least one.
For quadratic evaluation words, `m=2^r` and `r=O(sqrt(lambda))`,
giving the explicit lower width bound
`2lambda-2r-O(log r)`. This closes the quadratic coefficient candidate
without any condition on its weight distribution or piece count.

A sharper direct criterion is retained as a useful corollary. If a linear
code containing the constant word has every nonconstant weight strictly
between `m/3` and `2m/3`, then some piece has both ranks at least
`D^2/4-D`. Quadratic prefix codes also contain such balanced subcodes
of dimension at least `lambda-2r`.

The theorem concerns this coefficient construction, not factoring or arbitrary
effective crossings. Tensors created outside the stated span remain possible
escape routes. The factoring exponent stays `1/2`. These are ordinary
mathematical proofs with exact replay, without new Lean formalization or
a claim of literature priority.

## 1. Sources and the precise span

The local 23-term tensor is pinned to
[`hamming_tensor_power.py` at PR #24](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/tools/hamming_tensor_power.py).
The earlier
[`affine-code-barrier.md`](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/affine-code-barrier.md)
proves a barrier for affine evaluation codes and identifies quadratic
Boolean evaluation codes as a remaining candidate. The exact quadratic
scalar oracle is recorded separately in `quadratic-code-oracle.md`.

For a label set `C` contained in `F_2^m`, define

```math
P_{xyuv}(z)=\sum_{j=0}^{22}L_j(z)A_j(0000;z)
                 \prod_{t=1}^{m}A_j(x(t),y(t),u(t),v(t);z).             (4)
```

The scale-one identity gives `P(1)=E`, where
`E(x,y,u,v)=[x=v][y=u]`. The eight allowed leg permutations are precisely
the permutations preserving the partition `{{x,v},{y,u}}`.

The coefficient symmetry span contains every coefficient tensor of (4),
each of those eight images of it, and all linear combinations of these
tensors. An individual piece may use arbitrarily many coefficient images.
There is also no bound on the number of pieces in (2). General label
bijections, different code embeddings, or different local tensors are
outside this defined span.

## 2. A source-table criterion using only pair projections

For a local support `S_j={p in F_2^4:A_j(p;z) != 0}`, project onto a
pair of coordinates joining the two crossing pairs. These coordinate
pairs are `(x,y),(x,u),(v,y),(v,u)`.

Every support except `S_22` omits at least one pattern on one of these
pairs. The following witnesses are sufficient; bits in the final column
are ordered as in the stated pair.

| Terms | Projected pair | A missing pattern |
| --- | --- | --- |
| 0, 5 | `(v,y)` | `01` |
| 1, 6 | `(v,u)` | `11` |
| 2, 7 | `(v,y)` | `10` |
| 3, 8 | `(v,y)` | `11` |
| 4, 10, 19, 20 | `(v,u)` | `11` |
| 9, 21 | `(v,u)` | `01` |
| 11, 12, 15, 16, 17, 18 | `(v,u)` | `10` |
| 13, 14 | `(v,u)` | `01` |

These statements are exact checks of the fixed 16-entry local tensors.
There are no assumptions about the number or uniformity of evaluation
histogram values.

The remaining term is exactly

```math
A_{22}(x,y,u,v;z)=z^2[x=v][y=u],\qquad L_{22}(z)=z^4.   (5)
```

The fixed leading factor `A_22(0000;z)=z^2` accounts for the remaining
two powers in (3).

**Rectangular restriction lemma.** Let `X,Y` be sets of binary words of
length `m`. Suppose every `f in X,g in Y` realizes all four ordered
patterns `00,01,10,11` at some word coordinates. Then (3) holds.
The same identity holds after any of the eight allowed leg symmetries.

To prove it, take `x,v in X` and `y,u in Y`. For every noncrossing
local term, its missing pair pattern occurs at some coordinate of the
chosen quartet, so the product of its local entries in (4) is zero.
Only (5) remains. Every crossing-preserving leg symmetry still puts
one word from each set on every projected cross pair. Exchanging the
order of the two sets preserves the four-pattern property. This proves
the symmetric version as well.

This lemma applies to arbitrary sets of binary words, without linearity
or balance. Those assumptions below provide large suitable sets explicitly.

## 3. Every large binary-word set has a suitable rectangle

Remove the at most two constant words of `C`. Let `D_w` denote the
number of remaining words of Hamming weight `w`, for
`1<=w<=m-1`. For each weight and ordered pair of distinct coordinates
`(a,b)`, consider the fiber

```math
C_{w,a,b}=\{f\in C:\operatorname{wt}(f)=w,\ f(a)=1,\ f(b)=0\}.
```

Let `L_max` be the maximum size of these fibers. Every weight-`w` word
belongs to exactly `w(m-w)` ordered-coordinate fibers. Averaging over
the `m(m-1)` ordered coordinate pairs gives

```math
D_w\le \frac{L_{\max}m(m-1)}{w(m-w)}.                  (6)
```

Since

```math
\sum_{w=1}^{m-1}\frac1{w(m-w)}
  =\frac2m H_{m-1},
\qquad
D_*\le2L_{\max}(m-1)H_{m-1},                          (7)
```

some fiber has size at least the integer `L` in (1).

Any two distinct words in this fiber realize every pair pattern.
The fixed coordinates `a,b` supply `11` and `00`. Distinct words
of equal weight cannot contain one another's supports, so both `01`
and `10` also occur. Split the fiber into disjoint sets `X,Y` whose
sizes differ by at most one. The rectangular restriction lemma gives
(3), with `|X||Y|>=floor(L^2/4)`. Section 5 proves the rank conclusion.

All quantities in (1) can be computed with exact rational arithmetic.
A simpler but weaker bound is
`L>=ceil(D_*/(m(m-1)))`: choose a largest interior weight layer and use
`w(m-w)>=m-1`. The harmonic bound is used in the theorem.

The rectangle can be found explicitly by listing weight layers and
counting their ordered-coordinate fibers. This takes
`O(Dm^2)` elementary word-coordinate operations. For subexponential
word length and `D=2^lambda`, that cost is `2^O(lambda)`.

The combinatorial argument may also use the number of *distinct coordinate
columns* in the `D`-by-`m` word table. Removing duplicate columns
preserves every pair's set of realized patterns and preserves distinct
word labels. A rectangle found after this compression remains a valid
rectangle in the original table. Repeating coordinates cannot evade
the obstruction by artificially increasing the stated word length.

## 4. Balanced codes provide larger rectangles

In this section, assume `C` is a linear binary code, `1 in C`, and
every `f in C` other than zero and one satisfies

```math
m/3<\operatorname{wt}(f)<2m/3.                         (8)
```

**Disjoint-support lemma.** If `f,g in C` have pointwise product zero,
then `f=0`, `g=0`, or `g=1+f`.

Constant cases are immediate. Otherwise disjoint support gives
`wt(f+g)=wt(f)+wt(g)>2m/3`. By (8), `f+g` must be constant.
It cannot be zero, since `f=g` and `fg=0` would force `f=0`.
Therefore it is one.

Suppose now that `f,g,f+g` are all nonconstant. Every pair pattern
occurs. If one were missing, one of

```math
fg,\quad f(1+g),\quad(1+f)g,\quad(1+f)(1+g)
```

would vanish identically. Applying the disjoint-support lemma would
force `f`, `g`, or `f+g` to be constant, a contradiction.

Let `B=C/<1>`, so `|B|=D/2`, and choose a nonzero linear functional
`h:B -> F_2`. With bars denoting quotient classes, set

```math
X=\{x\in C:h(\bar x)=0,\ \bar x\ne0\},\qquad
Y=\{y\in C:h(\bar y)=1\}.                            (9)
```

For every `x in X,y in Y`, the three words `x,y,x+y` are nonconstant.
Thus `X,Y` satisfy the rectangular restriction lemma. Their sizes are

```math
|X|=D/2-2,\qquad |Y|=D/2,\qquad
|X||Y|=D^2/4-D.                                       (10)
```

Both sets include the two constant lifts of each of their quotient
classes. The quotient-zero words zero and one are excluded from `X`.

## 5. A one-dimensional restriction forces one common bad piece

On `X x Y x Y x X`, each coefficient image is zero except at degree
`2m+6`, where it equals `E_(X,Y)`. The same is true for all eight
leg symmetries. Every tensor in their span therefore restricts to a
scalar multiple of `E_(X,Y)`.

If `E=sum_i F_i`, write these restrictions as `a_i E_(X,Y)`.
When the lower bound is nontrivial, the rectangle is nonempty, so
`sum_i a_i=1`. At least one scalar
`a_i` is nonzero.

In the row flattening `(x,y)|(u,v)`, the restricted crossing pairs
`(x,y)` bijectively with `(u,v)=(y,x)`. In the column flattening
`(x,u)|(y,v)`, it similarly pairs `(x,u)` with `(y,v)=(u,x)`.
Both matrices are permutation matrices of order `|X||Y|`. A nonzero
scalar multiple has the same rank. Restriction cannot increase rank,
so the same chosen piece has both ranks at least `|X||Y|`. This proves
(2) for the arbitrary-code fiber and the sharper bound (10) for balanced
codes.

This argument works over every characteristic-two extension field,
allows arbitrary collective recombinations, and is independent of the
number of final pieces.

## 6. Large balanced subcodes of quadratic evaluation codes

Use the standard prefix labels in the monomial order

```math
1,t_1,\ldots,t_r,t_1t_2,\ldots,t_{r-1}t_r
```

at the least `r` with `1+r+binom(r,2)>=lambda`. For sufficiently
large `lambda`, all affine monomials are present. Let `Q` be the
pure-quadratic coefficient subspace, of dimension `n=lambda-r-1`.

Every nonzero alternating form of rank two is decomposable:
`a b^T+b a^T` for two independent vectors. Including zero, there are
at most `2^(2r)` alternating forms of rank at most two, including
those in any prescribed subspace `Q`.

A greedy construction therefore gives `Q_0<=Q` with
`dim Q_0>=n-2r` and no nonzero form of rank at most two. If the
current dimension is `k<n-2r`, the translates of the bad set by the
current subspace contain at most `2^k 2^(2r)<2^n` points. A vector
outside that union extends the subspace without introducing a bad
nonzero member. Enumerating coefficient vectors and checking ranks
gives an explicit uniform construction in `2^O(lambda)` time.

Adjoin all `r` linear forms and the constant one. The resulting
subcode has dimension at least

```math
(n-2r)+(r+1)=\lambda-2r.                              (11)
```

For a quadratic of alternating rank `2s`, variable-pair elimination
gives a Gauss sum either zero or of absolute value `2^(r-s)`.
At rank at least four, this is at most `2^(r-2)`, so its Hamming
weight is in

```math
[3\cdot2^{r-3},\ 5\cdot2^{r-3}].
```

This interval is strictly inside `(2^r/3,2^(r+1)/3)`. Adding affine
terms leaves the alternating rank unchanged. Pure nonconstant affine
forms have weight `2^(r-1)`. The constructed subcode consequently
satisfies (8), proving the balanced-subcode corollary.

The argument concerns a loss of only `O(r)` coefficient dimensions,
which is the quantity needed for the exponent calculation.

## 7. Consequences for the effective coefficient candidates

For `D=2^lambda` and `m>=2`, the harmonic bound yields, whenever the
right side is nontrivial, one piece with both logarithmic widths at least

```math
2\lambda-2\log_2 m-2\log_2 H_{m-1}-O(1).              (12)
```

The harmless constant covers the at most two removed words, rounding,
and the split into two sets. Since `H_(m-1)<=1+ln(m-1)`,
subexponential word length gives widths `2lambda-o(lambda)`.
No restriction on the number of final pieces is required.

For quadratic evaluation words `m=2^r`, this becomes
`2lambda-2r-O(log r)`. It is stronger asymptotically than applying
only the balanced-subcode construction, which gives
`2lambda-4r-3`. Both arguments are retained because the balanced
criterion directly gives the larger bound (10) when the whole label
code is balanced.

The result covers raw coefficient pieces, growing-support coefficient
searches, and collective refactorizations inside the defined span. For
the full affine evaluation code, every nonconstant word has weight
exactly half the length; (10) gives
`r_i,c_i>=D^2/4-D`, improving and simplifying the earlier affine
construction obstruction.

For a family with `m=D^(gamma+o(1))`, (12) gives the necessary
contraction bound `rho>=1-gamma`. Thus this construction cannot
reach `rho<1/sqrt(2)` unless its effective word length has exponent
at least `1-1/sqrt(2)` relative to the number of labels. This is a
necessary condition, not a construction attaining that parameter.

**Coordinate weights do not remove the restriction.** If each local
coordinate uses the supported substitution `z -> z^(w_t)`, with a
fixed leading weight `w_*`, all omitted local patterns remain zero.
The surviving term on the same rectangle becomes

```math
z^{\,4+2w_*+2\sum_t w_t}E_{X,Y}.
```

Hence the entire span of coefficients from any number of such weight
choices still restricts to a single line of scalar multiples of
`E_(X,Y)`. This statement concerns substitutions into the fixed
23 local tensors; assigning independent variables before the source's
local cancellations is a different construction.

**Several label embeddings can share one rectangle.** Suppose a common
set of `D` labels has `J` injective word embeddings of lengths `m_j>=2`,
including arbitrary label permutations. Put

```math
C_j=2(m_j-1)H_{m_j-1},\qquad P_J=\prod_{j=1}^J C_j.
```

Apply the fiber selection successively, always to the currently retained
labels. Its size satisfies `L_j >= (L_(j-1)-2)/C_j`. Since every `C_j>=2`,
induction gives `L_J>=D/P_J-2`: if the preceding loss is at most two,
the next loss is at most `4/C_j<=2`. All pair projections established
at earlier steps remain complete after further restriction. Split the
final common fiber into two sets. Every coefficient profile from every
embedding, its permitted leg symmetries, and its coordinate weight
variants restricts to a scalar multiple of the same crossing rectangle.

Consequently the span mixing all these templates still forces leading
contraction at least one whenever
`sum_j log_2(m_j+2)=o(log_2 D)`. In particular, for quadratic word
lengths `2^O(sqrt(lambda))`, a union of `o(sqrt(lambda))` arbitrary label
relabelings remains obstructed. This is a sufficient uniform scope
condition; the exact product `P_J` gives the quantitative finite bound.

## 8. Exact replay and a necessary negative control

`check_word_restriction_barrier.py` checks every missing-pattern
witness and the complete term-22 identity from the pinned local table.
It checks all eight leg symmetries and all cross pairs in an explicit
balanced quadratic code, including its full affine extension. It then
compares the source expansion and the scalar oracle on every histogram
in the finite rectangle and its eight images. Histograms need not be
uniform. The harmonic fiber bound is also checked exhaustively on every
subset of binary words through length four, and the exact crossing
restriction is checked on a fiber from an unbalanced full word set.
Direct weighted source profiles verify the coordinate-weight extension.

A deliberate unbalanced example shows why the larger balanced-code
rectangle (9) requires its balance hypothesis. On three Boolean variables take

```math
x=v=t_1,\qquad y=t_1t_3,\qquad u=t_2+t_1t_2.
```

A linear functional can be zero on `t_1` and one on the other two
displayed independent words. The quartet then fits (9), but its
profile is nonzero although it is not a crossing tuple:

```math
P(z)=z^{11}(1+z)^{13}+z^{15}(1+z)^{13}.
```

The words `y,u` have weight one quarter of the length and violate
(8). Their cross pair projections omit patterns, as required for this
negative control.

Run:

```sh
python tools/check_word_restriction_barrier.py \
  --output evidence/word-restriction-barrier.json
```

The JSON records exact coverage and source hashes. The arbitrary-size
conclusion uses the proof above, not extrapolation from the replay.

Substantial OpenAI assistant support was used in the proof and exact replay.
