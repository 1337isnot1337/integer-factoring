# Crossing pieces from noncommuting group words

## Result and scope

An ordered group-word deformation gives an exact characteristic-two crossing
family on the Heisenberg groups `H(F_(2^k))`. The label width is `lambda=3k`,
the unweighted family has `3k+1` pieces, and its entries and exact rank
factorizations have uniform `2^O(lambda)` construction algorithms. Its two
cuts are ordinary scalar matrices. No matrix-valued coefficient is being
treated as a scalar, and no cross-half commutation is assumed.

This is a new explicit construction to test, not a stronger factoring bound.
The canonical coefficient family is excluded by the all-`k` calculation below:
its constant coefficient alone has normalized width tending to

```math
\rho_0=\frac{7+\log_2 3}{12}>\frac1{\sqrt2}.
\tag{1}
```

The exact constant-coefficient rank calculation reduces a `q^6`-row matrix
to a `q^3`-row matrix, where `q=2^k`. This also permits exhaustive tests of
many different noncommuting orders that would be expensive as dense cuts.
Noncanonical orders and arbitrary recombinations of whole coefficients are
not excluded asymptotically by (1). An exhaustive finite central-shear search
reduces the constant-cut rank from 31 to 28 on 64 labels, but the complete
profile and its binary recombinations do not improve the factoring score.

The point of using a nonabelian group is precise: the row cut remains a sum
of group-convolution blocks, but the other cut is no longer the second
abelian convolution cut used by the binary conservation-kernel obstruction
in [PR27, at `69ee5a3`](https://github.com/1337isnot1337/integer-factoring/tree/69ee5a3aa8f575846a94be5c92a51bd8af8cb3a5).
That earlier theorem cannot simply be applied by replacing the label group.
The present calculation supplies a separate obstruction to the canonical
coefficient family, with a narrower scope than the entire orbit-kernel space.

## 1. An orbit kernel and its exact row blocks

Let `G` be a finite group of order `d`, with identity `e`. For any scalar
kernel `K:G×G→F`, define

```math
T_K(x,y,u,v)=
 [xy^{-1}=vu^{-1}]\,K(xy^{-1},y^{-1}u).
\tag{2}
```

Equivalently, for `S,a,y in G`, the entry with

```math
x=Sy,\qquad u=ya,\qquad v=Sya
\tag{3}
```

is `K(S,a)`, and all other entries vanish. Equation (3) is a bijective
parameterization of the support. The kernel

```math
K_E(S,a)=[a=e]
\tag{4}
```

gives exactly `E_d(x,y,u,v)=[x=v][y=u]`.

In the row cut `(x,y)|(u,v)`, arrange rows and columns by
`S=xy^{-1}=vu^{-1}`. In the block indexed by `S`, use `y,u` as free
coordinates. This proves the exact identity

```math
\operatorname{rank}_{row}T_K
 =\sum_{S\in G}\operatorname{rank}
       \bigl(K(S,y^{-1}u)\bigr)_{y,u\in G}.
\tag{5}
```

Every block is an ordinary group-convolution matrix. The column cut must
still be accounted for separately. Noncommutativity does not supply a cheap
second factorization by itself.

## 2. The noncommuting family

Put `q=2^k`, and define

```math
G=H(F_q)=F_q^3,\qquad
(a,b,c)(A,B,C)=(a+A,b+B,c+C+aB).
\tag{6}
```

Associativity follows by expanding the bilinear term. The inverse is
`(a,b,c+ab)`, and the center contains all `(0,0,c)`. In particular,

```math
(a,b,c)^2=(0,0,ab),\qquad
[(a,b,c),(A,B,C)]=(0,0,aB+Ab).
\tag{7}
```

Fix a binary basis `e_1,...,e_k` of `F_q`. The `3k` generators are

```math
X_i=(e_i,0,0),\quad Y_i=(0,e_i,0),\quad Z_i=(0,0,e_i).
```

Each has order two. In the canonical order `X_1,...,X_k,Y_1,...,Y_k,
Z_1,...,Z_k`, every group element has one binary word

```math
S=g_1^{s_1}\cdots g_{3k}^{s_{3k}},\qquad s_j\in\{0,1\}.
\tag{8}
```

For example, the final central coordinate of (8) is the binary sum of
the selected `Z` coordinates plus the product of the selected `X` and
`Y` coordinates. Thus finding this unique word is effective and does
not require searching the group.

Work in the group algebra `F_2[G][z]`, where `z` is central. Choose a
permutation `sigma` of the generators and positive integer weights `w_j`.
Define

```math
F_S(z)=\prod_{j\text{ in order }\sigma}
       \bigl(g_j+z^{w_j}(1+g_j)\bigr)^{s_j},
\qquad K_z(S,a)=[a]F_S(z).
\tag{9}
```

The products retain their specified order. In particular, factors involving
`X_i,Y_j` cannot be interchanged without proof. The word bits `s_j` always
come from (8), even when the factors in (9) use a different order.

At `z=1`, each selected factor is the group-algebra identity. Therefore

```math
F_S(1)=e,\qquad T_{K_1}=E_d,\qquad
\sum_{i=0}^{D} A_i=E_d,
\quad A_i=[z^i]T_{K_z},\quad D=\sum_jw_j.
\tag{10}
```

This is an exact scale-one crossing, including in characteristic two.
There is no trace normalization, group-order division, or interpolation
at elements of the coefficient field. When all weights are one, there
are at most `3k+1` pieces.

Each augmentation `sum_a K_z(S,a)` is the constant polynomial one:
every factor in (9) has augmentation one. Consequently any coefficient
functional `L` with `L(1) != 0` has a unit group-algebra element in every
row block. For a finite 2-group over characteristic two, augmentation-one
elements are units. One proof is induction on the group order: quotient
by a central involution, lift the inverse from the quotient, and use that
the kernel ideal generated by `1+g` has square zero. Hence

```math
L(1)\ne0\quad\Longrightarrow\quad
\operatorname{rank}_{row}L(T_{K_z})=d^2.
\tag{11}
```

This explains why reducing the corresponding column rank matters. It does
not give a lower bound for that column rank of an arbitrary recombination.

### Effectiveness and uniformity

For a fixed computable order rule and a fixed polynomial bound on `D`,
the family is effective. An irreducible degree-`k` polynomial over `F_2`
can be found by enumerating monic polynomials and possible proper factors,
in `2^O(k)` bit time. Field multiplication and group operations then use
polynomial time in `k`.

Compute (9) as an array of at most `d` group elements with polynomial
coefficients of degree at most `D`. Multiplying the `3k` binomial group
factors takes `2^O(lambda)` bit time. This computes a scalar piece entry
without assuming commuting generators. Complete cut matrices have `d^2`
rows and columns, so binary elimination and exact factors also take
`2^O(lambda)` bit time. Their ranks are at most `d^2`; each cut can
therefore be capped at `2lambda` bits in the usual crossing interface.

For a binary width not divisible by three, tensor the nearest smaller
family with a crossing on at most two tail bits. This changes the width
sum by at most eight and preserves polynomial piece count. Small widths
are handled directly. A fixed order rule is part of the construction;
an unbounded search through all `(3k)!` orders is not being charged as
`2^O(lambda)`.

These bounds establish the effective interface, not an improved width
inequality. The constructed matrices can consume exponential local space;
a claim of linear local workspace needs a separate blocking argument.

## 3. Exact rank of the canonical constant coefficient

Use the canonical order in (9), with any positive weights. The constant
coefficient has `K_0(S,a)=[a=S]`. Each row block in (5) is a permutation,
so its row rank is `d^2=q^6`.

The column rank has the exact formula

```math
\boxed{\quad
\operatorname{rank}_{col} A_0
 =q+(q-1)3^k.
\quad}
\tag{12}
```

Here the coefficient field is `F_2`. Extending it to `F_q` does not change
the rank, and will simplify the proof.

From `a=S=xy^{-1}`, the support conditions are

```math
u=yxy^{-1},\qquad v=x^2y^{-1}.
```

Because `x^2` is central, these are equivalently

```math
y=x^2v^{-1},\qquad u=v^{-1}xv.
\tag{13}
```

Write `x=(a,b,c)` and `v=(A,B,C)`. Equation (7) shows that `u` has
the same first two coordinates as `x`, and its central displacement
is `aB+bA`. Also `yv=x^2=(0,0,ab)`. Group the column matrix by the
central value `z=ab=yv`. Remove zero rows and columns, then remove
identical central translates. These operations leave its rank equal
to the sum of the ranks of the following `q` incidence matrices:

```math
M_z((a,b,t),(A,B))=[t=aB+bA],\qquad ab=z.
\tag{14}
```

This is an exact compression of the scalar cut.

### The block `z=0`

Here `a=0` or `b=0`. The row space is all functions of `A` plus all
functions of `B`. Their intersection is the constant functions. Thus

```math
\operatorname{rank}M_0=2q-1.
\tag{15}
```

### Each block `z != 0`

All `a,b` are nonzero and satisfy `ab=z`. Their slopes `a/b=a^2/z`
run through `F_q^*`, since squaring is a bijection of `F_q`. The row
space of (14) is therefore the span of all functions of `A+sB`, over
all nonzero slopes `s`.

Every function on `F_q` is represented by a polynomial of degree less
than `q`. It suffices to span

```math
(A+sB)^j,\qquad 0\le j<q,\quad s\in F_q^*.
```

In characteristic two, the coefficient of `A^(j-b) B^b` is nonzero
exactly when the set of one-bits of `b` is contained in that of `j`.
This follows by expanding
`(1+Z)^j=product_{i:j_i=1}(1+Z^(2^i))`. Thus there are
`2^(popcount(j))` possible monomials at degree `j`.

For `j<q-1`, their distinct exponents of `s` lie between zero and
`q-2`. Evaluation on `F_q^*` is injective on this polynomial span:
a nonzero polynomial of degree at most `q-2` cannot have `q-1`
distinct roots. Consequently all these monomials are spanned. At
`j=q-1`, all `q` exponents occur, but the functions `s^0` and
`s^(q-1)` coincide. The other `q-2` exponents are independent, so
this degree contributes exactly `q-1` dimensions.

Different total degrees are independent as functions on `F_q^2`,
because every individual variable degree is below `q`. Hence

```math
\operatorname{rank}M_z
 =\sum_{j=0}^{q-1}2^{\operatorname{popcount}(j)}-1
 =3^k-1\qquad(z\ne0).
\tag{16}
```

Combining (15) and (16) gives
`2q-1+(q-1)(3^k-1)=q+(q-1)3^k`, proving (12).

Since `lambda=3k`, equations (12) and `rank_row A_0=q^6` imply

```math
\frac{\log_2\operatorname{rank}_{row}A_0+
      \log_2\operatorname{rank}_{col}A_0}{4\lambda}
\longrightarrow \frac{7+\log_2 3}{12}.
```

This exceeds `1/sqrt(2)`. For example, `log_2 3>3/2`, so the limit
exceeds `17/24`, while `(17/24)^2=289/576>1/2`. Thus the canonical
coefficient family cannot satisfy the stronger half-rate width target,
even before examining its other pieces. The same conclusion covers
every positive weight schedule with canonical factor order, because
its constant coefficient is unchanged.

## 4. A small exact matrix for every factor order

For any permutation in (9), its constant word at an element
`S=(a,b,c)` has the form

```math
\tau(S)=(a,b,c+f(a,b)),
```

where `f:F_q×F_q→F_q` is an explicitly computed bilinear map over
`F_2`. The central generators commute with every factor, so only the
relative order of the `2k` noncentral generators matters.

For the constant coefficient `K_0(S,a)=[a=tau(S)]`, the column cut
has the same rank as a `q^3`-by-`q^3` binary matrix. Its rows are
`(a,b,delta)` and columns are `(A,B,epsilon)`, with entry one exactly
when

```math
\begin{aligned}
\delta&=aB+Ab+f(a+A,b+B),\\
\epsilon&=ab+AB+f(a+A,b+B).
\end{aligned}
\tag{17}
```

To verify this, write `x=(a,b,c)` and `y=(A,B,C)`, substitute the
formula for `tau(xy^{-1})` in (3), and put
`u=(a,b,c+delta)`, `v=(A,B,C+epsilon)`. Expanding (6) gives (17).
The remaining `c,C` coordinates duplicate rows and columns. Every
other row or column is zero, which proves the rank statement.

Formula (17) is useful for a broader candidate: **any** function
`f:F_q^2->F_q` gives an invertible central shear `tau`. The derivation
of (17) did not use bilinearity. One can take
the canonical binary word of `tau(S)` as the control bits in (9).
The endpoint (10) remains valid; if `f` has a uniform `2^O(lambda)`
entry algorithm, the effectivity bounds remain valid too. A fixed
finite table is also a legitimate finite construction. No assumption
that `tau` is a group automorphism or that `f(0,0)=0` is needed. The new
constant coefficient has exactly the matrix (17). This extends the
order family; it does not assert that one of these shears improves
the full crossing ranks.

### Exact symmetries of the compressed matrix

Write `C_f` for (17), and put `Q(a,b)=ab`. For all functions `f`,
all `t=(t_1,t_2) in F_q^2`, and all `c in F_q`,

```math
\operatorname{rank}C_{f+t_1b+t_2a+c}
 =\operatorname{rank}C_f,\qquad C_{f+Q}=C_f^T.
\tag{18}
```

For the first equality, translate both horizontal row and column
coordinates by `t`, and both central coordinates by `c`. The two
expressions `aB+Ab` and `ab+AB` each gain
`t_1(b+B)+t_2(a+A)`. This is exactly the added function at
`(a+A,b+B)`. The second equality follows because
`Q(a+A,b+B)=ab+AB+aB+Ab`: adding it interchanges the two prescribed
central outputs in (17). These are identities of actual matrices
up to the displayed permutations, not just coincident finite ranks.

In particular, every search can normalize
`f(0,0)=f(1,0)=f(0,1)=0`, then identify `f` with `f+Q`.
The normalization uses `F_q`-affine additions; it does not permit
discarding arbitrary `F_2`-linear functions without proof.

For `q=4`, a Boolean quadratic function from four input bits to two
output bits has 22 coefficients. The three normalized function
values remove six bits, and `f` versus `f+Q` removes one more.
Thus only `32768` normalized cases remain. One explicit enumeration
uses the two free `F_4` linear coefficients of `a_1,b_1` and the six
`F_4` quadratic coefficients; set the low bit of the coefficient of
`a_0b_0` to zero using `f` versus `f+Q`. Section 6 uses this exact
finite reduction to find and independently replay the rank-28 cases.

## 5. Exact profiles and bounded recombination tests

The checker constructs the group algebra polynomial by dynamic
programming and separately by its subset expansion. It verifies the
group operations, both the row-block formula and actual dense cuts,
and the complete crossing sum. The compressed constant matrices are
checked against dense cuts at the bounded widths where both are used.

For canonical unit weights, the exact profiles are:

| `d` | Row ranks in increasing coefficient order | Column ranks |
| ---: | --- | --- |
| 8 | `64,32,10,1` | `5,19,32,64` |
| 64 | `4096,2592,1388,518,133,20,1` | `31,226,706,1379,2164,2592,4096` |

All nonzero `F_2` recombinations of these coefficient spans were
ranked. The spans have dimensions four and seven. The smallest
maximum rank product in an exact crossing decomposition within each
span is respectively `608` and `979928`, attained by the original
coefficient decomposition. This finite optimality is checked by
sorting the exact integer rank products and testing when their
coefficient masks first span the crossing mask. It concerns the
maximum product, not the minimum sum of square roots of products.
Extension-field recombinations are not included in this exhaustive
binary search.

All noncentral factor orders were also screened through `k=3`:

| `k` | Labels | Orders | Distinct central corrections | Minimum constant-cut rank |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 8 | 2 | 2 | 5 |
| 2 | 64 | 24 | 14 | 31 |
| 3 | 512 | 720 | 230 | 197 |

In these cases the minimum is attained by the canonical and reverse
orders. These finite facts do not prove that they minimize the rank
at all widths or among all bilinear central shears.

## 6. A quadratic shear lowers the constant-cut rank from 31 to 28

The broader construction in Section 4 does improve one exact finite
rank. Let `F_4=F_2[t]/(t^2+t+1)`, and write
`a=a_0+a_1t`, `b=b_0+b_1t`, with binary coordinates. Take

```math
f(a,b)=a_1+a_0a_1+b_0b_1
       +t(a_0b_0+a_0b_1+a_1b_0).
\tag{19}
```

Then `rank C_f=28`, compared with the canonical rank 31. Among
**all functions** `F_4^2->F_4`, the minimum is exactly 28 in the
exhaustive finite screen. The gauge normalization reduces the
`2^32` functions to `2^25=33554432` representatives. The complete
native run used exact elimination for every representative and
found exactly 12 minimizers, all Boolean quadratic. Larger fields
and other crossing constructions are outside this finite result.

The 32768 normalized Boolean quadratics were also completely
replayed by independent Python matrix construction and elimination,
with the same histogram and all 12 minimizers. For the larger
arbitrary-function screen, the independent Python checker replays
384 widely separated lookup cases and every minimizing function.
It does not independently replay all 33554432 cases. The full
native record, its exact histogram, command, source hash, and
the twelve minimizing tables are saved in
[`noncommutative-crossing-lookup-complete.json`](../evidence/noncommutative-crossing-lookup-complete.json).
Its ordered rank digest is `1415671051223566181`.

For reproduction, the a-major table of (19), with field elements
encoded as `0,1,t,1+t`, is

```text
0,0,0,1, 0,2,2,1, 1,3,1,2, 0,0,2,3.
```

The complete crossing family obtained by controlling the canonical
words with `tau(S)` has the exact profile

```text
row:    4096,2592,1388,518,133,20,1
column: 28,288,813,1519,2168,2592,4096.
```

Thus its other coefficients become more expensive: its maximum
rank product is `1128444`, greater than the canonical `979928`.
The rank-28 improvement is a local construction result, not an
improvement of the whole decomposition or a new factoring score.

All 127 nonzero binary recombinations of these seven coefficients
were also ranked. Their minimum possible bottleneck remains
`1128444`. Among recombinations whose constant coefficient is
selected, the minimum column rank is still 28. Their minimum
rank sum lies in the certified interval
`[2696.508822194921,2696.508822194923]`; the upper witness selects
coefficient masks `62,1,64`. The interval follows from shortest
paths in the seven-dimensional binary coefficient space, with
each square-root edge weight bounded above and below by integer
multiples of `10^-12`. Thus arbitrary binary regrouping does not
turn the improved constant piece into a successful finite seed.

The construction also cannot be extrapolated by retaining the same
short field polynomial in a larger field. Rewriting (19) using
`t^2+t+1=0` gives

```math
f(a,b)=a^3+b^3+t^2ab+(ab)^2
       +t^2a+ta^2+tb+t^2b^2.
```

Over `F_16`, the same polynomial has constant-cut rank 3572 for
either embedded root `t`, rather than `28^2`. The canonical rank
there is 1231. This is a negative finite result for that particular
extension, and it prevents assuming a tensor-product rank law
between different field sizes.

The native implementation is
[`explore_group_correction_rank.cpp`](../tools/explore_group_correction_rank.cpp).
The independent checker is
[`check_noncommutative_crossing_native.py`](../tools/check_noncommutative_crossing_native.py).
Its complete quadratic replay has digest `7695949591803970661`;
[`noncommutative-crossing-native-audit.json`](../evidence/noncommutative-crossing-native-audit.json)
records the full histogram, all minimizing indices, and additional
widely separated checks of the larger native searches. Capped
native ranks are lower bounds once their stopping threshold is
reached; only runs with `--exact` provide full rank histograms.

## 7. Reproduction

Run the complete recorded suite with:

```sh
python3 tools/check_noncommutative_crossing_orbit.py \
  --degree-two --span-degree-two --order-degree 3 \
  --output evidence/noncommutative-crossing-orbit.json
```

The default command runs the smaller construction checks. The
ordinary arguments above establish the all-size endpoint, effective
construction, rank compression, and canonical obstruction. Finite
tests support only their stated cases. The work used substantial
OpenAI assistant support; the factoring evaluator and transfer
framework retain the repository's existing source attribution.

For the independent full quadratic replay, compile the native
screen and run:

```sh
c++ -O3 -std=c++17 tools/explore_group_correction_rank.cpp \
  -o /tmp/intfact-group-correction-screen
python3 tools/check_noncommutative_crossing_native.py \
  --binary /tmp/intfact-group-correction-screen --full \
  --lookup-summary evidence/noncommutative-crossing-lookup-complete.json \
  --output evidence/noncommutative-crossing-native-audit.json
```

To repeat the full arbitrary-function enumeration:

```sh
/tmp/intfact-group-correction-screen --k 2 --family lookup \
  --normalized --count 33554432 --exact
```
