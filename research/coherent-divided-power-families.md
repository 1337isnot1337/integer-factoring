# Graded dimensions and protected jets of divided-power inverse systems

There is a polynomial-size formula for the complete physical height and
occupation distribution of the derivative space of a paired divided power.
It applies over every characteristic-two field and does not enumerate its
potentially exponential basis. For the divided square it gives

\[
(1,m,\tbinom m2-[m/2\text{ odd}],m,1),\qquad m\ge4\text{ even}.
\]

These are physical-space dimensions, not cut ranks. The actual two-cut
profiles still require a separate computation. A protected low-height
subspace does, however, determine an initial segment of both profiles
exactly. The present construction and finite controls establish no better
factoring exponent; the score remains **1/2**.

## 1. The family and its disjointness blocks

Put `m=2s`, pair mode `i` with `m-1-i`, and let

\[
g=\sum_{i=0}^{s-1}x_i x_{m-1-i},\qquad
f_r=\gamma_r(g)=\sum_{|R|=r}\prod_{i\in R}x_i x_{m-1-i},
\quad 0\le r\le s.
\]

Each pair has height `2m`, so `f_r` has occupation `2r` and height
`H=2mr`. Let `W_r` be the span of all its bit derivatives. Its compatible
coordinate retraction is the existing
[inverse-system construction](coherent-hot-compound-retractions.md).

Fix a set of `b` broken pairs and choose one surviving mode in each.
A derivative producing this pattern deletes the other mode in those
pairs and deletes both modes in `a` further pairs. Its output has
`k=r-b-a` complete surviving pairs. On the `n=s-b` unbroken pairs, the
coefficient matrix is exactly

\[
M_{A,C}=[A\cap C=\varnothing],\qquad |A|=a,\ |C|=k.
\]

Different broken patterns have disjoint monomial supports. For one fixed
pattern, distinct `k` have distinct output occupation. The resulting rank
contributions therefore add directly.

For `u=min(a,k)` and `v=max(a,k)`, complementation and transposition turn
this matrix into the inclusion matrix of `u`-subsets in `(n-v)`-subsets.
Wilson's formula gives its rank

\[
R(n;a,k)=\sum_{i=0}^{u}
 \left[\binom{n-v-i}{u-i}\text{ is odd}\right]
 \left(\binom ni-\binom n{i-1}\right).                 \tag{1}
\]

The imported rank theorem and primary attribution are already stated in
[the binary coefficient note](binary-coefficient-entropy.md).
A primary proof is Frankl,
[*Intersection Theorems and Mod p Rank of Inclusion Matrices*, Theorem 4.1](https://www.renyi.hu/~pfrankl/j62.pdf).
The hypothesis holds here: `a+k=r-b<=n`, hence
`u<=min(n-v,v)`. All matrices have binary entries, so their rank is the
same over every characteristic-two extension. The application to the
derivative space and the height distribution below is separate from the
earlier crossing-coefficient application.

## 2. Exact physical type costs without a basis

Define

\[
E_b(z)=[t^b]\prod_{i=0}^{s-1}
 (1+t(z^{2i+1}+z^{2(m-1-i)+1})).
\]

It counts choices of `b` broken pairs with their surviving heights. The
complete occupation-height enumerator is

\[
\boxed{F_r(y,z)=\sum_{b=0}^r E_b(z)y^b
 \sum_{a=0}^{r-b}R(s-b;a,r-b-a)
                 y^{2(r-b-a)}z^{2m(r-b-a)}.}         \tag{2}
\]

The [implementation](../tools/coherent_divided_powers.py) uses an
elementary-symmetric recurrence for the `E_b`, followed by integer
binomial coefficients and parity tests in (1). Its state count and bit
cost are polynomial in `m,r,H`: at most `O(rH)` coefficient cells,
`O(srH+r^3)` integer operations, and polynomial-size integers. This
computes the type partition, not a complete crossing table or a Smith
profile. The native probe explicitly enumerates its finite top form and
basis and pays those separate costs.

Specializing `y=z=1` gives the physical dimension. For `r=1` the Hilbert
function is `(1,m,1)`. For `r=2`, the middle disjointness block on the
complete pairs is `I+J` of order `s`. It has rank `s` for even `s` and
`s-1` for odd `s`. All remaining quadratic blocks are permutations.
This proves the divided-square formula at the start. At `r=s`, the top
form is one monomial and `W_r` is the full Fock space of dimension `2^m`.

More generally a homogeneous top form has transposed catalecticants at
complementary height and occupation. Thus

\[
[y^d z^h]F_r=[y^{2r-d}z^{H-h}]F_r.                 \tag{3}
\]

In particular the mean physical height is `H/2`. These identities are
useful checks, and do not identify the two crossing profiles.

## 3. A protected coordinate space and exact low Smith orders

Set `K=min(2m,(r+1)^2)`. Every occupation word of height below `K`
contains at most `r` modes and contains no complete mirror pair. Such a
word `A` is itself one derivative of `f_r`: delete the partners of its
occupied modes, together with all modes in any `r-|A|` additional pairs.
Exactly one term of `f_r` survives, and it becomes `A`.

Consequently `W_r` contains every coordinate state below height `K`.
Reduced elimination fixes these coordinates, while all other embedding
columns and retraction rows have height at least `K`.

For the regular one-particle matrix `B(t)`, each nonzero entry has
valuation at least both its row height and its column height. Every
compound entry therefore has valuation at least the maximum of its
total row and column heights. This persists under homogeneous physical
compression. Modulo `t^K`, any row or column of pair height at least `K`
vanishes. The remaining matrix is exactly the corresponding block of
the full Gaussian, at either cut. Adding the zero rows or columns does
not change finite Smith orders below `K`.

**Protected-jet theorem.** Both natural regular cut profiles of this
coordinate-retraction code agree with the full Gaussian below order `K`.
For dyadic `m`, the known one-particle Gaussian orders are
`1,3,...,4m-1`; hence the multiplicity at any order below `K` is the
number of partitions of that order into distinct positive odd parts.
The odd-profile input is the existing
[reflection reduction](reflection-dual-rail-compression.md), with its
full interval extension in [Pascal complements](schur-pascal-complement.md).
This theorem does not determine the unprotected orders or exclude a
better pressure elsewhere.

## 4. An explicit regular symplectic identity

The reflection family already has a
[symplectic structure](gaussian-complement-symplectic-relations.md).
Here is its direct form for the regular block used by the native engine.
For dyadic `m`, let `S` be the mirror permutation on `m` modes and
`Omega=diag(S,S)`. With `B=[[D,J],[J,D]]`,

\[
\boxed{B(t)\Omega B(t)^{\mathsf T}=t^{4m}\Omega.}    \tag{4}
\]

Indeed `J S J^T=t^(2m)S`. Write
`D=(1+t)t tensor_j [[1,z_j],[z_j,z_j]]`,
where `z_j=t^(2^(j+1))`. Each two-by-two factor preserves its alternating
form up to determinant `z_j(1+z_j)`. The product identity gives
`D S D^T=t^(2m)(1+t^(2m))S`. The cross terms
`D S J^T` and `J S D^T` agree entrywise and cancel in characteristic two.
Adding the diagonal terms proves (4).

This explains the complementary one-particle Smith orders. It does not
make a rectangular compression invariant, and does not justify discarding
the opposite cut. The supplied unequal-alphabet interface computes the
two distinct compound restrictions explicitly.

## 5. Exact controls and exploratory interfaces

The [checker](../tools/check_coherent_divided_powers.py) and
[receipt](../evidence/coherent-divided-powers-controls.json) verify 220
disjointness ranks by direct elimination; all twenty divided powers on
two through ten modes by literal derivatives, including 2,957 nonzero
derivative vectors; six larger enumerators through 2,048 modes; all 5,456
polynomial entries of (4) through 32 modes; and 28 one-particle Smith
orders through eight modes. The protected jets are checked on two saved
complete profile pairs.

The [unequal-alphabet interface](../tools/coherent_unequal_compounds.py)
is additionally checked against 780 literal ambient polynomial minors
and 76 separate fraction-free Smith orders. One test is a genuinely
nonorthogonal five-by-six alphabet, so the two cut orientations differ.
The equal-alphabet specialization reproduces both existing profiles.

```sh
make verify-coherent-divided-powers
python3 tools/search_coherent_divided_powers.py --modes 10 --power 2 --output /tmp/divided-power-m10.json
python3 tools/search_coherent_unequal_divided.py --modes 8 --left-power 2 --right-power 3 --output /tmp/divided-power-unequal.json
```

The discovery drivers retain complete profiles, maps, source hashes and
numerical fits. Unequal alphabets use both physical partitions and both
supporting weights; six retained rational points have exact integer
comparisons. A negative finite test does not establish a global optimum.
Larger occupation blocks can exceed the native implementation's limit
of 4,096 and are not silently treated as complete. These author controls
do not constitute independent review or a new Lean theorem.

## 6. Completed larger probes

The [retained native results](../evidence/coherent-divided-power-probes.json)
now include the ten-mode divided square with 66 physical letters and
both complete 4,356-order profiles. Its numerical balanced fit is
`rho approximately 0.7155596412`, compared with `0.7155813177` for the
eight-mode member. This small within-family gain remains above
`1/sqrt(2)` and does not improve the existing 96-letter construction.

Eight unequal paired codes also completed both profiles and continuous
supporting-weight searches. Their support parameter `R` must be below
one to pass the asymmetric criterion; it is not the balanced `rho`.

| Modes | Divided powers | Alphabet sizes | Best numerical R |
|---|---|---|---:|
| 4 | 1, 2 | 6, 16 | 1.034130641 |
| 6 | 1, 2 | 8, 28 | 1.034444645 |
| 6 | 1, 3 | 8, 64 | 1.034320701 |
| 8 | 1, 2 | 10, 46 | 1.034358972 |
| 8 | 1, 3 | 10, 120 | 1.034317806 |
| 10 | 1, 2 | 12, 66 | 1.034345401 |
| 10 | 1, 3 | 12, 232 | 1.034316851 |
| 8 | 2, 3 | 46, 120 | 1.009658920 |

Each retains six rational Pythagorean supporting weights and its exact
integer comparison at the declared rational scale. All 48 retained
points fail the strict threshold. The numerical minimizations do not
prove global exclusions, and those six rational points need not coincide
with the best continuous weight.

The [receipt checker](../tools/check_coherent_divided_probes.py) rebuilds
all maps, checks the complete native receipts, evaluates the fixed
numerical equations, and recomputes all 48 integer comparisons. Its
[replay receipt](../evidence/coherent-divided-probes-replay.json) covers
31,888 recorded complete orders across nine pairs. It explicitly does
not rerun the large native profiles. Reproduce that receipt with
`make verify-coherent-divided-probes`; use the discovery drivers above
to repeat a full native construction.
