# Quotient symbols for the modular partition crossing

## Result and scope

The normalized modular construction can be computed entirely in the
regular-partition quotient, in both its state and operator-parameter
arguments. This fact was already proved in
[PR #22, Section 3](https://github.com/1337isnot1337/integer-factoring/blob/d94e4a99fec453768d69ae083038ba65cf31ae9e/research/past-half-pieces.md#3-reducing-the-operator-labels-is-valid-as-well).
Section 1 restates it with a normal-symbol proof as a foundation.

In characteristic two, the original source dual becomes an exterior
Boolean-zeta substitution. Its full angle polynomial has an explicit
determinant formula with monomial entries. This gives a direct exact
constructor, a transpose identity, and the vanishing of every odd angle
and odd crossing shift. It avoids expanding the full integral dual
polynomial before evaluating a coefficient.

These statements do not improve the factoring exponent. They give an
explicit description of the individual angle factors and an independent
way to inspect their algebra. The earlier
[complete-profile Pfaffian formula in PR #22](https://github.com/1337isnot1337/integer-factoring/blob/d94e4a99fec453768d69ae083038ba65cf31ae9e/research/modular-pfaffian.md)
already bypasses the integral expansion and evaluates all complete pieces
for four labels in `O(h^4)` bit operations. Its cost is stronger than the
subset determinant algorithm used here. The new formula is useful for
factor identities and independent verification; it is not a new fastest
complete-piece algorithm. In particular, the filtered source dual used here is
different from the homogeneous inverse-Gram dual in Section 6 of
[modular-crossing.md](modular-crossing.md).

The starting integral operators and dual polynomial are Anthropic's
Appendix B construction at revision
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The quotient and normalization are established in the existing
[modular note](modular-crossing.md). The derivations and checker here
used OpenAI assistant support; no new Lean work is claimed.

## 1. The operator parameter also factors through the quotient

Fix a prime `p` and use

```math
W_p=\mathbf F_p[t_1,t_2,\ldots]/
       (t_{pk}:k\ge1, t_k^p:k\ge1).
```

The source coefficient operator for a linear parameter is

```math
t_m[r]=\binom{r-1}{m-1}\nabla_{-r}.
```

For a positive mode `k`, the creation and annihilation coefficients are
respectively

```math
\binom{k-1}{m-1}M_{t_k},\qquad
(-1)^{m-1}\binom{k+m-1}{m-1}k\partial_{t_k}.       \tag{1}
```

If `p|m` and `p∤k`, both binomial coefficients vanish modulo `p`.
The last base-`p` digit of `m-1` is `p-1`, whereas the corresponding
upper digit is `k-1 mod p`, at most `p-2`. Lucas' identity proves the
claim. If `p|k`, creation is zero in the quotient and the weighted
derivative is zero. Thus every coefficient of the linear field with
parameter `t_m`, for `p|m`, acts as zero on `W_p`.

For a monomial parameter, the source operator is the normally ordered
product of the linear fields (1). Represent a normal symbol using
commuting creation and annihilation symbols, putting all actual
multiplications before all actual derivatives when it is evaluated.
In characteristic `p`, the `p`th power of a linear symbol is the sum
of the `p`th powers of its individual terms. Each such term contains
either `M_(t_k)^p=0` or `(k partial_(t_k))^p=0` on `W_p`.
Consequently a parameter factor `t_m^p` also acts as zero.

These remain zero after multiplication by any other normal symbol:
the resulting normal monomial still contains the vanishing creation
or derivative factor. Every fixed state-to-state coefficient involves
only finitely many modes, so no multiplication of undefined two-sided
series is needed. This proves

```math
a\in I_p\quad\Longrightarrow\quad a[r]|_{W_p}=0
\quad\text{for every mode }r.                         \tag{2}
```

The angle operation preserves this conclusion. The ideal is homogeneous,
and an angle multiplies each homogeneous component by a scalar
binomial coefficient. Thus both ordinary and angle operator parameters
may be reduced modulo `I_p` before computation.

## 2. The source dual is an involution of the quotient

The integral source substitution is

```math
\Theta(t_m)=(-1)^{m-1}\sum_{j=1}^m\binom mjt_j,
\qquad a^\mu=\Theta(t^\mu).                           \tag{3}
```

It preserves `I_p` modulo `p`. For `m` divisible by `p`, only terms
with `j` divisible by `p` survive, and a `p`th power maps to a sum of
`p`th powers. The binomial inversion identity gives `Theta^2=Id`.
Hence `Theta` descends to an automorphism of `W_p`, preserving ordinary
monomial length and the filtration by weighted degree.

For `p=2`, every surviving part is odd. Write it as `m=2a+1`.
Lucas' identity turns (3) into the particularly simple formula

```math
\Theta(t_{2a+1})=\sum_{b\subseteq a}t_{2b+1},          \tag{4}
```

where `b subseteq a` means that every set binary digit of `b` is also
set in `a`. The product is taken in the square-zero commutative algebra
on the odd modes. Repeated variables vanish.

Equation (4) is a Boolean-zeta transform. Applying it twice leaves the
coefficient of `t_(2c+1)` equal to the number of intermediate subsets
between `c` and `a`, modulo two. This is one precisely when `c=a`,
giving an independent proof of the involution in characteristic two.

## 3. One determinant gives every angle coefficient

Let `mu` be a regular characteristic-two partition of weight `h` and
length `L`. Its parts are `m_j=2a_j+1`. Let `sigma` and `pi` be regular
input and output partitions. Define the polynomial

```math
\mathcal A_{\mu;\pi,\sigma}(Z)
 =\sum_{i=0}^{2h}
    [t^\pi](a^\mu)^{\langle i\rangle}t^\sigma\ Z^i,   \tag{5}
```

where the mode is the one fixed by the input and output weights.
The notation in (5) suppresses that uniquely determined mode; it does
not sum different modes. The angle convention is the source convention
`a^(angle i)=sum_w binom(h+w,i) a_w`.

Fix an unaffected subset `eta` of the common parts of `pi` and `sigma`,
of cardinality

```math
|\eta|=(|\pi|+|\sigma|-L)/2.                          \tag{6}
```

If this is not a feasible nonnegative integer, the coefficient is zero.
Make one column for each output part outside `eta` (a creation), and
one for each input part outside `eta` (an annihilation). For a column
with mode `k=2c+1`, put

```math
e(a,c)=
\begin{cases}
 2(a\mathbin{\&}c),&\text{creation},\\
 2(a\mathbin{\&}\mathord{\sim}c),&\text{annihilation}.
\end{cases}                                         \tag{7}
```

The complement in (7) is masked by the finite integer `a`; it requires
no infinite binary word. Let `M_eta(Z)` be the `L` by `L` matrix with
entries `Z^(e(a_j,c))`. Then

```math
\boxed{
\mathcal A_{\mu;\pi,\sigma}(Z)
  =(1+Z)^{h+L}\sum_{\eta}\det M_\eta(Z)
  \quad\text{in }\mathbf F_2[Z].}                    \tag{8}
```

To prove (8), first sum all angle indices in (5). A homogeneous
parameter term of weight `w` receives the multiplier `(1+Z)^(h+w)`.
For one odd part, combine (4) with the raw creation and annihilation
coefficients. On odd modes these coefficients are respectively
`[b subseteq c]` and `[b & c = 0]`. Therefore the two resulting sums are

```math
\begin{aligned}
\sum_{b\subseteq a\cap c}(1+Z)^{2b+1}
   &=(1+Z)Z^{2(a\mathbin{\&}c)},\\
\sum_{\substack{b\subseteq a\\b\cap c=\varnothing}}
        (1+Z)^{2b+1}
   &=(1+Z)Z^{2(a\mathbin{\&}\mathord{\sim}c)}.
\end{aligned}                                        \tag{9}
```

Each identity follows by multiplying `1+(1+Z)^(2^j)=Z^(2^j)` over
the available binary digits. The outside angle contribution is
`(1+Z)^h`, and the `L` rows contribute another `(1+Z)^L`.
The normal-order coefficient sums all assignments of operator parts to
the selected creation and annihilation modes. This is a permanent;
over `F_2` it equals the determinant. The unaffected parts are exactly
the choices in (6), completing the proof.

No square root of a matrix, numerical rank, or unproved field inversion
enters this construction.

## 4. Degree, transpose, and parity consequences

For every entry in row `j` of (7), `0<=e<=2a_j=m_j-1`.
Thus the determinant sum in (8) has degree at most `h-L`, and the
complete polynomial has degree at most `2h`, as the source requires.

Interchanging input and output exchanges creation and annihilation in
every column. Every exponent `e` becomes `2a_j-e`. Factoring a monomial
from each row and using the reciprocal symmetry of `(1+Z)^(h+L)` gives

```math
\mathcal A_{\mu;\sigma,\pi}(Z)
 =Z^{2h}\mathcal A_{\mu;\pi,\sigma}(Z^{-1}).          \tag{10}
```

In particular, if `O_mu(pi,sigma)` denotes the raw monomial-parameter
operator coefficient, then

```math
[Z^{2h}]\mathcal A_{\mu;\pi,\sigma}
   =O_\mu(\pi,\sigma),\qquad
[Z^0]\mathcal A_{\mu;\pi,\sigma}
   =O_\mu(\sigma,\pi).                               \tag{11}
```

The first endpoint also follows because the highest-weight part of
`Theta(t^mu)` is `t^mu`. The second is its exact transpose. All regular
characteristic-two Fock norms are one, so no diagonal norm factors
are missing from (10)--(11).

Every exponent in (7) is even, and `h+L` is even because all parts
are odd. Thus (8) has only even powers of `Z`: every odd angle is zero.

For a crossing whose four labels all have weight `h`, a nonzero first
operator application has intermediate monomial length congruent to
`L+h=0 mod 2`. Its energy is therefore even. In the row order that
energy is `s+i`. Since `i` is even, a nonzero piece must have `s` even.
**Every odd shift piece in this normalized characteristic-two family
is zero**, for arbitrary `h`, not only powers of two. This reduces
constant factors in the explicit construction; it does not change its
leading contraction rate.

## 5. Polynomial-time queries for the individual factors

The determinant sum in (8) is itself one Pfaffian. This gives polynomial
time for an individual angle factor, with no unaffected-subset enumeration.
The earlier PR22 formula cited above computes the **complete** four-label
piece profile; the matrix below computes a single normal-ordered factor
between two states, which may have different weights.

Form a symmetric zero-diagonal matrix `B_mu,pi,sigma(Z)` with three blocks:
one vertex `O_j` for each operator part `m_j`, one `U_k` for each output
part `k`, and one `Y_k` for each input part. Set

```math
B(O_j,U_{2c+1})=Z^{2(a_j\mathbin{\&}c)},\qquad
B(O_j,Y_{2c+1})=Z^{2(a_j\mathbin{\&}\mathord{\sim}c)},
\qquad B(U_k,Y_l)=[k=l].                              \tag{12}
```

All other entries, including every within-block entry, are zero. If the
matrix order is odd its Pfaffian is zero. For even order, its perfect
matchings pair each operator vertex with one state vertex. Every remaining
output vertex must pair with its equal input part. These equal pairs are
exactly an unaffected subset `eta` from (6). After choosing them, the
remaining matching sum is `det M_eta`, because the signs disappear in
characteristic two. Thus

```math
\mathcal A_{\mu;\pi,\sigma}(Z)
  =(1+Z)^{h+L}\operatorname{Pf}B_{\mu,\pi,\sigma}(Z). \tag{13}
```

For states of energy at most `4h`, this matrix has at most `5sqrt(h)`
vertices. A nonzero matching uses each operator vertex exactly once, so
the Pfaffian degree is at most `sum_j(m_j-1)=h-L`. More generally, a
determinant minor has degree at most `2(h-L)`: each operator vertex can
occur once as a row and once as a column, and every nonconstant entry is
incident to such a vertex.

Compute the determinant by fraction-free Bareiss elimination over
`F_2[Z]`, with row pivots. Its exact divisions follow from the determinant
minor identity, and the stored entries are minors of the pivoted matrix.
The degree bound just proved applies to every such minor. Each numerator
has degree at most `4h`; schoolbook polynomial multiplication and exact
division therefore take `O(h^2)` bit operations. There are
`O(h^(3/2))` updates. In particular, `O(h^4)` bounds the complete
determinant calculation and matrix construction.

For an alternating matrix, `det B=(Pf B)^2`. In characteristic two,
squaring a polynomial just doubles its exponents. Check that all odd
determinant coefficients are zero and read its even coefficients to
obtain `Pf B`. Finally multiply by `(1+Z)^(h+L)`. This proves the same
`O(h^4)` bit bound for the full angle polynomial (13).

The raw monomial-parameter factors are its highest coefficients by (11).
Consequently both displayed factor tables of the normalized construction
have polynomial-time **entry queries** when their state indices are
supplied. This does not make the number of channels polynomial or give
minimum-rank factorizations; the width bounds and their unresolved
improvement remain unchanged.

## 6. Independent checks and complete-piece comparisons

The determinant of monomials in (8) can be computed by subset dynamic
programming in `poly(h) L 2^L` bit operations, representing a polynomial
by its coefficient bits. The number of possible unaffected subsets is
at most exponential in the number of parts of the two states. A
regular partition of energy `E` has at most `sqrt(E)` parts, since
the sum of the first `k` odd positive integers is `k^2`.

Consequently, for labels of weight `h` and states through energy `4h`,
one complete angle polynomial takes `2^O(sqrt(h)) poly(h)` bit
operations. Under the existing quadratic-energy label schedule this
is `2^O(lambda)`. This direct algorithm keeps the original effective
family bound while bypassing the ordinary integral-partition expansion.
It supplies an independent finite-checking path for the polynomial-time
Pfaffian method, rather than the best bound on these local queries.

Run:

```sh
python3 tools/check_homogeneous_central_symbols.py --central 16 \
  --output evidence/homogeneous-central-symbols.json
python3 tools/check_angle_pfaffian.py --output evidence/angle-pfaffian.json
```

The independent reference route uses the original integral
`dual_polynomial` and `operator` functions, including the monomials that
the new constructor eliminates. It checks the operator-parameter ideal
in characteristics `2,3,5`, the complete angle polynomial, both endpoints,
the reciprocal identity, and the original filtered dual itself. Complete
central pieces at `h=6,8` are also compared entry by entry with the prior
piece constructor.

The angle-Pfaffian checker separately compares fraction-free determinant
elimination with a perfect-matching recursion on 70 alternating polynomial
matrices, and compares 1,200 state/operator triples with the unaffected-set
determinant formula. Its small-weight cases also use the original integral
operator expansion. Exact-division, alternating-matrix, and omitted
unchanged-state-edge controls reject invalid substitutions. Its evidence is
[angle-pfaffian.json](../evidence/angle-pfaffian.json).

The filtered-source central ranks at `h=16` are `(15,14)` on its five
labels, reproducing the earlier PR #22 result. The homogeneous
inverse-Gram family has ranks `(12,12)` there,
as recorded in the modular note. They are different valid families;
the direct symbolic formula does not turn one into the other and does
not supply a smaller rank by itself. The finite output is
[homogeneous-central-symbols.json](../evidence/homogeneous-central-symbols.json).

The remaining constructive question is whether a different choice of
duals or recombination of complete pieces uses this explicit binary
structure to prove a better two-cut rank profile. Formula (8), or a low
rank in one finite example, does not answer that question.
