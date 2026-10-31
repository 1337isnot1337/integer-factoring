# Growing quartic words and removable endpoint poles in reflection coordinates

## Result and scope

This note gives explicit nonquadratic perturbations of the reflection tensor.
The number of perturbed modes can grow with its interval M. It includes both
literal crossing endpoints and endpoints that become the crossing after the
paired physical word types are imposed. Some latter perturbations have local
poles at w=1; the poles cancel in every retained coefficient before endpoint
evaluation.

The construction and all degree costs are exact. No improved factoring bound
is claimed. Full-profile numerical screens are separate finite experiments.
The earlier [positive-energy endpoint argument](positive-energy-quartic-endpoints.md)
is adapted here to the compressed reflection coordinates and to ordered
creation/annihilation words.

## 1. Literal occupation gates

Let M be a power of two, M at least 4, and put m=M/2. The four physical
letters have m binary occupation modes. On a mode i, C_i adds an occupation
when it is absent and A_i removes it when it is present. Take four distinct
physical modes and an ordered product Q of their chosen C or A operators.
Then Q^2=0. For a scalar rational function g, the gate

\[
 U=1+gQ                                                     \tag{1}
\]

has inverse itself in characteristic two. A word of such gates acts directly
on the complete coefficient tensor. Different gates may share modes; their
order is retained.

If every g is a Laurent polynomial with g(1)=0, the word equals the identity
at the endpoint. The resulting tensor is therefore a literal strict crossing
on the complete physical alphabet. The supported colored two-power example is

\[
 g(w)=\gamma(w^a+w^b),\qquad a\ne b,
 \qquad \gamma\in\mathbb F_4^\times.                         \tag{2}
\]

These are genuinely quartic families. For example, for four modes on one
physical letter, U=1+g C_iC_jC_kC_l satisfies

\[
 U A_i U^{-1}=A_i+g C_jC_kC_l.                              \tag{3}
\]

The second term sends the vacuum to a three-particle state. It cannot be a
linear combination of creation and annihilation operators. Thus U does not
normalize the Clifford vector space. Multiplying an invertible Gaussian cut
by U leaves that class. Disjoint quartets can be repeated on a number of modes
proportional to M, so this conclusion is not restricted to a fixed small group
of modes.

There is also a direct coefficient check. An all-creation quartet leaves the
vacuum and every two-particle coefficient unchanged, while adding g times the
vacuum coefficient to its four-particle coefficient. It violates the quadratic
Gaussian Pfaffian relation whenever g is nonzero.

## 2. A paired occupation charge

For a physical occupation state define

\[
 E=\sum_{i=0}^{m-1}i(x_i+y_i-u_i-v_i).                      \tag{4}
\]

In a tensor power with the same full-letter word type on X and V, and the
same full-letter word type on Y and U, the total charge is zero. The two paired
types can be different. No independent-bit assumption is needed.

Every occupation gate Q has a fixed charge increment Delta: add the mode's
coefficient from (4) for a creation and subtract it for an annihilation. Each
nonzero application of Q changes the charge by exactly Delta. This remains
true for ordered words whose gates share modes.

If every gate with nonzero endpoint coefficient has Delta>0 and all
coefficients are Laurent polynomials, a nonempty endpoint perturbation has
positive charge. It vanishes after the paired word-type restriction. Thus the
typed endpoint remains the crossing even when g(1) is nonzero. In particular,
monomial coefficients gamma w^a are permitted for positive-charge gates.

The local endpoint in this case has the form crossing plus terms of strictly
positive charge. Both cost flattenings are invertible at w=1. At the zero cut
sort by the potentials sum i(x_i+y_i) and sum i(v_i+u_i); at the infinity cut
sort by sum i(x_i-u_i) and sum i(v_i-y_i). The endpoint is triangular with
identity diagonal blocks after the crossing permutation.

One growing example consists of creation quartets

\[
 x_i y_i u_{m-1-i}v_{m-1-i},\qquad m/2\le i<m.             \tag{5}
\]

Their increments are 4i-2m+2>0. The quartets have disjoint mode supports and
their number is proportional to M.

## 3. An exact charge budget for the reflection Gaussian

Put chi=1+w and L=m-1. Every quadratic edge in the reflection kernel obeys

\[
 E(\text{edge})\ge-L\operatorname{ord}_{\chi}(\text{edge}). \tag{6}
\]

Here the valuation is at w=1, and Laurent monomials in w are units. The
kernel in [`reflection_dual_rail_gaussian.py`](../tools/reflection_dual_rail_gaussian.py)
has the following edge classes:

| Edge class | Most negative charge | chi-order |
|---|---:|---:|
| Constant matched edges v_i x_i and y_i u_i | 0 | 0 |
| External nonconstant edges | -2(m-1) | 2 |
| Internal reversal edges on a physical letter | -(m-1) | 2M |

Every external nonconstant coefficient has the factor 1+w^2=chi^2.
An internal coefficient is 1+w^(2M) or 1+w^(-2M). Since 2M is a power of
two, its chi-order is exactly 2M. The original common scalar is a monomial
in w, so it has charge and chi-order zero.

A Gaussian coefficient is a sum of matching terms. Charges and chi-orders
add within each term. Therefore each individual matching term satisfies (6).
Cancellation between matching terms can only raise its final chi-order.

## 4. Removable poles for mixed quartic words

Consider an ordered word of gates with

\[
 g_j(w)=\frac{\gamma_j w^{a_j}}{\chi^{s_j}},\qquad
 s_j\ge0,qquad
 \Delta_j>L s_j.                                           \tag{7}
\]

**Typed polynomiality and endpoint theorem.** Every tensor power of this
word applied to the reflection tensor, restricted to the paired full-letter
types, is a Laurent polynomial in w. Its value at w=1 is the original
crossing. The statement permits arbitrary creation/annihilation kinds and
overlap between different gates.

To prove the theorem, expand a retained coefficient into selections of gates
and original Gaussian matching terms. A nonzero selected operator word has
charge increment sum Delta_j, counted with multiplicity across the tensor
copies. Since the target charge is zero, the original matching terms have
total charge minus this sum. Equation (6) forces their chi-order V to satisfy

\[
 LV\ge\sum_j\Delta_j>L\sum_j s_j                         \tag{8}
\]

for every nonempty selection. Since L>0 and the orders are integers,
V is at least 1+sum s_j. All endpoint denominators cancel term by term, and
every term using a gate still vanishes at the endpoint. The empty selection
is the original crossing. This proves the theorem.

The condition is sufficient and deliberately strict. Equality Delta=L s
does not guarantee endpoint vanishing. A numerator with additional factors
of chi can permit broader conditions, but those are not used in (7).

The same charge proof works for occupation monomials of other fixed degrees.
The executable family and controls in this note use degree four.

## 5. Degree, scalar, and construction costs

The original polynomial reflection tensor has degree at most D_base=3M^2/2.
For a word (7), let S=sum s_j and use the common denominator chi^S. One
step multiplies the previous numerator by chi^(s_j) and adds gamma_j w^(a_j)
times an occupation transition. Its minimum exponent drops by at most
max(0,-a_j), and its maximum exponent rises by at most max(s_j,a_j).

After multiplication by w^C, where

\[
 C\le\sum_j\max(0,-a_j),
\]

the local numerator is polynomial of degree at most

\[
 D_{\rm num}\le D_{\rm base}
       +\sum_j\bigl(\max(0,-a_j)+\max(s_j,a_j)\bigr).       \tag{9}
\]

In an r-fold tensor power, the retained numerator is divisible by chi^(Sr)
by Section 4. After exact division, its degree is at most r(D_num-S). Thus
the effective cleared degree is bounded by

\[
 \boxed{D_{\rm eff}\le\frac{3M^2}{2}
   +\sum_j\bigl(\max(0,-a_j)+\max(0,a_j-s_j)\bigr).}       \tag{10}
\]

A numerator with powers ranging from a_j to b_j has the corresponding bound
max(0,-a_j)+max(0,b_j-s_j). For a purely strict polynomial word, set s_j=0.

For O(M) gates with powers O(M) and endpoint-pole orders O(1), all these
degrees are O(M^2). The complete local coefficient table has 2^(2M) entries.
Direct occupation transitions and exact polynomial division retain a
2^O(M) construction bound. These costs do not treat a growing set of gates
or an endpoint pole as free.

If N is the common polynomial numerator after the w clearing, the original
local tensor is w^(-C)N/chi^S. Consequently its full cut orders are obtained
from N by the shifts

\[
 \sigma_0=\sigma_0(N)-C,\qquad
 \sigma_\infty=\sigma_\infty(N)+C+S.                       \tag{11}
\]

Both shifts are included in the profile interface. The opposite C terms
cancel in the product of the two partition functions. The S term remains
and is justified by the exact typed division, not by discarding a pole.

## 6. Rank and exact controls

Strict words and polynomial positive-charge words have full local rank by
their endpoint arguments. Removable-pole words need not have a regular local
endpoint. The construction theorem does not assert full local rank for every
such word. The checker supplies a nonzero determinant specialization over
F16 before using a full-rank Smith oracle. An unsuccessful finite-field
specialization is an unresolved rank control, not a proof of rank deficiency.

The direct constructor is
[`reflection_quartic_tiles.py`](../tools/reflection_quartic_tiles.py).
Its exact checker verifies the kernel edge budget, literal occupation words,
strict and triangular endpoints, restored original rational Fock profiles,
two-copy paired-type denominator cancellation and crossing, and a nonzero
quartic Pfaffian residual. Small profile controls use independent rational
DVR elimination with no series cutoff.

```sh
python tools/check_reflection_quartic_tiles.py \
  --output evidence/reflection-quartic-tiles.json
```

These controls establish the construction and its accounting. Searching its
complete physical types and heights still requires a strict pressure witness
before the factoring score changes.
