# Forward and reverse quartics in one growing paired word

## Result and scope

The [paired-quartic construction](reflection-paired-quartic-poles.md) can
use both occupation directions in one word, provided different gates use
disjoint logical pairs from the same fixed-bit pairing. The two directions
have different removable-pole budgets under a common charge.

Write h=4b for the selected logical bit b. A forward gate permits pole
h+7, while a reverse gate permits pole h-1. For the highest logical bit,
these bounds are M+7 and M-1. After the displayed common monomial clearing,
every paired-type tensor power is polynomial and has the crossing endpoint.
The degree and denominator accounting remain explicit. No improved
factoring pressure is claimed.

## 1. The two directions and one common charge

Use E=2(|x|-|v|) and the chi-gauged Gaussian kernel from the preceding
construction, with chi=1+w. On a pair i,j=i xor b define

\[
 Q^+_{ij}=C_{X,i}C_{X,j}A_{V,i}A_{V,j},\qquad
 Q^-_{ij}=A_{X,i}A_{X,j}C_{V,i}C_{V,j}.                    \tag{1}
\]

Their charge increments are eight and minus eight. The forward pullback
has its designated V modes occupied and their X partners absent. The
reverse pullback instead has its designated X modes occupied and their
V partners absent.

The same statement applies after replacing X,V by Y,U or reversing the
chosen common charge. One unordered physical pair is used throughout.

## 2. Extra row factors for a reverse pullback

Every designated X row in the reverse source has its only constant
matched V column absent. Its remaining external entries have a factor
chi^4. Internal XX entries have order 2M+4. Thus chi^4 can be factored
from each designated X row and column.

For k_r reverse pairs there are 2k_r such rows. Their initial contribution
to the Pfaffian is chi^{8k_r}. This simultaneous extraction is regular:
an internal XX entry between two extracted rows retains order 2M-4,
and an XX entry with only one extracted endpoint retains order 2M.

After these factors, the two X rows of each reverse pair agree modulo
chi^h, just as the V rows of a forward pair do. Their external monomial
exponents differ by zero or plus or minus h. Since h is at most M,
the remaining internal orders are sufficient as well.

## 3. A source bound for every mixed subword

Select k_f forward pairs and k_r reverse pairs from a disjoint fixed-bit
pairing. First extract the reverse X row factors described above. Then
perform the row and column addition on each selected logical pair.
Each designated difference row is divisible by chi^h.

All intersections between designated rows vanish exactly. For two V
differences or two scaled X differences, the remaining internal block is
one common scalar times reversal; the fixed-bit pair vectors are mutually
isotropic for that reversal form. A mixed X,V intersection is also zero:
the only X,V kernel edges are the constant matched edges, and the two
selected logical pairs have disjoint indices.

Therefore every selected source coefficient satisfies

\[
 v_\chi(T_{\rm source})+E({\rm source})
 \ge8k_r+h(k_f+k_r).                                      \tag{2}
\]

The selected operators add charge 8k_f-8k_r. Before their scalar poles
are inserted, their total correction margin is at least

\[
 (h+8)k_f+h k_r.                                          \tag{3}
\]

Consequently the growing product

\[
 \prod_{ij}\left(1+\frac{P_{ij}(w)}{\chi^{s_{ij}}}
                          Q^{\pm}_{ij}\right)
 \tag{4}
\]

has strictly positive weighted correction valuation for every nonempty
selected subword whenever

\[
 s_{ij}\le h+7\quad\text{in the forward direction},
 \qquad
 s_{ij}\le h-1\quad\text{in the reverse direction}.       \tag{5}
\]

Each numerator P is a Laurent polynomial in w and is regular at w=1.
Summing subwords preserves the
strict margin. Paired physical word types cancel the total charge, so
every nonempty perturbation term in every tensor power is regular and
vanishes at the endpoint. The original crossing is unchanged.

Disjointness is essential to this proof. It does not permit a forward
and reverse product on the same logical pair: their cross term has not
been paid for by (2).

## 4. Additive corrections on the same pair

There is a separate valid construction that uses both directions on
one pair without that cross term:

\[
 T'=T+\frac{P_+(w)}{\chi^{s_+}}Q^+T
       +\frac{P_-(w)}{\chi^{s_-}}Q^-T,
 \quad s_+\le h+7,\quad s_-\le h-1.                       \tag{6}
\]

Each summand is a strict correction under the same charge E. Their sum
therefore lies in the complete-entry typed charge lattice and preserves
all paired-type endpoints. Formula (6) is an additive tensor family;
expanding a product of the two gates would introduce an additional term
and is a different construction.

## 5. Costs and exact controls

The raw common denominator is the sum of the selected gate pole orders.
The strict charge condition again forces the actual local pole to be at
most M-1. The constructor removes its exact common chi content before
profile evaluation. If C is the common w clearing, N the reduced
numerator, and S its remaining pole, the original profiles are

\[
 \sigma_0(N)-C,\qquad\sigma_\infty(N)+C+S.                \tag{7}
\]

For numerator powers in [a_j,b_j], the effective degree remains bounded
by

\[
 \frac{3M^2}{2}+
 \sum_j\bigl(\max(0,-a_j)+\max(0,b_j-s_j)\bigr).          \tag{8}
\]

In particular, numerator supports in [0,s_j] keep the original degree
bound even for a growing mixed word. Full generic local rank follows
from the all-length paired-type endpoint theorem.

[`reflection_mixed_paired_words.py`](../tools/reflection_mixed_paired_words.py)
implements the disjoint mixed-word family and enforces the directional
budgets. Its finite constructor uses the complete M=4 and M=8 tensors;
the theorem and validator are uniform in M.

The checker verifies the initial X row divisions, the scaled row
congruences, the exact zero mixed intersections, complete original
subword margins, paired-type polynomial division and endpoints, actual
Fock profiles, and independent rational profiles. It rejects both
overlapping opposite-direction gates and an excessive reverse pole.

```sh
python tools/check_reflection_mixed_paired_words.py \
  --output evidence/reflection-mixed-paired-words.json
```

These constructions enlarge the explicit non-Gaussian search family.
Their endpoint proofs do not replace the strict full-pressure comparison
needed for a new factoring bound.
