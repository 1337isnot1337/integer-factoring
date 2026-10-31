# Growing paired quartic words with removable poles of order M+7

## Result and scope

The reflection tensor admits a growing family of nonquadratic occupation
gates whose coefficients have endpoint poles of order M+7 per gate.
The gates act on disjoint pairs of logical modes related by one fixed
binary bit flip. Their number can be proportional to M. Every paired
full-letter-type restriction of every tensor power is polynomial and has
the crossing endpoint.

The proof uses exact divisibility of Gaussian Pfaffians beyond the earlier
edgewise charge budget. For numerator powers between zero and the gate's
pole order, the effective polynomial degree is at most the original
3M^2/2. The raw local denominator can have order proportional to M^2, but
a forced common factor reduces the actual local pole to at most M-1.

This is an explicit uniform construction and a stronger pole certificate.
It does not provide a negative full-pressure witness or a better factoring
bound. Its local coefficients belong to the
[complete-entry typed charge lattice](typed-charge-lattice-synthesis.md).

## 1. The occupation gauge and its local kernel

Let M be a power of two at least four and m=M/2. Work first with the X,V
matched physical pair. Set

\[
 E(x,y,u,v)=2(|x|-|v|),\qquad \chi=1+w.
 \tag{1}
\]

The [charge theorem](nonlinear-symbol-charge-lattices.md) gives
v_chi(T_lambda)+E(lambda) at least zero. Equivalently, the entrywise
gauged tensor G_lambda=chi^{E(lambda)}T_lambda is regular at chi=0.
It remains Gaussian: multiply each kernel edge by chi to the sum of its
two vertex charges.

The gauged kernel K is regular over the chi-adic local ring. Modulo chi^4
it is

\[
 K_0=\sum_i(x_i v_i+y_i u_i)
      +\left(\sum_i v_i\right)\left(\sum_j(y_j+u_j)\right).
 \tag{2}
\]

Here the display specifies the alternating quadratic kernel. The VY and
VU external edges become powers w^{4a}, all congruent to one modulo
chi^4. The XU and XY edges acquire a factor chi^4. Internal VV edges
have order 2M-4, at least four; every other internal edge has still higher
order. The original common scalar is a chi-adic unit.

## 2. A stronger row congruence for each pair

Choose a logical bit b=2^r with b<m and pair i with j=i xor b. Put

\[
 h=4b\le M.
 \tag{3}
\]

In a Pfaffian coefficient with V_i,V_j occupied and X_i,X_j absent, the
two V rows agree modulo chi^h. For the Y and U columns, their exponents
as powers of w differ by zero or plus or minus h. Thus their row sums
are divisible by chi^h. Internal V entries have order 2M-4, at least h.
The only constant matched X entries would be X_i and X_j, and both are
absent from this coefficient.

Adding row V_j to row V_i and performing the same column addition
preserves the Pfaffian and makes the resulting V_i row divisible by
chi^h. For one pair, this already proves a factor chi^h in the gauged
coefficient.

For a single arbitrary pair i,j, the same reasoning gives
h=4 times the largest power of two dividing i-j. The growing-word
interface below uses the fixed-bit pairing so all simultaneous row
factors are controlled exactly.

## 3. Simultaneous factors for every selected subword

Take any k disjoint pairs from the same fixed-bit pairing. Suppose all
their 2k V modes are occupied and all matched X modes are absent. Apply
the row and column addition above for every pair.

Each designated difference row is divisible by chi^h. To extract a
factor from each row and column simultaneously, intersections of two
designated rows must be divisible by chi^{2h}. In this family they are
exactly zero.

Indeed, the internal VV block is a common scalar times the reversal
matrix R, with R(i)=m-1-i. A row-difference vector is e_i+e_(i xor b).
Bitwise complement commutes with the fixed xor flip. Therefore reversal
maps each complete flip-pair to another complete flip-pair, and

\[
 (e_i+e_{i\mathbin{\mathrm{xor}}b})^T
 R(e_j+e_{j\mathbin{\mathrm{xor}}b})=0
 \tag{4}
\]

in characteristic two: there are either zero or two incidences.

After the simultaneous congruences, factor chi^h from each of the k
designated rows and columns. The remaining alternating matrix is regular.
The Pfaffian scales by the determinant of this diagonal factor, giving

\[
 v_\chi(T_{\rm source})+E({\rm source})\ge kh.
 \tag{5}
\]

The statement holds for every selected subset of the pairing, not only
for the complete set. This is the condition needed for a gate-word
expansion and for tensor powers of that word.

## 4. The growing word and its endpoint theorem

For each selected logical pair i,j define

\[
 Q_{ij}=C_{X,i}C_{X,j}A_{V,i}A_{V,j},\qquad
 U_{ij}=1+\frac{P_{ij}(w)}{\chi^{s_{ij}}}Q_{ij},
 \qquad 0\le s_{ij}\le h+7.
 \tag{6}
\]

Each P is a Laurent polynomial regular at w=1. The gates have disjoint
mode supports and charge increment eight under (1). A nonzero selected
product of k gates pulls its target coefficient back to a source of the
form used in (5). Its weighted valuation is at least

\[
 kh+8k-\sum s_{ij}\ge k>0.
 \tag{7}
\]

Numerator zeros at w=1 can only increase this margin. Summing the selected
terms preserves the strict inequality for the complete correction.
Consequently the full local word minus T satisfies the strict entrywise
charge-lattice condition.

In every retained paired-type tensor-power coefficient, total charge is
zero. Every nonempty perturbation selection therefore has positive total
chi-order. All local endpoint poles cancel, those terms vanish, and the
empty selection gives the crossing. This proves polynomiality and the
typed endpoint for every tensor-power length.

Choose the highest logical bit b=m/2. Then h=M, and the complete pairing
has m/2=M/4 gates. Each gate may have pole order

\[
 \boxed{s_{ij}=M+7.}                                       \tag{8}
\]

The same construction works for the Y,U pair, and with creation and
annihilation directions reversed. A word uses one fixed oriented physical
pair throughout. Mixing those different words requires an additional
joint charge analysis and is not asserted here.

## 5. Common denominator cancellation and degree costs

There is a general local-pole corollary for the complete-entry charge
spaces. For any symbol potentials with L_f+L_g at most two, their combined
oscillation is at most m(L_f+L_g), hence the maximum paired charge is at
most M. If a correction satisfies v_chi(delta T_lambda)+E(lambda)>0,
integrality of the valuation gives

\[
 v_\chi(\delta T_\lambda)\ge1-M.                           \tag{9a}
\]

Adding a polynomial baseline preserves this lower bound. Thus any raw
common denominator chi^S must cancel by at least max(0,S-(M-1)) in every
numerator entry. This conclusion does not depend on a quartic-word
representation. Dividing all entries and the denominator by chi^c lowers
both the actual numerator degree and S by exactly c, so deg(N)-S is
unchanged.

Let S_raw be the sum of the gate pole orders, and clear all negative
powers of w by a common w^C. The literal local numerator is N_raw, so the
word is w^{-C}N_raw/chi^{S_raw}. The maximum charge in (1) is 2m=M.
Since the entire correction has weighted valuation strictly positive,
its local pole is at most M-1. The baseline has no chi pole. Therefore

\[
 \chi^{\max(0,S_{\rm raw}-(M-1))}\mid N_{\rm raw,\lambda}
 \quad\text{for every entry.}
 \tag{9}
\]

The constructor removes the exact common chi content. If c is removed,
both numerator degree and denominator order fall by c. The effective
degree and both original rational cut profiles are unchanged.

If the numerator of gate j has powers from a_j to b_j, the complete
effective cleared degree obeys

\[
 D_{\rm eff}\le\frac{3M^2}{2}
  +\sum_j\left(\max(0,-a_j)+\max(0,b_j-s_j)\right).
 \tag{10}
\]

This is the same ledger as for the earlier quartic words, now applied to
the larger proved poles. In particular, if each numerator has powers in
[0,s_j], the degree bound is exactly the original 3M^2/2. With O(M)
gates, powers O(M), and poles O(M), the raw degrees and the number of
coefficient operations remain O(M^2) and 2^O(M), respectively.

After the exact common cancellation, let the remaining numerator be N
and its denominator order S. The original weighted profiles are

\[
 \sigma_0=\sigma_0(N)-C,\qquad
 \sigma_\infty=\sigma_\infty(N)+C+S.
 \tag{11}
\]

The remaining S is included in every profile evaluation. Full generic
rank follows from the all-length paired-type endpoint theorem; finite
specialization checks remain independent controls.

## 6. Exact interface and checks

[`reflection_paired_quartic_words.py`](../tools/reflection_paired_quartic_words.py)
validates the common xor bit, the oriented matched physical pair,
disjointness, and each pole bound. Its literal constructor supports the
existing M=4 and M=8 complete coefficient controls, while the theorem and
validator apply to every allowed interval. The implemented numerators
are arbitrary distinct powers with one nonzero F4 color; the proof also
permits general characteristic-two polynomial numerators.

The checker verifies row congruences and exact zero intersections through
M=64 in all four orientations, full original coefficient bounds for every
subword in small complete intervals, paired-type polynomial division,
restored Fock profiles, independent rational DVR profiles, and nonzero
quartic Pfaffian residuals. It also checks the forced common denominator
cancellation.

The strict endpoint margin is necessary in general. For M=4, the single
X,V gate on logical pair (0,1) at the first excluded pole s=12 has a
noncrossing balanced three-copy coefficient with endpoint one. The word
uses entries (3,0,0,0), (0,0,1,1), and (1,1,0,3). Its X,V word multisets
and Y,U word multisets agree, but the pointwise crossing conditions fail.
This literal negative control prevents interpreting polynomiality alone
as a valid endpoint.

```sh
python tools/check_reflection_paired_quartic_words.py \
  --output evidence/reflection-paired-quartic-words.json
```

The family supplies larger exact rational perturbations for subsequent
full-type and full-height searches. A change in the factoring score still
requires a strict pressure crossing and complete recursive accounting.
