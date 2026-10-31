# The exact nonlinear symbol-charge domain of the reflection tensor

## Result and scope

The paired charges used for typed rational synthesis can be arbitrary
functions of each physical letter. They need not be sums of occupation-bit
weights. For the complete reflection alphabet, their admissible domain has
an exact description by Hamming-cube Lipschitz constants.

At four symbols the domain is a six-dimensional cross-polytope with twelve
vertices. Four of those vertices are nonlinear parity charges, and are
absent from the additive-mode family. This supplies an exact parameter
domain for the complete-entry rational synthesis interface. It is not a
pressure-floor theorem, and does not assert a factoring improvement.

## 1. The exact charge criterion at every interval

Let M be a power of two at least four, m=M/2, and q=2^m. Let T be the
polynomial reflection tensor, chi=1+w, and v the valuation at chi=0.
Choose arbitrary real functions f,g on the m-dimensional binary cube and
put

\[
 E(x,y,u,v)=f(x)+g(y)-g(u)-f(v).
 \tag{1}
\]

Their constants are irrelevant. Define

\[
 L_f=\max_{d_H(x,v)=1}|f(x)-f(v)|,
 \qquad
 L_g=\max_{d_H(y,u)=1}|g(y)-g(u)|.
 \tag{2}
\]

**Theorem.** The inequalities

\[
 v(T_{xyuv})+E(x,y,u,v)\ge0
 \quad\text{for every nonzero coefficient}
 \tag{3}
\]

hold if and only if

\[
 \boxed{L_f+L_g\le2.}
 \tag{4}
\]

The theorem concerns the complete physical alphabet. Restricting the
alphabet can remove inequalities and give a larger admissible domain.

### Sufficiency

Consider one nonzero Gaussian matching term for a physical entry. Let r be
the number of external nonconstant edges, i_X the number of internal edges
within X or V, and i_Y the number within Y or U. The other edges are the
constant matched X,V and Y,U edges.

Every occupied bit that differs between X and V must be incident to a
nonconstant edge. An external edge consumes one vertex from each paired
species, and an internal edge consumes two from its own species. Hence

\[
 d_H(x,v)\le r+2i_X,\qquad d_H(y,u)\le r+2i_Y.
 \tag{5}
\]

The charge of the entry is therefore at least

\[
 -(L_f+L_g)r-2L_f i_X-2L_g i_Y.
 \tag{6}
\]

The matching term has chi-order 2r+2M(i_X+i_Y). Under (4), both L_f and
L_g are at most two, and 2M is at least eight. Its valuation pays (6).
Summing matching terms can only increase the final valuation. The original
common power of w is a unit at chi=0, so it changes no conclusion.

### Necessity

Take any directed Hamming-neighbor pair x,v and any directed neighbor pair
y,u. Their physical entry has

\[
 v(T_{xyuv})=2.                                             \tag{7}
\]

There is exactly one matching with one external nonconstant edge and all
other edges constant and matched. The two differing bits specify the
external edge, and every remaining occupied pair specifies its constant
edge. Its leading chi-squared coefficient is one. Every other matching
uses at least two external nonconstant edges or an internal edge, so has
higher valuation and cannot cancel it.

Choose an edge attaining L_f and orient it so f(x)-f(v)=-L_f. Independently
choose and orient an edge attaining L_g so g(y)-g(u)=-L_g. Applying (3) to
their entry gives 2-L_f-L_g at least zero. This proves (4).

## 2. Complete typed rational families from nonlinear charges

For rational f,g satisfying (4), use the entrywise construction from
[typed charge lattices](typed-charge-lattice-synthesis.md):

\[
 k_\lambda=\lfloor-E(\lambda)\rfloor+1,
 \quad S=\max(0,-\min k_\lambda),
 \quad
 N_\lambda=\chi^S T_\lambda+\chi^{S+k_\lambda}P_\lambda(w).
 \tag{8}
\]

At negative integer charge, the optional boundary extension can replace
k_lambda by -E(lambda). All P_lambda are independent polynomials with the
explicit numerator degree caps. Paired full-letter types still force the
total charge (1) to be zero, so the original cancellation and endpoint
proof applies without change. This uses no independent-bit assumption on
the physical types.

The cube diameter is m. Thus the oscillations of f and g have sum at most
m(L_f+L_g), at most 2m. Consequently S is at most 2m-1, exactly as in the
mode-weight family. Effective degrees, the complete infinity scalar shift,
full generic local rank, and the 2^O(M) construction accounting remain as
stated in the preceding interface.

## 3. The exact four-symbol cross-polytope

Normalize f(0)=g(0)=0 and identify a symbol with its two occupation bits.
For f define

\[
 p_1=\frac{f(3)}4,\qquad
 p_2=\frac{f(1)-f(2)}4,\qquad
 p_3=\frac{f(1)+f(2)-f(3)}4,
 \tag{9}
\]

and define q_1,q_2,q_3 similarly from g. Equivalently,

\[
 (f(1),f(2),f(3))
 =p_1(2,2,4)+p_2(2,-2,0)+p_3(2,2,0).
 \tag{10}
\]

The directed differences on the four square edges are precisely twice all
signed sums of the three p coordinates. Therefore

\[
 L_f=2(|p_1|+|p_2|+|p_3|),
 \qquad
 L_g=2(|q_1|+|q_2|+|q_3|).
 \tag{11}
\]

The entire admissible charge polytope is

\[
 \boxed{|p_1|+|p_2|+|p_3|+|q_1|+|q_2|+|q_3|\le1.}
 \tag{12}
\]

It has twelve signed coordinate vertices and sixty-four sign facets. The
vertices on either physical species are the signed potentials

\[
 (0,2,2,4),\qquad(0,2,-2,0),\qquad(0,2,2,0),             \tag{13}
\]

with zero potential on the other species. The last potential is twice the
occupation parity. It is nonlinear in the occupation bits. The first two
are additive count and signed-count potentials.

The exact coefficient certificate independently verifies every original
inequality at all twelve vertices. It also identifies all sixty-four facet
inequalities among the literal coefficients of chi-order two. This is a
complete polytope certificate; no numerical vertex enumeration is needed
to establish (12).

## 4. Finite exact synthesis interface

[`nonlinear_symbol_charge_synthesis.py`](../tools/nonlinear_symbol_charge_synthesis.py)
implements the exact Lipschitz test, the four-symbol coordinates, and the
complete M=4 affine constructor. It shares coefficient decoding, complete
profile evaluation, and witness checks with the frozen typed-charge
interface. Its constructor also verifies (3) directly on every nonzero
baseline coefficient before building the two regular frames.

The nonlinear coefficients are physical-symbol charges used to certify the
endpoint. They are separate from the independent physical-symbol heights
used to optimize the two cost cuts. Both choices remain explicit in a
candidate. Finite target walks provide exploratory witnesses and do not
establish global optimality over the charge polytope or the row frames.

```sh
python tools/check_nonlinear_symbol_charges.py \
  --output evidence/nonlinear-symbol-charges.json
```

The checker verifies the full finite charge polytope, literal neighbor
coefficients, nonlinear affine tables, paired-type division and endpoints,
and complete restored local profiles. The theorem is useful even when a
particular finite synthesis run yields no pressure improvement.
