# Cancellation-aware finite modules for Smith channel words

## Result and scope

The degree-sensitive Smith bound can retain cancellation between auxiliary
paths. A finite deterministic automaton tracks a module of possible auxiliary
row vectors over a truncated power-series ring. Its weighted word count gives
the same effective crossing criterion as the earlier min-plus automaton.

For the same precision cap, this word bound is never weaker than capped
entry-order counting. Increasing the precision gives a monotone hierarchy.
For each fixed word, a finite exact order is eventually recovered. The state
space is finite because the coefficient field, auxiliary dimension, and
precision are fixed.

This is a sufficient-certificate interface, not an improved factoring
exponent. The checked elementary, ordered-pairing, and nilpotent-loop controls
do not pass the threshold. A separate dense nilpotent example proves a
strict improvement over min-plus counting.

The predecessor is PR #30 at
[28105165edb8e9a0427468dc833f3c09bc0a7fb7](https://github.com/1337isnot1337/integer-factoring/tree/28105165edb8e9a0427468dc833f3c09bc0a7fb7),
especially its
[local Smith profile](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md),
[degree-sensitive paths](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/degree-sensitive-smith-paths.md),
and [min-plus determinization](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/tropical-smith-channel-words.md).
Their ordered whole-word coefficient split is retained. The refinement
changes the test for which channel words can contribute.

## 1. Exact local data

Let \(F\) be one fixed finite field, \(O=F[[\tau]]\), and \(D\ge1\).
For one cut, suppose the fixed local crossing has a regularized factorization

\[
\sum_{i=1}^{b}\tau^{-e_i}U_i(\tau;a)v_i(\tau;b),                 \tag{1}
\]

where \(U_i\) is a \(D\)-by-\(D\) matrix and \(v_i\) is scalar, both with
regular rational entries. The physical half-labels range over fixed finite
sets. The auxiliary boundaries \(f^T,e\) are regular rational vectors.
Fixed poles from a rational gauge are charged to the predecessor theorem's
boundary shifts.

Discard an identically zero channel. Divide each scalar factor by its common
power of \(\tau\), and subtract that order from \(e_i\). Thus its minimum
order is zero, and (1) is unchanged. Each remaining matrix family is nonzero.
Put

\[
\delta_i=\min_{a,r,s}\operatorname{ord}_\tau(U_i(\tau;a)_{r,s}),
                                                                  \tag{2}
\]

ignoring zero entries. This is a finite nonnegative integer.

For a word \(w=(i_1,\ldots,i_t)\), define

\[
N_w=\operatorname{span}_{O}
 \{f^TU_{i_1}(a_1)\cdots U_{i_t}(a_t):a_1,\ldots,a_t\},
\qquad
v(w)=\min_{x\in N_w}\operatorname{ord}_\tau(xe).                 \tag{3}
\]

The zero ideal has infinite order. Taking the module span does not lower
the minimum order of its literal generators. Thus \(v(w)\) measures
whole-word cancellation that entrywise minimum-path costs can miss.

## 2. A finite module state

Fix an integer \(P\ge1\). A state is a submodule

\[
M\subseteq (O/\tau^P O)^D.
\]

Interpret it as its full inverse image \(\widetilde M\) in \(O^D\). In
particular, \(\tau^P O^D\subseteq\widetilde M\). This full tail is essential:
unrecorded high coefficients can enter the recorded window after division
by a common power.

Initially reduce \(Of^T\) modulo \(\tau^P\); the full lift is
\(Of^T+\tau^PO^D\). On channel \(i\), form

\[
R_i(M)=\sum_a\widetilde M U_i(a),\qquad
\gamma_i(M)=\min_{x\in R_i(M)}\operatorname{ord}_\tau x,
\qquad
T_i(M)=\tau^{-\gamma_i(M)}R_i(M)\bmod\tau^P.                    \tag{4}
\]

The order of a vector is the minimum coordinate order. The full tail
supplies \(\tau^P\) times every coordinate vector, so some generated entry
has order \(P+\delta_i\). Therefore

\[
0\le\gamma_i(M)\le P+\delta_i.                                \tag{5}
\]

The next state is uniquely determined. There are finitely many states:
each is an \(F\)-linear subspace of a \(PD\)-dimensional vector space,
invariant under multiplication by \(\tau\). The number of all such linear
subspaces is an explicit finite upper bound. No polynomial bound in \(P\)
or \(D\) is claimed; both are fixed certificate parameters.

### Exact finite computation

Represent a state by a field basis of polynomial row vectors of degree
below \(P\), and append the \(D\) generators \(\tau^P e_j\). These generate
its full lift over \(O\). Multiply them by every \(U_i(a)\) modulo

\[
\tau^{\,2P+\delta_i}.                                        \tag{6}
\]

The coefficient determining (5) is visible, and after removing that power
at least \(P\) coefficients remain. Reduction and module span compute
\(T_i(M)\) exactly. Rational denominators are units and have unique
computable jets. No numerical approximation or unproved zero test is used.

Keeping only \(P\) coefficients before dividing is insufficient. For
example, \((1,\tau^P)\) reduces to \((1,0)\), but multiplication by the
matrix with sole nonzero entry \(U_{21}=1\) gives
\((\tau^P,0)\). Declaring the truncated zero to be the exact product would
lose this vector. The appended full tail retains it.

## 3. Word bounds and completeness at fixed words

Let \(G_P(w)\) be the sum of transition increments and \(M_P(w)\) the
final state. For nonzero terminal vector \(e\), put

\[
\kappa(M)=\min_{x\in\widetilde M}\operatorname{ord}_\tau(xe),
\qquad m_P(w)=G_P(w)+\kappa(M_P(w)).                            \tag{7}
\]

If \(\delta_e\) is the minimum coordinate order of \(e\), the full tail
gives \(0\le\kappa(M)\le P+\delta_e\). This terminal cost is also a finite
exact jet calculation. A zero terminal vector makes the whole contraction
zero and can be handled separately.

The containment invariant is

\[
N_w\subseteq A_P(w):=\tau^{G_P(w)}\widetilde M_P(w).             \tag{8}
\]

It holds initially. Taking the image under every next matrix and summing
preserves containment; replacing a reduction by its full lift only
enlarges the result. Hence \(m_P(w)\le v(w)\).

An unnormalized version of the update makes two stronger statements
transparent. If \(R=\sum_a A_P(w)U_i(a)\) and \(g=\operatorname{ord}R\),
then

\[
A_P(wi)=R+\tau^{g+P}O^D.                                    \tag{9}
\]

**Monotonicity.** For \(P_2\ge P_1\),

\[
A_{P_2}(w)\subseteq A_{P_1}(w),\qquad
m_{P_2}(w)\ge m_{P_1}(w).                                   \tag{10}
\]

The initial tails are nested. Image modules remain nested, and their
orders satisfy \(g_2\ge g_1\). Thus the new tails in (9) are nested in
the same direction, proving the induction.

**Exactness below the precision.** All matrices are regular, so they
preserve \(\tau^PO^D\). Every removed increment is nonnegative. The
same induction proves

\[
N_w\subseteq A_P(w)\subseteq N_w+\tau^PO^D.
\]

The terminal vector is also regular. Consequently

\[
\min(m_P(w),P)=\min(v(w),P).                                 \tag{11}
\]

A fixed word of finite exact order is recovered once \(P>v(w)\).
A zero word receives arbitrarily large lower bounds. This does not
assert uniform convergence of exponential word-count rates or existence
of a finite improving certificate.

## 4. Comparison with capped min-plus counting

Assign auxiliary edge \((r,s)\) of symbol \(i\) the cost

\[
k_{i,r,s}=\min\!\left(P,\min_a
             \operatorname{ord}_\tau(U_i(a)_{r,s})\right),     \tag{12}
\]

with an identically zero entry assigned \(P\). Initialize coordinate
\(r\) with \(\min(P,\operatorname{ord}f_r)\), assigning \(P\) to a
zero coordinate. Update by min-plus multiplication, obtaining an
unnormalized vector \(c_w\). Use exact coordinate orders of \(e\) for
the terminal cost and omit its zero coordinates.

Then

\[
A_P(w)\subseteq\bigoplus_r\tau^{(c_w)_r}O.                    \tag{13}
\]

It holds initially. Matrix multiplication proves the inclusion for
the image \(R\). Since all edge costs lie in \([0,P]\), the next
coordinate costs have spread at most \(P\). Also
\(g=\operatorname{ord}R\ge\min_r(c_{wi})_r\). Every coordinate of the
added tail \(\tau^{g+P}O^D\) therefore has at least its required order.
Equation (9) proves the induction.

Applying the terminal vector yields

\[
m_P(w)\ge\text{the capped min-plus terminal cost of }w.        \tag{14}
\]

The boundary initialization and cap must agree in this comparison.
An uncapped exact-zero graph is a separate option; no infinite edge
cost is silently replaced by a finite one.

## 5. Weighted count and the factoring criterion

Fix rational \(s\ge1\). On reachable module states define

\[
\mathcal T(s)_{M,M'}=
 \sum_{i:T_i(M)=M'}s^{e_i-\gamma_i(M)},\qquad
r(s)_M=s^{-\kappa(M)}.                                      \tag{15}
\]

If \(\iota\) indicates the initial state, then

\[
Q_P(t,s)=\iota^T\mathcal T(s)^t r(s)
       =\sum_{w\in[b]^t}s^{\sum_j e_{i_j}-m_P(w)}.             \tag{16}
\]

Each channel word is counted once. Reversing the entire word, transposing
each auxiliary matrix, and swapping the boundaries gives another valid
certificate.

At infinity, a word can contribute to coefficient \(k\) only if
\(\sum e_i-m_P(w)\ge k-H_\infty\); at zero the condition is
\(\sum e_i-m_P(w)\ge-k-H_0\). The predecessor theorem's single
whole-word coefficient split uses at most \(C(t+1)\) rank-one terms
per admitted word. Thus

\[
\begin{aligned}
r(A_{t,k})&\le C_\infty(t+1)s^{H_\infty-k}Q_{P,\infty}(t,s),\\
c(A_{t,k})&\le C_0(t+1)s^{H_0+k}Q_{P,0}(t,s).
\end{aligned}                                               \tag{17}
\]

Supply strictly positive rational vectors and rational positive bounds
satisfying the finite entrywise inequalities

\[
\mathcal T_\infty(s)h_\infty\le\beta_\infty h_\infty,\qquad
\mathcal T_0(s)h_0\le\beta_0h_0.                              \tag{18}
\]

Each gives
\(Q_P(t,s)\le h_{\rm initial}\max_M(r_M/h_M)\beta^t\).
The coefficient index cancels between the two cuts, giving

\[
r(A_{t,k})c(A_{t,k})\le C(t+1)^2B_P^t,\qquad
B_P=\beta_\infty\beta_0.                                    \tag{19}
\]

There are \(O(t)\) coefficient pieces, so their total rank sum is at
most \(C(t+1)^2 B_P^{t/2}\). For physical alphabet size \(q\), the
effective crossing transfer therefore has width parameter

\[
\rho=\frac{\log B_P}{4\log q}.                               \tag{20}
\]

The exact sufficient test \(B_P^5<q^{14}\) would imply
\(\rho<7/10<1/\sqrt2\), hence a factoring exponent above one half.
No certificate satisfying it is supplied here.

All fields, rational matrices, automata, and precisions are fixed data.
Full coefficient tables of the original polynomial matrix product and
their rank factorizations still take \(2^{O(t)}\) bit operations.
The automaton is not an input-dependent oracle. The predecessor's
uniform crossing-to-factoring transfer and bit-cost obligations remain.

## 6. Completeness for the full auxiliary-matrix word rate

A scoped completeness theorem is available when both auxiliary boundaries
are left free. Start with \(N_\varnothing=O^D\), and define \(v_*(w)\) as
the minimum order of every entry of every literal matrix product for word
\(w\). Use the initial full module and terminate with minimum coordinate
order. The normalized terminal cost is zero.

This gives a valid upper bound for any fixed regular scalar boundaries:
multiplying by those vectors cannot lower the matrix-entry order. It may
lose cancellations special to their contraction.

For the fixed rational \(s>1\), define \(s^{-\infty}=0\) and

\[
Q_*(t)=\sum_{|w|=t}s^{e(w)-v_*(w)},\qquad
Q_{P,*}(t)=\sum_{|w|=t}s^{e(w)-G_P(w)}.                       \tag{21}
\]

The exact order satisfies \(v_*(uv)\ge v_*(u)+v_*(v)\), by ordinary
matrix multiplication with regular coefficients. Hence
\(Q_*(t+u)\le Q_*(t)Q_*(u)\).

The approximate update (9) is monotone under module inclusion and
homogeneous under multiplication by a scalar power of \(\tau\). Write
\(A_P(w)=\tau^{G_P(w)}\widetilde M_P(w)\), with
\(\widetilde M_P(w)\subseteq O^D\). Normalization gives
\(\operatorname{ord}\widetilde M_P(w)=0\), including the initial full
module, so \(G_P(w)=\operatorname{ord}A_P(w)\).
Continuing with a word \(v\) therefore
gives a submodule of \(\tau^{G_P(w)} A_P(v)\). Thus

\[
G_P(wv)\ge G_P(w)+G_P(v),
\qquad Q_{P,*}(t+u)\le Q_{P,*}(t)Q_{P,*}(u).                 \tag{22}
\]

For these nonnegative submultiplicative sequences, the exponential rates are

\[
\beta_*=\lim_{t\to\infty}Q_*(t)^{1/t}
       =\inf_{t\ge1}Q_*(t)^{1/t},
\qquad
\beta_P=\inf_{t\ge1}Q_{P,*}(t)^{1/t}.                        \tag{23}
\]

If an exact count is zero, all sufficiently long counts are zero, and its
rate is zero. The formulas include that case. For every fixed \(t\),
there are finitely many words. Monotonicity and (11) imply
\(Q_{P,*}(t)\downarrow Q_*(t)\). It follows that

\[
\boxed{\quad
\inf_{P\ge1}\beta_P
 =\inf_{P\ge1}\inf_{t\ge1}Q_{P,*}(t)^{1/t}
 =\inf_{t\ge1}\inf_{P\ge1}Q_{P,*}(t)^{1/t}
 =\beta_*.\quad}                                            \tag{24}
\]

The two infima commute; no uniform convergence is assumed.

Each \(\beta_P\) is the spectral radius of its finite reachable
nonnegative transfer matrix, because the initial state reaches every
retained state and the terminal vector is strictly positive. For every
rational \(\beta>\beta_P\), the rational vector

\[
h=(I-\mathcal T/\beta)^{-1}{\bf1}
\]

is strictly positive and satisfies
\(\mathcal T h=\beta(h-{\bf1})\le\beta h\). Thus positive rational
potential certificates approximate the rate with arbitrary strict slack.

Consequently, if the exact full-matrix word rates of two fixed local
factorizations satisfy a strict crossing threshold at one fixed rational
\(s>1\), some finite precisions and rational potentials certify it.
There is no guarantee about how large those precisions or state spaces
must be. This theorem does not assert completeness for cancellations
that exist only after fixed final boundary contraction.


## 7. Strict cancellation and actual crossing controls

Over \(F_2\), take one exponent-zero channel and

\[
U=\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad f=e=(1,0)^T.
\]

Min-plus gives zero cost at every length, but \(U^2=0\). The first
module transition reaches the span of \((1,1)\) with increment zero.
Subsequent transitions annihilate that exact row; the full tail
contributes order \(P\), and the normalized state is unchanged. Hence

\[
m_P(0^t)=P(t-1),\qquad Q_P(t,s)=s^{-P(t-1)}.                   \tag{25}
\]

The old word weight is one. This is an auxiliary-word control, not a
crossing decomposition.

Run the standard-library checker from the repository root:

    python3 tools/check_cancellation_module_automata.py \
      --output evidence/cancellation-module-automata.json

It reconstructs exact rational Schur factors, canonical module states,
literal ordered rational products, and both cuts of actual polynomial
crossings. It checks min-plus domination, exact weighted identities,
and rational potential bounds.

The default run uses precisions 1, 2, 3, word lengths through two, both
cuts and auxiliary-side choices, three integral gauges in dimension two,
and the rational points 1, 3/2, 2, 3, 5.

| Family | Automata checked | Largest reachable state count | Best sampled \(B_P\) |
| --- | ---: | ---: | ---: |
| Scalar elementary | 24 | 1 | \(256/27\) |
| Auxiliary elementary | 72 | 1 | \(256/27\) |
| Ordered pairing | 72 | 12 | \(16\) |
| Mixed nilpotent loop | 72 | 7 | \(16\) |

The first value recovers the exact elementary Smith result. These are
bounded certificate searches, not lower bounds on the crossing families.

The next constructive question is whether a larger auxiliary crossing
has recurring cancellations on enough words to make (19) competitive.
The finite module interface supplies a rigorous test. The local crossing
and a passing potential certificate remain needed.

The construction and exact checks used substantial OpenAI assistant
support. No Lean formalization or new leaderboard score is claimed.
