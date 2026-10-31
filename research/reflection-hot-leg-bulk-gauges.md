# Separate physical-leg bulk heights for reflection hot cores

## Result and scope

The [cold-compensated hot-core reduction](reflection-hot-compensated-gauges.md)
can assign a separate nonnegative bulk height to each of the four
physical legs. A simple ordering cone preserves the same cross-only
Gaussian basis and the same \(\min(M,4k)\) active modes. No additional
matrix coefficients need to be compiled: only explicit coordinate,
spectator, and physical type heights change.

The construction applies to an even exterior polynomial defect on the
selected \(k\) hottest reflection bits, where \(k\) is a power of
two. It retains every internal Gaussian factor involving two selected
hot variables. Arbitrary coefficient masks are outside this statement.
This is an evaluator and gauge extension, with no improved factoring
exponent asserted.

## 1. A cone for separate leg heights

Use the physical leg order \((x,y,u,v)\). At zero the two boundaries
are \((x,y)\mid(u,v)\); at infinity use \((u,x)\mid(v,y)\).
At either cut denote the boundaries by
\((\ell_1,\ell_2)\mid(r_1,r_2)\). Give every bit of leg \(a\)
a bulk height \(b_a\ge0\), subject to

\[
b_{\ell_1}\ge b_{\ell_2},\qquad
b_{r_2}\ge b_{r_1}.
\tag{1}
\]

Thus the zero inequalities are \(b_x\ge b_y\), \(b_v\ge b_u\),
and the infinity inequalities are \(b_u\ge b_x\), \(b_y\ge b_v\).
Each cut has its own four bulk heights.

Write \(A,B\) for the lower and upper halves of the natural diagonal
Gaussian modes. The natural regular cross reduction remains regular
with the following boundary mode heights:

| Boundary | Block \(A\) | Block \(B\) |
|---|---:|---:|
| Left | \(b_{\ell_2}\) | \(b_{\ell_1}\) |
| Right | \(b_{r_1}\) | \(b_{r_2}\) |

Indeed, its two weighted cross shears acquire the extra nonnegative
orders \(b_{\ell_1}-b_{\ell_2}\) and
\(b_{r_2}-b_{r_1}\). The same is true of their inverses.
The residue concentration shears act within each block, so their
orders are unchanged. All concentration companions stay regular.

## 2. Which internal factors may be removed

Let \(h_a(s)\) be an additional height on the selected \(k\)-bit
state of physical leg \(a\). The full added physical height is
\(b_a\) times total bit occupation, plus \(h_a(s)\).

Consider an internal Gaussian factor joining a selected hot bit on
leg \(a\) to an unselected cold bit on leg \(d\). Let \(\nu\ge0\)
be its coefficient order in the natural gauged kernel, after the
ramification substitution. It is enough to require

\[
\nu+b_a+b_d+h_a(s+e_i)-h_a(s)\ge0
\tag{2}
\]

whenever the selected bit \(i\) is absent from \(s\). Condition
(2) is imposed for every such factor. Cold-cold factors are regular
because all bulk heights are nonnegative. Factors joining two selected
hot variables remain in the active tensor and require no inequality.

An elementary internal factor and its inverse have the same coefficient
orders in characteristic two. Hence all removed factors remain regular
units under (2), on both boundaries and at both cuts. The exterior
defect commutes with them, exactly as in the compensated construction.

The implementation checks (2) directly against the exact natural
kernel. This also covers small intervals where a potential cold partner
is absent. For one hot bit with cold partners present, the zero-order
same-leg reflection and cross-leg factors give the sufficient bound

\[
h_a(1)-h_a(0)\ge-b_a-\min(b_a,b_d),
\tag{3}
\]

where \(d\) is the other leg of its boundary.

## 3. Why the existing hot-state coordinates still work

Let \(D_b\) be the diagonal matrix of the mode heights in the table
of Section 1. The cross-only transported physical creation and
annihilation operators have a triangular form. On the left, the
\(\ell_1\) creator is supported in block \(B\), while the
\(\ell_2\) creator has its principal support in \(A\) and a
regular tail in \(B\). The dual annihilators have the opposite
triangular form. Conjugating by \(D_b\), and normalizing a creator
by \(z^{-b_a}\) and its annihilator by \(z^{b_a}\), multiplies
the spreading coefficients by
\(z^{b_{\ell_1}-b_{\ell_2}}\). They remain regular by (1).

For the column boundary the conjugation is reversed. The regular
normalized right creator has its additional factor
\(z^{b_{r_2}-b_{r_1}}\); its dual obeys the same condition.
All canonical anticommutation relations are preserved. The residue
concentration occurs within blocks, so the same argument holds after
the reduction to \(\min(M,4k)\) modes.

Choose a principal unit minor of the physical vacuum projector. Such a
minor always exists: an idempotent of residue rank \(v\) has
\(\det(I+tQ)=(1+t)^v\), so the sum of its principal minors of
order \(v\) is one in the residue field. At least one is a unit.
If that projector is \(Q\),
the weighted projector is \(D_bQD_b^{-1}\), or its column-side
counterpart. Its selected principal minor has exactly the same
determinant: the row and column monomial factors cancel. Thus the
same vacuum labels give a regular invertible coordinate basis after
the leg weights are introduced. The existing cached cross-only bases
use such a principal minor in the replayed families; the implementation
verifies that their saved vacuum rows and columns coincide and rejects
an unsupported cached basis before using this rule.

A coordinate with selected physical hot state \(s\) and vacuum
Fock label \(t\) therefore receives the height

\[
\sum_a b_a|s_a|+h_a(s_a)
 +\sum_{i\in t} b_{\mathrm{block}(i)}.
\tag{4}
\]

The block height in the last sum is the appropriate left or right
entry from the table. In particular, it is not in general a common
multiple of total occupation. Conjugating the coordinate bases gives
exactly these row and column scalars on the previously compiled
coefficient matrix. This proves that the old complete coefficients
can be reused without a new truncation or source expansion.

## 4. Spectators and the paired physical partitions

The original natural diagonal mode indexed by \(i\),
\(0\le i<M\), has order \(P(2+4i)\). Every unselected
spectator gets the additional order

\[
\begin{cases}
b_{\ell_2}+b_{r_1},&i<M/2,\\
b_{\ell_1}+b_{r_2},&i\ge M/2.
\end{cases}
\tag{5}
\]

Write \(b_a^{(0)},b_a^{(\infty)}\) for the two-cut bulks.
The added energy of every cold bit in the paired \(x,v\) type is

\[
\Delta_x=b_x^{(0)}+b_v^{(0)}
          +b_x^{(\infty)}+b_v^{(\infty)},
\tag{6}
\]

while the paired \(y,u\) type gets the corresponding sum
\(\Delta_y\). These can differ. A hot state receives the same
per-bit contribution plus the paired sums of its relative heights
\(h\). Every spectator and both physical type partitions are
included in the pressure evaluator.

Setting the two bulks equal on each boundary recovers the previous
compensated cone exactly. The original polynomial tensor, endpoint,
and degree are unchanged by this choice of profile gauges.

## Exact replay

The evaluator is
[`reflection_hot_leg_bulk_gauges.py`](../tools/reflection_hot_leg_bulk_gauges.py).
The [independent checker](../tools/check_reflection_hot_leg_bulk_gauges.py)
and [saved extended replay](../evidence/reflection-hot-leg-bulk-gauges.json)
verify twelve complete original weighted cuts and 132,096 endpoint
coefficients, including nonlinear two-hot-bit states and a ramification
of three. Weighted creation and annihilation controls verify the
principal unit minors through \(M=64\), and common-boundary gauges
agree with the preceding compensated evaluator. All seventeen analytic
derivatives agree with centered finite differences, with maximum error
below \(1.5\times10^{-10}\).

Profile calculations retain the complete source coefficients and use
the existing adaptive Smith workspace; a larger workspace never
replaces missing source coefficients. The
[search driver](../tools/search_reflection_hot_leg_bulk_gauges.py)
uses continuous optimization within the stated cone. Its finite runs
are exploratory evidence, separate from the exact reduction above.
