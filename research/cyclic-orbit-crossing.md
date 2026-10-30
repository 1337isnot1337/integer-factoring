# Odd cyclic orbit sums of a finite crossing seed

## Construction and finite result

Let \(A_0,\ldots,A_{19}\) be the twenty three-label tensors from the
[pinned PR #24 certificate](https://github.com/1337isnot1337/integer-factoring/tree/e96ee50560284301428d72421330d8773155d25b).
Their exact matrices are included in
[finite-seed-amplification.md](finite-seed-amplification.md).
The coefficient tensors of

\[
\left(\sum_{i=0}^{19}X_iA_i\right)^{\otimes t}
\tag{1}
\]

sum to the crossing on \(3^t\) word labels.

Suppose \(t\) is odd. On each external leg, replace the word basis by
the sums of its cyclic rotation orbits. These sums are self-dual:
different orbits are disjoint, and every orbit has odd size.
Consequently restricting all four legs of (1) to this basis preserves
the scale-one crossing identity.

There are

\[
D_t=\frac1t\sum_{j=0}^{t-1}3^{\gcd(t,j)}
\tag{2}
\]

labels. In particular \(D_3=11\).
The exact three-copy computation gives **1,522 nonzero coefficient
pieces, with maximum rank product \(210\)**. The four maximizing types
are

\[
(10,15,16),\quad(11,15,16),\quad
(15,16,17),\quad(15,16,19),
\]

with ranks \(14,15\) or \(15,14\).
The coefficient tables and all ranks are reproduced by
[check_cyclic_orbit_crossing.py](../tools/check_cyclic_orbit_crossing.py).
This is a finite exact construction, not an asymptotic factoring
improvement.

On the **full orbit-sum label domain** for one growing odd cyclic
group, the literal types have an asymptotic obstruction:

\[
\max_i\log_2 r_i,\ \max_i\log_2 c_i
\ \ge\ t\log_2 6-O(\log(t+2)),
\tag{3}
\]

where the same explicit type achieves both inequalities.
Thus its joint-width coefficient is at least

\[
\rho\ge\frac{\log 6}{2\log3}=0.81546\ldots.
\tag{4}
\]

The proof below identifies blocks that survive cyclic averaging.
This lower bound does not automatically survive the binary prefix
restriction of Section 1. A rank lower bound for a selected prefix
would require a separate argument.
It leaves open taking the fixed eleven-label polynomial seed and
amplifying that seed by another construction. Such a restriction uses
many separate three-cycles and is different from one growing cycle.
The result also leaves arbitrary coefficient recombinations open.

## 1. Uniform effectivity of the cyclic construction

The number of piece-count vectors in (1) is
\(\binom{t+19}{19}\), which is polynomial in \(t\).
For an orbit-sum entry, sum over the four external orbits; each has at
most \(t\) words. A fixed word quadruple contributes a coefficient
of a product of \(t\) linear forms in twenty variables.
A coefficient dynamic program has at most \((t+1)^{20}\) states.
Thus an entry can be computed in polynomial time in \(t\).

The orbit labels can be generated canonically by comparing the cyclic
rotations of each word. Complete tables and both rank factorizations
therefore take \(2^{O(t)}\) time.
For a binary label length \(\lambda\), take the least odd \(t\) with
\(D_t\ge2^\lambda\), and retain the first \(2^\lambda\) orbit sums.
Their Gram matrix remains the identity.
Since \(3^t/t\le D_t\le3^t\), minimality gives

\[
t=\Theta(\lambda),\qquad
\log_2D_t=\lambda+O(\log(\lambda+2)).
\]

This verifies a uniform, scale-one, polynomial-piece construction
with \(2^{O(\lambda)}\) generation time. No uncomputed factorization
oracle is used.

## 2. A literal type with forced position patterns

Write

\[
t=3q+r,\qquad 0\le r\le2,
\]

and consider the type

\[
12^r\,4^q\,5^q\,17^q.
\tag{5}
\]

Restrict further to the orbit sums of words containing exactly
\(m=q+r\) copies of \(w=e_2\).
This gives a subspace of the full cyclic label space, so its ranks
lower-bound the full ranks.

Put \(U=\operatorname{span}(e_0,e_1)\), \(v=e_0+e_1\), and let \(J\)
swap \(e_0,e_1\). The original source matrices give

\[
A_4=w_xv_vJ_{yu},\qquad
A_5=v_xw_vJ_{yu},\qquad
A_{17}=w_u(w_yI_{xv}+J_{xy}v_v).
\tag{6}
\]

The entry of \(A_{12}\) at \(wwww\) is one.
In a surviving term of (5), all \(q\) copies of \(A_4\) supply an
\(x\)-occurrence of \(w\), and the other core pieces never supply one.
The \(x\)-count \(q+r\) therefore forces all \(r\) copies of \(A_{12}\)
to have \(x=w\). The same argument for \(v\), using \(A_5\), and
for \(u\), using \(A_{17}\), forces their corresponding coordinates
to be \(w\). Finally the \(y\)-count forces every copy of \(A_{17}\)
to use its first term and every \(A_{12}\) to have \(y=w\).

The positions therefore partition uniquely into

\[
K,\ A,\ B,\ C,\qquad
|K|=r,\quad |A|=|B|=|C|=q,
\tag{7}
\]

carrying pieces \(12,4,5,17\), respectively.
The four sets of \(w\)-positions are

\[
S_x=K\cup A,\qquad S_y=S_u=K\cup C,\qquad S_v=K\cup B.
\tag{8}
\]

Before orbit averaging, each allocation (7) gives a block of rank
\(2^{3q}\) in either flattening: each of the \(3q\) core positions
is a rank-two map, and each common position has rank one.

## 3. Unique alignment preserves an entire color block

Call (7) good if, for every nonzero cyclic shift \(\delta\),

\[
|S_x\cap(S_y+\delta)|>r,\qquad
|S_v\cap(S_y+\delta)|>r.
\tag{9}
\]

At shift zero both intersections are \(K\), of size \(r\).
Thus (9) makes zero the unique relative shift at which either pair
can occur in a term of (5).
It also makes each of \(S_x,S_y,S_v\) aperiodic. For example, a
nonzero period of \(S_x\) would give another shift with intersection
size \(r\); the other cases follow in the same way.

In the row cut, a row pair consists of orbit labels for \(x,y\).
Its unique relative alignment recovers \(K,A,C\), and then \(B\),
up to a common cyclic rotation. The column pair \(u,v\) recovers the
same allocation. Different global rotation orbits of good allocations
therefore give disjoint row and column blocks, with no entries joining
two such blocks.

The column cut uses pairs \(x,u\) and \(y,v\). Because \(S_u=S_y\),
the same two uniqueness conditions prove the identical block statement.

All choices of the \(U\)-colors are retained within each block.
Indeed aperiodicity of the \(w\)-position sets gives a unique
representative for every colored-word orbit once its position pattern
is fixed. For a consistent row and column block, the four orbit sums
have exactly \(t\) valid simultaneous alignments. Their entries agree
by position invariance of the type coefficient.
Since \(t\) is odd, their sum over \(\mathbb F_2\) is one copy of the
original color matrix, not zero. Its rank is \(2^{3q}\).

If \(G_t\) is the number of good allocations (7), this proves in both
cuts

\[
\operatorname{rank}\ge\frac{G_t}{t}\,2^{3q}.
\tag{10}
\]

The factor \(t\) is exact: a good allocation is aperiodic, so its
global cyclic orbit has size \(t\).

## 4. Almost every balanced allocation is good

Fix \(K\) and a nonzero shift \(\delta\). The directed shift graph has
edges \(i\longrightarrow i-\delta\).
Because \(t\) is odd, its underlying graph is a union of cycles of
length at least three, with maximum degree two.
Deleting the \(r\) vertices of \(K\) removes at most \(2r\) edges.
A greedy matching removes at most three edges at each step, so the
remaining graph has at least \(q-1\) disjoint edges.

Color the \(3q\) remaining positions with the symbols \(A,B,C\).
On any selected directed edge, the ordered color pair \(A,C\)
contributes to \(S_x\cap(S_y+\delta)\).
If this intersection has size at most \(r\), at most \(r\) of the
selected \(M=q-1\) disjoint edges have that pair.
Among unrestricted three-colorings, the number satisfying this
necessary condition is at most

\[
3^{3q}\left(\frac89\right)^M
\sum_{j=0}^{\min(r,M)}\binom Mj8^{-j}.
\tag{11}
\]

The summation is at most \((M+1)^r\).
There are
\(T_q=(3q)!/(q!)^3\)
balanced colorings, and

\[
T_q\ge\frac{3^{3q}}{(3q+1)^2}.
\tag{12}
\]

For (12), the balanced multinomial coefficient is the largest of the
at most \((3q+1)^2\) possible count coefficients whose sum is \(3^{3q}\).
The same argument applies to the pair \(B,C\).
Union-bound over the two tests and \(t-1\) nonzero shifts.
Uniformly in the choice of \(K\), the fraction of bad balanced
allocations is at most

\[
\varepsilon_t=
2t(t+1)^4\left(\frac89\right)^{q-1}.
\tag{13}
\]

This tends to zero exponentially.
An optional explicit threshold is \(q\ge512\), when
\(\varepsilon_t<1/2\).
The checker proves the two integer inequalities

\[
4\cdot1538\cdot1539^4\cdot8^{511}<9^{511},
\qquad
8\cdot1539^5<9\cdot1536^5.
\tag{14}
\]

The first verifies the worst base value with \(q=512,r\le2\).
For each fixed \(r\), increasing \(q\) by one multiplies the polynomial
part of (13) by at most \((1539/1536)^5\), and its exponential part
by \(8/9\). The second inequality therefore proves the threshold
for all larger \(q\).

There are \(\binom trT_q\) allocations in total. Equations
(10)--(13) give the explicit lower bound

\[
r_{\rm row},r_{\rm col}
\ \ge\
\frac{1-\varepsilon_t}{t}
\binom tr\frac{(3q)!}{(q!)^3}\,2^{3q}.
\tag{15}
\]

For \(q\ge512\), replacing \(1-\varepsilon_t\) by \(1/2\) gives a
positive fully explicit bound. Since \(r\le2\),
the logarithm of (15) is \(t\log_2 6-O(\log(t+2))\).
Together with \(\log_2D_t=t\log_2 3-O(\log(t+2))\), this proves
(3)--(4).

## 5. Reproducible checks and remaining distinction

Run:

    python3 tools/check_cyclic_orbit_crossing.py \
      --output evidence/cyclic-orbit-crossing.json

The checker verifies the three local support identities in (6) on
all \(81\) index quadruples; constructs the complete one-copy and
three-copy coefficient profiles; checks the summed crossings; and
verifies four explicit good allocations at lengths \(27,31,39,41\).
For those allocations it independently enumerates all relative
alignments after fixing the first leg, finding exactly one.
It also checks the two exact threshold inequalities in (14).

The fixed eleven-label family remains an independently usable finite
seed. Its tensor powers, coefficient substitutions, and basis changes
are separate questions. Neither the finite maximum \(210\) nor the
large-cycle obstruction proves their asymptotic rank rate.
