# Two-state formal Ward queries for the fan

The [duplicated fan family](ward-fan-family.md) admits a division-free
physical entry evaluator with only two accumulated ring elements. It
computes the entire formal polynomial, using \(O(P)\) operations in
the coefficient ring. It also gives a physical operator representation
of bond dimension two across duplicated-edge groups.

This is a stronger entry representation for this structured family.
It does not give the row or column ranks of its restricted coefficient
tensors. Section 5 specifies a small exact rank computation that the
representation makes accessible.

## 1. Collapse parallel columns

Use the scalar problem from the
[universal Pfaffian query theorem](ward-pfaffian-entry-queries.md):
\[
 [z_1\cdots z_k]\prod_i(m_i+b_i\ell_i(z)),
 \qquad \ell_i(z)=\sum_a A_{ai}z_a.
 \tag{1}
\]
All calculations in this note work over a commutative
characteristic-two ring. In the squarefree form-variable ring,
every linear form satisfies \(\ell^2=0\).

If several sites have the same form column, their factors therefore
collapse to
\[
 \prod_{i\in e}(m_i+b_i\ell_e)=M_e+B_e\ell_e,
 \qquad
 M_e=\prod_{i\in e}m_i,\quad
 B_e=\sum_{i\in e}b_i\prod_{\substack{j\in e\\j\ne i}}m_j.
 \tag{2}
\]
No inverses are used. For the two copies in the fan,
\[
 M_e=m_1m_2,\qquad B_e=b_1m_2+m_1b_2.
 \tag{3}
\]
For any number of copies, retain a pair \((M,B)\), initialized to
\((1,0)\), and update it by
\[
 (M,B)\longmapsto(Mm,\;Bm+Mb).
\]
This remains valid when some or all of the factors are zero or
nilpotent.

## 2. The selected type edges form a spanning tree

For the fan, the distinct form columns are the reduced incidence
columns of its spokes and path edges. A set of exactly \(k\) such
columns has determinant one precisely when its edges form a
spanning tree on vertices \(0,1,\ldots,k\). If the edges are
disconnected, the incidence rank is smaller than \(k\). A connected
set of \(k\) edges is a tree, and deleting a leaf proves its reduced
incidence determinant is one.

Consequently (1) is the weighted spanning-tree polynomial
\[
 \boxed{\quad
 \sum_{T\ {\rm spanning\ tree}}
       \prod_{e\in T}B_e\prod_{e\notin T}M_e.
 \quad}
 \tag{4}
\]
The same identity applies to any graph whose reduced incidence
columns give the form weights, with the expected full-rank component
condition.

There are two different graph enumerations in this construction.
The admissible **code sectors** in the dimension formula use all
connected spanning subgraphs, with any number of selected types.
The **Ward operator entry** in (4) selects exactly \(k\) types, so
only spanning trees occur.

## 3. Two exact transfer matrices

Process the fan in the order
\[
 s_1,\ p_1,\ s_2,\ p_2,\ldots,p_{k-1},s_k,
\]
where \(s_j\) is the spoke \(0j\) and \(p_j\) is the path edge
\(j(j+1)\). Retain two partial-forest states. In state \(C\), all
processed vertices lie in the component of \(0\). In state \(U\),
the only component not connected to \(0\) contains the current
vertex. Components that can no longer meet a later edge are
discarded.

A path edge introduces the next vertex. If its predecessor is
pending, the path edge is required. A spoke can close a pending
component, but selecting it when the current vertex is already
connected would create a cycle. Thus the transfer matrices are
\[
 T_p=
 \begin{pmatrix}B_p&0\\M_p&B_p\end{pmatrix},
 \qquad
 T_s=
 \begin{pmatrix}M_s&B_s\\0&M_s\end{pmatrix}.
 \tag{5}
\]
Start with the vector \((C,U)=(0,1)\), apply the transfers in the
displayed edge order, and return \(C\). Equivalently, after the
first spoke the state is \((B_{s_1},M_{s_1})\), and each new
path/spoke pair updates it by
\[
 \begin{aligned}
 C'&=(B_sM_p+M_sB_p)C+B_sB_pU,\\
 U'&=M_sM_pC+M_sB_pU.
 \end{aligned}
 \tag{6}
\]
Every retained final object is a connected acyclic spanning
subgraph, and every spanning tree has exactly one such sequence
of choices. This proves (4)--(6).

The evaluator performs a constant number of ring operations at
each physical site, and retains a constant number of ring
elements. Polynomial ring elements can grow with \(P\); this is
not a claim of constant bit storage. For a fixed number of formal
parameters and a supplied polynomial degree bound, ordinary
polynomial representations give polynomial bit cost.

There is also an operator version. On each duplicated-edge pair,
replace the scalar \(M_e,B_e\) by
\[
 \begin{aligned}
 M_e&=S_1\otimes S_2,\\
 B_e&=(T_1S_1)\otimes S_2+S_1\otimes(T_2S_2),
 \qquad T_i=R_i+\eta O_i.
 \end{aligned}
 \tag{7}
\]
Operators for distinct edge groups act on disjoint physical
spaces. The same transfers therefore represent the full physical
operator \(W(\eta)M_P\), with bond dimension two between
duplicated-edge groups. Inside a two-site group, retaining an
additional selected/not-selected flag gives bond dimension at
most four. This is a representation along the site order; the
cut ranks used in the factoring objective separate tensor legs
at every site and are different quantities.

## 4. An explicit sparse basis for the first new fan code

Consider \(k=3,P=L=10\), with its 84 orientation-selected labels.
There are eight admissible type sets of size three and five of
size four; the full size-five sector has zero syndrome and is
omitted by the selector.

Write \(x_{e,0},x_{e,1}\) for the even-position variables of the
two copies of a selected edge, and set
\[
 l_e=x_{e,0}+x_{e,1}.
\]
The creation forms are the three independent linear combinations
\(\Omega_a=\sum_e A_{ae}l_e\).

For a size-three type set, the graph is a tree. The code basis
vector is \(\Omega_1\Omega_2\Omega_3\), its homogeneous Ward
preimage is \(1\), and its support has eight occupation monomials.

For a size-four type set, the graph has a unique cycle \(C\) of
length three or four. Choose one cycle edge \(e_*\). Two residual
linear forms are
\[
 v_1=l_{e_*},\qquad v_2=\sum_{e\in C}x_{e,1}.
 \tag{8}
\]
They belong to the common kernel of the annihilation forms.
Modulo the span of the \(\Omega_a\), they are independent: the
first supplies the one missing direction in the span of the
\(l_e\), and the second has a nonzero component in the cycle
direction outside that span. The two code vectors are
\[
 I_j=\Omega_1\Omega_2\Omega_3v_j,\qquad j=1,2.
 \tag{9}
\]
Their homogeneous Ward preimages are \(v_j\). The zero-Gram
condition lets the annihilation operators commute through the
creation forms, and each annihilation operator kills \(v_j\).
Thus the product of the three \(Q_a\) acts on \(v_j\) exactly as
the creation product in (9).

The preimage support sizes are \(2,3\) on a triangle and \(2,4\)
on a four-cycle. The input support sizes are \(16,32\) and
\(16,56\), respectively. In every case each input vector has
an occupation coordinate outside the other's support. Selecting
one such coordinate from each gives a coordinate output map
\(\Pi I=I_2\), without a further basis change. The checker stores
these exact pivots.

After the single-mode orientation multiplicities are included,
the 84 labels have a total of 116 preimage monomials and 1,040
input monomials. Their preimages all have physical grade four,
as required by \(L-2k=4\).

## 5. A bounded route to restricted coefficient ranks

Let \(\sigma\) record which of the five type groups are singly
occupied. The code has thirteen nonzero syndromes: eight sectors
of dimension eight and five of dimension four. For a pair
syndrome \(\tau\), the dimension of a logical cut block is
\[
 n_\tau=\sum_\sigma d_\sigma d_{\sigma+\tau}.
 \tag{10}
\]
The exact convolution gives 31 nonempty blocks, with largest
dimension 592. The number of coefficient entries across the
blocks is
\[
 \sum_\tau n_\tau^2=1,917,184,
 \qquad 84^4=49,787,136.
 \tag{11}
\]
These counts specify the restriction before any rank computation.

A small starting block is \(\tau=7\), in the type order
\((e_1,e_1+e_2,e_2,e_2+e_3,e_3)\). Its dimension is 128. The
explicit sparse input/output maps above reduce each restricted
entry to a finite sum of the physical polynomial queries
(4)--(6), through the exact Ward factorization
\[
 G=y^3\Gamma,\qquad p=1+x,\quad q=1+y,\quad r=1.
\]
In particular,
\[
 [x^3y^3]G=[x^3]\Gamma(x,0).
 \tag{12}
\]
For this coefficient set \(q=r=1,\eta=0\) in the physical
oracle. The remaining polynomials are univariate; no
multivariate evaluation grid is needed. The homogeneous
factor \(\Gamma\) has degree seven in \(p,q,r\). Its entries
after this substitution therefore have degree at most seven.

For \(\tau=7\), the sparse contractions use at most 458,752
physical polynomial queries before any caching. Equal local
entry patterns and exchanges of the two copies of an edge can
be cached exactly. This is a bounded \(128\times128\) rank
target, rather than an assembly of the full tensor.

If its exact rank in one cut is \(r_\tau<128\), the complete
coefficient rank in that cut has the valid upper bound
\[
 84^2-128+r_\tau=6928+r_\tau.
 \tag{13}
\]
All uncomputed syndrome blocks are bounded by their dimensions.
The current complete-rectangle bound at coefficient \((3,3)\)
is the cap \(84^2=7056\), so any verified deficiency improves
that coefficient bound. The same procedure applies independently
to the other cut. More generally, for any collection of computed
blocks \(\mathcal T\),
\[
 r_\gamma\le
 \min\left\{r_\gamma^{\rm old},
 d^2-\sum_{\tau\in\mathcal T}(n_\tau-r_{\gamma,\tau})\right\}.
 \tag{14}
\]
An exact coefficient oracle is essential here: a rank measured
after substituting numeric values for the formal parameters is
not an upper bound for the requested coefficient rank.

The package records the sparse interface and this computation
target. It does not claim that the target block is deficient or
that an improved rank profile has already been obtained.

## 6. Exact controls

The evaluator is
[ward_fan_entry_queries.py](../tools/ward_fan_entry_queries.py).
Run

~~~sh
python3 tools/check_ward_fan_entry_queries.py --output /tmp/ward-fan-entry-queries.json
~~~

The [checker](../tools/check_ward_fan_entry_queries.py) verifies
4,112 exhaustive binary cases; 600 three-way comparisons over
\(\mathbf F_4,\mathbf F_8,\mathbf F_{16}\) against direct
squarefree multiplication and the universal Pfaffian engine;
90 full bivariate polynomial cases; and 160 cases in the
nonreduced ring \(\mathbf F_2[\epsilon]/(\epsilon^2)\).

It also verifies every sparse preimage and coordinate pivot
needed in the thirteen sectors above, computes the complete
syndrome convolution, and records the uncached contraction
counts. The [evidence](../evidence/ward-fan-entry-queries.json)
pins all source hashes. No dense restricted coefficient matrix
is constructed by this checker.
