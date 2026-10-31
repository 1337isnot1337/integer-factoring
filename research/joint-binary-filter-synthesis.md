# Residue-compressed synthesis for joint prime-state identity filters

## Result and scope

The [prime-power identity theorem](prime-power-joint-identity-filters.md)
turns an arbitrary joint diagonal filter on \(r\) selected binary
identity directions into \(2^r\) scalar finite-difference polynomials.
Its degree/order feasibility problem can be solved using at most
\(2^{r+1}+1\) residue patterns, regardless of the value of \(nP\).
The local systems depend on degree and order **quotients** by \(nP\).

A fixed coefficient extension of order \(2^{r+1}\) permits arbitrary
nonzero endpoint colors. Exact feasibility then reduces to checking
that no endpoint coordinate is forced to zero. The tool constructs a
literal polynomial witness and checks every degree and zero-order
bound. It does not infer feasibility from a numerical pressure value.

This provides a faster constructive search interface for correlated
three- and four-mode filters, including nonconstant full-occupation
weights. No improved factoring exponent is claimed.

## 1. The joint transform

Let \(n\) be a power of two, \(P\) a positive odd integer,
\(m=nP\), and \(\mathcal A=\{0,1\}^r\). For each
\(\mathbf d\in\mathcal A\), let \(\alpha_{\mathbf d}(x)\)
be a polynomial over a fixed binary extension field. Put
\(C=1+x^m\). The exact transform is

\[
\beta_{\mathbf k}(x)
 =\sum_{\mathbf h\le\mathbf k}
 C^{|\mathbf h|}
 \alpha_{\mathbf1-\mathbf k+\mathbf h}(x).
\tag{1}
\]

Thus the coefficient matrix is the \(r\)-fold tensor power of

\[
\begin{pmatrix}0&1\\1&1+x^m\end{pmatrix}.
\tag{2}
\]

For prescribed nonnegative bounds \(D_{\mathbf d},B_{\mathbf k}\),
the synthesis problem is

\[
\deg\alpha_{\mathbf d}\le D_{\mathbf d},\qquad
\mathrm{ord}_x\beta_{\mathbf k}\ge B_{\mathbf k},\qquad
\alpha_{\mathbf d}(1)\ne0.
\tag{3}
\]

A common Laurent monomial clearing reduces any finite Laurent instance
to this form. The clearing shifts every \(D\) and every \(B\) equally;
the two physical profile shifts cancel in their product cost.

## 2. At most \(2^{r+1}+1\) residue patterns

Write \(X=x^m\). Each polynomial has the unique decomposition

\[
\alpha_{\mathbf d}(x)
=\sum_{\eta=0}^{m-1}x^\eta A_{\mathbf d,\eta}(X).
\tag{4}
\]

Because (1) has coefficients in \(\mathbf F_2[X]\), distinct residues
do not interact. For residue \(\eta\), the local bounds are

\[
d_{\mathbf d}(\eta)
 =\left\lfloor\frac{D_{\mathbf d}-\eta}{m}\right\rfloor,
\qquad
b_{\mathbf k}(\eta)
 =\max\left(0,\left\lceil\frac{B_{\mathbf k}-\eta}{m}\right\rceil\right).
\tag{5}
\]

Negative degree bounds mean the corresponding local polynomial is zero.
The tuple in (5) changes only at

\[
D_{\mathbf d}\bmod m+1,
\qquad B_{\mathbf k}\bmod m,
\tag{6}
\]

besides the endpoints zero and \(m\). There are \(2^r\) coordinates
of each kind, so there are at most \(2^{r+1}+1\) nonempty intervals
with distinct local bound patterns.

For each pattern, form its exact homogeneous coefficient system over
\(\mathbf F_2\), solve its nullspace, and map that nullspace to
the endpoint vector \((A_{\mathbf d,\eta}(1))_{\mathbf d}\).
The endpoint image of the full problem is the sum of these image spaces.
All residues in one interval have the same image, so keeping one
representative residue preserves the image exactly. Repeated image
spaces do not enlarge their sum.

This proves the compression used in
[`joint_binary_filter_synthesis.py`](../tools/joint_binary_filter_synthesis.py).
The coefficient lifts retained during row elimination give explicit
polynomials after feasibility is established.

## 3. A fixed extension gives an exact nonzero-endpoint criterion

Let \(E\subseteq\mathbf F_2^{2^r}\) be the computed endpoint image.
Extend coefficients to any fixed field
\(K=\mathbf F_{2^{r+1}}\). Then

\[
E\otimes_{\mathbf F_2}K
\text{ contains a vector with every coordinate nonzero}
\quad\Longleftrightarrow\quad
\text{every coordinate projection of }E\text{ is nonzero}.
\tag{7}
\]

Necessity is immediate. For sufficiency, the coordinate-zero conditions
are at most \(2^r\) proper hyperplanes in the extended image. Their union
has fewer than \(|K|^{\dim E}\) points because \(2^r<|K|\).

The implementation also supplies a direct construction. It adjoins
binary endpoint basis vectors one at a time, choosing a scalar that
preserves the current nonzero coordinates and makes the newly reached
coordinates nonzero. Fewer than \(|K|\) scalars are forbidden at each
step. The coefficient lifts are combined using the same scalars.

Only addition and multiplication by prime-field coefficients occur in
(1), the feasibility equations, and the beta-order replay. Hence the
tool represents \(K\) in any fixed additive binary basis. A physical
endpoint repair additionally uses the inverse of each resulting nonzero
constant; this exists in the fixed finite field and contributes no
polynomial degree or Smith order.

The [independent checker](../tools/check_joint_binary_filter_synthesis.py)
constructs the full, uncompressed coefficient systems from the tensor
matrix (2). Its [record](../evidence/joint-binary-filter-synthesis.json)
contains 288 exact feasibility comparisons for two, three, and four
modes, at residue scales 1, 2, 7, and 32, plus 32 complete polynomial
witness replays. The full and compressed systems agree in every case.

## 4. The index dependence is an explicit finite expectation

This identity also applies to \(p\)-state filters in characteristic
\(p\), using the general transform from the prime-power theorem.
For selected indices \(j_1,\ldots,j_r\), define

\[
\lambda_\ell=2n-1-2j_\ell,\qquad t=P\log s,
\qquad
\mu_{\lambda}(\mathbf k)
=\frac{\exp(-t\lambda\cdot\mathbf k)}
 {\prod_{\ell=1}^r\sum_{d=0}^{p-1}\exp(-t\lambda_\ell d)}.
\tag{8}
\]

Let \(B_{\mathbf k}\) be the exact \(x\)-order of the corresponding
beta polynomial and \(D_{\mathbf d}\) the full degree of the
corresponding alpha polynomial. Set

\[
U_{\mathbf k}
=\exp\!\left(2t\left[n|\mathbf k|-B_{\mathbf k}/P\right]\right),
\qquad
V_{\mathbf k}
=\exp\!\left(2tD_{\mathbf p-\mathbf1-\mathbf k}/P\right).
\tag{9}
\]

The exact ratio of the filtered two-cut cost to the Gaussian cost at
the same charge \(J\) is

\[
\boxed{
R=\mathbb E_{\mu_{\lambda+J\mathbf1}}U\;
   \mathbb E_{\mu_\lambda}V.
}
\tag{10}
\]

To derive (10), factor the unchanged selected-Y digit partitions from
the zero cut. Its remaining Gaussian energies are
\(\lambda_\ell+J\). At infinity, replace each digit \(d\) by
\(p-1-k\); the common ground factor cancels between the filtered and
Gaussian partitions. What remains is precisely (8)--(10).

The coefficients and their ordered \((B,D)\) profile are independent
of the selected indices. Equation (10) therefore permits a single
profile to be evaluated over every selected subset. Permuting the
filter modes gives every assignment to that subset. For a physical
witness, the \(\lambda_\ell\) must be distinct positive odd integers
in \([1,2n-1]\), arising from actual indices. Continuous optimization
is an exploratory selection aid only.

## 5. The same construction in every prime characteristic

Let \(p\) be prime, let \(r\ge1\), and replace the binary state set
by \(\{0,\ldots,p-1\}^r\). The prime-power identity theorem gives

\[
\beta_{\mathbf k}(x)
=\sum_{\mathbf0\le\mathbf h\le\mathbf k}
(-1)^{|\mathbf k|-|\mathbf h|}
\left(\prod_{\ell=1}^r\binom{k_\ell}{h_\ell}\right)
(1-x^m)^{|\mathbf h|}
\alpha_{\mathbf p-\mathbf1-\mathbf k+\mathbf h}(x).
\tag{11}
\]

Every coefficient lies in \(\mathbf F_p[x^m]\). Therefore the residue
decomposition and local bounds (4)--(6) remain valid verbatim. There are
now at most \(2p^r+1\) patterns. Row elimination takes place over
\(\mathbf F_p\). The endpoint image has \(p^r\) coordinates, so the
fixed field \(K=\mathbf F_{p^{r+1}}\) is strictly larger than that
coordinate count. The proper-hyperplane argument and constructive
support extension prove the same exact feasibility criterion.

At \(x=1\), (11) gives
\(\beta_{\mathbf k}(1)=(-1)^{|\mathbf k|}
\alpha_{\mathbf p-\mathbf1-\mathbf k}(1)\).
Thus every beta polynomial in an admissible witness is nonzero. Its
order is finite, and the exact physical profile formulas apply without
a missing-state convention.

[`joint_prime_filter_synthesis.py`](../tools/joint_prime_filter_synthesis.py)
implements this construction. Its coefficient representation stores
the \(r+1\) prime-field coordinates of an element of \(K\).
Equation (11), elimination, and witness reconstruction use only
addition and multiplication by elements of \(\mathbf F_p\); they do
not require a choice of irreducible polynomial. Any fixed field model
with this additive basis supplies the endpoint inverses for the
physical diagonal repair.

The [independent prime checker](../tools/check_joint_prime_filter_synthesis.py)
builds the full coefficient system by tensoring the one-mode transform.
Its [record](../evidence/joint-prime-filter-synthesis.json) contains
100 exact feasibility comparisons and 23 componentwise tensor-transform
replays over additive models of \(\mathbf F_9,\mathbf F_{27},\mathbf F_{25}\).
These are checks of the generic construction. The finite search
catalogues currently use binary three- and four-mode filters and
qutrit two-mode filters.

## 6. Search records and limits

[`search_joint_binary_filter_synthesis.py`](../tools/search_joint_binary_filter_synthesis.py)
maximizes the beta-order targets in several fixed coordinate priorities.
Each reported target is feasible and comes with a deterministic
polynomial lift. Every actual degree and beta order is replayed before
the profile is retained. The complete ordered profile catalogue contains
degree bounds, order bounds, endpoint values, and a hash of the literal
witness, so each profile can be reconstructed directly.

This finite Pareto search does not enumerate all possible order targets
or all polynomial coefficients. A missing gain does not rule out other
targets, other degree caps, or other gauges. The all-index pressure
driver records its complete finite subset and scale/charge grids, exact
profile counts, catalogue hashes, and representative profiles. Its
floating comparisons remain search evidence rather than an all-scale
obstruction or an improved factoring theorem.

The [search record note](joint-filter-search-records.md) lists the
completed catalogues, their replay commands, the all-index pressure
screens, and interrupted calculations that are excluded from the
reported results. The
[catalogue checker](../tools/check_joint_filter_catalogues.py)
reconstructs every retained literal witness from its stored bounds and
verifies its exact orders, full degrees, endpoints, and polynomial hash.
