# Ward entry evaluation controlled by support overlap

The scalar form memory can be bounded by the number of forms that
remain active across the current site. Combined with the
[eight-dimensional ancillary algebra](ward-stable-operator-interface.md),
this gives operator width \(8\cdot2^w\). A row-basis algorithm makes
the best possible \(w\) for any fixed site order explicit.

These are entry-representation and entry-evaluation results. They do
not bound the ranks of the formal coefficient matrices.

## 1. Active intervals and the exact transition

Let \(A=(w_{ai})\) be a \(k\times P\) form matrix over a
characteristic-two field. Suppose each row is nonzero. For row \(a\),
write \(f_a,\ell_a\) for its first and last nonzero sites. At site
\(i\), the active forms are
\[
 I_i=\{a:f_a\le i\le\ell_a\},
 \qquad w=\max_i|I_i|.
 \tag{1}
\]
Both endpoints are included. A gap inside a row's support does not
end that row's activity.

Use the [form-variable recurrence](ward-form-operator-transfer.md)
with local transition
\[
 B_i+\left(\sum_a w_{ai}z_a\right)(D_i^R+\eta D_i^O).
 \tag{2}
\]
At the first support of a form, introduce its selection bit with
value zero. Retain only bits for active forms. The ordinary term in
(2) keeps the current subset. An inserted term selects **one**
currently unselected form \(a\), with weight \(w_{ai}\).

All insertions at a site must use the old state vector. A transition
must not select one form and then select another at the same site.
This restriction is exact: the local generating factor is linear in
the form variables, and the physical operator satisfies
\((R_i+\eta O_i)^2=0\).

After processing the last support of a form, retain only states in
which its bit is one, then drop that bit. If several forms end at the
same site, all their bits must be one. If this cannot be achieved,
the contribution is zero. The eight ancillary channels are retained
through every introduction and removal of a form bit.

This proves the all-size width bound
\[
 \boxed{\qquad \text{operator width}\le8\cdot2^w.\qquad}
 \tag{3}
\]
Zero columns use the ordinary transition. If a form row is zero,
the entire form product is zero and needs no recurrence.

## 2. Entry evaluation and grade recovery

For a fixed stabilized source and target, the grade argument from
[the disjoint-support evaluator](ward-disjoint-operator-transfer.md)
still applies. The product of \(k\) forms has at most one power of
\(\eta\), namely
\[
 j=R_A^{\rm target}-R_A^{\rm source}-k.
 \tag{4}
\]
It is zero unless \(0\le j\le k\); the state-grade difference
\(t\) must also satisfy \(|t|\le k-j\) and
\(t\equiv k-j\pmod2\). Evaluation at \(\eta=1\) therefore
recovers its exact coefficient.

For supplied read-only first/last-support metadata and an ordered
stream of weights and basis indices, the active-interval algorithm
uses
\[
 O\!\left(P(1+w)2^w\right)
 \quad\text{field operations},
 \qquad O(2^w)\quad\text{field elements of workspace}.
 \tag{5}
\]
The constant ancillary factor is eight. Active row identifiers require
\(O(w\log k)\) additional bits; position and grade counters require
\(O(\log P)\) bits. This separates field-element storage from the
bit cost of representing those elements. Storage or preprocessing of
the input matrix and its endpoint metadata is not included in the
streamed workspace assertion.

At fixed \(w\), the entry computation is linear in \(P\) in field
operations. The statement concerns one operator entry, not a rank
calculation over all entries.

Coefficient extraction in form variables commutes with all aggregate
grade projections. The stable isomorphism preserves those grades,
and no ancillary reset is performed. Thus (3)--(5) retain the global
grade constraints of the Ward construction.

## 3. An optimal row basis for a fixed site order

Assume \(A\) has full row rank \(k\). Define inclusive prefix and
suffix column matrices
\[
 A_{\le i}=A[:,1,\ldots,i],\qquad
 A_{\ge i}=A[:,i,\ldots,P].
\]
Their union contains every column. Set
\[
 \omega_i=\operatorname{rank}A_{\le i}
            +\operatorname{rank}A_{\ge i}-k.
 \tag{6}
\]
Equivalently, \(\omega_i\) is the dimension of the intersection of
the two column spans. The current column appears on both sides.

**Theorem.** Among all invertible changes of the form-row basis, the
smallest possible active interval count at every site is
\(\omega_i\), simultaneously. In particular, the optimal width for
the fixed site order is
\[
 \boxed{\qquad w_*=\max_i\omega_i.\qquad}
 \tag{7}
\]

**Construction.** First perform row elimination from the left, so
that all nonzero rows have distinct first positions. Whenever two
rows share a last position, replace the row with the earlier first
position by a suitable multiple of the other row plus itself, chosen
to cancel their last entry.

The earlier first entry is unchanged, since the other row is zero
there and at every preceding position. The changed last position
strictly decreases. Thus the sum of last positions decreases at every
step, and the procedure terminates with distinct first positions and
distinct last positions. Every operation is invertible.

**Optimality.** In this final basis, the number of rows beginning by
site \(i\) equals \(\operatorname{rank}A_{\le i}\): their
restricted prefixes have distinct leading positions, and the other
rows have zero prefixes. Similarly, the number of rows ending at or
after \(i\) equals \(\operatorname{rank}A_{\ge i}\).
Every row belongs to at least one of these two sets, and their
intersection is the active set. Its size is exactly (6).

For any other row basis, the two row counts are at least the
corresponding ranks. Their intersection therefore has size at least
\(\omega_i\). The constructed basis attains this lower bound at
every site, proving (7). \(\square\)

The elimination uses ordinary finite-field operations and terminates
after at most \(kP\) trailing-end reductions. It may store the
transformed matrix during preprocessing. The later entry evaluator
uses that matrix and its endpoint metadata through read-only access.

This theorem optimizes the row basis for a **fixed column order**.
Choosing a better ordering of physical sites is a separate problem;
no optimal ordering algorithm is claimed here.

## 4. Row changes preserve the Ward product

Let \(Z_a=L_a+\eta O_a\). These operators commute and satisfy
\(Z_a^2=0\). If \(A'=CA\) for an invertible row matrix \(C\),
then expansion in this square-zero algebra gives
\[
 \prod_{a=1}^k Z'_a=\det(C)\prod_{a=1}^k Z_a.
 \tag{8}
\]
The original product is recovered by dividing by the known nonzero
determinant. For binary row changes over \(\mathbf F_2\), the
factor is one. The same argument applies to the annihilating form
operators \(Q_a\), so their product image, including each
homogeneous intersection, is unchanged.

If the form rows are linearly dependent, their product is zero:
after a row change, one row is zero. Equivalently, every full-size
minor in the site expansion vanishes. Thus the full-rank hypothesis
in the basis-optimization theorem loses no nonzero Ward product.

## 5. Exact controls

The checker
[check_ward_interval_transport.py](../tools/check_ward_interval_transport.py)
compares the full form ring, the active-memory recurrence, and scalar
evaluation at \(\eta=1\) on 800 entries. The cases include shared
endpoints, supports of one site, internal gaps, zero columns, and a
rank-deficient collision of two forms at one site. All 444 nonzero
entries satisfy the exact grade conditions. Seven full-rank cases
also check 700 entries after changing the row basis.

The row-basis checker verifies 1,854 full-rank matrices, including
complete small binary ranges. It checks all 35,409 invertible row
changes in those complete ranges against the claimed pointwise lower
bound. One recorded example starts with width three and reduces to
width one, while leaving all operator entries unchanged.

The evidence is
[ward-interval-transport.json](../evidence/ward-interval-transport.json).
The general proof is the elimination and coefficient-extraction
argument above; the finite tests validate the indexing, same-site
selection rule, boundary ordering, and explicit transport matrices.

