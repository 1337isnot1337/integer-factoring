# Sixteen scalar channels for disjoint Ward forms

For forms with disjoint supports, the form-subset memory in the
[eight-state generating recurrence](ward-form-operator-transfer.md)
can be removed. Order the physical sites by support and extract each
form variable when its block ends. This gives a scalar matrix-product
representation of width at most 16, independent of the number of forms.
The eight ancillary channels continue through every boundary.

The representation also evaluates an individual stabilized operator
entry with constant field-element workspace. It does not provide an
algorithm for the ranks of the coefficient matrices.

## 1. Process one form support at a time

Suppose the nonzero supports \(B_1,\ldots,B_k\) of the form rows are
pairwise disjoint. Weights within a support may be arbitrary elements
of a characteristic-two field. A simultaneous permutation of physical
sites makes these supports consecutive. This permutation preserves
all aggregate exterior grades.

At a site \(i\in B_a\), the form polynomial is
\(\ell_i(z)=w_{ai}z_a\). No other site block uses \(z_a\).
Therefore coefficient extraction factors as a sequence of linear
operations: process \(B_a\), retain its coefficient of \(z_a\),
then continue to \(B_{a+1}\). A site outside every support has an
ordinary transport transition and needs no selection flag.

For fixed source and target basis indices, let \(B_i,D_i\) be the
eight-by-eight scalar transition matrices read from the local
transport certificate, where
\[
 D_i=D_i^R+\eta D_i^O.
\]
Keep two eight-component vectors, \(v_0,v_1\), for degree zero or
one in the current form variable. The exact site update is
\[
 \boxed{\quad
 v_0'=B_i v_0,\qquad
 v_1'=B_i v_1+w_{ai}D_i v_0.
 \quad}
 \tag{1}
\]
At the start of a support block, set \((v_0,v_1)=(v,0)\). At its
end, set \(v=v_1\), and pass this entire eight-component vector to
the next block.

Only the form-selection flag is reset. In particular, the ancillary
state is not reset to the identity or contracted at a block boundary.
The continuing vector contains all eight channels produced by the
previous blocks. Formula (1) is consequently the same ordered stable
isomorphism as the full generating recurrence.

The update has width at most
\[
 \boxed{\qquad 2\dim\mathcal A=16.\qquad}
 \tag{2}
\]
This holds over \(F[\eta]\) before fixing \(\eta\). At the matrix
level, the same statement uses the sixteen-by-sixteen local emitted
operators in place of their individual scalar entries.

For binary forms of minimal quadratic-charge rank, the
[charge-rank classification](ward-code-charge-entropy.md) supplies an
invertible row change to this disjoint-support situation. The product
of the form operators changes by the determinant of that row change,
as follows from the site-minor expansion. Thus the construction also
applies after such a change of form basis.

## 2. Grade constraints are preserved

The local stable isomorphism preserves both aggregate exterior grades.
Coefficient extraction in a form variable commutes with every grade
projection, since that variable is scalar. Ordering the sites and
performing the block extractions therefore does not discard any
global grade restriction.

One may keep grade-generating variables through the whole recurrence,
or select the required matrix rows and columns afterward. This is
not a claim that each site's grade is preserved separately by the
stable isomorphism. No such stronger property is needed here.

## 3. A streamed entry evaluator

Consider the product of \(k\) form operators
\[
 W(\eta)=\prod_{a=1}^k(L_a+\eta O_a)
\]
in the stabilized basis. Each incidence factor raises the external
grade by one; each creation factor raises it by two and supplies one
power of \(\eta\). If a term uses \(j\) creation factors, its
external-grade increase is exactly \(k+j\).

Thus for a fixed source and target basis vector, an entry of
\(W(\eta)\) has at most one possible power:
\[
 \boxed{\quad
 W(\eta)_{\rm target,source}
   =c\eta^j,\qquad
 j=R_A^{\rm target}-R_A^{\rm source}-k.
 \quad}
 \tag{3}
\]
It is zero unless \(0\le j\le k\). If the state-grade change is
\(t\), it also vanishes unless
\[
 |t|\le k-j,\qquad t\equiv k-j\pmod2.
 \tag{4}
\]
These statements hold in the stabilized basis because the basis
change preserves the total grades, including the ancillary factor.
The physical operator acts trivially on that ancillary factor before
conjugation.

Consequently evaluation at \(\eta=1\) returns the exact coefficient
\(c\) in (3). No array of powers of \(\eta\) is needed. The
sixteen scalar channels in (1), with \(D_i=D_i^R+D_i^O\), evaluate
the entry directly. Grade counters accumulated while reading the
basis indices determine \(j\) and the zero conditions (4).

For an ordered stream of the \(P\) site weights, support boundaries,
and source/target basis digits, this takes \(O(P)\) field operations,
\(O(1)\) field elements of working storage, and \(O(\log P)\) bits
for position and grade counters. The fixed transition tables have
constant size. This accounting treats the input stream as read-only;
it does not count storage of the input weights or basis-index strings
as working memory. Field-operation cost and the bit length needed to
represent a field element are separate from the constant-element
workspace assertion.

The outside Ward factor \(\eta^k\) changes the final exponent to
\(k+j\) and does not change the coefficient evaluator.

At fixed numerical parameters, the same representation can be
composed locally with
\(G\otimes F^{\otimes2P}\) by multiplying each emitted operator
matrix by the corresponding local kernel matrix. This retains a
constant number of scalar channels. Extracting formal bivariate
coefficients and calculating their large matrix ranks remain separate
tasks; the entry evaluator does not settle them.

## 4. Verification

The checker
[check_ward_disjoint_transport.py](../tools/check_ward_disjoint_transport.py)
compares the full square-free form-ring calculation with (1) on 800
entries over five support profiles, with up to four forms and nine
sites. It then compares the streamed scalar evaluation at
\(\eta=1\) with every full polynomial entry.

There are 537 nonzero controls. Each has exactly one power of
\(\eta\), and all satisfy the exact grade conditions (3)--(4).
The remaining controls verify cancellation to zero. The evidence is
[ward-disjoint-transport.json](../evidence/ward-disjoint-transport.json).

The proof of (1) is coefficient extraction at successive support
boundaries; the finite controls verify the indexing and the absence
of an unintended ancillary reset. Equations (3)--(4) then give the
entry-complexity statement for all support sizes and numbers of forms.

