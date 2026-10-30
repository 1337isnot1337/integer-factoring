# A PBW degeneration closes the graded Hopf crossing lane

## Main result

Every coefficient of the restricted-enveloping Hopf construction has both
scalar cut ranks at least those of an ordinary weighted abelian mixing
coefficient. The same statement holds for every scalar recombination,
including parameter-dependent local gauges whose leading PBW map is the
identity.

This has a stronger consequence than the earlier coefficient calculations.
For every exact crossing decomposition within this coefficient span,

\[
\boxed{\qquad
\sum_j\sqrt{r(T_j)c(T_j)}\ \ge\ d^{3/2}.
\qquad}
\tag{1}
\]

The bound holds at every size, over every characteristic-two extension
field, with arbitrary piece count and arbitrary positive primitive weights.
It follows from a four-leg degeneration and the
[full conservation-kernel theorem, pinned in PR27](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/full-conservation-rank-sum-barrier.md).

Because \(3/2>\sqrt2\), these families cannot satisfy the finite-seed
criterion \(\sum_j\sqrt{r_jc_j}<d^{\sqrt2}\). Polynomial-size families
also require \(\rho\ge3/4\) in the leading width interface.
This is a restriction on the specified construction, not a lower bound
on factoring or on arbitrary crossing tensors.

The earlier [graded Hopf construction](noncommutative-crossing-graded-hopf.md)
and its finite ranks remain correct. The present theorem explains why
improving those ranks within this lane cannot improve the target score.

## 1. The filtration is ordinary PBW length

Let \(H\) be one of the finite-dimensional restricted enveloping algebras
under discussion, with primitive PBW generators \(g_1,\ldots,g_\lambda\),
and square-free ordered PBW basis \(g^a\), \(a\in\{0,1\}^{\lambda}\).
The argument also permits square relations whose right side has smaller
PBW length.

Define

\[
\ell(g^a)=|a|,\qquad
F_bH=\operatorname{span}\{g^a:\ell(g^a)\le b\}.
\tag{2}
\]

This length gives every primitive generator weight one. It is separate
from the positive weights \(w_i\) used in the deformation parameter.

Commuting two generators replaces their commutator by a linear Lie
generator. A square relation likewise lowers PBW length. Hence

\[
F_aH\,F_bH\subseteq F_{a+b}H.
\tag{3}
\]

The primitive coproduct preserves total PBW length exactly; the antipode
reverses words and never increases length. The associated graded Hopf
algebra is

\[
\operatorname{gr}_F H
 \cong \mathbb F[g_1,\ldots,g_\lambda]/(g_1^2,\ldots,g_\lambda^2),
\tag{4}
\]

with commuting primitive generators.

The standard dilation
\(\phi_\tau(g^a)=\tau^{\sum_iw_i a_i}g^a\) preserves PBW length.
The argument needs only this diagonal linear map and its endpoint
\(\phi_0=\eta\epsilon\). Compatibility of the \(w_i\) with the Lie
brackets is needed when calling it a Hopf endomorphism, but is not needed
for the filtered rank comparison below.

## 2. A single degeneration works for both cuts

Use the Hopf formula

\[
R_\tau(u\otimes v)
=\sum v_{(1)}a_{(1)}\otimes u_{(1)}a_{(2)},\qquad
a=\phi_\tau(u_{(2)}S(v_{(2)})).
\tag{5}
\]

Every operation in (5) preserves or lowers total PBW length. Thus a
nonzero tensor entry satisfies

\[
\ell(x)+\ell(y)\le\ell(u)+\ell(v).
\tag{6}
\]

Introduce a new formal variable \(s\), independent of \(\tau\), and let

\[
D_s(g^a)=s^{\ell(g^a)}g^a.
\]

Over \(\mathbb F(s)\), define

\[
R_{\tau,s}=
(D_s^{-1}\otimes D_s^{-1})R_\tau(D_s\otimes D_s).
\tag{7}
\]

Its \((x,y,u,v)\) entry is multiplied by

\[
s^{\,\ell(u)+\ell(v)-\ell(x)-\ell(y)}.
\tag{8}
\]

The exponent is nonnegative by (6), so (7) is polynomial in \(s\).
At \(s=0\), only the length-preserving terms remain. Applying the
associated-graded operations in (4) gives the ordinary abelian mixing
operator \(M_\tau\).

For \(s\ne0\), (7) consists of one invertible linear transformation on
each of the four tensor legs. It therefore preserves rank in both
\((x,y)\mid(u,v)\) and \((x,u)\mid(y,v)\), despite their different
arrangements. Specialization at \(s=0\) cannot increase matrix rank:
any nonzero minor of the specialization is a nonzero polynomial minor
before specialization. Scalar extension to \(\mathbb F(s)\) also
preserves matrix rank.

Consequently, for every polynomial coefficient or scalar recombination
\(T\), with corresponding leading tensor \(\overline T\),

\[
r(T)\ge r(\overline T),\qquad
c(T)\ge c(\overline T).
\tag{9}
\]

The degeneration is common to all pieces. No different summand is
selected for the two cuts.

Equivalently, the row matrix is block triangular by total PBW length.
For the column matrix, rearrange (6) as

\[
\ell(x)-\ell(u)\le\ell(v)-\ell(y).
\]

Ordering rows and columns by these differences gives the corresponding
block-triangular explanation. The four-leg formulation (7) avoids
having to coordinate separate minor selections.

## 3. The leading family lies in the full conservation space

On one primitive coordinate of weight \(w\), the abelian mixing map is

\[
g_L\mapsto\tau^w g_L+(1+\tau^w)g_R,\qquad
g_R\mapsto(1+\tau^w)g_L+\tau^w g_R.
\tag{10}
\]

For the four binary basis labels, put
\(S=x+y=u+v\) and \(b=x+u\). The local kernel is

\[
\begin{array}{c|cc}
 &b=0&b=1\\ \hline
S=0&1&0\\
S=1&\tau^w&1+\tau^w
\end{array}.
\tag{11}
\]

The complete leading tensor is the tensor product of these local
kernels. Therefore every coefficient of \(M_{1+z}\), and every scalar
combination of those coefficients, has the form

\[
T_K(x,y,u,v)=[x+y=u+v]K(x+y,x+u)
\tag{12}
\]

on the group \(\mathbb F_2^\lambda\).

If \(\sum_jT_j=E_d\), applying the common degeneration gives
\(\sum_j\overline T_j=E_d\). The conservation theorem applies to that
entire decomposition, without position symmetry or a degree restriction:

\[
\sum_j\sqrt{r(\overline T_j)c(\overline T_j)}\ge d^{3/2}.
\tag{13}
\]

Combining (9) and (13) proves (1).

The pinned theorem proves (13) using actual independent row and column
functions of the kernels. Double Cauchy--Binet selection gives matching
per-piece counts in two bases. Boolean-degree convolution bounds and
the minimum total degree of a Boolean function basis then give the
rank-sum estimate. In particular, it does not require that the leading
pieces be the original monomial coefficients.

This proves an obstruction even when all primitive weights are distinct.
A large number of weight classes does not evade (1).

## 4. Scalar changes, gauges, and formal coefficients

All scalar recombinations use the same degeneration variable \(s\).
Their scalars may lie in any characteristic-two extension field.
Linearity commutes with taking the leading tensor, so no coefficient-basis
change can avoid (9)--(13).

Consider a local dual substitution

\[
z_j\mapsto z_j+a(\tau)Q_j(x),\qquad a(0)=0,
\tag{14}
\]

where \(Q_j\) is quadratic in lower primitive coordinates.
In the dual PBW basis it raises ordinary monomial length whenever the
correction is used. Its primal map therefore lowers PBW length, with
identity associated graded. The four independent local gauges used in
the bounded search are all of this form.

After any combination of these gauges, the same filtered degeneration
still has leading family \(M_\tau\). The endpoint remains crossing
because each gauge is the identity at \(\tau=0\). Thus the theorem
covers these gauge families as well.

More generally, the argument covers filtered local modifications whose
associated-graded maps are fixed legwise invertible transformations of
the same conservation family. Arbitrary parameter-dependent leading
maps require a new argument and are outside this statement.

One may also use scalar coefficients in a Laurent-series field
\(\mathbb F((t))\), with \(s\) independent of \(t\). The same rank
comparison holds over that field. If the summed tensors are integral
and specialize to crossing, so do their leading conservation kernels.
Their kernel sum is invertible over \(\mathbb F((t))\), and Section 4
of the pinned conservation theorem gives the same generic rank-sum
bound. This does not cover arbitrary tensors outside the filtered
coefficient space.

## 5. The earlier product-class consequence

There is also a weaker route through position symmetry, which helps
connect the result to prior work. Suppose \(m_1,\ldots,m_b\) are the
multiplicities of the distinct primitive weights. The leading pieces
are invariant under independent permutations within those classes.

Section 6 of the
[pinned binary-word theorem](https://github.com/1337isnot1337/integer-factoring/blob/023e4fc4482d72f25738f6b1d0f9cd31583baed1/research/modular-binary-word-barrier.md)
uses a product of common functionals on dyadic position blocks. Apply
that construction independently inside each weight class and take
one product over all resulting blocks. The same detected summand
satisfies both rank lower bounds. With \(\lambda=\sum_jm_j\), an explicit
loss is

\[
E=\frac92 b+
8\sum_{j=1}^b(1+\log_2m_j)\log_2(m_j+4),
\tag{15}
\]

giving ranks at least \(2^{3\lambda/2-E}\) in both cuts.
The omitted at most three bits per class account for the first term.

Thus \(\sum_j\log^2(m_j+4)=o(\lambda)\) implies \(\rho\ge3/4\).
All two-step examples have two weight classes. For the strictly upper
triangular \(r\)-by-\(r\) Lie algebra, the weight-\(j\) class has size
\(r-j\), and the loss is \(O(r\log^2r)=o(r^2)\).

This source was fetched at the stated commit and its product-group
argument was checked. The stronger conservation argument above removes
the class-size condition entirely.

## 6. Scope and exact verification

The theorem closes the graded two-step Hopf family, its general
alternating commutator forms, the local-ring examples, and the graded
strictly upper triangular examples. It closes their stated PBW-filtered
gauges and scalar coefficient spans. It does not, by itself, apply to
nonabelian group-word families or arbitrary nonlinear control
permutations. A separate compatible filtration is required there.

The label domain here is the full PBW basis. Arbitrary monomial subsets
or non-graded linear retracts are not automatically covered: the common
scaling \(D_s\) need not preserve the selected label spaces, and their
leading tensors need not satisfy the full conservation theorem's hypotheses.
A compatible filtered restriction requires its own verification of those
hypotheses. The present result is not a prohibition on every restricted
or retracted Hopf crossing.

The checker constructs the actual coefficient matrices for:

- One Heisenberg pair and two pairs with shared center.
- The irreducible field pencil and the nilpotent Jordan pencil.
- A parameter-dependent local gauge.
- The restricted upper triangular algebra of order four.

It checks every nonzero entry for (6), verifies that the length-preserving
part equals the abelian coefficient tensor entrywise, and checks both
rank inequalities. It also checks all binary recombinations in the
smallest example and selected cancellations in the larger examples.
The all-size conclusion comes from the ordinary proof, not extrapolation.

Run:

    python3 tools/check_noncommutative_crossing_pbw.py \
      --output evidence/noncommutative-crossing-pbw.json

The evidence records the pinned source commits and returned Git blob
identifiers. The complete checker is
[check_noncommutative_crossing_pbw.py](../tools/check_noncommutative_crossing_pbw.py),
with results in
[noncommutative-crossing-pbw.json](../evidence/noncommutative-crossing-pbw.json).
