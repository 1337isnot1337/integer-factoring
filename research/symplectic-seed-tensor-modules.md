# Exact small modules for the 26-label radical seed

## Result and scope

The two cuts of the 26-label symplectic radical seed admit **constant direct-sum
decompositions into small matrix polynomials** on the curve

\[
R(t,0)=E+tA+t^2B+t^3C.
\]

The row cut has ten module types, with largest dimension 22. The column cut has
nine types, with largest dimension 6. Each cut acts on a space of dimension
\(26^2=676\). These are actual direct summands with explicit invertible basis
changes, rather than composition factors or a filtration estimate.

Consequently, the complete tensor-square coefficient span can be evaluated
using 55 row type blocks of dimension at most 484 and 45 column type blocks of
dimension at most 36. The original tensor-square cuts would have dimension
456,976. Exact ranks of all 127 nonzero binary coefficient combinations take
less than a second in the reference environment.

This finite test does **not** improve the factoring score. For two copies, the
minimum binary decomposition cost lies in the certified interval

\[
250528\leq S_{\min}\leq250530,
\]

whereas the required threshold for \(d=676\) is approximately 10049.61.
The later [common-character theorem](symplectic-seed-character-barrier.md)
excludes the complete six-generator span at all powers over every
characteristic-two field. The module decomposition remains useful for
exact rank calculations without constructing the large ambient matrices.

## Seed conventions

The frozen seed is
`evidence/symplectic-radical-P4-L4-basis.json`. Its four relevant basis
operators are

\[
E=B_0,\qquad A=B_1,\qquad B=B_2,\qquad C=B_3.
\]

The complete six-operator normalization from
`research/symplectic-radical-seed-algebra.md` is

\[
T=u^4E+u^3vA+u^2v^2B+uv^3C+su^2D+suvF,
\quad
u=q+r,\ v=p+q,\ s=qr.
\]

Setting \(u=1,v=t,s=0\) gives \(R(t,0)\). Its constant coefficient is the
crossing tensor \(E\). The same convention is used in both cuts and in every
rank table here. A mask \(m\) denotes the sum of coefficient matrices whose
degree bits occur in \(m\); it does not denote evaluation at the integer
\(m\).

The scalar specialization \(t\mapsto1+t\) is convenient for inspecting
singularities, but the complete binary span is invariant under that invertible
change of coefficient coordinates. The finite searches below use the original
\(t\)-coordinate, so their target mask is 1.

## A five-dimensional column algebra

Let the operators in this section denote the column cut. Define

\[
e=A^3,\qquad f=I+e,\qquad N=A+B.
\]

The exact matrix relations are

\[
e^2=e,\quad eA=eB=e,\quad
N^2=A^2N=C^2=0,\quad fC=AN.
\]

The projector \(e\) has rank 8. On its 668-dimensional complement, \(A^3=0\)
and

\[
B=A+N,\qquad C=AN.
\]

Thus this part is a module over

\[
\mathbb F_2[a,n]/(a^3,n^2,a^2n),
\]

whose spanning monomials are \(1,a,n,a^2,an\). The polynomial action is

\[
R(t,0)=I+(t+t^2)A+t^2N+t^3AN.
\]

On the remaining eight dimensions, \(A=B=I\), while \(C^2=0\) and
\(\operatorname{rank}C=4\). This is four copies of the two-dimensional
polynomial block

\[
\begin{pmatrix}
1+t+t^2&t^3\\0&1+t+t^2
\end{pmatrix}.
\]

The checker also verifies the useful global identity

\[
R(t+1,0)R(1,0)=R(t,0),\qquad R(1,0)^2=I
\]

on the column cut. Hence multiplication by the constant invertible matrix
\(R(1,0)\) converts every centered coefficient matrix into its corresponding
uncentered coefficient matrix. No coefficient rank changes under this
normalization.

For additional independent numerical checks, the complement ranks are

\[
\operatorname{rank}(fA)=179,\quad
\operatorname{rank}(fA^2)=16,\quad
\operatorname{rank}N=321,\quad
\operatorname{rank}(AN)=22.
\]

These relations are verified directly on the original column blocks, before
using the small-module certificate.

## Constant direct-sum certificate

The explicit certificate is
`evidence/symplectic-seed-tensor-modules-certificate.json`. For each original
block it lists the small summands and the columns of their embeddings. The
concatenation of those columns has full rank. If \(P_i\) is an embedding and
\(M_i\) is the corresponding small matrix, the checker verifies

\[
MP_i=P_iM_i
\]

for each of \(E,A,B,C\). Since all embeddings together form a basis, these
identities prove the direct sum and preserve every polynomial coefficient rank.

The row source is the already verified adapted representation in
`evidence/symplectic-radical-seed-algebra.json`; the column source is the
original seed export. Both source hashes are pinned in the certificate. There
are 5,408 verified column identities, namely four operators on 676 basis
columns in each of the two cuts.

The following table groups summands by size. Equal sizes need not denote equal
module types; the certificate distinguishes all ten row types and nine column
types.

| Summand dimension | Row multiplicity | Column multiplicity |
|---:|---:|---:|
| 1 | 94 | 26 |
| 2 | 176 | 191 |
| 4 | 24 | 55 |
| 6 | 0 | 8 |
| 8 | 4 | 0 |
| 10 | 8 | 0 |
| 22 | 1 | 0 |

The decomposition was found through three elementary operations. First,
undirected support components give immediate invariant direct sums. Second,
cyclic submodules are tested for a module projection that is the identity on
the submodule. Such a projection certifies a complementary kernel. Third, an
endomorphism commuting with all generators gives a Fitting decomposition into
the stable kernel and stable image of a sufficiently high power. These spaces
are complementary invariant subspaces.

The search can use random commuting endomorphisms, but the final verification
is deterministic. It checks the delivered bases and all intertwining
identities. No assertion that the delivered summands are indecomposable is
needed.

## Exact common-parameter tensor powers

Suppose a cut has the certified decomposition

\[
R(t,0)\simeq\bigoplus_{i=1}^h m_i R_i(t).
\]

For \(n\) copies, group the tensor factors by their type counts
\(\nu_1+\cdots+\nu_h=n\). The corresponding polynomial matrix is

\[
\bigotimes_{i=1}^h R_i(t)^{\otimes\nu_i}
\]

and its multiplicity is the ordinary integer

\[
\frac{n!}{\prod_i\nu_i!}\prod_i m_i^{\nu_i}.
\]

All factors use the same parameter. Coefficients are merged inside each tensor
block before rank is computed. This differs from multiplying costs of
independent one-copy decompositions.

Every binary linear functional on the \(3n+1\) coefficient coordinates can be
evaluated by the same type calculation. The row and column type lists may
differ. One applies the identical coefficient functional to both lists and
then pairs their exact ranks.

The number of type-count vectors is polynomial in \(n\) for these fixed type
lists, but their matrix dimensions can still grow exponentially. Thus this
reduction does not itself establish a polynomial-time all-size rank algorithm
or a better asymptotic exponent. It is already a large reduction for the
finite common-parameter tests.

## Complete binary recombination

For each nonzero coefficient mask \(m\), write its exact ranks as
\(r_m,c_m\). Its cost is \(\sqrt{r_mc_m}\). A decomposition of the crossing
tensor is a path in the XOR graph of coefficient masks, starting at 0 and
ending at 1. An edge labeled \(m\) has that cost. A shortest path therefore
searches every binary recombination, including decompositions that use sums of
many coefficient matrices.

The checker runs this calculation three ways. Floating-point costs provide a
numerical minimizer. Integer floor costs provide a rigorous lower bound for
the real optimum. Integer ceiling costs provide a rigorous upper bound from a
concrete decomposition. The finite failure statements use only the integer
bounds.

| Copies | Labels \(d\) | Nonzero masks | Certified cost interval | Numerical \(\rho\) |
|---:|---:|---:|---:|---:|
| 1 | 26 | 15 | \([494,495]\) | 0.9520889330 |
| 2 | 676 | 127 | \([250528,250530]\) | 0.9538798482 |

Here \(\rho=\log S/(2\log d)\). The desired threshold is
\(\rho<1/\sqrt2\). In fact, both certified lower bounds exceed
\(d^{3/2}\), which is a stronger integer-checkable exclusion because
\(\sqrt2<3/2\).

For two copies, one numerical minimizer uses masks 24, 102, and 127. Their XOR
is 1, and their rank pairs are

\[
(456176,102453),\qquad(5002,226598),\qquad(1,456976).
\]

Its real cost is approximately 250529.14106148. The complete rank tables are
included in the evidence. An independent evaluator reproduced all 127 rank
pairs. The one-copy table is also checked against the independently computed
16-entry table from the support-component calculation.

These exclusions concern the binary coefficient span of the stated one-copy
and two-copy \(z=0\) family. They do not exclude different codes, larger tensor
powers, or extra directions involving \(D,F\). Those latter two extensions
are covered by the separate common-character theorem; different codes
still require their own analysis.

## Reproduction

Run from the repository root:

```sh
python tools/check_symplectic_seed_tensor_modules.py \
  --output /tmp/symplectic-seed-tensor-modules.json
cmp /tmp/symplectic-seed-tensor-modules.json \
  evidence/symplectic-seed-tensor-modules.json
```

The checker verifies the constant module certificate, the column relations,
all one-copy controls, both complete tensor-square rank tables, and the
integer shortest-path bounds. It does not require replaying the exploratory
module-splitting search.
