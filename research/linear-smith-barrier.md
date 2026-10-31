# A uniform obstruction for degree-one Smith certificates

## Result and scope

Every degree-one polynomial crossing with alphabet size \(q\ge2\), over
any field and with any fixed auxiliary dimension, has the bare Smith
profile bound

\[
B(s)\ge q^3\qquad(s\ge1).                                  \tag{1}
\]

Consequently this certificate has \(\rho\ge3/4>1/\sqrt2\). Increasing
the alphabet or auxiliary dimension cannot make a degree-one bare
profile pass the strict half-rate threshold. The proof combines the
two endpoint ranks with the integral nature of Smith orders.

For a binary alphabet the stronger universal bound is \(B\ge256/27\).
It is sharp in characteristic two. For **scalar** binary crossings in
odd characteristic the sharp bound is \(B\ge12\). An exact enumeration
of all 65,536 binary scalar pencils finds 60 paired Smith profiles and
exactly six optimal pencils, all changes of basis of the elementary one.

These are restrictions on the specified certificate, not lower bounds
for factoring or for actual crossing ranks. Higher degree, rational
auxiliary gauges, [fixed blocks](block-smith-crossings.md), and
[cancellation modules](cancellation-module-smith-automata.md) are outside
the asserted obstruction. A fixed block of a linear local crossing has
higher polynomial degree and is not covered by (1).

The polynomial crossing interface and the earlier single-cut degree
budget are from [PR #30 at `28105165`](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md).
The coupled endpoint argument and exact catalogue below were developed
with substantial OpenAI assistant support. No new Lean formalization or
improved factoring exponent is claimed.

## 1. Polynomial and rank conventions

Let \(F\) be any field, and let

\[
K(z;x,y,u,v)=A(x,y,u,v)+zB(x,y,u,v)
\]

have \(D\)-by-\(D\) matrix entries. The physical labels range over
\(q\) values. Suppose fixed vectors satisfy

\[
f^Te=1,\qquad (A+B)e=[x=v][y=u]e.                            \tag{2}
\]

Use the cut \((x,y)\mid(u,v)\) at infinity and the cut
\((x,u)\mid(y,v)\) at zero. In each cut put both auxiliary indices on
one physical side; either side is allowed. Transpose the flattening if
necessary so it has \(n=q^2\) columns.

The regularized matrix at infinity is \(B+\tau A\), with Smith orders
\(\mu_1,\ldots,\mu_n\). The zero matrix is \(A+\tau B\), in the other
cut, with orders \(\nu_1,\ldots,\nu_n\). Contracting by \(f,e\) at
\(\tau=1\) gives an invertible crossing matrix. Each rectangular matrix
therefore has a nonzero \(n\)-by-\(n\) minor of degree at most \(n\).
Its determinant bounds the valuation of the determinantal divisor, so

\[
\mu_i,\nu_i\in\mathbb Z_{\ge0},\qquad
\sum_i\mu_i\le n,\quad\sum_i\nu_i\le n.                    \tag{3}
\]

Let \(r\) and \(c\) be the numbers of zero orders in the two lists.
Equivalently, \(r\) is the first-cut rank of the auxiliary flattening
of \(B\), and \(c\) is the second-cut rank of that of \(A\).

For \(s\ge1\), put \(x=s^{-1}\) and

\[
P_r(x)=\sum_i x^{\mu_i},\qquad
P_c(x)=\sum_i x^{\nu_i},\qquad
B(s)=\frac{P_r(x)P_c(x)}x.                                  \tag{4}
\]

The letter \(x\) in (4) is a positive real variable, independent of the
physical label used in (2). Formula (4) is the degree-one specialization
of the existing bare Smith transfer.

## 2. The two endpoint ranks satisfy \(r+c\ge q\)

First consider a scalar decomposition of the crossing into two tensors:

\[
E=U+V,\qquad
\mathrm{rank}_{xy\mid uv}U=r',\quad
\mathrm{rank}_{xu\mid yv}V=c'.
\]

Write each tensor as the corresponding sum of rank-one matrices. Fix
the physical labels \(x=v=a\), for any one label \(a\). The remaining
\(y\)-by-\(u\) slice of \(E\) is the identity matrix. Each of the
\(r'+c'\) summands restricts to an outer product in \(y,u\). Hence

\[
q\le r'+c'.                                                 \tag{5}
\]

In the matrix-valued case, contract (2) by \(f^T\), obtaining
\(E=f^TBe+f^TAe\). Because both auxiliary indices are on one side of
each flattening, this contraction is a linear map on that side and
cannot increase its matrix rank. Applying (5) proves

\[
r+c\ge q.                                                   \tag{6}
\]

No positivity, symmetry, field-size restriction, or full matrix
endpoint \(K(1)=EI_D\) is needed.

## 3. Integral orders give the uniform bound

For \(0<x\le1\) and every integer \(k\ge1\), discrete convexity gives

\[
x^k\ge x+(k-1)(x^2-x).                                     \tag{7}
\]

For example, for \(k\ge2\) this follows by summing the increasing
successive differences of \(1,x,x^2,\ldots\). The case \(k=1\) is
equality. Sum (7) over the \(n-r\) positive orders and use (3). The
coefficient \(x^2-x\) is nonpositive, so the upper bound on their sum
has the required direction. Including the \(r\) zero orders gives

\[
P_r(x)\ge nx+r(1-x)^2,
\qquad P_c(x)\ge nx+c(1-x)^2.                               \tag{8}
\]

Both right sides are nonnegative. Multiplication and (6) imply

\[
\begin{aligned}
B(s)
&\ge n^2x+n(r+c)(1-x)^2+\frac{rc(1-x)^4}{x}\\
&\ge nq\bigl(qx+(1-x)^2\bigr)\\
&=q^3\bigl(1+(q-2)x+x^2\bigr)\ge q^3.
\end{aligned}                                               \tag{9}
\]

This proves (1). The consequence \(\rho\ge3/4\) uses
\(\rho=\log B/(4\log q)\) from the existing transfer; it is a lower
bound on this sufficient width parameter only.

## 4. The sharp binary bound with arbitrary auxiliaries

Now \(n=4\). If both \(A\) and \(B\) are nonzero, then \(r,c\ge1\),
so (8) yields

\[
P_r(x),P_c(x)\ge(1+x)^2,
\qquad B(s)\ge\frac{(1+x)^4}{x}\ge\frac{256}{27}.             \tag{10}
\]

The last function has its minimum at \(x=1/3\). If \(A=0\) or
\(B=0\), the nonzero coefficient has a full-rank scalar contraction
by (2); the two order lists instead give \(B(s)=16\).

In characteristic two, the scalar pencil

\[
K(z)=I+z(E+I),\qquad I=[x=u][y=v],
\]

has both order lists \((0,1,1,2)\). Thus it attains (10). The bound is
sharp even without auxiliary matrices.

## 5. Scalar binary pencils in odd characteristic

The following small rank fact identifies the exceptional case.

**Lemma.** For any two \(2\)-by-\(2\) matrices \(X,Y\), the matrix
\(E-X\otimes Y\) never has rank two. It has rank one exactly when
\(X,Y\) are invertible and \(YX=I\).

**Proof.** Independent changes of basis preserving the swap reduce a
singular \(X\) to zero or \(\mathrm{diag}(1,0)\). In the latter
case three rows have unavoidable pivots in distinct columns, so the
rank is at least three. If \(X\) is invertible, left multiplication by
\(X^{-1}\otimes I\) and right multiplication by \(I\otimes X\) give
\(E-I\otimes W\), where \(W=YX\). Extend the field to an algebraic
closure, which does not change ranks, and triangularize \(W\).
For \(W=\left(\begin{smallmatrix}a&b\\0&d\end{smallmatrix}\right)\),
the matrix is

\[
\begin{pmatrix}
1-a&-b&0&0\\
0&-d&1&0\\
0&1&-a&-b\\
0&0&0&1-d
\end{pmatrix}.
\]

If either diagonal entry of \(W\) differs from one, direct elimination
gives rank at least three. If \(a=d=1\), the rank is one when \(b=0\)
and three otherwise. This also proves the rank-one characterization.
\(\square\)

For a scalar pencil, \(r=1\) implies that \(B\) is a rank-one matrix
in the first cut. Reshuffling it into the second cut gives \(X\otimes Y\).
Since \(A=E-B\), the lemma shows that \(c\ne2\). The same argument
with the cuts exchanged rules out \((r,c)=(2,1)\).

If \(r=c=1\), the lemma gives a one-leg change of basis reducing the
pencil to \(I+z(E-I)\). Such changes preserve both Smith profiles.
In odd characteristic the two order lists are \((0,1,1,1)\), giving

\[
B(s)=\frac{(1+3x)^2}{x}\ge12,
\qquad x=1/3\text{ at equality}.                            \tag{11}
\]

If \(r,c\ge2\), (8) gives

\[
B(s)\ge\frac{4(1+x^2)^2}{x}
\ge\frac{64\sqrt3}{9}>12.                                  \tag{12}
\]

In the remaining nonzero cases one count is one and the other is at
least three. Its order polynomial is at least \(3+x^4\), since its
only possible positive order is at most four. Consequently

\[
B(s)\ge\frac{(1+x)^2(3+x^4)}x\ge12,                         \tag{13}
\]

where the last numerator minus \(12x\) is
\(3(1-x)^2+x^4(1+x)^2\). The zero-coefficient cases again give 16.
This proves the claimed sharp odd-characteristic scalar bound.

## 6. Exact binary scalar catalogue

Over \(\mathbb F_2\), every scalar degree-one pencil with endpoint
\(E\) is uniquely \(A+z(E+A)\), with 16 binary entries of \(A\).
The enumerator visits all \(2^{16}=65,536\) masks. It computes every
determinantal divisor from the exact polynomial determinants of all
minors, then takes consecutive differences of their valuations.

There are 60 distinct paired order lists. Exactly six pencils have both
lists \((0,1,1,2)\), hence attain the universal binary minimum. They are
precisely

\[
A(x,y,u,v)=S_{x,u}(S^{-1})_{y,v},
\qquad S\in\mathrm{GL}_2(\mathbb F_2).                       \tag{14}
\]

For a proof independent of enumeration, equality in (10) forces
\(r=c=1\). Write the rank-one second-cut coefficient as \(X\otimes Y\).
The lemma gives \(YX=I\), and the six invertible binary matrices give
(14). Conversely each is a one-leg basis change of the elementary
pencil. Thus the enumeration does not supply an unsupported numerical
optimization claim.

Run:

```sh
python3 tools/check_linear_smith_barrier.py \
  --output evidence/linear-smith-barrier.json
```

The driver compiles the C++17 enumerator, independently reconstructs
the reported profiles from ranks of truncated block-Toeplitz matrices,
checks the six equality cases, checks the rank lemma over \(\mathbb F_2\)
and \(\mathbb F_3\), and verifies (8) on all order partitions of total
size at most \(q^2\), for \(q=2,3,4\), at exact rational points.
The all-field and all-alphabet claims have the proofs above.

Files: [exact evidence](../evidence/linear-smith-barrier.json),
[verification driver](../tools/check_linear_smith_barrier.py), and
[complete enumerator](../tools/enumerate_binary_linear_crossings.cpp).
