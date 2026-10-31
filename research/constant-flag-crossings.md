# Unimodular crossings from two constant flags

The one-species construction extends beyond the Boolean Pascal matrices.
Two constant invertible matrices, a list of exponents, and any unitriangular
endpoint matrix give a polynomial Gaussian crossing with determinant one.
This provides a way to change the internal matrix geometry while preserving
the endpoint and the determinant cancellation.

The first exact finite screen is negative: a family of 4096 four-mode flags
over \(\mathbb F_4\), including all 64 binary flags, cannot improve the
displayed Pascal gauge. The statement concerns this finite flag family and
this gauge. It does not exclude different gauges or larger constructions.

## 1. The constant-flag constructor

Work over any field of characteristic two. Let \(M\) be a power of two,
choose \(n\) integers

\[
0<d_i<2M,
\]

and let \(U,V\) be invertible constant \(n\times n\) matrices. Let \(P\)
be a constant unitriangular matrix. Put \(c=1+w\) and define

\[
\begin{aligned}
A(w)&=U\mathrm{diag}(c^{d_i})V,\\
B(w)&=V^{-1}\mathrm{diag}(c^{2M-d_i})U^{-1}P,\\
\mathcal D(w)&=
\begin{pmatrix}
A(w)&w^{2M}P\\
I&B(w)
\end{pmatrix}.
\end{aligned} \tag{1}
\]

All entries are polynomial. Direct multiplication gives

\[
A(w)B(w)=c^{2M}P.
\]

Taking the block determinant using the identity block,

\[
\det\mathcal D(w)
=\det\bigl(w^{2M}P+A(w)B(w)\bigr)
=\det P=1, \tag{2}
\]

because \(c^{2M}=1+w^{2M}\) in characteristic two.

The strict bounds on the \(d_i\) ensure \(A(1)=B(1)=0\), so

\[
\mathcal D(1)=
\begin{pmatrix}0&P\\I&0\end{pmatrix}. \tag{3}
\]

Interpret its rows as the plus modes on boundaries \(V,Y\), and its
columns as the minus modes on \(B,X\). Its minors define the one-species
Gaussian tensor.

### Why a triangular endpoint is sufficient

Take tensor powers and prescribe a common occupation count for every mode
on all four physical boundaries. The identity block in (3) forces the
corresponding two subset labels to agree. Every nonzero term in a minor
of a lower-unitriangular \(P\) matches each column index to a row index at
least as large. Fixed mode counts equalize the total index sums across
the tensor factors. Hence every such matching must use equality, and the
remaining endpoint factor is also the identity.

For an upper-unitriangular \(P\), reverse the inequalities. The diagonal
entries are one, so the surviving endpoint coefficient is one. Thus the
fixed-mode restriction is an exact crossing. Binary digit inclusion is
not required for this argument.

## 2. A family preserving the useful zero gauge

Now set \(n=M=4\), let \(P\) be the Boolean Pascal matrix, and choose

\[
d_i=2i+1,\quad U=I,\quad
V=P^{\mathsf T}C^{-1},
\]

where \(C\) is lower unitriangular. The two polynomial factors are

\[
A=\mathrm{diag}(c^{2i+1})P^{\mathsf T}C^{-1},
\qquad
B=C P^{\mathsf T}\mathrm{diag}(c^{2M-2i-1})P.
\tag{4}
\]

At \(C=I\), this is the original full-interval one-species kernel from
[the odd-spectrum construction](schur-one-species-interval.md).

Keep that construction's physical heights, in boundary order
\((B_-,X_-,V_+,Y_+)\):

\[
h_0(i)=(0,-2M+1+2i,0,1+2i),
\]

\[
h_\infty(i)=(2i-n+1,2M-n,2i-n+1,-n).
\tag{5}
\]

The zero row cut has cross blocks \(w^{2M}P\) and \(I\), independent
of \(C\). Each therefore retains orders \(1,3,5,7\).

The internal block \(A\) is regular because its two zero-height
contributions vanish. At \(C=I\), the other internal block has entry
orders at least

\[
2M-2-2i-2j.
\]

Multiplying it by lower-unitriangular \(C\) mixes row \(i\) only with rows
\(k\le i\), whose required divisibility is at least as strong. Hence the
same zero gauge remains regular for every \(C\) in this family.

The physical mean weights stay fixed:

\[
G=(-12,-4,4,12).
\tag{6}
\]

Only the infinity Smith profile can change.

## 3. Exact finite result over two fields

The checker enumerates all \(4^6=4096\) lower-unitriangular \(4\times4\)
matrices over \(\mathbb F_4=\mathbb F_2[\omega]/(\omega^2+\omega+1)\).
It includes all \(2^6=64\) binary matrices.

There are exactly nine distinct infinity profiles. The baseline is

\[
\lambda=(-7,-5,-3,-1,1,3,5,7).
\]

For every enumerated profile \(\kappa\),

\[
\sum_{i=1}^{j}\kappa_i
\le \sum_{i=1}^{j}\lambda_i
\quad(1\le j<8),
\qquad
\sum_i\kappa_i=\sum_i\lambda_i=0. \tag{7}
\]

Equality of the full profile occurs only at \(C=I\). Thus every other
profile is more spread out than the baseline. Since
\(x\mapsto\log(1+s^{-x})\) is strictly convex for \(s>1\),

\[
\sum_i\log(1+s^{-\kappa_i})
\ge
\sum_i\log(1+s^{-\lambda_i}). \tag{8}
\]

The zero profiles and physical weights are unchanged, so this finite flag
family cannot lower the pressure at the gauge (5), for any \(s>1\) or any
independent mode probabilities.

The complete profile multiplicities are recorded in
[constant-flag-crossings.json](../evidence/constant-flag-crossings.json).
This exact result does not extend the exclusion to arbitrary constant
flags, other physical gauges, or higher dimensions.

## 4. Reproduction and arithmetic controls

Run:

~~~sh
python tools/check_constant_flag_crossings.py \
  --output evidence/constant-flag-crossings.json
~~~

The checker verifies the polynomial identity \(AB=c^{2M}P\), the zero
regularity conditions, determinant one, and every majorization comparison.
Its additional controls include:

- 65,536 comparisons of packed \(\mathbb F_4\) multiplication with an
  independent coefficient-convolution implementation, together with exact
  division controls;
- independent enumeration of all nonzero minors of sizes one through
  three for four selected \(8\times8\) weighted matrices;
- agreement with the binary Smith oracle on every one of the 64 binary
  flags.

The supporting
[four-element-field oracle](../tools/weighted_smith_oracle_f4.py) uses two
bits per polynomial coefficient. It applies the same fraction-free
minimum-valuation argument as the
[binary oracle](weighted-smith-separation-oracle.md); field multiplication
and exact polynomial division are the only arithmetic changes.
