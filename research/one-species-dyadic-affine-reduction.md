# A dyadic affine reduction for one-species mode sets

## Result and scope

Appending the same string of binary ones to every mode does not give a
new one-species crossing. The complete typed polynomial tensor is exactly
the smaller-mode tensor after a change of variable. Its nonzero coefficient
tensors, including both flattening ranks, are unchanged.

More precisely, let \(L\) be a power of two, let \(N\) be a power of two,
and let \(T\subseteq\{0,\ldots,N-1\}\) be nonempty. Set

~~~math
S=\{La+L-1:a\in T\},\qquad M=LN.                       \tag{1}
~~~

For every tensor power and every common product type on the four physical
boundaries, the one-species tensors satisfy

~~~math
\boxed{\quad \mathcal T_{S,M,\mathbf k}(w)
       =\mathcal T_{T,N,\mathbf k}(w^L).\quad}          \tag{2}
~~~

The label sets are identified by \(La+L-1\leftrightarrow a\). Thus an
improvement obtained from the left side is already available from the
right side. This is an exact equivalence, not an asymptotic entropy
comparison or a numerical optimization result.

The restriction to one species matters. The full two-species kernel has
an additional scalar defect described in Section 5. This note does not
exclude arbitrary mode sets, different coefficient decompositions, or a
new factoring exponent. The derivation and checks used substantial
OpenAI assistant support.

## 1. The polynomial kernel

Use the one-species construction in
[the interval note](schur-one-species-interval.md), with arbitrary finite
mode sets as justified by [fixed-mode crossings](fixed-mode-crossings.md).
Write

~~~math
z=w^2,\qquad q=1+w^2,\qquad c=1+w,
\qquad P_{ab}=\mathbf 1_{b\subseteq a},
~~~

where \(b\subseteq a\) means inclusion of binary digits. The rows are
\((V,Y)\), and the columns are \((B,X)\). The kernel is

~~~math
\mathcal D_{S,M}(w)=
\begin{pmatrix}
c\,\mathrm{diag}(q^a)P^{\mathsf T}&z^MP\\
I&c\,(z^{M-1-(a\mathbin{|}b)})_{a,b\in S}
\end{pmatrix}.                                        \tag{3}
~~~

A local tensor entry is the selected minor of (3), or zero when the
number of selected rows differs from the number of selected columns.
The occupied modes of physical labels \((x,y,u,v)\) select rows
\((v,y)\) and columns \((u,x)\), respectively.

## 2. A rational boundary identity

Put \(\lambda=(1+w)^{L-1}\). Under the index identification in (1),

~~~math
\boxed{\quad
\mathcal D_{S,M}(w)=
\begin{pmatrix}\lambda I&0\\0&I\end{pmatrix}
\mathcal D_{T,N}(w^L)
\begin{pmatrix}I&0\\0&\lambda^{-1}I\end{pmatrix}.
\quad}                                                \tag{4}
~~~

This is an identity over \(\mathbb F_2(w)\). Both kernels themselves
have polynomial entries.

**Proof.** Binary digit inclusion is preserved by the common affine
map, so the two Pascal matrices agree. Also

~~~math
(La+L-1)\mathbin{|}(Lb+L-1)=L(a\mathbin{|}b)+L-1.
~~~

Because \(L\) is a power of two, \(1+w^L=(1+w)^L\) and
\(1+w^{2L}=(1+w^2)^L\). The upper-left block of (3) therefore differs
from the smaller kernel at \(w^L\) by the factor

~~~math
\frac{(1+w)(1+w^2)^{L-1}}{1+w^L}
=(1+w)^{L-1}=\lambda.
~~~

The lower-right block differs by \((1+w)/(1+w^L)=\lambda^{-1}\).
The two off-diagonal blocks agree. These are exactly the four block
identities in (4).

In particular, if a selected minor uses \(r_V\) rows from \(V\) and
\(r_X\) columns from \(X\), then

~~~math
\det\mathcal D_{S,M}[I,J]
=\lambda^{r_V-r_X}
 \det\mathcal D_{T,N}(w^L)[I,J].                       \tag{5}
~~~

One can verify (5) entirely in the polynomial ring by multiplying its
two sides by \(\lambda^{r_X}\) and \(\lambda^{r_V}\), respectively.

## 3. Common type cancels the boundary factors exactly

Take \(m\) tensor copies. For each mode \(a\in T\), select words with
exactly \(k_a\) occupied positions. Use this same word set at all four
physical boundaries. Its size is

~~~math
D_m=\prod_{a\in T}\binom m{k_a}.                       \tag{6}
~~~

The product of the factors in (5) is

~~~math
\lambda^{\sum_{j=1}^m(r_{V,j}-r_{X,j})}
=\lambda^{\sum_a k_a-\sum_a k_a}=1.                   \tag{7}
~~~

Multiplying (5) across the copies proves (2) over the rational-function
field. The inclusion \(\mathbb F_2[w]\hookrightarrow\mathbb F_2(w)\)
is injective, so (2) is a polynomial identity. It can therefore be
evaluated at \(w=1\), even though the separate boundary factors in (4)
are singular there. No singular boundary map is applied at that point.

The proof needs only equal total occupation counts on \(V\) and \(X\)
over the copies. The common product type supplies that equality. The
usual fixed-mode endpoint theorem supplies the crossing endpoint for
these labels. Thus neither an approximation by a small second-species
density nor an extra projection is used.

## 4. Coefficients, ranks, and search normalization

Write the smaller typed tensor as

~~~math
\mathcal T_{T,N,\mathbf k}(w)=\sum_{j=0}^{d} A_jw^j.
~~~

Equation (2) states exactly that

~~~math
[w^r]\mathcal T_{S,M,\mathbf k}(w)
=\begin{cases}A_{r/L},&L\mid r,\\0,&L\nmid r.\end{cases} \tag{8}
~~~

Consequently the multiset of nonzero coefficient tensors is the same.
Both flattening ranks of every corresponding coefficient are equal,
as is any coefficient cost depending only on these tensors and ranks.
An effective decomposition or evaluation query can be transported by
reindexing the exponent. Zero coefficients need not be enumerated.

This gives a canonical reduction for mode searches. Let

~~~math
e=\min_{a\in S}v_2(a+1),\qquad L=2^e,
\qquad T=\{(a+1)/L-1:a\in S\}.                        \tag{9}
~~~

When \(e>0\), remove the common trailing ones by (9), and divide the
ambient power of two by \(L\). At least one reduced mode is even, so
the reduction is maximal. An initial ambient power of two containing
\(S\) is necessarily divisible by \(L\).

Examples include every singleton \(S=\{2^e-1\}\), which reduces to
\(\{0\}\), and every odd progression
\(S=\{1,3,\ldots,2n-1\}\), which reduces to
\(\{0,1,\ldots,n-1\}\). Arbitrary noninterval sets can be reduced as
well: \(\{3,11,23\}\) reduces to \(\{0,2,5\}\).

The common-height pressure formula in
[the Pascal complement note](schur-pascal-complement.md) is one specific
gauge choice. It need not assign the same numerical value to two
equivalent typed tensors. Such a discrepancy measures the chosen
certificate, not a difference between the polynomial tensors in (2).

There is also an explicit correspondence for independent monomial gauges.
List the smaller one-particle heights in order \((V,Y;B,X)\), and let
\(h_0,h_\infty\) denote their zero and infinity values. Set

~~~math
\widetilde h_0=Lh_0,\qquad
\widetilde h_\infty=Lh_\infty+(L-1)(1,0;0,-1),          \tag{10}
~~~

with the displayed corrections on every mode. Then every cross-block
Smith order is multiplied by \(L\), and internal regularity is
preserved. Indeed, \(\lambda\) is a unit at zero and has infinity
valuation \(-(L-1)\); after the corrections in (10), the remaining
boundary factors are units at both places. Each restored physical mode
height is also multiplied by \(L\), since its two corrections cancel.

Taking the pressure parameter to be \(\widetilde s=s^{1/L}\) gives
exactly the original pressure and occupation probabilities. Conversely,
the inverse of (10) gives rational heights, which are allowed after a
common ramification. Thus the family of independent monomial-gauge
pressure certificates is the same on both sides of (2). A numerical
search that returns different optima has not found an intrinsic
distinction between these families.

## 5. The remaining defect with two occupied species

For completeness, let \(K_{S,M}\) denote the full polynomial kernel in
boundary order \((B,X,V,Y)\), as defined in
[the down-set checker](../tools/check_schur_downset_crossing.py). Its
formula applies to arbitrary mode sets. Define

~~~math
R=\mathrm{diag}(I,\lambda I,\lambda I,I),\qquad
C=\mathrm{diag}(I,\lambda^{-1}I,\lambda^{-1}I,I),
\qquad
\zeta=\left(\frac{w}{1+w}\right)^{2(L-1)}.             \tag{11}
~~~

Then \(K_{S,M}=R\widetilde K_{T,N}(w^L)C\), where
\(\widetilde K\) is the smaller kernel except that its two \(G\)
blocks, from row boundary \(X\) to column boundaries \(B\) and \(Y\),
are multiplied by \(\zeta\). This follows by the same digit
calculation; \(G_{S,ab}=z^{L-1}G_{T,ab}(z^L)\), whereas the
\(H\) block has no extra monomial factor.

The row boundary \(X\) is absent from the one-species selected minors,
so these two blocks do not occur in (3). With both species occupied they
can occur, and their factors do not reduce to a common boundary scalar.
Equation (2) is therefore not asserted for that larger family.

## 6. Exact reproduction

Run

```sh
python3 tools/check_one_species_dyadic_affine.py \
  --output evidence/one-species-dyadic-affine.json
```

The checker verifies the polynomial block identities, all minors for
small gapped mode sets, complete typed tensor entries, and the reindexing
of coefficient flattening ranks. It also checks examples where the
occupation counts vary between copies, and a negative control showing
that an unbalanced minor does retain a boundary factor. The general
proof is (4)--(8); the finite checks do not replace it.
