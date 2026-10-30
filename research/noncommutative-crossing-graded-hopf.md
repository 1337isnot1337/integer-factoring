# Graded Hopf crossing families and scalar rank factorizations

## Result and scope

A graded, noncommutative Hopf algebra gives an exact characteristic-two
crossing decomposition with a linear number of coefficients. This note proves
the construction for every size, derives its generator substitution, and gives
quantitative bounds and explicit factorizations for both scalar cuts.

The main reusable results are:

- A change of coordinates that makes every row coefficient a map of a fixed
  bidegree, giving an explicit Hilbert-function rank bound.
- A factorization into an abelian mixing operator and a central translation.
- An exact compression of every translation coefficient's column rank from
  a matrix of order \(d^2\) to a matrix of order \(d\).
- An independent-set bound for those compressed ranks and a convolution
  bound for every complete column coefficient.
- Parameter-dependent local changes that preserve the crossing endpoint.

The finite examples improve some tested noncommutative profiles. A subsequent
[PBW degeneration theorem](noncommutative-crossing-pbw-degeneration.md) proves
that every complete decomposition on the full PBW label domain has rank sum
at least \(d^{3/2}\), including scalar recombinations and the stated gauges.
These full-domain families cannot meet the improved finite-seed criterion.
The explicit constructions and factorizations below remain valid; no
stronger factoring score is claimed.

## 1. The graded algebra

Let \(L\) be an \(r\)-by-\(m\) binary matrix with columns \(v_i\). Work over
any field of characteristic two. Define \(H_L\) by generators

\[
x_1,\ldots,x_m,\quad y_1,\ldots,y_m,\quad z_1,\ldots,z_r
\]

with square zero, central \(z_j\), and relations

\[
y_i x_i+x_i y_i=\sum_j L_{ji}z_j,\qquad
[x_i,y_k]=0\quad(i\ne k).
\tag{1}
\]

All other distinct generator pairs commute. Give \(x_i,y_i\) degree one
and \(z_j\) degree two. Put

\[
\lambda=2m+r,\qquad d=2^\lambda,\qquad D=2m+2r.
\tag{2}
\]

The ordered square-free monomials \(x^a y^b z^c\) form a basis. One
elementary construction proves this without an assumed normal-form theorem.
Start with the commutative square-zero algebra on the \(x_i,z_j\), and
adjoin \(y_i\) using the derivation
\(\delta_i=(\sum_jL_{ji}z_j)\partial_{x_i}\). The derivations commute,
their squares vanish, and the Ore-extension relations permit \(y_i^2=0\).
This gives an associative algebra with exactly the stated basis.

Declare every generator primitive. The square and commutator relations
generate a Hopf ideal in characteristic two. Thus

\[
\Delta(g)=g\otimes1+1\otimes g,\quad
\epsilon(g)=0,\quad S(g)=g
\tag{3}
\]

on generators, with \(S\) reversing products. The coproduct of a PBW
monomial is the sum over all partitions of its set of generators.
In particular, the coalgebra is the ordinary square-free subset coalgebra.
The grading gives a Hopf endomorphism
\(\phi_\tau(h)=\tau^{\deg h}h\) on homogeneous elements.

## 2. An exact crossing endpoint

Define \(R_\tau:H_L\otimes H_L\longrightarrow H_L\otimes H_L\) by

\[
R_\tau(u\otimes v)=
\sum v_{(1)}a_{(1)}\otimes u_{(1)}a_{(2)},\qquad
a=\phi_\tau\!\left(u_{(2)}S(v_{(2)})\right).
\tag{4}
\]

The coproducts of \(u,v,a\) are expanded linearly. At \(\tau=0\),
\(\phi_0\) is the unit followed by the counit. The two counit identities give

\[
R_0(u\otimes v)=v\otimes u.
\tag{5}
\]

Write

\[
R_{1+z}=\sum_{i=0}^{D}z^i A_i.
\tag{6}
\]

The degree bound \(D\) follows directly from (4): a nonzero homogeneous
product \(u_{(2)}S(v_{(2)})\) has degree at most the largest PBW degree
\(D\). At \(z=1\), characteristic two turns (6) into \(R_0\). Hence

\[
\sum_i A_i=E_d,\qquad E_d(x,y,u,v)=[x=v][y=u].
\tag{7}
\]

There are at most \(D+1=O(\lambda)\) pieces. All cuts are ordinary scalar
matrices. PBW expansion and exact elimination construct every entry and
rank factorization in \(2^{O(\lambda)}\) time, provided the sequence of
code matrices \(L\) is uniformly constructible within that bound.
This does not hide a search for an effective code.

Every map used in (4) is a coalgebra map; cocommutativity permits the
coproduct and antipode maps in this assertion. Thus \(R_\tau\) is a
coalgebra map. Its dual is an algebra homomorphism of the commutative
square-zero algebra on the dual variables in two copies of \(H_L\).

## 3. The complete dual generator substitution

Use subscripts \(L,R\) for the two tensor factors. On horizontal variables,

\[
\begin{aligned}
x_{i,L}&\longmapsto \tau x_{i,L}+(1+\tau)x_{i,R},\\
x_{i,R}&\longmapsto (1+\tau)x_{i,L}+\tau x_{i,R},
\end{aligned}
\tag{8}
\]

and the same formulas hold for \(y_i\). On central variables,

\[
\begin{aligned}
z_{j,L}\longmapsto {}&
\tau^2z_{j,L}+(1+\tau^2)z_{j,R}\\
&+\sum_iL_{ji}\bigl(
\tau x_{i,L}y_{i,R}+\tau^2y_{i,L}x_{i,R}
+(\tau+\tau^2)x_{i,R}y_{i,R}\bigr),\\
z_{j,R}\longmapsto {}&
(1+\tau^2)z_{j,L}+\tau^2z_{j,R}\\
&+\sum_iL_{ji}\bigl(
\tau x_{i,L}y_{i,L}+(\tau+\tau^2)y_{i,L}x_{i,R}
+\tau^2x_{i,R}y_{i,R}\bigr).
\end{aligned}
\tag{9}
\]

These identities follow from (4) on inputs of weighted degree at most two.
They determine the full map: the dual algebra is generated by these
variables, and a central coordinate of weight two cannot have terms of
higher weighted degree. Inputs from distinct horizontal pairs commute.

The horizontal linear transformation has determinant one. The central
linear transformation also has determinant one; its added terms use only
horizontal variables. Thus \(R_\tau^*\) is an invertible triangular
substitution over the polynomial parameter ring.

The construction is natural under graded Hopf algebra isomorphisms.
In particular, invertible changes of central basis and permutations of
the horizontal pairs preserve every complete coefficient rank in both
cuts. The isomorphism commutes with \(\phi_\tau\), so (4) conjugates by
one linear map on each tensor leg. Here the central basis change is
compatible with the entire dilation.

## 4. Coordinates that expose every row coefficient

Set

\[
\begin{gathered}
p_{x,i}=x_{i,L}+x_{i,R},\quad q_{x,i}=x_{i,R},\\
p_{y,i}=y_{i,L}+y_{i,R},\quad q_{y,i}=y_{i,R},\\
P_j=z_{j,L}+z_{j,R},\quad Q_j=z_{j,R},\\
\mathcal A_j=\sum_iL_{ji}p_{x,i}p_{y,i},\quad
\mathcal B_j=\sum_iL_{ji}p_{x,i}q_{y,i},\quad
\mathcal C_j=\sum_iL_{ji}q_{x,i}p_{y,i}.
\end{gathered}
\tag{10}
\]

Then (8)--(9) become

\[
\begin{aligned}
p&\longmapsto p,&q&\longmapsto q+(1+\tau)p,\\
P&\longmapsto P+\tau\mathcal A,&
Q&\longmapsto Q+(1+\tau^2)P+
\tau\mathcal A+\tau\mathcal B+\tau^2\mathcal C .
\end{aligned}
\tag{11}
\]

In particular \(R_1^*\) is an involution. Multiplying the family by this
fixed invertible map preserves every row-coefficient rank. For
\(F_z=R_{1+z}^*\circ R_1^*\), the formulas are

\[
p\mapsto p,\quad q\mapsto q+zp,\quad
P\mapsto P+z\mathcal A,\quad
Q\mapsto Q+z(\mathcal A+\mathcal B)+z^2(P+\mathcal C).
\tag{12}
\]

Make the invertible polynomial coordinate change

\[
Z=P+\mathcal C,\qquad W=Q+P.
\tag{13}
\]

Its inverse is \(P=Z+\mathcal C,\ Q=W+Z+\mathcal C\). Since
\(\mathcal C\mapsto\mathcal C+z\mathcal A\), (12) simplifies to

\[
\boxed{\quad
p\mapsto p,\quad Z\mapsto Z,\quad
q\mapsto q+zp,\quad W_j\mapsto W_j+z\mathcal B_j+z^2Z_j.
\quad}
\tag{14}
\]

This is a global basis change for the row cut. It is not a product of
four separate leg maps and cannot be used to preserve column ranks.

Give invariant variables \(p\) weight one and \(Z\) weight two.
Give active variables \(q\) weight one and \(W\) weight two.
The coefficient of \(z^i\) raises invariant weighted degree by exactly
\(i\), and lowers active weighted degree by exactly \(i\).
Both graded spaces have Hilbert coefficients

\[
h_a=[x^a]\,(1+x)^{2m}(1+x^2)^r.
\tag{15}
\]

Distinct bidegree blocks have disjoint input and output spaces. Hence

\[
\boxed{\quad
r_i\le B_i:=
\sum_{a=0}^{D-i}\sum_{b=i}^{D}
\min\{h_a h_b,\ h_{a+i}h_{b-i}\}.
\quad}
\tag{16}
\]

The bound is independent of the code. It is exact at the ends:
\(r_0=d^2\) and \(r_D=1\).

There is a real obstruction to treating (12) as an additive-group action.
Direct composition gives

\[
F_sF_t=F_{s+t}C_{st},\qquad C_v:\ Q\mapsto Q+v\mathcal A,
\tag{17}
\]

with other coordinates fixed by \(C_v\). Thus \(F_t^2=C_{t^2}\).
Equivalently, for the one-pair algebra, if
\(F=\operatorname{Swap}\circ R_1\), then

\[
F^2(x\otimes y)=x\otimes y+z\otimes1+1\otimes z.
\tag{18}
\]

An invariant coordinate such as \(P+\mathcal C\) does not remove this carry.

## 5. The last coefficient has full column rank

Let \(\Omega\) be the top PBW monomial, and let \(\ell\) select its
coefficient. Multiplication by every positive-degree generator kills
\(\Omega\) on either side. Thus \(\Omega\) is a two-sided integral.
The subset coproduct gives

\[
(\operatorname{id}\otimes\ell)\Delta(a)=\ell(a)1.
\tag{19}
\]

The primitive relations also give

\[
(v\otimes1)\Delta\Omega=(1\otimes S(v))\Delta\Omega.
\tag{20}
\]

Taking the top parameter coefficient in (4), then applying (19)--(20),
yields

\[
A_D(u\otimes v)=\ell(uS(v))\,\Delta\Omega.
\tag{21}
\]

The coefficient matrix of \(\Delta\Omega\) is the complement-permutation
matrix of the PBW subsets, so it has rank \(d\). The form \(\ell(uv)\)
is nondegenerate. Choose a nonzero term of an element with maximum
horizontal PBW degree and multiply by its complement. Terms of smaller
horizontal degree cannot reach \(\Omega\); a different horizontal support
of equal degree cannot reach it either. The remaining central supports
distinguish the selected term. The same argument works on the other side.

Since \(S\) is invertible, \(\ell(uS(v))\) also has rank \(d\).
The rank-one row tensor in (21) therefore has column rank

\[
c_D=d^2.
\tag{22}
\]

## 6. An exact factorization for the column cut

The column rank of an operator on \(H_L\otimes H_L\) is its operator
Schmidt rank: the least number of terms \(U_k\otimes V_k\) in such a
sum. This is exactly scalar matrix rank in the \((x,u)\mid(y,v)\) cut.

Let \(M_\tau\) use only the linear parts of (8)--(9). Let \(K_\tau\)
fix the horizontal variables and, in the coordinates (10), translate

\[
P\mapsto P+\tau\mathcal A,\qquad
Q\mapsto Q+\tau\mathcal B+\tau^2\mathcal C+\tau^3\mathcal A.
\tag{23}
\]

Substitution in the generators proves the commuting factorization

\[
R_\tau^*=M_\tau K_\tau=K_\tau M_\tau.
\tag{24}
\]

### The abelian factor

After \(\tau=1+z\), each horizontal local operator is \(I+zN\),
where \(N\) has rank one and square zero on a four-dimensional space.
Each central local operator is \(I+z^2N\). Let \(b_i\) be the row ranks
of their tensor-product coefficients.

In the local column cut, the coefficient of \(z\) is an invertible
permutation matrix. Multiplying by its inverse turns the constant
coefficient into another rank-one square-zero matrix, similar to \(N\).
Reversing every local polynomial therefore gives exactly

\[
c_i(M_{1+z})=b_{D-i}.
\tag{25}
\]

The \(b_i\) are constructible in \(2^{O(\lambda)}\) time and satisfy
\(b_i\le B_i\) with the Hilbert coefficients (15).

### Exact compression of the translation coefficients

On each leg, let \(D_j=\partial_{z_j}\) act on the square-zero dual
algebra. Horizontal multiplication operators and central derivatives
commute and square to zero. Their monomials

\[
\{x^a y^b D^c:\ a,b,c\text{ binary}\}
\tag{26}
\]

are linearly independent operators: apply them to the product of all
central variables, with no horizontal factor, to distinguish their
horizontal monomials and central complements.

Thus (26) is a \(d\)-dimensional local operator space. Every coefficient
of \(K_{1+z}\) belongs to the tensor product of these two spaces.
Its Schmidt rank is exactly the rank of its \(d\)-by-\(d\) coefficient
matrix in the bases (26). This compresses the full
\(d^2\)-by-\(d^2\) cut without changing rank. Denote these exact ranks
by \(k_a\).

The degree of \(K_\tau\) is at most \(3r\). Expand a product of \(p\)
translations of type \(P\) and \(q\) of type \(Q\). If \(a\) of the
latter supply their cubic \(\mathcal A\) term, then \(p+a\le r\),
because the \(\mathcal A_j\) span at most \(r\) square-zero elements.
The parameter degree is at most \(p+3a+2(q-a)\le3r\).

### A code-dependent symbolic upper bound

For pair \(i\), abbreviate its four horizontal variables by \(X,Y,X',Y'\).
Let \(P_{i,L},P_{i,R}\) be the translations of the original left and
right central coordinates in (23). Put
\(D_{v_i,L}=\sum_jL_{ji}D_{j,L}\), and similarly on the right. Then

\[
K_\tau=
\prod_{i=1}^{m}
\bigl(1+P_{i,L}(\tau)D_{v_i,L}
       +P_{i,R}(\tau)D_{v_i,R}\bigr).
\tag{27}
\]

There is no double-selection term for a pair:
\(P_{i,L}P_{i,R}=0\). In the \(p,q\) variables, this follows from
\((p_xp_y)(p_xq_y)=(p_xp_y)(q_xp_y)=0\).

After \(\tau=1+z\), the coefficient polynomials are:

| Power of \(z\) | Left translation | Right translation | Schmidt-rank bound |
| ---: | --- | --- | ---: |
| 0 | \(XY'+YX'\) | \(XY+X'Y'\) | 2 |
| 1 | \((X+X')Y'\) | \(Y(X+X')\) | 2 |
| 2 | \(X(Y+Y')\) | \(X(Y+Y')\) | 2 |
| 3 | \((X+X')(Y+Y')\) | \((X+X')(Y+Y')\) | 4 |

Products of central derivative forms vanish if their selected code
columns are linearly dependent. Define

\[
S_L(y)=\sum y^{|J_L|+|J_R|},
\tag{28}
\]

where the sum is over disjoint sets \(J_L,J_R\), such that the columns
in each set are linearly independent. Their union need not be independent.
Set

\[
\kappa(z)=S_L(2+2z+2z^2+4z^3).
\tag{29}
\]

Expanding (27) and the rank factorizations in the table proves
\(k_a\le[z^a]\kappa(z)\). Also impose \(k_a=0\) for \(a>3r\).

### Every complete column coefficient

Operator Schmidt rank is submultiplicative under composition and
subadditive under sums. Taking coefficients in (24), and using (25), gives

\[
\boxed{\quad
c_i\le
\min\left\{d^2,\
\sum_{\substack{a+b=i\\0\le a\le3r\\0\le b\le D}}
k_a\,b_{D-b}\right\}.
\quad}
\tag{30}
\]

One may replace \(k_a\) by (29) and \(b_j\) by (16). The proof supplies
actual scalar rank-one factorizations: compose the local factors of the
two operators. This is an upper bound; it does not assert that coefficient
sums have independent images.

## 7. The first coefficient for a single central bit

At \(\tau=1\), the mixing factor is the identity. Remove the invertible
local right-leg translation by its own quadratic form. For one central
bit, the remaining translation polynomial is

\[
(1+\mathcal B_{\rm cross}D_L)(1+Q_LD_R),
\quad
\mathcal B_{\rm cross}=\sum_i(x_{i,L}y_{i,R}+y_{i,L}x_{i,R}),
\quad Q_L=\sum_ix_{i,L}y_{i,L}.
\tag{31}
\]

Its right-hand monomials separate the left spaces spanned by \(1,Q_L\),
by the \(2m\) horizontal linear forms times \(D_L\), and by those forms
times \(Q_LD_L\). Horizontal degrees distinguish the spaces.
Multiplication by \(Q_L\) is injective on degree one when \(m\ge2\):
pair any selected coordinate with a different pair to read its
coefficient in the degree-three image. For \(m=1\), that image is zero.
Therefore

\[
c_0=4m+2\quad(m\ge2),\qquad c_0=4\quad(m=1).
\tag{32}
\]

For arbitrary \(L\), the same expansion gives \(c_0\le2^r I(L)\), where
\(I(L)\) counts independent sets in the two-copy parallel extension from
[coded-central-crossings.md](coded-central-crossings.md).

## 8. Endpoint-preserving local changes

For a polynomial \(a(\tau)\) with \(a(0)=0\), the single-leg substitution

\[
z_j\longmapsto z_j+a(\tau)\sum_iL_{ji}x_i y_i
\tag{33}
\]

is invertible and is the identity at \(\tau=0\). The added expression
is independent of central variables and has square zero. Apply four
such transformations independently, to the two input and two output
legs of \(R_\tau\). The resulting family still has endpoint
\(R_0=\operatorname{Swap}\), hence an exact crossing decomposition
after \(\tau=1+z\).

For bounded-degree choices of \(a\), the degree remains \(O(\lambda)\)
and the construction remains uniformly effective. The choices
\(0,\tau,\tau^2,\tau+\tau^2\) give 256 four-leg cases.
This lemma justifies that search but makes no numerical claim about
its outcome. Parameter-dependent local changes can alter individual
coefficient ranks, despite their invertibility as polynomial families.

## 9. Exact finite checks and reproduction

The independent symbolic checker verifies the complete generator
identities, both orders in (24), the nonlinear change (13), the row
bound, abelian reversal (25), and column convolution bound.

| Horizontal pairs | Central bits | Row ranks | Column ranks |
| ---: | ---: | --- | --- |
| 1 | 1 | \(64,28,18,6,1\) | \(4,12,24,28,64\) |
| 2 | 1 | \(1024,560,348,142,48,10,1\) | \(10,45,136,258,468,560,1024\) |

The maximum products are \(432\) and \(47328\). For the second profile,
(16) gives \(1024,652,349,150,48,10,1\); its middle bound \(349\) is
close to the exact rank \(348\). The compressed translation coefficient
ranks are \(4,4,4,8\) and \(10,10,10,14\), respectively.
These finite computations do not prove asymptotic profile quality.

Run:

    python3 tools/check_noncommutative_crossing_hopf.py \
      --output evidence/noncommutative-crossing-hopf.json

The separate literal Hopf and native substitution implementations are
[explore_graded_hopf_crossings.py](../tools/explore_graded_hopf_crossings.py)
and [explore_graded_hopf_crossings.cpp](../tools/explore_graded_hopf_crossings.cpp).
Their independent comparison is handled by
[check_graded_hopf_crossings.py](../tools/check_graded_hopf_crossings.py).
The symbolic evidence is
[noncommutative-crossing-hopf.json](../evidence/noncommutative-crossing-hopf.json).

### Larger native profiles and complete binary recombinations

The independently audited native implementation also computes these complete
profiles. The entries in the last two columns run from degree zero to degree
\(D\), with no omitted middle coefficients.

| Code columns | Labels | Row ranks | Column ranks |
| --- | ---: | --- | --- |
| \(1,1,1\) | 128 | \(16384,9660,6280,2924,1278,394,94,14,1\) | \(14,98,442,1285,2898,4744,8352,9660,16384\) |
| \(1,2,3\) | 256 | \(65536,40872,29324,16074,9024,3768,1466,422,96,14,1\) | \(52,294,1039,2701,5786,10566,17186,25042,38248,40872,65536\) |

Their maximum products are \(3757340\) and \(52212864\). The evidence
also includes the code \((1,1,2)\), whose maximum product is \(48567168\).
These are finite comparisons, not estimates of an asymptotic exponent.

For the one-center profiles on 32 and 128 labels, every nonzero binary
combination of coefficients was computed: 127 and 511 combinations,
respectively. A linear-span certificate determines the least possible
maximum rank product among binary decompositions of the crossing. It is
\(47328\) and \(3757340\), respectively, so recombination does not
improve the bottleneck in these two cases.

The minimum total cost \(\sum\sqrt{r_i c_i}\) can improve. A shortest-path
calculation on the coefficient vector space, with integer lower and upper
bounds for every square root, gives

\[
\begin{aligned}
662.123715363140&\le S_{32}\le662.123715363142,\\
7356.777950072195&\le S_{128}\le7356.777950072198.
\end{aligned}
\]

Both lower bounds exceed \(d^{3/2}\), hence also \(d^{\sqrt2}\).
Thus no binary recombination within either complete coefficient span
satisfies the finite seed criterion. The argument does not exclude
other fields or other crossing tensors.

All 256 gauges in Section 8 were also evaluated on both profiles. The
ungauged family is the unique minimum of the maximum rank product in
each sweep. Every candidate's maximum product exceeds \(d^3\), which
already rules out \(S<d^{\sqrt2}\) without a floating-point comparison.
Generic polynomial composition independently checks all 256 eight-label
gauges and nine selected 32-label gauges, including complete tensor
digests and both cut ranks.

The records are in
[graded-hopf-complete-profiles.json](../evidence/graded-hopf-complete-profiles.json),
[graded-hopf-native-audit.json](../evidence/graded-hopf-native-audit.json), and
[graded-hopf-gauge-audit.json](../evidence/graded-hopf-gauge-audit.json).
To reproduce the native comparisons, run:

    c++ -O3 -std=c++17 tools/explore_graded_hopf_crossings.cpp -o /tmp/hopf-screen
    python3 tools/check_graded_hopf_crossings.py --binary /tmp/hopf-screen
    /tmp/hopf-screen --horizontal 2 --center 1 --columns 1,1 --span > /tmp/hopf-span.jsonl
    python3 tools/summarize_hopf_coefficient_span.py /tmp/hopf-span.jsonl --output /tmp/hopf-span.json
    c++ -O3 -std=c++17 tools/explore_graded_hopf_gauges.cpp -o /tmp/hopf-gauges
    python3 tools/check_graded_hopf_gauges.py --binary /tmp/hopf-gauges --all-small

The native screen checks coefficient independence, the complete crossing
endpoint, and conservation of the exact block gradings. The span summarizer
rejects incomplete or overlapping chunks and requires matching full-tensor
digests before combining them.
