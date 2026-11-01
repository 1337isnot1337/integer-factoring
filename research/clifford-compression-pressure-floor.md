# Low-order Smith constraints for linear elementary restrictions

## Result and scope

Let `F` be any field of characteristic two and let `n=2^t`. Restrict
the tensor product of `t` binary elementary crossings by constant
maps `P:F^q -> F^n` and `L:F^n -> F^q` with `LP=I`. The
[linear restriction constructor](linear-physical-restriction-crossings.md)
allows arbitrary left inverses, including nonorthogonal ones.

If `q>n/2`, the column cut at zero has exactly one Smith order zero
and exactly `2t` Smith orders one. If `q>3n/4`, its order-two
multiplicity is also forced to be `binom(2t,2)`. These statements
hold for every such restriction, without enumerating subspaces or
left inverses. They exclude a half-rate improvement for all four-site
restrictions with `q=9,...,16`, including cases left open by the
general integer shared-minor relaxation.

A further commutator argument supplies higher-order constraints in
both cuts when `q` is closer to `n`. Exact all-scale certificates cover
57 site/dimension pairs through seven sites. There is also an
asymptotic consequence: restrictions of codimension
`n^(1/2+o(1))` retain a pressure exponent of at least
`h_2(1/4)=0.811278...`. The full-dimensional cases are already tensor
powers of the elementary crossing up to constant basis changes.

These are restrictions on the bare local Smith pressure of this family
in characteristic two. They do not cover nonconstant physical frames, additional
gauges, adaptive channels, other scalar seeds, or general crossings.
The complete supported factoring exponent is unchanged.

## 1. The first Clifford layer cannot disappear

On the ambient binary tensor space `V=(F^2)^(tensor t)`, let `e_i`
and `f_i` be the single-site matrix units `E_01` and `E_10`. They obey

$$e_i^2=f_i^2=0,\qquad e_i f_i+f_i e_i=I,$$

and operators at different sites commute. Every operator

$$A=cI+\sum_{i=1}^t(a_i e_i+b_i f_i) \tag{1}$$

therefore has `A^2=lambda I` for a scalar `lambda` in `F`.
If `lambda!=0`, it is invertible. If `lambda=0` and `A!=0`, choose
a scalar multiple `B` of one of the `e_i,f_i` for which `AB+BA=I`. Such a partner
exists unless `A` is scalar, and a nonzero scalar does not have square
zero. For `v` in `ker A`,

$$v=(AB+BA)v=A(Bv).$$

Hence `ker A=im A`, and both spaces have dimension `n/2`.

Write `W=im P` and `U=ker L`, so `V=W direct-sum U`. Suppose
`LAP=0`. Then `A(W)` is contained in `U` and is disjoint from `W`.
An invertible `A` would give `q=dim A(W)<=dim U=n-q`.
For a nonzero square-zero `A`, the two spaces

$$W\cap\ker A,\qquad A(W)$$

are disjoint subspaces of `ker A`. Their dimensions sum to `q` by
rank-nullity on `A|_W`. Thus this case also gives `q<=n/2`.

Consequently, when `q>n/2`, compression `A -> LAP` is injective on
the entire `(2t+1)`-dimensional space in (1). In particular the
compressed operators

$$I_q,\ Le_1P,Lf_1P,\ldots,Le_tP,Lf_tP \tag{2}$$

are linearly independent. This proof uses no inner product.

## 2. The column tangent is a nondegenerate alternating form

For one binary site write `v=vec(I)`, `h=vec(E_00)`,
`e=vec(E_01)`, and `f=vec(E_10)`. The column matrix of
`T(z)=(1-z)C+zE` is, in characteristic two,

$$K(z)=vv^{\mathsf T}
 +z(vh^{\mathsf T}+hv^{\mathsf T}+ef^{\mathsf T}+fe^{\mathsf T}). \tag{3}$$

Vectorization of `A -> LAP` is the same constant map on both sides
of this Gram matrix. After tensoring (3), its constant term is
`vec(I_q) vec(I_q)^T`. Modulo that one-dimensional unit direction,
the coefficient of `z` is the alternating Gram form pairing the
compressed `e_i` with the compressed `f_i`.

The vectors (2) are independent, so this quotient Gram form has
rank exactly `2t`. Eliminating the single unit pivot therefore
leaves a matrix whose first nonzero coefficient has rank `2t`.
This proves the exact multiplicities

$$\boxed{m_0=1,\qquad m_1=2t\quad\text{if }q>2^{t-1}.} \tag{4}$$

The dimension threshold cannot simply be removed. A coordinate
subcube of dimension `q=2^(t-1)` retains only `t-1` active sites;
its order-one multiplicity is `2(t-1)`.

The actual polynomial degree is `t` in the regime of (4). The leading
ambient row coefficient is the outer product obtained from
`J=[[0,1],[1,0]]^(tensor t)`. Its compressed factors are
`vec(L J L^T)` and `vec(P^T J P)`. A subspace of dimension greater
than `n/2` cannot be totally isotropic for this nondegenerate
alternating form. Both factors are nonzero, so the degree does not drop.

## 3. A higher-layer rank bound

Let `F_d` denote the filtration space spanned by ordered squarefree monomials of degree at most
`d` in `e_1,f_1,...,e_t,f_t`. This is a filtration of the full matrix
algebra, with

$$\dim F_d=\sum_{j=0}^d\binom{2t}{j}\qquad(0<=d<=2t).$$

For a nonzero `A` in `F_d`, choose a nonzero monomial coefficient of
largest degree `k<=d`. Commuting with the partner of each generator
in that monomial removes that generator. The resulting `k` iterated
commutators equal a nonzero scalar multiple of `I`: no other monomial
of degree at most `k` contains the chosen full set. Since

$$\mathrm{rank}(AB+BA)\le2\,\mathrm{rank}(A),$$

we obtain the elementary rank bound

$$\mathrm{rank}(A)\ge n/2^k\ge n/2^d. \tag{5}$$

If `LAP=0`, decompose the domain as `W direct-sum U`. Its image
from `W` has dimension at most `n-q`; the image from `U` also has
dimension at most `n-q`. Thus `rank A<=2(n-q)`. Combining with
(5) shows that compression is injective on `F_d` whenever

$$2(n-q)<n/2^d. \tag{6}$$

This is a sufficient condition. The stronger first-layer conclusion
(4) uses the additional square-zero structure and does not follow
from (6) alone.

To translate injection into Smith orders, rewrite (3) exactly as

$$K(z)=(v+zh)(v+zh)^{\mathsf T}
       +z(ef^{\mathsf T}+fe^{\mathsf T})+z^2hh^{\mathsf T}. \tag{7}$$

This is a Gram factorization with local weights `0,1,1,2`.
Tensoring gives a polynomial frame whose constant vectors are the
ordered Clifford monomials. At weight `j` there are `binom(2t,j)`
vectors; the Gram pairing is a permutation within each equal-weight
block. If compression is injective on `F_d`, those frame columns
of weights at most `d` remain independent modulo `z` and can be
completed to a basis over `F[[z]]`. All remaining terms have a
factor `z^(d+1)` and do not affect the lower orders. Therefore

$$\boxed{m_j=\binom{2t}{j}\quad(0<=j<=d)\quad
\text{under condition (6).}} \tag{8}$$

The implementation keeps the exact Gram expansion (7), including
every coefficient cancellation; it does not replace the entries by
generic noncommutative symbols.

### A stronger second-layer theorem

For an operator `A` of filtration degree at most two, every first
commutator `[A,v]`, `v` in the first filtration layer, belongs
to `F_1`. We claim

$$A^2=0\quad\text{or}\quad\mathrm{rank}(A)\ge n/2. \tag{8a}$$

If some `[A,v]` is invertible, the commutator rank bound proves
(8a). Otherwise every first commutator has square zero by Section 1.
Their image is a linear space, so any two such commutators commute:
expand the square of their sum. The Jacobi identity then gives

$$0=[[A,v],[A,w]]=[v,[A^2,w]].$$

Here the discarded Jacobi term is `[A,[v,[A,w]]]`, which vanishes
because commutators of two first-layer elements are scalar.
Thus `[A^2,w]` commutes with every first-layer element. These
elements generate the full matrix algebra, whose center is scalar.
But `[A^2,w]=[A,[A,w]]` is itself a first commutator and has square
zero. Being scalar, it is zero. It follows that `A^2` is scalar.
A nonzero scalar square makes `A` invertible; a zero square is the
other branch of (8a).

Now suppose `LAP=0`, `A!=0`, and `q>3n/4`. The rank-at-least-`n/2`
branch contradicts `rank A<=2(n-q)`. In the square-zero branch,
(5) gives `rank A>=n/4`. The disjoint spaces `W intersect ker A`
and `A(W)` lie in `ker A` and have total dimension `q`, giving
`q<=n-rank A<=3n/4`, again a contradiction. Consequently

$$\boxed{q>3n/4\quad\Longrightarrow\quad
 m_0=1,\ m_1=2t,\ m_2=\binom{2t}{2}.} \tag{8b}$$

This threshold is sharp for injection on the second layer. If
`t>=2`, take `W=ker(e_1e_2)`, which has dimension `3n/4`, and any
retraction onto it. The nonzero degree-two operator `e_1e_2` has
zero compression. A coordinate version also loses two second-order
Smith factors, as checked in the evidence.

The quadratic argument uses the fact that a double first-layer
commutator is scalar. For degree greater than two, its Jacobi
expansion has an additional term that this proof does not eliminate.
No higher-degree version of (8a) is asserted. Higher layers use the
proved threshold (6).

### The same generic threshold resolves the row reversal

There is a separate regular frame for the reversed row cut. Put
`j=vec(J)`, `h=vec(E_01)`, `a=vec(E_00)`, and `b=vec(E_11)`, with
`J=[[0,1],[1,0]]`. The one-site row reversal is `jj^T+zI`, and
the exact identity in `F[[z]]` is

$$\begin{aligned}
jj^{\mathsf T}+zI
={}&(1+z)\left(j+\frac z{1+z}h\right)
             \left(j+\frac z{1+z}h\right)^{\mathsf T}\\
 &+zaa^{\mathsf T}+zbb^{\mathsf T}
   +\frac{z^2}{1+z}hh^{\mathsf T}. \tag{8c}
\end{aligned}$$

The four weights are again `0,1,1,2`; their scalar factors are
units. At zero, the four frame vectors are `vec(AJ)` for
`A=I,e,f,ef`. Tensoring (8c) gives the same Clifford filtration.
Write `J_t=J^(tensor t)`. On the two sides of the restricted row
matrix, the constant frame maps are

$$A\longmapsto L A J_t L^{\mathsf T},\qquad
  A\longmapsto P^{\mathsf T} A J_t P.$$

Each map has the following property: if its value is zero, `A`
maps a `q`-dimensional subspace into an `(n-q)`-dimensional one.
Hence `rank A<=2(n-q)`, exactly as in Section 3. Under (6), both
maps are injective on `F_d`. Complete their low-weight columns to
separate regular bases on the two sides. The high-weight terms
again carry `z^(d+1)`. Thus both cuts have the forced multiplicities

$$\boxed{m_j^{\rm row}=m_j^{\rm column}=\binom{2t}{j}
\quad(0<=j<=d)\quad\text{under (6).}} \tag{8d}$$

The unit `1+z` in (8c) need not have a square root. The proof uses
the displayed unit factors and regular left and right bases directly.

## 4. All-scale pressure consequences

Put `N0=q^2`, `r=1/s` for `s>1`, and let `mu` be the row-infinity
orders and `nu` the column-zero orders. Both total sums are at most
`N0*t`, and the shared minor gives

$$\sum_{i=1}^q(\mu_i+\nu_i)\le qt. \tag{9}$$

Use any proved initial multiplicities `c_0,...,c_d` from (4), (8), or (8b).
Set `A=sum c_j` and `B=sum j*c_j`. The remaining column orders have
total budget at most `N0*t-B`. The balanced integer polynomial

$$J_m(a;r)=(m-v)r^u+v r^{u+1},\qquad a=mu+v,\quad0<=v<m,$$

therefore gives

$$C(r)=\sum_{j=0}^d c_jr^j+J_{N0-A}(N0t-B;r)
       \le\sum_i r^{\nu_i}. \tag{10}$$

If no orders remain, the last term is zero. Lower-bound the first
`q` column orders by filling from these proved multiplicities and
putting every remaining required order at `d+1`; call the resulting
integer `c`. Then the row prefix budget in (9) is at most `b=qt-c`.
For the cases certified here, `c>=t`, so `b<=(q-1)t`. On this range
the balanced row polynomial is nonincreasing as its prefix budget
increases. Consequently

$$R(r)=1+J_{q-1}(b;r)+J_{N0-q}(N0t-b;r)
       \le\sum_i r^{\mu_i}. \tag{11}$$

This monotonicity follows directly from the difference formula in
the [integer shared-minor proof](integer-shared-minor-pressure-floor.md):
the prefix average is at most `t`, and the remaining average is at
least `t`.

The final lower bound using the shared minor is

$$B(s)\ge r^{-t}R(r)C(r). \tag{12}$$

When (8d) supplies initial row multiplicities, applying the same
balanced-sum formula (10) to the row profile gives another valid
choice of `R(r)`. Each exact certificate records which row lower
bound it uses. No assumption about a generic completion of the
unresolved profile is made.

Exact rational Bernstein certificates show that (12) is strictly
greater than `q^(577/204)` in the following cases. Hence
`rho>577/816>1/sqrt(2)` throughout the full scale range.

| Sites `t` | Physical dimensions `q` | Available information used |
| --- | --- | --- |
| 1 | 2 | Orders zero and one |
| 2 | 3,4 | First column layer; complete profiles at full dimension |
| 3 | 5,6,7,8 | First/second column layers and generic row layers |
| 4 | 9 through 16 | First/second column layers and generic row layers |
| 5 | 25 through 32 | Second column layer and generic layers in both cuts |
| 6 | 49 through 64 | Second column layer and generic layers in both cuts |
| 7 | 97,98; 113 through 128 | Second column layer and generic layers in both cuts |

The stated relaxation gives no exclusion for the omitted dimensions.
The theorem does not assert that any omitted case has an improving
crossing.

## 5. An asymptotic barrier near full dimension

Consider a sequence with `t` tending to infinity, `n=2^t`, and
positive integral codimension

$$c=n-q=n^{\beta+o(1)},\qquad0<=\beta<1.$$

Take `d=floor(t-log_2 c)-2`, which satisfies the strict inequality
(6) for all sufficiently large `t`, and put
`j=min(d,floor(t/2))`. One forced level from (8d) in each cut gives,
uniformly over every scale `s>=1`,

$$B(s)\ge\binom{2t}{j}^{\!2}s^{t-2j}
       \ge\binom{2t}{j}^{\!2}. \tag{13}$$

Since `log_2 q=t+o(t)`, the elementary binomial asymptotic yields

$$\boxed{\liminf\rho\ge
 h_2\!\left(\min\left\{\frac{1-\beta}{2},\frac14\right\}\right),}
 \tag{14}$$

where `h_2(a)=-a log_2 a-(1-a)log_2(1-a)` is binary entropy.
In particular, codimension at most `n^(1/2+o(1))` retains the floor

$$h_2(1/4)=2-\frac34\log_2 3=0.811278\ldots.$$

This is an asymptotic statement for the specified sequence of
constant restrictions, uniform in the optimizing scale. It does
not exclude a finite exceptional seed outside the certified table.
For `c=0`, the physical maps are invertible and both complete
profiles equal the elementary tensor-power profile; the same
floor holds exactly.

## 6. Further exact invariants

The column matrix is symmetric and has constant diagonal
`K_(x,u),(x,u)=delta_xu`. To see the diagonal identity, pair the
distinct ambient input indices in the contraction of any partial
swap against equal output covectors and equal input vectors. Those
terms cancel in characteristic two. The surviving equal indices
give `(LP)_(x,u)^2`, independently of the partial swap.

Thus every diagonal entry remains a square under congruence over
`F[[z]]`. In a minimum-valuation elimination step, a diagonal
pivot has even valuation. An odd minimum must occur off the
diagonal; its two-by-two symmetric block has two equal Smith
orders. Removing that block preserves square diagonals. Induction
proves that every odd column Smith order has even multiplicity.

The row matrix commutes with `N=swap+I` and induces the identity on
`ker N/im N`, by the same diagonal identity. In a basis ordered as
`im N`, its diagonal quotient, and preimages of `im N`, it has form

$$\begin{pmatrix}A&B&C\\0&I&E\\0&0&A\end{pmatrix}.$$

Its determinant is therefore `(det A)^2`. In particular the sum of
its degree-`t` reversed Smith orders is congruent to `tq^2` modulo
two. These invariants are useful for checking candidates; they are
not needed by the pressure certificates above.

## 7. Replay

Run:

```sh
python3 tools/check_clifford_compression.py
```

The [checker](../tools/check_clifford_compression.py) checks the exact
Clifford identities and rank bound, tests arbitrary graph retractions,
compares the Gram jets with the full tensor constructor, and resolves
the proved orders by independent block-Toeplitz ranks. The exact code
works over `F_2`; the characteristic-two generality follows from the
proofs above. The checker also records the half- and three-quarter-
dimension countercontrols and all-scale rational certificates.
The [evidence](../evidence/clifford-compression.json) includes the
relevant source hashes. No floating optimization enters certification.
